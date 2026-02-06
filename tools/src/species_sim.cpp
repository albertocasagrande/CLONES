/**
 * @file mutants_sim.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the mutants simulator
 * @version 1.1
 * @date 2026-02-06
 *
 * @copyright Copyright (c) 2023-2026
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>
#include <vector>
#include <chrono>
#include <csignal>
#include <fstream>
#include <regex>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

#include "common.hpp"
#include "simulation.hpp"
#include "ending_conditions.hpp"
#include "progress_bar.hpp"

#ifdef WITH_SDL2
#include "SDL_plot.hpp"
#endif


CLONES::Mutants::Evolutions::Simulation simulation;
CLONES::UI::ProgressBar *bar;

void termination_handling(int signal_num)
{
    simulation.make_snapshot(bar);

    if (bar != nullptr) {
        bar->set_message("Aborted");

        delete bar;
    }

    CLONES::UI::ProgressBar::show_console_cursor();
    exit(signal_num);
}

class DriverSimulator : public BasicExecutable
{
    std::filesystem::path simulation_filename;
    std::filesystem::path logging_dir;
    std::string snapshot_path;
    long double time_horizon;
#ifdef WITH_SDL2
    unsigned int frames_per_second;
#endif
    bool plot;
    bool quiet;
    bool disable_storage;
    bool duplicate_internal_cells;

    static double get_rate(const nlohmann::json& rate_json)
    {
        double rate = rate_json.template get<double>();

        if (rate < 0) {
            throw std::domain_error("the rate must be non-negative");
        }

        return rate;
    }

    static double get_epigenetic_rate(const nlohmann::json& rate_json,
                                      const std::vector<std::string>& field_aliases)
    {
        for (const auto& field_name : field_aliases) {
            if (rate_json.contains(field_name)) {
                return get_rate(rate_json[field_name]);
            }
        }

        std::ostringstream oss;
        auto alias_it = std::begin(field_aliases);

        oss << "Every \"epigenetic rates\" element must contain a \""
            << *alias_it;
        while (++alias_it != std::end(field_aliases)) {
            oss << "/" << *alias_it;
        }
        oss << "\" field";

        throw std::domain_error(oss.str());
    }

    static CLONES::Mutants::MutantProperties create_mutant(const nlohmann::json& mutant_json)
    {
        using namespace CLONES::Mutants;
        std::vector<EpigeneticRates> epigenetic_rates;

        if (!mutant_json.is_object()) {
            throw std::domain_error("The mutant specification must be an object");
        }

        if (!mutant_json.contains("epigenetic rates")) {
            throw std::domain_error("The mutant specification must contain a \"epigenetic rates\" field");
        }

        if (!mutant_json["epigenetic rates"].is_array()) {
            throw std::domain_error("The \"epigenetic rates\" field must be an array of epigenetic rates");
        }

        if (!mutant_json.contains("epigenetic types")) {
            throw std::domain_error("The mutant specification must contain a \"epigenetic types\" field");
        }

        if (!mutant_json["epigenetic types"].is_array()) {
            throw std::domain_error("The \"epigenetic types\" field must be an array of epigenetic types");
        }

        if (static_cast<size_t>(1<<mutant_json["epigenetic rates"].size())!=mutant_json["epigenetic types"].size()) {
            throw std::domain_error("The epigenetic types must be exponential in number with respect "
                                    "to the number of promotor epigenetic rates");
        }

        for (const auto& rate_json : mutant_json["epigenetic rates"]) {

            if (!rate_json.is_object()) {
                throw std::domain_error("Every \"epigenetic rates\" must be an object");
            }

            const auto methylation_rate = get_epigenetic_rate(rate_json,
                                                              {"methylation", "on"});
            const auto demethylation_rate = get_epigenetic_rate(rate_json,
                                                                {"demethylation", "off"});
            epigenetic_rates.push_back({methylation_rate, demethylation_rate});
        }

        if (!mutant_json.contains("name")) {
            throw std::domain_error("The mutant specification must contain a \"name\" field");
        }

        MutantProperties mutant(mutant_json["name"].template get<std::string>(), epigenetic_rates);

        for (const auto& type_json : mutant_json["epigenetic types"]) {
            auto& type = mutant[type_json["status"].template get<std::string>()];
            type.set_rate(CellEventType::DEATH, get_rate(type_json["death rate"]));
            type.set_rate(CellEventType::DUPLICATION, get_rate(type_json["duplication rate"]));
        }

        return mutant;
    }

    static void configure_tissue(const nlohmann::json& tissue_json)
    {
        using namespace CLONES::Mutants::Evolutions;

        if (!tissue_json.contains("size")) {
            throw std::domain_error("The tissue specification must contain a \"size\" field");
        }

        if (!tissue_json["size"].is_array()) {
            throw std::domain_error("The tissue specification size must be an array");
        }

        std::vector<AxisSize> sizes;
        for (const auto& value : tissue_json["size"]) {
            sizes.push_back(value.template get<AxisSize>());
        }
        if (sizes.size()==2 || sizes.size()==3) {
            std::string name;
            if (tissue_json.contains("name")) {
                name = tissue_json["name"].template get<std::string>();
            }
            simulation.set_tissue(name, sizes);

            return;
        }

        throw std::domain_error("tissue is either a 2D or a 3D space");
    }

    static CLONES::Mutants::Evolutions::PositionInTissue
    get_position(const nlohmann::json& position_json)
    {
        using namespace CLONES::Mutants::Evolutions;

        if (!position_json.is_array()) {
            throw std::domain_error("The \"position\" field must be an array of Natural values");
        }
        std::vector<AxisPosition> position;
        for (const auto& value : position_json) {
            position.push_back(value.template get<AxisPosition>());
        }

        if (position.size()==2) {
            return {position[0], position[1]};
        }

        if (position.size()==3) {
            return {position[0], position[1], position[2]};
        }

        throw std::domain_error("tissue may be either a 2D or a 3D space");
    }

    static void configure_initial_cells(const nlohmann::json& initial_cells_json,
                                        const std::map<std::string, CLONES::Mutants::MutantProperties>& mutants)
    {
        if (!initial_cells_json.is_array()) {
            throw std::domain_error("The \"initial cells\" field must contain an array");
        }
        for (const auto& initial_cell_json : initial_cells_json) {

            if (!initial_cell_json.is_object()) {
                throw std::domain_error("Every element in the \"initial cells\" field must be an object");
            }

            if (!initial_cell_json.contains("mutant")) {
                throw std::domain_error("Every initial cell must contain a \"mutant\" field");
            }

            if (!initial_cell_json.contains("mutant")) {
                throw std::domain_error("Every initial cell must contain a \"mutant\" field");
            }
            auto& mutant_json = initial_cell_json["mutant"];

            if (!mutant_json.contains("name")) {
                throw std::domain_error("Every initial cell \"genotoype\" field must contain "
                                         "a \"name\" field");
            }

            std::string mutant_name = mutant_json["name"].template get<std::string>();

            if (!mutant_json.contains("epistate")) {
                throw std::domain_error("Every initial cell \"genotoype\" field must contain "
                                         "a \"epigenetic status\" field");
            }

            std::string epigenetic_status = mutant_json["epistate"].template get<std::string>();

            auto& species = mutants.at(mutant_name)[epigenetic_status];

            if (!initial_cell_json.contains("mutant")) {
                throw std::domain_error("Every initial cell must contain a \"position\" field");
            }
            simulation.place_cell(species, get_position(initial_cell_json["position"]));
        }
    }

    static void configure_timed_events(const nlohmann::json& timed_events_json,
                                       const std::map<std::string, CLONES::Mutants::MutantProperties>& mutants)
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;

        if (!timed_events_json.is_array()) {
            throw std::domain_error("The \"timed events\" field must contain an array");
        }

        for (const auto& timed_event_json : timed_events_json) {
            if (!timed_event_json.is_object()) {
                throw std::domain_error("Every element in the \"timed events\" field must be an object");
            }

            Evolutions::TimedEvent timed_event = ConfigReader::get_timed_event(simulation, mutants,
                                                                               timed_event_json);

            simulation.schedule_timed_event(timed_event);
        }
    }

    static void configure_simulation(const std::string& simulation_filename)
    {
        using namespace CLONES::Mutants;

        std::ifstream f(simulation_filename);
        nlohmann::json simulation_cfg = nlohmann::json::parse(f);

        std::map<std::string, MutantProperties> mutants;

        if (!simulation_cfg.contains("tissue")) {
            throw std::domain_error("The simulation configuration file must contain a \"tissue\" field");
        }
        configure_tissue(simulation_cfg["tissue"]);

        if (!simulation_cfg.contains("mutants")) {
            throw std::domain_error("The simulation configuration file must contain a \"mutants\" field");
        }
        if (!simulation_cfg["mutants"].is_array()) {
            throw std::domain_error("The \"mutants\" field must be an array of mutant specifications");
        }
        for (const auto& mutant_json: simulation_cfg["mutants"]) {
            auto mutant = create_mutant(mutant_json);
            simulation.add_mutant(mutant);

            std::string mutant_name = mutant.get_name();
            mutants.emplace(mutant_name, std::move(mutant));
        }

        if (!simulation_cfg.contains("initial cells")) {
            throw std::domain_error("The simulation configuration file must contain a \"initial cells\" field");
        }
        configure_initial_cells(simulation_cfg["initial cells"], mutants);

        if (simulation_cfg.contains("timed events")) {
            configure_timed_events(simulation_cfg["timed events"], mutants);
        }

        if (simulation_cfg.contains("death activation level")) {
            simulation.death_activation_level = simulation_cfg["death activation level"].template get<size_t>();
        }

        using namespace std::chrono_literals;

        simulation.set_interval_between_snapshots(5min);
    }

    void init_simulation()
    {
        using namespace CLONES;

        configure_simulation(simulation_filename);

        simulation.duplicate_internal_cells = duplicate_internal_cells;
    }

public:

    DriverSimulator(int argc, char* argv[]):
        BasicExecutable(argv[0], {{"options", "Options"}})
    {
        namespace po = boost::program_options;

        visible_options.at("options").add_options()
            ("recover-simulation,r",
             po::value<std::string>(&snapshot_path)->default_value(""),
             "recover a simulation")
            ("uniform-grow-model,u", "admit duplications of any cell in the tissue")
            ("log directory,o", po::value<std::filesystem::path>(&logging_dir),
             "the directory in which all the simulation data are saved")
            ("disable-storage,d", "disable result storage")
            ("seed,s", po::value<int>()->default_value(0), "random generator seed")
    #if WITH_INDICATORS
            ("quiet,q", "disable progress bar")
    #endif // WITH_INDICATORS
    #ifdef WITH_SDL2
            ("plot,p",
             "plot a graphical representation of the simulation")
            ("frames-per-second,f", po::value<unsigned int>(&frames_per_second)->default_value(1),
             "the number of frames per second")
    #endif
            ("help,h", "get the help");

        hidden_options.add_options()
            ("simulation file", po::value<std::filesystem::path>(&simulation_filename),
             "the name of the file describing the simulation")
            ("time horizon", po::value<long double>(&time_horizon),
             "the simulation time horizon");

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("simulation file", 1);
        positional_options.add("time horizon", 1);

        po::variables_map vm;
        try {
            po::command_line_parser parser{argc, argv};
            po::store(parser.options(program_options).positional(positional_options).run(), vm);
            po::notify(vm);
        } catch (boost::wrapexcept<boost::program_options::unknown_option> &except) {
            print_help_and_exit(except.what(), 1);
        }

        if (vm.count("help")) {
            print_help_and_exit(0);
        }

        if (!vm.count("simulation file")) {
            print_help_and_exit("The simulation file is mandatory", 1);
        }

        if (!vm.count("time horizon")) {
            print_help_and_exit("The simulation time horizon is mandatory", 1);
        }

        if (vm.count("log directory")) {
            simulation.rename_log_directory(logging_dir);
        }

        quiet = vm.count("quiet")>0;
        plot = vm.count("plot")>0;
        disable_storage = vm.count("disable-storage")>0;
        duplicate_internal_cells = (vm.count("border-duplications-only")==0);
    }

    void run()
    {
        using namespace CLONES;
        using namespace CLONES::Mutants::Evolutions;

        if (!quiet) {
            UI::ProgressBar::hide_console_cursor();

            bar = new UI::ProgressBar(std::cout);
            bar->set_message("Creating tissue");
        }

        if (snapshot_path != "") {
            snapshot_path = get_last_snapshot_path(snapshot_path, "simulation");

            simulation = load_species_simulation(snapshot_path, quiet);
        } else {
            init_simulation();
        }

        simulation.storage_enabled = !disable_storage;

        TimeTest time_test(time_horizon);

#ifdef WITH_SDL2
        if (plot) {

            UI::TissuePlotter<UI::SDLWindow> plotter(simulation.tissue());

            plotter.set_frames_per_second(frames_per_second);

            if (bar != nullptr) {
                simulation.run(time_test, plotter, *bar);
            } else {
                simulation.run(time_test, plotter);
            }

            while (plotter.waiting_end()) {
                plotter.plot(simulation.get_statistics());
            }
        } else {
#endif // WITH_SDL2

            if (bar != nullptr) {
                simulation.run(time_test, *bar);
            } else {
                simulation.run(time_test);
            }

#ifdef WITH_SDL2
        }
#endif // WITH_SDL2

        UI::ProgressBar::show_console_cursor();

        if (bar != nullptr) {
            delete bar;
        }
    }

};

int main(int argc, char* argv[])
{
    std::signal(SIGINT, termination_handling);

    DriverSimulator simulator(argc, argv);

    simulator.run();

    return 0;
}
