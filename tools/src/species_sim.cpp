/**
 * @file mutants_sim.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the mutants simulator
 * @version 1.4
 * @date 2026-07-14
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
#include "tissue_simulation.hpp"
#include "ending_conditions.hpp"
#include "progress_bar.hpp"
#include "error.hpp"

#ifdef WITH_SDL2
#include "SDL_plot.hpp"
#endif


CLONES::Mutants::Evolutions::TissueSimulation simulation;
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
    double time_horizon;
#ifdef WITH_SDL2
    uint8_t frames_per_second;
#endif
    bool plot;
    bool quiet;
    bool disable_storage;
    bool duplicate_internal_cells;

    struct Epistate
    {
        std::string name;
        double death_rate;
        double duplication_rate;

        Epistate():
            name{""}, death_rate{0.0}, duplication_rate{0.0}
        {}

        explicit Epistate(const nlohmann::json& epistate_json):
            death_rate{0.0}, duplication_rate{0.0}
        {
            using namespace CLONES;

            if (!epistate_json.is_object()) {
                throw Error<std::runtime_error>("Every epigenetic state must be an object.");
            }

            if (!epistate_json.contains("name")) {
                throw Error<std::runtime_error>("Every epigenetic state must contain a "
                                               "\"name\" field.");
            }

            name = epistate_json["name"].template get<std::string>();

            if (epistate_json.contains("death rate")) {
                death_rate = epistate_json["death rate"].template get<double>();
            }

            if (epistate_json.contains("duplication rate")) {
                duplication_rate = epistate_json["duplication rate"].template get<double>();
            }
        }
    };

    static double get_rate(const nlohmann::json& rate_json)
    {
        using namespace CLONES;

        double rate = rate_json.template get<double>();

        if (rate < 0) {
            throw Error<std::runtime_error>("The rate must be non-negative. Got "
                                           + std::to_string(rate) + ".");
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
        oss << "\" field.";

        using namespace CLONES;
        throw Error<std::runtime_error>(oss.str());
    }

    static std::list<Epistate> get_epistate_names(const nlohmann::json& epistates_json)
    {
        if (!epistates_json.is_array()) {
            using namespace CLONES;

            throw Error<std::runtime_error>("The \"epigenetic states\" field must be "
                                           "an array of epigenetic states");
        }

        std::list<Epistate> epistates;
        for (const auto& state_json : epistates_json) {
            epistates.push_back(Epistate(state_json));
        }

        return epistates;
    }

    static void set_switch_rate(CLONES::Mutants::Evolutions::Species& species,
                                const nlohmann::json& switch_json)
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;

        if (!switch_json.is_object()) {
            throw Error<std::runtime_error>(("Any epigenetic switch specification must be an object. "
                                             "This is not always the case for species \"")
                                            + species.get_name() + "\".");
        }

        if (!switch_json.contains("source")) {
            throw Error<std::runtime_error>(("Any epigenetic switch specification must contain a "
                                            "\"source\" field. This is not always the case "
                                            "for species \"") + species.get_name() + "\".");
        }

        const std::string src_epistate_name = switch_json["source"].template get<std::string>();

        if (!switch_json.contains("destination")) {
            throw Error<std::runtime_error>(("Any epigenetic switch specification must contain a "
                                            "\"destination\" field. This is not always the case "
                                            "for species \"") + species.get_name() + "\".");
        }

        const std::string dst_epistate_name = switch_json["destination"].template get<std::string>();

        if (!switch_json.contains("rate")) {
            throw Error<std::runtime_error>(("Any epigenetic switch specification must contain a "
                                            "\"rate\" field. This is not always the case for "
                                            "species \"") + species.get_name() + "\".");
        }

        const double rate = switch_json["rate"].template get<double>();

        SpeciesName dst_name{species.get_mutant_name(), dst_epistate_name};

        auto& dst_species = simulation.tissue().get_species(dst_name);

        species.set_rate(CellEventType::DUP_AND_EPI_SWITCH, dst_species, rate);
    }

    static void set_species_rates(const std::string& mutant_name, const nlohmann::json& species_json)
    {
        using namespace CLONES;

        if (!species_json.contains("epigenetic states")) {
            throw Error<std::runtime_error>(("Any mutant specification must contain a "
                                            "\"epigenetic states\" field. This is not the case "
                                            "for mutant \"") + mutant_name +"\".");
        }

        std::list<Epistate> epistates = get_epistate_names(species_json["epigenetic states"]);
        std::set<std::string> epistate_names;

        auto& tissue = simulation.tissue();

        for (const Epistate& epistate : epistates) {
            if (epistate_names.count(epistate.name)>0) {
                throw Error<std::runtime_error>("The epistate \"" + epistate.name + "\" has "
                                               + "multiple definitions in mutant \""
                                               + mutant_name + "\".");
            }
            epistate_names.insert(epistate.name);

            if (tissue.get_epigenetic_state_names().count(epistate.name)==0) {
                tissue.add_epigenetic_state(epistate.name);
            }

            const CLONES::Mutants::SpeciesName species_name(mutant_name, epistate.name);

            auto& species = tissue.get_species(species_name);

            using namespace CLONES::Mutants;
            species.set_rate(CellEventType::DEATH, epistate.death_rate)
                   .set_rate(CellEventType::DUPLICATION, epistate.duplication_rate);

            if (species_json.contains("epigenetic switches")) {
                if (!species_json["epigenetic switches"].is_array()) {
                    throw Error<std::runtime_error>(("Any \"epigenetic switches\" field must be an "
                                                     "array of epigenetic switches. This is not the "
                                                     "case for mutant \"") + mutant_name + "\".");
                }

                for (const auto& switch_json : species_json["epigenetic switches"]) {
                    set_switch_rate(species, switch_json);
                }
            }

        }
    }

    static void add_mutant(const nlohmann::json& mutant_json)
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;

        if (!mutant_json.is_object()) {
            throw Error<std::runtime_error>("Any mutant specification must be an object.");
        }

        if (!mutant_json.contains("name")) {
            throw Error<std::runtime_error>("Any mutant specification must contain a \"name\" field.");
        }

        const std::string mutant_name = mutant_json["name"].template get<std::string>();

        simulation.add_mutant(mutant_name);

        if (mutant_json.contains("epigenetic states")) {
            set_species_rates(mutant_name, mutant_json["epigenetic states"]);
        } else {
            auto& species = simulation.tissue().get_species(mutant_name);

            if (mutant_json.contains("death rate")) {
                const double rate = mutant_json["death rate"].template get<double>();

                species.set_rate(CellEventType::DEATH, rate);
            }

            if (mutant_json.contains("duplication rate")) {
                const double rate = mutant_json["duplication rate"].template get<double>();

                species.set_rate(CellEventType::DUPLICATION, rate);
            }
        }
    }

    static void configure_tissue(const nlohmann::json& tissue_json)
    {
        using namespace CLONES;
        using namespace CLONES::Mutants::Evolutions;

        if (!tissue_json.contains("size")) {
            throw Error<std::runtime_error>("Any tissue specification must contain "
                                           "a \"size\" field.");
        }

        if (!tissue_json["size"].is_array()) {
            throw Error<std::runtime_error>("The tissue specification size must be "
                                           "an array.");
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

        using namespace CLONES;
        throw Error<std::runtime_error>(("Any tissue must be either a 2D "
                                        "or a 3D space. Got a ")
                                       + std::to_string(sizes.size())
                                       + " specification.");
    }

    static CLONES::Mutants::Evolutions::PositionInTissue
    get_position(const nlohmann::json& position_json)
    {
        using namespace CLONES;
        using namespace CLONES::Mutants::Evolutions;

        if (!position_json.is_array()) {
            throw Error<std::runtime_error>("The \"position\" field must be an "
                                           "array of Natural values.");
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

        throw Error<std::runtime_error>(("Any tissue may be either a 2D or "
                                        "a 3D space. Got ")
                                       + std::to_string(position.size())
                                       + " dimension specification.");
    }

    static void configure_initial_cells(const nlohmann::json& initial_cells_json)
    {
        using namespace CLONES;

        if (!initial_cells_json.is_array()) {
            throw Error<std::runtime_error>("The \"initial cells\" field must contain an array.");
        }
        for (const auto& initial_cell_json : initial_cells_json) {

            if (!initial_cell_json.is_object()) {
                throw Error<std::runtime_error>("Every element in the \"initial cells\" "
                                               "field must be an object.");
            }

            if (!initial_cell_json.contains("species")) {
                throw Error<std::runtime_error>("Every initial cell must contain a \"species\" field.");
            }

            const std::string species_name = initial_cell_json["species"].template get<std::string>();

            auto& species = simulation.tissue().get_species(species_name);

            if (!initial_cell_json.contains("position")) {
                throw Error<std::runtime_error>("Every initial cell must contain "
                                               "a \"position\" field.");
            }
            simulation.place_cell(species, get_position(initial_cell_json["position"]));
        }
    }

    static void configure_timed_events(const nlohmann::json& timed_events_json)
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;

        if (!timed_events_json.is_array()) {
            throw Error<std::runtime_error>("The \"timed events\" field must contain an array.");
        }

        for (const auto& timed_event_json : timed_events_json) {
            if (!timed_event_json.is_object()) {
                throw Error<std::runtime_error>("Every element in the \"timed events\" field "
                                               "must be an object.");
            }

            Evolutions::TimedEvent timed_event = ConfigReader::get_timed_event(simulation,
                                                                               timed_event_json);

            simulation.schedule_timed_event(timed_event);
        }
    }

    static std::set<std::string> collect_epistate_names(const nlohmann::json& simulation_cfg)
    {
        std::set<std::string> epistate_names;
        if (!simulation_cfg.contains("epigenetic states")) {
            return epistate_names;
        }

        if (!simulation_cfg["epigenetic states"].is_array()) {
            using namespace CLONES;

            throw Error<std::runtime_error>("The \"epigenetic states\" field must be an "
                                            "array of epigenetic state names.");
        }

        for (const auto& epi_json: simulation_cfg["epigenetic states"]) {
            std::string name = epi_json.template get<std::string>();

            if (epistate_names.count(name)>0) {
                std::cerr << "The epigenetic state \"" << name
                          << "\" has been added multiple times." << std::endl;
            } else {
                epistate_names.insert(name);
            }
        }

        return epistate_names;
    }

    static void configure_simulation(const std::string& simulation_filename)
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;

        std::ifstream f(simulation_filename);
        nlohmann::json simulation_cfg = nlohmann::json::parse(f);

        std::map<std::string, MutantProperties> mutants;

        if (!simulation_cfg.contains("tissue")) {
            throw Error<std::runtime_error>("The simulation configuration file "
                                           "must contain a \"tissue\" field.");
        }
        configure_tissue(simulation_cfg["tissue"]);

        if (!simulation_cfg.contains("mutants")) {
            throw Error<std::runtime_error>("The simulation configuration file "
                                           "must contain a \"mutants\" field.");
        }
        if (!simulation_cfg["mutants"].is_array()) {
            throw Error<std::runtime_error>("The \"mutants\" field must be an "
                                           "array of mutant specifications.");
        }

        for (const auto& mutant_json: simulation_cfg["mutants"]) {
            add_mutant(mutant_json);
        }

        if (!simulation_cfg.contains("initial cells")) {
            throw Error<std::runtime_error>("The simulation configuration file must "
                                           "contain a \"initial cells\" field.");
        }
        configure_initial_cells(simulation_cfg["initial cells"]);

        if (simulation_cfg.contains("timed events")) {
            configure_timed_events(simulation_cfg["timed events"]);
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
            ("frames-per-second,f", po::value<uint8_t>(&frames_per_second)->default_value(1),
             "the number of frames per second")
    #endif
            ("help,h", "get the help");

        hidden_options.add_options()
            ("simulation file", po::value<std::filesystem::path>(&simulation_filename),
             "the name of the file describing the simulation")
            ("time horizon", po::value<double>(&time_horizon),
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
