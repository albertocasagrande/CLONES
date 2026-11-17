/**
 * @file tissue_sampler.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the CLONES tool that sample tissues
 * @version 1.4
 * @date 2026-11-17
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
#include <fstream>
#include <string>
#include <filesystem>
#include <list>
#include <set>

#include <boost/program_options.hpp>

#include "common.hpp"

#include "tissue_simulation.hpp"

class TissueSampler : public BasicExecutable
{
    std::filesystem::path snapshot_path;
    std::filesystem::path species_directory;
    std::filesystem::path config;

    bool quiet;

    nlohmann::json load_sampling_config() const
    {
        using namespace CLONES;

        std::ifstream config_stream(config);
        if (!config_stream) {
            throw Error<std::runtime_error>("Cannot open \"" + to_string(config)
                                            + "\" for reading.");
        }

        auto sampling_cfg = nlohmann::json::parse(config_stream);

        if (!sampling_cfg.is_array()) {
            throw Error<std::runtime_error>("The sampling configuration does not "
                                            "consist in an array.");
        }

        return sampling_cfg;
    }

    void perform_sampling(const nlohmann::json& sampling_cfg) const
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;
        using namespace CLONES::Mutants::Evolutions;

        // parse all the sample specifications before loading the
        // (possibly large) simulation, so that errors are reported early
        std::list<SampleSpecification> sample_specs;
        for (const auto& sample_specs_json : sampling_cfg) {
            sample_specs.push_back(ConfigReader::get_sample_specification(sample_specs_json));
        }

        auto simulation = load_species_simulation(snapshot_path, quiet);

        // `sample_tissue()` saves a snapshot per sample: reject any
        // already used name before sampling, so that a failure does not
        // leave a partially sampled simulation
        std::set<std::string> used_names;
        for (const auto& sample : simulation.get_tissue_samples()) {
            used_names.insert(sample.get_name());
        }
        for (const auto& sample_spec : sample_specs) {
            if (!used_names.insert(sample_spec.get_name()).second) {
                throw Error<std::runtime_error>("Sample name \"" + sample_spec.get_name()
                                                + "\" has been already used.");
            }
        }

        for (const auto& sample_spec : sample_specs) {
            simulation.sample_tissue(sample_spec);
        }

        simulation.make_snapshot<CLONES::UI::ProgressBar>(nullptr);
    }
public:

    TissueSampler(int argc, char* argv[]):
        BasicExecutable(argv[0], {{"options", "Options"}})
    {
        namespace po = boost::program_options;
        namespace fs = std::filesystem;

        visible_options.at("options").add_options()
            ("quiet,q", "avoid output messages")
            ("help,h", "get the help")
        ;

        hidden_options.add_options()
            ("species simulation", po::value<std::filesystem::path>(&species_directory),
            "the species simulation directory")
            ("sampling configuration", po::value<std::filesystem::path>(&config),
            "the sampling configuration file")
        ;

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("species simulation", 1);
        positional_options.add("sampling configuration", 1);

        po::variables_map vm;
        try {
            po::command_line_parser parser{argc, argv};
            po::store(parser.options(program_options).positional(positional_options).run(), vm);
            po::notify(vm);
        } catch (const po::error& except) {
            print_help_and_exit(except.what(), 1);
        }

        if (vm.count("help")) {
            print_help_and_exit(0);
        }

        quiet = (vm.count("quiet")>0);

        if (!vm.count("species simulation")) {
            print_help_and_exit("The species simulation directory is mandatory. "
                                "You can produce it by using `species_sim`.", 1);
        }

        if (!vm.count("sampling configuration")) {
            print_help_and_exit("The sampling configuration file is mandatory.", 1);
        }

        if (!fs::is_regular_file(config)) {
            print_help_and_exit("\"" + to_string(config) + "\" is not an existing "
                                "regular file.", 1);
        }

        try {
            snapshot_path = get_last_snapshot_path(species_directory, "species simulation");
        } catch (const std::exception& except) {
            print_help_and_exit(except.what(), 1);
        }
    }

    void run() const
    {
        try {
            const auto sampling_cfg = load_sampling_config();

            perform_sampling(sampling_cfg);
        } catch (const std::exception& ex) {
            print_help_and_exit(ex.what(), 1);
        }
    }
};

int main(int argc, char* argv[])
{
    TissueSampler sampler(argc, argv);

    sampler.run();

    return 0;
}
