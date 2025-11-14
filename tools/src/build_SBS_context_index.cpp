/**
 * @file build_SBS_context_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Builds the context index
 * @version 1.3
 * @date 2025-10-31
 *
 * @copyright Copyright (c) 2023-2025
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
#include <string>
#include <filesystem>
#include <algorithm>

#include <boost/program_options.hpp>

#include "common.hpp"

#include "genome_mutations.hpp"
#include "sbs_context_index.hpp"
#include "driver_storage.hpp"

#include "progress_bar.hpp"

template<typename INDEX_TYPE>
class ContextIndexBuilder : public BasicExecutable
{
    std::string index_directory;
    std::string genome_fasta_filename;
    std::string driver_mutations_filename;
    size_t cache_size;
    bool quiet;

    void build_and_save_context_index() const
    {
        using namespace RACES;
        using namespace RACES::Mutations;

        std::set<GenomicRegion> regions_to_avoid;

        if (driver_mutations_filename!="") {
            auto driver_storage = DriverStorage::load(driver_mutations_filename);

            for (const auto& [name, mutation_entry] : driver_storage.get_code2mutation_map()) {
                regions_to_avoid.emplace(mutation_entry.mutation,
                                         std::max(static_cast<size_t>(1),
                                                  mutation_entry.mutation.ref.size()));
            }
        }

        {
            INDEX_TYPE context_index;

            std::mt19937_64 random_generator(0);
            if (quiet) {
                INDEX_TYPE::build(random_generator, index_directory, genome_fasta_filename,
                                  regions_to_avoid, cache_size);
            } else {
                UI::ProgressBar progress_bar(std::cout);

                INDEX_TYPE::build(random_generator, index_directory, genome_fasta_filename,
                                  regions_to_avoid, cache_size, progress_bar);
            }

            if (!quiet) {
                std::cout << " Cleaning memory..." << std::flush;
            }
        }

        if (!quiet) {
            UI::ProgressBar::show_console_cursor();

            std::cout << "done" << std::endl;
        }
    }

public:

    ContextIndexBuilder(const int argc, const char* argv[]):
        BasicExecutable(argv[0], {{"generic", "Options"}})
    {
        namespace po = boost::program_options;

        visible_options.at("generic").add_options()
            ("driver-mutations,d", po::value<std::string>(&driver_mutations_filename),
             "the driver mutations file")
            ("index-directory,o", po::value<std::string>(&index_directory),
             "index directory")
            ("cache-size,c", po::value<size_t>(&cache_size)->default_value(1000),
             "cache size in Mbytes")
#if WITH_INDICATORS
            ("quiet,q", "disable output messages")
#endif // WITH_INDICATORS
            ("help,h", "get the help");

        hidden_options.add_options()
            ("genome file", po::value<std::string>(&genome_fasta_filename),
            "the path to the genome in FASTA file format")
        ;

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("genome file", 1);

        po::variables_map vm;
        try {
            po::command_line_parser parser{argc, argv};
            po::store(parser.options(program_options).positional(positional_options).run(), vm);
            po::notify(vm);
        } catch (boost::wrapexcept<boost::program_options::unknown_option> &except) {
            std::cerr << except.what() << std::endl<< std::endl;
            print_help_and_exit(1);
        }

        if (vm.count("help")) {
            print_help_and_exit(0);
        }

        quiet = vm.count("quit")>0;

        if (!vm.count("genome file")) {
            print_help_and_exit("Missing genome FASTA filename.", 1);
        }

        if (!vm.count("index-directory")) {
            index_directory = "context_index";
        }

        if (std::filesystem::exists(index_directory)) {
            const auto msg = "The input directory \"" + index_directory + "\" already exists.";

            print_help_and_exit(msg, 1);
        }
    }

    void run() const
    {
        using namespace RACES::IO::FASTA;

        build_and_save_context_index();
    }
};

int main(const int argc, const char* argv[])
{
    using namespace RACES::Mutations;
    using Index = SBSContextIndex<std::mt19937_64>;

    ContextIndexBuilder<Index> builder(argc, argv);

    builder.run();

    return 0;
}
