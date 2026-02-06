/**
 * @file build_context_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Builds the context index
 * @version 1.2
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
#include <string>
#include <filesystem>
#include <algorithm>

#include <boost/program_options.hpp>

#include "common.hpp"

#include "genome_mutations.hpp"
#include "context_index.hpp"
#include "driver_storage.hpp"

#include "progress_bar.hpp"

class ContextIndexBuilder : public BasicExecutable
{
    std::string output_filename;
    std::string genome_fasta_filename;
    std::string driver_mutations_filename;
    size_t sampling_rate;
    unsigned int byte_per_abs_position;
    bool quiet;

    template<typename GENOME_WIDE_POSITION>
    std::list<CLONES::Mutations::GenomicRegion> build_and_save_context_index() const
    {
        using namespace CLONES;
        using namespace CLONES::Mutations;

        std::list<GenomicRegion> chr_regions;
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
            using Index = ContextIndex<GENOME_WIDE_POSITION>;

            Index context_index;

            if (quiet) {
                context_index = Index::build_index(genome_fasta_filename, regions_to_avoid, sampling_rate);
            } else {
                UI::ProgressBar::hide_console_cursor();

                UI::ProgressBar progress_bar(std::cout);

                context_index = Index::build_index(genome_fasta_filename, regions_to_avoid, sampling_rate, &progress_bar);
            }

            chr_regions = context_index.get_chromosome_regions();

            if (chr_regions.size() > 0) {
                Archive::Binary::Out archive(output_filename);

                if (quiet) {
                    archive & context_index;
                } else {
                    archive.save(context_index, "context index");
                }
            }

            if (!quiet) {
                std::cout << " Cleaning memory..." << std::flush;
            }
        }

        if (!quiet) {
            UI::ProgressBar::show_console_cursor();

            std::cout << "done" << std::endl;

            if (chr_regions.size()==0) {
                std::cout << " No chromosome processed. Try to select the proper "
                          << "sequence name decoder by using \"-s\"" << std::endl;
            } else {
                std::cout << " Processed:" << std::endl;
                for (const auto& chr_region: chr_regions) {
                    std::cout << "  - Chromosome "
                                << GenomicPosition::chrtos(chr_region.get_chromosome_id())
                                << " (size: " << chr_region.size() << " bps)" << std::endl;
                }
            }
        }

        return chr_regions;
    }

public:

    ContextIndexBuilder(const int argc, const char* argv[]):
        BasicExecutable(argv[0], {{"generic", "Options"}})
    {
        namespace po = boost::program_options;

        visible_options.at("generic").add_options()
            ("driver-mutations,d", po::value<std::string>(&driver_mutations_filename),
             "the driver mutations file")
            ("output-filename,o", po::value<std::string>(&output_filename),
             "output filename")
            ("sampling-rate,s", po::value<size_t>(&sampling_rate)->default_value(1),
             "context sampling rate (to sample contexts and reduce index memory and space requirements)")
            ("force-overwrite,f", "force overwriting output file")
            ("bytes-per-pos,b", po::value<unsigned int>(&byte_per_abs_position)->default_value(4),
             "bytes per absolute position (2, 4, or 8)")
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

        if (!vm.count("output-filename")) {
            output_filename = std::filesystem::path(genome_fasta_filename + ".cif").filename();
        }

        if (vm.count("sampling-rate") && sampling_rate==0) {
            print_help_and_exit("The sampling rate must be positive.", 1);
        }

        if (std::ifstream(output_filename).good() && !vm.count("force-overwrite")) {
            const auto msg = "The output file \"" + output_filename + "\" already exists. "
                             + "Use \"--force-overwrite\" to overwrite it.";

            print_help_and_exit(msg, 1);
        }
    }

    std::list<CLONES::Mutations::GenomicRegion> run() const
    {
        using namespace CLONES::IO::FASTA;

        switch (byte_per_abs_position) {
            case 2:
                return build_and_save_context_index<uint16_t>();
            case 4:
                return build_and_save_context_index<uint32_t>();
            case 8:
                return build_and_save_context_index<uint64_t>();
            default:
                std::cerr << "Unsupported bits per absolute position."
                        << " Supported values are 1, 2, 4, or 8." << std::endl<< std::endl;

                print_help_and_exit(1);
                exit(1);
        }
    }
};

int main(const int argc, const char* argv[])
{
    ContextIndexBuilder builder(argc, argv);

    builder.run();

    return 0;
}
