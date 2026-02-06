/**
 * @file build_repetition_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Builds the repetition index
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

class RepetitionIndexBuilder: public BasicExecutable
{
    std::string output_filename;
    std::string genome_fasta_filename;
    std::string driver_mutations_filename;
    size_t max_rep_per_type;
    size_t max_unit_size;
    unsigned int rep_type_size;
    bool quiet;

    template<typename GENOME_WIDE_POSITION>
    void build_and_save_rs_index() const
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
            RSIndex rs_index;

            if (quiet) {
                rs_index = RSIndex::build_index(genome_fasta_filename, regions_to_avoid,
                                                max_unit_size, max_rep_per_type, 0);
            } else {
                UI::ProgressBar::hide_console_cursor();

                UI::ProgressBar progress_bar(std::cout);

                rs_index = RSIndex::build_index(genome_fasta_filename, regions_to_avoid, max_unit_size,
                                                max_rep_per_type, 0, &progress_bar);
            }

            Archive::Binary::Out archive(output_filename);

            if (quiet) {
                archive & rs_index;
            } else {
                archive.save(rs_index, "repetition index");
                std::cout << " Cleaning memory..." << std::flush;
            }
        }
        std::cout << "done" << std::endl;
    }

public:

    RepetitionIndexBuilder(const int argc, const char* argv[]):
        BasicExecutable(argv[0], {{"generic", "Options"}})
    {
        namespace po = boost::program_options;

        visible_options.at("generic").add_options()
            ("driver-mutations,d", po::value<std::string>(&driver_mutations_filename),
             "the driver mutations file")
            ("output-filename,o", po::value<std::string>(&output_filename),
             "output filename")
            ("reps-per-type,r", po::value<size_t>(&max_rep_per_type)->default_value(500000),
             "maximum number of stored repetitions per type")
            ("max-motif-size,s", po::value<size_t>(&max_unit_size)->default_value(50),
             "maximum size of the repeated motif")
            ("force-overwrite,f", "force overwriting output file")
            ("rep-type-bytes,b", po::value<unsigned int>(&rep_type_size)->default_value(4),
             "size of the repetition type (one among 1, 2, 4, or 8)")
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
            output_filename = std::filesystem::path(genome_fasta_filename + ".rsif").filename();
        }

        std::set<size_t> admitted_rep_type_size{1,2,4,8};

        if (admitted_rep_type_size.count(rep_type_size)==0) {
            print_help_and_exit("The size of repetition type must be among 1, 2, 4, and 8.", 1);
        }

        if (max_unit_size<1 || max_unit_size > 255) {
            print_help_and_exit("The maximum motif size must lay in the interval [1, 255].", 1);
        }

        if (max_rep_per_type<=0) {
            print_help_and_exit("The maximum stored repetitions per type must be a positive value.", 1);
        }

        if (std::ifstream(output_filename).good() && !vm.count("force-overwrite")) {
            const auto msg = "The output file \"" + output_filename + "\" already exists. "
                             + "Use \"--force-overwrite\" to overwrite it.";

            print_help_and_exit(msg, 1);
        }
    }

    void run() const
    {
        switch (rep_type_size) {
            case 1:
                build_and_save_rs_index<uint8_t>();
                return;
            case 2:
                build_and_save_rs_index<uint16_t>();
                return;
            case 4:
                build_and_save_rs_index<uint32_t>();
                return;
            case 8:
                build_and_save_rs_index<uint64_t>();
                return;
            default:
                std::cerr << "Unsupported bits per absolute position."
                        << " Supported values are 1, 2, 4, or 8." << std::endl<< std::endl;

                print_help_and_exit(1);
        }
    }
};

int main(const int argc, const char* argv[])
{
    RepetitionIndexBuilder builder(argc, argv);

    builder.run();

    return 0;
}
