/**
 * @file build_context_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Builds the context index
 * @version 1.5
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
#include <string>
#include <filesystem>
#include <algorithm>
#include <system_error>

#include <boost/program_options.hpp>

#include "common.hpp"

#include "genome_mutations.hpp"
#include "context_index.hpp"
#include "driver_storage.hpp"
#include "utils.hpp"

#include "progress_bar.hpp"

class ContextIndexBuilder : public BasicExecutable
{
    std::filesystem::path output_filename;
    std::filesystem::path genome_fasta_filename;
    std::filesystem::path driver_mutations_filename;
    size_t sampling_rate;
    unsigned int bytes_per_abs_position;
    bool quiet;

    /**
     * @brief Save an index in the output file
     *
     * The index is written in a temporary file in the output directory,
     * which then replaces the output file. So, a failure neither leaves
     * a truncated output file nor overwrites an existing one.
     */
    template<typename INDEX>
    void save(const INDEX& index) const
    {
        auto output_dir = output_filename.parent_path();
        if (output_dir.empty()) {
            output_dir = ".";
        }

        const auto tmp_path = get_a_temporary_path(output_filename.filename().string() + ".tmp_",
                                                   output_dir);

        try {
            {
                CLONES::Archive::Binary::Out archive(tmp_path);

                if (quiet) {
                    archive & index;
                } else {
                    archive.save(index, "context index");
                }
            }

            std::filesystem::rename(tmp_path, output_filename);
        } catch (...) {
            std::error_code error;
            std::filesystem::remove(tmp_path, error);

            throw;
        }
    }

    template<typename GENOME_WIDE_POSITION>
    void build_and_save_context_index() const
    {
        using namespace CLONES;
        using namespace CLONES::Mutations;

        std::list<GenomicRegion> chr_regions;
        CLONES::set<GenomicRegion> regions_to_avoid;

        if (!driver_mutations_filename.empty()) {
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

            {
                // the progress bar hides the console cursor, when shown,
                // and restores it on destruction
                UI::ProgressBar progress_bar(std::cout, quiet);

                context_index = Index::build_index(genome_fasta_filename, regions_to_avoid,
                                                   sampling_rate, &progress_bar);
            }

            chr_regions = context_index.get_chromosome_regions();

            if (chr_regions.size()==0) {
                throw Error<std::runtime_error>("No chromosome has been processed: no "
                                                "sequence in \""
                                                + to_string(genome_fasta_filename)
                                                + "\" has been recognised as a chromosome.");
            }

            save(context_index);

            if (!quiet) {
                std::cout << " Cleaning memory..." << std::flush;
            }
        }

        if (!quiet) {
            std::cout << "done" << std::endl;

            std::cout << " Processed:" << std::endl;
            for (const auto& chr_region: chr_regions) {
                std::cout << "  - Chromosome "
                            << GenomicPosition::chrtos(chr_region.get_chromosome_id())
                            << " (size: " << chr_region.size() << " bps)" << std::endl;
            }
        }
    }

public:

    ContextIndexBuilder(const int argc, const char* argv[]):
        BasicExecutable(argv[0], {{"generic", "Options"}})
    {
        namespace po = boost::program_options;
        namespace fs = std::filesystem;

        visible_options.at("generic").add_options()
            ("driver-mutations,d", po::value<std::filesystem::path>(&driver_mutations_filename),
             "the driver mutations file")
            ("output-filename,o", po::value<std::filesystem::path>(&output_filename),
             "output filename")
            ("sampling-rate,s", po::value<size_t>(&sampling_rate)->default_value(1),
             "context sampling rate (to sample contexts and reduce index memory and space requirements)")
            ("force-overwrite,f", "force overwriting output file")
            ("bytes-per-pos,b", po::value<unsigned int>(&bytes_per_abs_position)->default_value(4),
             "bytes per absolute position (2, 4, or 8)")
#if WITH_INDICATORS
            ("quiet,q", "disable output messages")
#endif // WITH_INDICATORS
            ("help,h", "get the help");

        hidden_options.add_options()
            ("genome file", po::value<std::filesystem::path>(&genome_fasta_filename),
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
        } catch (const po::error& except) {
            print_help_and_exit(except.what(), 1);
        }

        if (vm.count("help")) {
            print_help_and_exit(0);
        }

        quiet = vm.count("quiet")>0;

        if (!vm.count("genome file")) {
            print_help_and_exit("Missing genome FASTA filename.", 1);
        }

        if (!fs::is_regular_file(genome_fasta_filename)) {
            print_help_and_exit("\"" + to_string(genome_fasta_filename)
                                + "\" is not an existing regular file.", 1);
        }

        if (!driver_mutations_filename.empty() && !fs::is_regular_file(driver_mutations_filename)) {
            print_help_and_exit("\"" + to_string(driver_mutations_filename)
                                + "\" is not an existing regular file.", 1);
        }

        if (!vm.count("output-filename")) {
            output_filename = fs::path(to_string(genome_fasta_filename) + ".cif").filename();
        }

        if (sampling_rate==0) {
            print_help_and_exit("The sampling rate must be positive.", 1);
        }

        if (bytes_per_abs_position!=2 && bytes_per_abs_position!=4
                && bytes_per_abs_position!=8) {
            print_help_and_exit("Unsupported bytes per absolute position ("
                                + std::to_string(bytes_per_abs_position)
                                + "). Supported values are 2, 4, and 8.", 1);
        }

        if (fs::exists(output_filename)) {
            if (!vm.count("force-overwrite")) {
                print_help_and_exit("The output file \"" + to_string(output_filename)
                                    + "\" already exists. Use \"--force-overwrite\" "
                                    "to overwrite it.", 1);
            }

            if (!fs::is_regular_file(output_filename)) {
                print_help_and_exit("\"" + to_string(output_filename)
                                    + "\" is not a regular file.", 1);
            }
        }
    }

    void run() const
    {
        try {
            switch (bytes_per_abs_position) {
                case 2:
                    build_and_save_context_index<uint16_t>();
                    break;
                case 4:
                    build_and_save_context_index<uint32_t>();
                    break;
                case 8:
                    build_and_save_context_index<uint64_t>();
                    break;
                default:
                    // already excluded by the constructor
                    print_help_and_exit("Unsupported bytes per absolute position.", 1);
            }
        } catch (const std::exception& except) {
            print_help_and_exit(except.what(), 1);
        }
    }
};

int main(const int argc, const char* argv[])
{
    ContextIndexBuilder builder(argc, argv);

    builder.run();

    return 0;
}
