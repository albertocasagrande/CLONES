/**
 * @file sample_context_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief The main source file for context index sampler
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
#include <string>
#include <filesystem>
#include <random>

#include <boost/program_options.hpp>

#include "context_index.hpp"


class ContextSampler
{
    const std::string program_name;
    boost::program_options::options_description visible_options;
    boost::program_options::positional_options_description positional_options;

    std::string context_index_filename;
    std::string output_filename;
    unsigned int sampling_ratio;
    size_t bytes_per_abs_position;
    bool quiet;

    template<typename GENOME_WIDE_POSITION>
    void sample_context_index() const
    {
        using namespace CLONES::Mutations;

        CLONES::UI::ProgressBar* bar{nullptr};

        {
            ContextIndex<GENOME_WIDE_POSITION> context_index;

            {
                CLONES::Archive::Binary::In archive(context_index_filename);
                if (quiet) {
                    archive & context_index;
                } else {
                    archive.load(context_index, "context index");
                }
            }

            if (!quiet) {
                CLONES::UI::ProgressBar::hide_console_cursor();

                bar = new CLONES::UI::ProgressBar(std::cout);
                bar->set_message("Sampling context index");
            }

            size_t total_removal{0}, removed{0};
            for (auto& [context, positions] : context_index.get_context_positions()) {
                total_removal += positions.size()*(sampling_ratio-1)/sampling_ratio;
            }

            std::mt19937_64 generator(0);
            for (auto& [context, positions] : context_index.get_context_positions()) {
                size_t aimed_size = positions.size()/sampling_ratio;
                while (positions.size() > aimed_size) {
                    std::uniform_int_distribution<> u_dist(0,positions.size()-1);

                    context_index.extract(context, u_dist(generator));
                    ++removed;

                    if (!quiet) {
                        size_t percentage = 100*removed/total_removal;
                        if (percentage<99 && percentage > bar->get_progress()) {
                            bar->set_progress(percentage);
                        }
                    }
                }
            }

            CLONES::Archive::Binary::Out archive(output_filename);

            if (quiet) {
                archive & context_index;
            } else {
                bar->set_progress(100, "Done");

                delete bar;

                CLONES::UI::ProgressBar::show_console_cursor();

                archive.save(context_index, "sampled index");

                std::cout << " Freeing memory..." << std::flush;
            }
        }

        if (!quiet) {
            std::cout << "done" << std::endl;
        }
    }

public:

    std::ostream& print_help(std::ostream& os) const
    {
        os << "Syntax: " << program_name;

        for (unsigned int i=0;i<positional_options.max_total_count(); ++i) {
            os << " <" << positional_options.name_for_position(i) << ">";
        }

        os << std::endl << visible_options << std::endl;

        return os;
    }

    ContextSampler(const int argc, const char* argv[]):
        program_name(argv[0]), visible_options("Options")
    {
        namespace po = boost::program_options;

        visible_options.add_options()
            ("force-overwrite,f", "force overwriting output file")
#if WITH_INDICATORS
            ("quiet,q", "disable output messages")
#endif // WITH_INDICATORS
            ("help,h", "get the help");

        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("context index filename", po::value<std::string>(&context_index_filename),
            "the genome context index filename")
            ("output filename", po::value<std::string>(&output_filename),
            "the resulting genome context index filename")
            ("sampling ratio", po::value<unsigned int>(&sampling_ratio),
            "the genome context index filename")
        ;

        po::options_description generic;
        generic.add(visible_options).add(hidden);

        positional_options.add("context index filename", 1);
        positional_options.add("output filename", 1);
        positional_options.add("sampling ratio", 1);

        po::variables_map vm;
        try {
            po::command_line_parser parser{argc, argv};
            po::store(parser.options(generic).positional(positional_options).run(), vm);
            po::notify(vm);
        } catch (boost::wrapexcept<boost::program_options::unknown_option> &except) {
            std::cerr << except.what() << std::endl << std::endl;
            print_help(std::cerr);
            exit(1);
        }

        if (vm.count("help")) {
            print_help(std::cout);
            exit(0);
        }

        if (!vm.count("context index filename")) {
            std::cerr << "The context index filename is mandatory. "
                      << "You can produce it by using `build_context_index`" << std::endl << std::endl;
            print_help(std::cerr);
            exit(1);
        }

        if (!std::filesystem::exists(context_index_filename)) {
            std::cerr << "\"" << context_index_filename << "\"  does not exists."
                      << std::endl << std::endl;
            print_help(std::cerr);
            exit(1);
        }

        if (!std::filesystem::is_regular_file(context_index_filename)) {
            std::cerr << "\"" << context_index_filename << "\"  is not a regular file."
                      << std::endl << std::endl;
            print_help(std::cerr);
            exit(1);
        }

        try {
            using namespace CLONES::Mutations;

            CLONES::Archive::Binary::In archive(context_index_filename);

            bytes_per_abs_position = ContextIndex<>::read_bytes_per_absolute_position(archive);

        } catch (std::exception& except) {
            std::cerr << except.what() << std::endl << std::endl;
            print_help(std::cerr);
            exit(1);
        }

        quiet = vm.count("quiet")>0;

        if (!vm.count("output filename")) {
            std::cerr << "The output context index filename is mandatory." << std::endl << std::endl;
            print_help(std::cerr);
            exit(1);
        }

        if (std::filesystem::exists(output_filename) && !vm.count("force-overwrite")) {
            std::cerr << "The output file \"" << output_filename << "\" already exists. "
                    << "Use \"--force-overwrite\" to overwrite it." << std::endl<< std::endl;
            print_help(std::cerr);
            exit(1);
        }

        if (!vm.count("sampling ratio")) {
            std::cerr << "The sampling ratio is mandatory." << std::endl << std::endl;
            print_help(std::cerr);
            exit(1);
        }
    }

    void run() const
    {
        switch (bytes_per_abs_position) {
            case 2:
                sample_context_index<uint16_t>();
                break;
            case 4:
                sample_context_index<uint32_t>();
                break;
            case 8:
                sample_context_index<uint64_t>();
                break;
            default:
                std::cerr << "Not supported bits per absolute position."
                          << " Supported values are 2, 4, or 8." << std::endl;
                exit(1);
        }
    }
};

int main(const int argc, const char* argv[])
{
    ContextSampler sampler(argc, argv);

    sampler.run();

    return 0;
}
