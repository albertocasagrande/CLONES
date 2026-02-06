/**
 * @file mutation_engine.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to place mutations on a descendants forest
 * @version 1.3
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

#include "mutation_engine.hpp"

namespace CLONES
{

namespace Mutations
{

MutationStatistics::SIDStatistics::SIDStatistics():
    mutated_alleles(0), num_of_cells(0)
{}

MutationStatistics::SampleMutationStatistics::SampleMutationStatistics():
    num_of_cells(0)
{}

MutationStatistics::MutationStatistics()
{}

MutationStatistics& MutationStatistics::record(const std::string& sample_name,
                                               const CellGenomeMutations& cell_mutations)
{
    std::set<std::shared_ptr<SID>> in_cell;

    auto& sample_stat = sample_statistics[sample_name];

    // for all chromosomes
    for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {

        for (const auto& cna : chromosome.get_CNAs()) {
            sample_stat.CNAs.push_back(cna);
            overall_statistics.CNAs.push_back(cna);
        }

        // for all chromosome alleles
        for (const auto& [allele_id, allele]: chromosome.get_alleles()) {

            // for all fragments in the allele
            for (const auto& [fragment_pos, fragment]: allele.get_fragments()) {

                // for all SIDs in the fragment
                for (const auto& [mutation_pos, mutation_ptr]: fragment.get_mutations()) {

                    ++(sample_stat.SIDs[mutation_ptr].mutated_alleles);
                    ++(overall_statistics.SIDs[mutation_ptr].mutated_alleles);

                    in_cell.insert(mutation_ptr);
                }
            }
        }
    }

    ++(sample_stat.num_of_cells);
    ++(overall_statistics.num_of_cells);
    for (const auto& mutation : in_cell) {
        ++(sample_stat.SIDs[mutation].num_of_cells);
        ++(overall_statistics.SIDs[mutation].num_of_cells);
    }

    return *this;
}

MutationStatistics& MutationStatistics::record(const std::list<CLONES::Mutations::SampleGenomeMutations>& mutations_list,
                                               UI::ProgressBar* progress_bar)
{
    using namespace CLONES;
    using namespace CLONES::Mutations;

    size_t total_steps{0};
    size_t recorded{0};

    if (progress_bar!=nullptr) {
        progress_bar->set_message("Collecting mutation data");

        for (const auto& sample_mutations : mutations_list) {
            total_steps += sample_mutations.mutations.size();
        }
    }

    for (const auto& sample_mutations : mutations_list) {
        const auto& sample_name = sample_mutations.name;
        for (const auto& mut_ptr: sample_mutations.mutations) {
            record(sample_name, *mut_ptr);

            if (progress_bar!=nullptr) {
                size_t percentage = 100*(++recorded)/total_steps;
                if (percentage>(progress_bar->get_progress())) {
                    progress_bar->set_progress(percentage);
                }
            }
        }
    }

    if (progress_bar!=nullptr) {
        progress_bar->set_message("Mutation data collected");
    }

    return *this;
}

std::ostream& MutationStatistics::write_SIDs_table(std::ostream& os, const char separator)
{
    os << "chr" << separator << "from" << separator << "to" << separator
       << "ref" << separator << "alt" << separator
       << "type" << separator << "cause" << separator << "class";

    if (sample_statistics.size()>1) {
        for (const auto& [sample_name, sample_stat]: sample_statistics) {
            os << separator << "# of mutated alleles in " << sample_name
               << separator << "# cells in " << sample_name;
        }
        os << separator << "# of mutated alleles in total"
           << separator << "# cells in total";
    } else {
        os << separator << "# of mutated alleles"
           << separator << "# of cells";
    }

    os << std::endl;

    for (const auto& [mutation_ptr, mutation_statistics] : overall_statistics.SIDs) {

        os << GenomicPosition::chrtos(mutation_ptr->chr_id) << separator
           << mutation_ptr->position
           << separator << mutation_ptr->position << separator
           << (mutation_ptr->ref==""?"-":mutation_ptr->ref) << separator
           << (mutation_ptr->alt==""?"-":mutation_ptr->alt) << separator
           << (mutation_ptr->is_SBS()?"SNV":"indel") << separator
           << mutation_ptr->cause << separator
           << mutation_ptr->get_nature_description();

        if (sample_statistics.size()>1) {
            for (const auto& [sample_name, sample_stat]: sample_statistics) {
                auto found = sample_stat.SIDs.find(mutation_ptr);
                if (found != sample_stat.SIDs.end()) {
                    const auto &mutation_stat = found->second;
                    os << separator << mutation_stat.mutated_alleles
                       << separator << mutation_stat.num_of_cells;
                } else {
                    os << separator << 0 << separator << 0;
                }
            }
        }

        os << separator << mutation_statistics.mutated_alleles
           << separator << mutation_statistics.num_of_cells
           << std::endl;
    }

    return os;
}

std::ostream& MutationStatistics::write_CNAs_table(std::ostream& os, const char separator)
{
    os << "chr" << separator << "from" << separator << "to" << separator
       << "type" << separator << "major" << separator << "minor" << std::endl;


    if (sample_statistics.size()>1) {
        for (const auto& cna: overall_statistics.CNAs) {
            os << GenomicPosition::chrtos(cna->chr_id) << separator
               << cna->begin() << separator
               << cna->end()  << separator
               << (cna->type == CNA::Type::AMPLIFICATION?"amplification":"deletion")
               << separator << "???" << separator << "???" << std::endl;
        }
    }

    return os;
}

}   // Mutations

}   // CLONES
