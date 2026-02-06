/**
 * @file id_context_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a ID context indices
 * @version 1.4
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

#include "id_context_index.hpp"

#include "utils.hpp"
#include "basic_IO.hpp"
#include "fasta_chr_reader.hpp"

#include "genomic_sequence.hpp"


namespace CLONES
{

std::list<Mutations::IDContext>
partition<Mutations::IDContext>::get_class_of(const Mutations::IDContext& context)
{
    using namespace Mutations;

    std::list<IDContext> class_contexts{context};
    switch (context.fragment_type()) {
        case IDContext::FragmentType::HOMOPOLYMER:
            {
                const auto complement_base = GenomicSequence::get_complement(context.unit_base());

                auto rev_context = IDContext::build_for_homopolymer(complement_base,
                                                                    context.num_of_repetitions());

                class_contexts.push_back(std::move(rev_context));
            }
            break;

        default:
            break;
    }

    return class_contexts;
}

namespace Mutations
{

RepetitionReference::RepetitionReference():
    position{}, unit_size{0}
{}

RepetitionReference::RepetitionReference(const ChromosomeId& chr_id, const ChrPosition begin,
                                         const RepetitionType unit_size):
    position{chr_id, begin}, unit_size{unit_size}
{
    if (unit_size==0) {
        throw std::domain_error("Unit size must be greater than 0.");
    }
}

size_t IDContextIndexBuilder::init_suffix_array(const char* s,
                                                std::vector<ChrPosition>& suffix_array,
                                                std::vector<ChrPosition>& classes)
{
    const size_t alphabet_size = 1<<8;

    std::vector<ChrPosition> counter(alphabet_size, 0);

    for (size_t i = 0; i < suffix_array.size(); ++i) {
        ++counter[s[i]];
    }
    for (size_t i = 1; i < alphabet_size; ++i) {
        counter[i] += counter[i-1];
    }
    for (size_t i = 1; i <= suffix_array.size(); ++i) {
        auto rev = suffix_array.size()-i;
        suffix_array[--counter[s[rev]]] = rev;
    }
    classes[suffix_array[0]] = 0;
    size_t num_of_classes = 1;
    for (size_t i = 1; i < suffix_array.size(); ++i) {
        if (s[suffix_array[i]] != s[suffix_array[i-1]]) {
            ++num_of_classes;
        }
        classes[suffix_array[i]] = num_of_classes - 1;
    }

    return num_of_classes;
}


void IDContextIndexBuilder::update_suffix_array(const size_t h,
                                                std::vector<ChrPosition>& h_suffix_array,
                                                std::vector<ChrPosition>& h_classes,
                                                size_t& num_of_classes,
                                                std::vector<ChrPosition>& tmp_a,
                                                std::vector<ChrPosition>& tmp_b)
{
    for (size_t i = 0; i < h_suffix_array.size(); ++i) {
        if (h_suffix_array[i] >= h) {
            tmp_a[i] = h_suffix_array[i] - h;
        } else {
            tmp_a[i] = h_suffix_array[i] + h_suffix_array.size() - h;
        }
    }
    auto& counter = tmp_b;
    fill(counter.begin(), counter.begin() + num_of_classes, 0);
    for (size_t i = 0; i < h_suffix_array.size(); ++i) {
        ++counter[h_classes[tmp_a[i]]];
    }
    for (size_t i = 1; i < num_of_classes; ++i) {
        counter[i] += counter[i-1];
    }
    for (size_t i = 1; i <= h_suffix_array.size(); ++i) {
        const auto& curr = tmp_a[h_suffix_array.size()-i];
        h_suffix_array[--counter[h_classes[curr]]] = curr;
    }

    auto& new_classes = tmp_b;
    new_classes[h_suffix_array[0]] = 0;
    num_of_classes = 1;
    for (size_t i = 1; i < h_suffix_array.size(); i++) {
        const auto& curr = h_suffix_array[i];
        const auto& prev = h_suffix_array[i-1];
        if ((h_classes[curr] != h_classes[prev])
                || (h_classes[(curr + h) % h_suffix_array.size()] !=
                    h_classes[(prev + h) % h_suffix_array.size()])) {
            ++num_of_classes;
        }
        new_classes[curr] = num_of_classes - 1;
    }
    std::swap(h_classes, new_classes);
}

void IDContextIndexBuilder::add_repetition(const ChromosomeId& chr_id, const char* seq,
                                           const ChrPosition& begin, const size_t unit_size,
                                           const ChrPosition& r_begin, const ChrPosition& r_end,
                                           std::vector<bool>& covered)
{
    const auto rep_begin = r_begin+begin;
    if (rep_begin>1) {
        size_t num_of_repetitions = 1+(r_end-r_begin)/unit_size;
        GenomicPosition g_pos{chr_id, rep_begin};
        add_polymer(g_pos, num_of_repetitions, seq+r_begin, unit_size);

        std::fill(covered.begin()+r_begin,
                  covered.begin()+r_end+unit_size, true);
    }
}

void IDContextIndexBuilder::add_null_heteropolymer(const ChromosomeId& chr_id, const size_t unit_size,
                                                   const ChrPosition& begin, const ChrPosition& r_begin)
{
    const auto rep_begin = r_begin+begin+1;
    const auto fl_code = get_unit_size_code(unit_size);
    const auto sl_code = get_num_of_repetitions_code(0);

    IDContext context{IDContext::FragmentType::HETEROPOLYMER,
                    fl_code, sl_code};

    RepetitionReference rep_ref{chr_id, rep_begin,
                                static_cast<RepetitionType>(unit_size)};
    this->insert(context, std::move(rep_ref));
}

void IDContextIndexBuilder::add_null_homopolymer(const size_t nucleotide_index, const char* seq,
                                                 const ChromosomeId& chr_id, const ChrPosition& begin,
                                                 const ChrPosition& r_begin)
{
    const GenomicPosition gpos{chr_id, r_begin+begin+1};

    add_polymer(gpos, static_cast<uint8_t>(0),
                seq+nucleotide_index, 1);
}

std::map<ChrPosition, std::map<size_t, ChrPosition>>
IDContextIndexBuilder::collect_candidates(const ChrPosition& begin, const size_t& h,
                                          std::vector<ChrPosition>& h_suffix_array,
                                          std::vector<ChrPosition>& classes)
{
    ChrPosition next_h = (h>std::numeric_limits<ChrPosition>::max()/2?
                            std::numeric_limits<ChrPosition>::max():2*h);

    std::map<ChrPosition, std::map<size_t, ChrPosition>> candidates;
    ChrPosition r_begin=0, r_end=0, curr_delta = next_h;
    for (size_t i = 1; i < h_suffix_array.size(); ++i) {
        const auto& curr = h_suffix_array[i];
        const auto& prev = h_suffix_array[i-1];
        const auto delta = curr-prev-h;

        if (classes[curr] == classes[prev] && curr >= h + prev
                && curr < next_h + prev && curr+delta < h_suffix_array.size()
                && classes[curr+delta] == classes[prev+delta]) {
            if (delta != curr_delta && curr_delta != next_h) {
                if (r_begin<r_end) {
                    if (begin+r_begin>1) {
                        candidates[r_begin][h+curr_delta] = r_end;
                    }

                    r_begin = curr;
                }
            }

            curr_delta = delta;
            r_end = curr;
        } else {
            if (r_begin<r_end && begin+r_begin>1) {
                candidates[r_begin][h+curr_delta] = r_end;
            }

            r_begin = curr;
            r_end = curr;
            curr_delta = next_h;
        }
    }
    if (r_begin<r_end && begin+r_begin>1) {
        candidates[r_begin][h+curr_delta] = r_end;
    }

    return candidates;
}

void IDContextIndexBuilder::add_repetitions(const ChromosomeId& chr_id, const char* seq,
                                            const ChrPosition& begin, const size_t& h,
                                            std::vector<ChrPosition>& h_suffix_array,
                                            std::vector<ChrPosition>& classes,
                                            std::vector<bool>& covered)
{
    auto candidates = collect_candidates(begin, h, h_suffix_array, classes);

    ChrPosition r_begin=0;
    std::map<size_t, ChrPosition> r_endings;
    for (auto c_it=candidates.begin(); c_it != candidates.end(); ++c_it, ++r_begin) {
        for (const auto& [unit_size, r_end]: c_it->second) {
            const auto& r_begin = c_it->first;
            auto e_it = r_endings.find(unit_size);

            if (e_it != r_endings.end()) {
                if (e_it->second < r_end) {
                    e_it->second = r_end;

                    add_repetition(chr_id, seq, begin, unit_size,
                                   r_begin, r_end, covered);
                }
            } else {
                r_endings.emplace(unit_size, r_end);

                add_repetition(chr_id, seq, begin, unit_size,
                               r_begin, r_end, covered);
            }
        }
    }
}

void IDContextIndexBuilder::add_microhomologies(const ChromosomeId& chr_id, const char* seq,
                                                const ChrPosition& begin,
                                                std::vector<bool>& covered)
{
    for (size_t i=1; i<covered.size()-2; ++i) {
        if (!covered[i]) {
            char const* head = seq+i;

            size_t j=i+2;
            auto cover_it = covered.begin()+j;
            while (j<std::min(covered.size()-1, i+50) && !(*cover_it)) {
                char const* head_z = head;
                char const* tail_z = seq+j;

                while (tail_z<seq+std::min(covered.size()-1, i+50)
                        && !(*cover_it) && *(head_z)==*(tail_z)
                        && head_z < seq+j) {
                    ++head_z;
                    ++tail_z;
                    ++cover_it;
                }

                if (head < head_z && head_z < seq+j) {
                    const size_t homology_distance = j-i;
                    const size_t homology_size = head_z-head;

                    const auto fl_code = get_unit_size_code(homology_distance);
                    const auto sl_code = get_homology_size_code(homology_size);

                    IDContext context{IDContext::FragmentType::MICROHOMOLOGY,
                                    fl_code, sl_code};

                    RepetitionReference rep_ref{chr_id, static_cast<ChrPosition>(begin+i),
                                                static_cast<RepetitionType>(homology_size)};

                    this->insert(std::move(context), std::move(rep_ref));
                }
                ++j;
                cover_it = covered.begin()+j;
            }
        }
    }
}

void IDContextIndexBuilder::add_non_repeated_seq(const ChromosomeId& chr_id, const char* seq,
                                                 const ChrPosition& begin,
                                                 std::vector<bool>& covered)
{
    ChrPosition begin_uncovered{0};
    std::vector<ChrPosition> last_char(1<<8, 0);
    for (size_t i=0; i<covered.size(); ++i) {
        if (covered[i]) {
            if (begin_uncovered != i) {
                for (size_t unit_size=2; unit_size<6; ++unit_size) {
                    for (size_t j=begin_uncovered; j+unit_size<i; ++j) {
                        add_repetition(chr_id, seq, begin, unit_size, j, j, covered);
                        add_null_heteropolymer(chr_id, unit_size, begin, j);
                    }
                }
            }
            begin_uncovered = i+1;
        } else {
            if (begin_uncovered == i) {
                last_char['A'] = last_char['C'] = last_char['G']
                                = last_char['T'] = i;
            }

            const char& curr_char = *(seq+i);
            if (last_char[curr_char]+4<i) {
                for (size_t j=last_char[curr_char]+2; j<i-2; ++j) {
                    add_null_homopolymer(i, seq, chr_id, begin, j);
                }
            }
            last_char[curr_char] = i;

            add_repetition(chr_id, seq, begin, 1, i, i, covered);
        }
    }
}

std::vector<bool>
IDContextIndexBuilder::add_repetitions(const ChromosomeId& chr_id, const char* seq,
                                       const ChrPosition begin, const size_t& length,
                                       UI::ProgressBar& progress_bar)
{
    std::vector<bool> covered(length, false);

    std::vector<ChrPosition> suffix_array(length), classes(length),
                            tmp_a(length), tmp_b(length);

    size_t num_of_classes = init_suffix_array(seq, suffix_array, classes);

    size_t h_max = ceil_div(max_unit_size, static_cast<RepetitionType>(2));
    if (h_max > length) {
        h_max = length;
    }

    size_t h;
    size_t next_h;
    for (h = 1; h < h_max; h=next_h) {
        next_h = (h>std::numeric_limits<size_t>::max()/2?
                    std::numeric_limits<size_t>::max():2*h);

        add_repetitions(chr_id, seq, begin, h, suffix_array, classes, covered);
        update_suffix_array(h, suffix_array, classes, num_of_classes, tmp_a, tmp_b);

        progress_bar.update_elapsed_time();
    }
    add_repetitions(chr_id, seq, begin, h, suffix_array, classes, covered);

    return covered;
}

void IDContextIndexBuilder::add_contexts_from(const ChromosomeId& chr_id,
                                              const std::string& sequence,
                                              const ChrPosition begin, size_t length,
                                              UI::ProgressBar& progress_bar)
{
    if (length < 2) {
        return;
    }

    const char *subseq = sequence.c_str()+begin-1;
    length = std::min(sequence.size()-begin+1, length);
    auto covered = add_repetitions(chr_id, subseq, begin, length, progress_bar);

    add_microhomologies(chr_id, subseq, begin, covered);

    add_non_repeated_seq(chr_id, subseq, begin, covered);
}

IDContext::FirstLevelType
IDContextIndexBuilder::get_unit_size_code(const size_t& unit_size)
{
    const auto unit_size_code{(unit_size>5?5:unit_size)};
    if (unit_size_code > std::numeric_limits<IDContext::FirstLevelType>::max()) {
        throw std::runtime_error("IDContextIndex::get_unit_size_code(): "
                                + std::to_string(unit_size_code) + " is not "
                                + "representable by IDContext::FirstLevelType.");
    }

    return static_cast<IDContext::FirstLevelType>(unit_size_code);
}

IDContext::SecondLevelType
IDContextIndexBuilder::get_num_of_repetitions_code(const size_t& num_of_repetitions)
{
    const auto num_of_rep_code{(num_of_repetitions>6?6:num_of_repetitions)};
    if (num_of_rep_code > std::numeric_limits<IDContext::SecondLevelType>::max()) {
        throw std::runtime_error("IDContextIndex::get_num_of_repetitions_code(): "
                                + std::to_string(num_of_rep_code) + " is not "
                                + "representable by IDContext::SecondLevelType.");
    }

    return static_cast<IDContext::SecondLevelType>(num_of_rep_code);
}

IDContext::SecondLevelType
IDContextIndexBuilder::get_homology_size_code(const size_t& homology_size)
{
    const auto size_code{(homology_size<5?homology_size:5)};
    if (size_code > std::numeric_limits<IDContext::SecondLevelType>::max()) {
        throw std::runtime_error("IDContextIndex::get_homology_size_code(): "
                                + std::to_string(size_code) + " is not "
                                + "representable by IDContext::SecondLevelType.");
    }

    return static_cast<IDContext::SecondLevelType>(size_code);
}

void IDContextIndexBuilder::add_polymer(const GenomicPosition& genomic_position,
                                        const size_t& num_of_repetitions,
                                        const char* unit, const size_t& unit_size)
{
    if (unit_size==0) {
        throw std::domain_error("Only initialized repetitions can be added.");
    }

    IDContext::FragmentType f_type;
    IDContext::FirstLevelType fl_code;

    if (unit_size == 1) {
        f_type = IDContext::FragmentType::HOMOPOLYMER;
        fl_code = unit[0];
    } else {
        f_type = IDContext::FragmentType::HETEROPOLYMER;
        fl_code = get_unit_size_code(unit_size);
    }

    const auto sl_code = get_num_of_repetitions_code(num_of_repetitions);

    IDContext context{f_type, fl_code, sl_code};
    RepetitionReference rep_ref{genomic_position.chr_id, genomic_position.position,
                                static_cast<RepetitionType>(unit_size)};

    this->insert(std::move(context), std::move(rep_ref));
}

void IDContextIndexBuilder::add_contexts_from(const ChromosomeId& chr_id,
                                              const std::string& chr_sequence,
                                              const std::set<GenomicRegion>& regions_to_avoid,
                                              UI::ProgressBar& progress_bar)
{
    ChrPosition begin=1, region_to_avoid_begin;
    size_t length=0;
    auto next_to_avoid = regions_to_avoid.begin();
    if (next_to_avoid != regions_to_avoid.end()) {
        region_to_avoid_begin = next_to_avoid->begin();
    } else {
        region_to_avoid_begin = chr_sequence.size()+1;
    }
    for (ChrPosition i=0; i<static_cast<ChrPosition>(chr_sequence.size()); ++i) {
        if (chr_sequence[i] != 'N' && chr_sequence[i] != 'n'
                && i < region_to_avoid_begin) {
            if (length == 0) {
                begin = i+1;
            }
            ++length;
        } else {
            if (length > 0) {
                add_contexts_from(chr_id, chr_sequence, begin, length,
                                  progress_bar);
                length = 0;
            }

            if (i>=region_to_avoid_begin
                    && next_to_avoid != regions_to_avoid.end()) {
                ++next_to_avoid;

                if (next_to_avoid != regions_to_avoid.end()) {
                    region_to_avoid_begin = next_to_avoid->begin();
                } else {
                    region_to_avoid_begin = chr_sequence.size()+1;
                }
            }
        }
    }

    add_contexts_from(chr_id, chr_sequence, begin, length, progress_bar);
}

IDContextIndexBuilder::IDContextIndexBuilder(const std::filesystem::path index_path,
                                             const RepetitionType max_unit_size,
                                             const size_t cache_size):
    Archive::IndexBuilder<IDContext, RepetitionReference>{index_path, cache_size},
    max_unit_size{max_unit_size}
{}

}   // Mutations

}   // CLONES

namespace std
{
    std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::RepetitionReference& rep_ref)
    {
        os << rep_ref.unit_size << " (" << rep_ref.position << ")";

        return os;
    }
}