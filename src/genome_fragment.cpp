/**
 * @file genome_fragment.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements genome DNA fragments
 * @version 1.3
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

#include <map>

#include "genome_fragment.hpp"
#include "union_map_proxy.hpp"

#include "error.hpp"

namespace CLONES
{

namespace Mutations
{

GenomeFragment::GenomeFragment()
{}

size_t GenomeFragment::Hamming_distance() const
{
    size_t total_mismatches{0};
    for (const auto& mismatch: alignment) {
        if (mismatch != MatchingType::MATCH) {
            ++total_mismatches;
        }
    }

    return total_mismatches;
}

GenomicRegion GenomeFragment::get_covered_reference_region() const
{
    if (alignment.size()==0) {
        return {this->genomic_position, 0};
    }

    size_t seq_size{0};
    for (const auto& matching : alignment) {
        // if part of the sequence have been deleted
        if (matching != MatchingType::INSERTION) { // i.e., MATCH || MISMATCH || DELETION
            ++seq_size;
        }
    }

    return {this->genomic_position,
            static_cast<GenomicRegion::Length>(seq_size)};
}

std::list<GenomicRegion> GenomeFragment::get_contained_reference_regions() const
{
    std::list<GenomicRegion> covered_list;

    if (alignment.size()==0) {
        return covered_list;
    }

    MatchingType last_seq_type = MatchingType::DELETION;
    GenomicPosition seq_begin(this->genomic_position);
    size_t seq_size{0};
    for (const auto& matching : alignment) {
        // if part of the sequence have been deleted
        if (matching == MatchingType::DELETION) {
            if (last_seq_type == MatchingType::MATCH) {
                covered_list.emplace_back(seq_begin, seq_size);
                seq_begin.position += seq_size;
            }
            ++seq_begin.position;
            seq_size = 0;
            last_seq_type = MatchingType::DELETION;
        } else if (matching != MatchingType::INSERTION) { // i.e., MATCH || MISMATCH
            ++seq_size;
            last_seq_type = MatchingType::MATCH;
        }
    }

    if (last_seq_type == MatchingType::MATCH) {
        covered_list.emplace_back(seq_begin, seq_size);
    }

    return covered_list;
}

std::string GenomeFragment::get_CIGAR() const
{
    std::ostringstream oss;

    std::map<MatchingType, std::string> matching_str{
        {MatchingType::MATCH, "M"},
        {MatchingType::MISMATCH, "X"},
        {MatchingType::INSERTION, "I"},
        {MatchingType::DELETION, "D"},
    };

    MatchingType last_seq_type = MatchingType::MATCH;
    size_t last_seq=0;
    for (const auto& matching : alignment) {
        if (last_seq_type != matching) {
            if (last_seq>0) {
                oss << last_seq
                    << matching_str.at(last_seq_type);
                last_seq = 0;
            }
            last_seq_type = matching;
        }
        ++last_seq;
    }

    if (last_seq>0) {
        oss << last_seq << matching_str.at(last_seq_type);
    }

    return oss.str();
}

inline
void validate_mutation(const std::string& reference, const SID& mutation)
{
    const auto length = std::min(mutation.ref.size(),
                                 reference.size()-mutation.position-1);

    std::string_view ref_substr(reference.c_str()+mutation.position-1, length);

    if (ref_substr != mutation.ref) {
        std::ostringstream oss;

        oss << "Wrong mutation ref sequence: expecting \""
            << mutation.ref << "\" got \"" << ref_substr << "\".";
        throw Error<std::runtime_error>(oss.str());
    }
}

void GenomeFragment::copy_reference(const std::string& reference,
                                     const size_t up_to_index,
                                     size_t& mutated_pos, size_t& reference_pos)
{
    for (; reference_pos <= up_to_index && mutated_pos < nucleotides.size();
            ++reference_pos,++mutated_pos) {
        nucleotides[mutated_pos] = reference[reference_pos-1];
        alignment_index.push_back(alignment.size());
        if (nucleotides[mutated_pos]=='N') {
            alignment.push_back(MatchingType::MISMATCH);
        } else {
            alignment.push_back(MatchingType::MATCH);
        }
    }
}

size_t get_common_prefix_size(const std::string& a, const std::string& b)
{
    const size_t max_size = std::min(a.size(), b.size());

    auto a_it = a.begin();
    auto b_it = b.begin();

    for (size_t i=0; i < max_size; ++i) {
        if (*a_it != *b_it) {
            return i;
        }
    }

    return max_size;
}

inline
void update_alignment(std::vector<MatchingType>& alignment,
                      std::vector<size_t>& alignment_index,
                      const SID& mutation, const size_t& to_be_copied,
                      size_t& last_base, const bool is_the_last)
{
    size_t common_size = get_common_prefix_size(mutation.ref,
                                                mutation.alt);

    common_size = std::min(common_size, to_be_copied);

    auto num_of_matches = std::min(mutation.ref.size(),
                                   to_be_copied);
    size_t i=0;
    for (; i<common_size; ++i) {
        alignment_index.push_back(alignment.size());
        alignment.push_back(MatchingType::MATCH);
    }
    for (; i<num_of_matches; ++i) {
        alignment_index.push_back(alignment.size());
        alignment.push_back(MatchingType::MISMATCH);
    }
    for (size_t i=num_of_matches; i<mutation.ref.size()
            && !is_the_last; ++i) {
        ++last_base;
        alignment.push_back(MatchingType::DELETION);
    }
    if (i==mutation.ref.size()) {
        for (size_t i=num_of_matches; i<to_be_copied; ++i) {
            alignment_index.push_back(alignment.size());
            alignment.push_back(MatchingType::INSERTION);
        }
    }
}

GenomeFragment::GenomeFragment(const std::string& reference,
                                 const std::map<GenomicPosition, std::shared_ptr<SID>>& germline,
                                 const std::map<GenomicPosition, std::shared_ptr<SID>>& somatic,
                                 const GenomicPosition& begin_pos,
                                 const size_t& size):
    GenomeFragment{reference, 0, germline, somatic, begin_pos, size}
{}

GenomeFragment::GenomeFragment(const std::string& reference_fragment,
                               const size_t& fragment_offset,
                               const std::map<GenomicPosition, std::shared_ptr<SID>>& germline,
                               const std::map<GenomicPosition, std::shared_ptr<SID>>& somatic,
                               const GenomicPosition& begin_pos,
                               const size_t& size):
    genomic_position{begin_pos}
{
    const auto mutation_union = union_map_proxy<GenomicPosition, std::shared_ptr<SID>>(germline, somatic);
    auto it = mutation_union.lower_bound(begin_pos);

    // remove initial deletions from the mutated fragment
    if (it != mutation_union.begin()) {
        auto prev = it;

        --prev;

        if (prev != mutation_union.begin()) {
            auto begin_ref = (prev->second)->get_region().end()+1;
            if (begin_ref>this->genomic_position.position) {
                this->genomic_position.position = begin_ref;
            }
        }
    }

    mutations = std::vector<SID>();
    nucleotides = std::string(size,'N');
    alignment.reserve(2*size);
    alignment_index.reserve(size);
    mutation_index = std::vector<MutationBaseIndex>(size);

    size_t frag_end=0, ref_end=this->genomic_position.position-fragment_offset;
    size_t last_base = std::min(this->genomic_position.position+size,
                                reference_fragment.size()+fragment_offset)-1;

    while (it != mutation_union.end() && frag_end < size
            && ref_end <= reference_fragment.size()) {

        // copy from the reference up to the mutation begin
        const auto ref_up_to = std::min(static_cast<size_t>(it->first.position-1),
                                        last_base)-fragment_offset;
        copy_reference(reference_fragment, ref_up_to, frag_end, ref_end);

        if (frag_end < size) {
            const SID& mutation = *(it->second);

#ifdef __DEBUG__
            validate_mutation(reference, mutation);
#endif // __DEBUG__

            // remove the reference from the mutated fragment
            ref_end += mutation.ref.size();

            // add the altered sequence to the mutated fragment
            auto to_be_copied = std::min(mutation.alt.size(),
                                         size-frag_end);
            nucleotides.replace(frag_end, to_be_copied, mutation.alt);
            for (size_t i=0; i<to_be_copied; ++i) {
                mutation_index[frag_end] = {mutations.size(), i};
                ++frag_end;
            }
            update_alignment(alignment, alignment_index, mutation,
                             to_be_copied, last_base, size<=frag_end);

            mutations.push_back(mutation);
        }

        ++it;
    }

    // copy from the reference until the last base
    last_base = size - frag_end - 1 + ref_end + fragment_offset;
    copy_reference(reference_fragment, last_base, frag_end, ref_end);

    nucleotides.resize(frag_end);
    mutation_index.resize(frag_end);
}

void GenomeFragment::remove_mutation(const size_t& pos)
{
    MutationBaseIndex& mb_index = mutation_index[pos];

    for (size_t i = 0; i < mutation_index.size(); ++i) {
        MutationBaseIndex& mb_index2 = mutation_index[i];
        if (mb_index2.is_valid()) {
            if (mb_index2.index == mutations.size()-1) {
                mb_index2.index = mb_index.index;
            }
        }
    }

    std::swap(mutations[mb_index.index], mutations[mutations.size()-1]);

    mutations.resize(mutations.size()-1);
}

void GenomeFragment::alter_base(const size_t pos, const char base)
{
    if (nucleotides[pos] != base) {
        nucleotides[pos] = base;

        auto& match_type = alignment[alignment_index[pos]];
        if (match_type == MatchingType::MATCH) {
            match_type = MatchingType::MISMATCH;
        }

        MutationBaseIndex mb_index = mutation_index[pos];
        if (mb_index.is_valid()) {
            auto& mutation = mutations[mb_index.index];

            mutation.alt[mb_index.base_pos] = base;
            if (!is_suffix_of(" + errors", mutation.cause)) {
                mutation.cause = mutation.cause + " + errors";
            }

            if (mb_index.base_pos == 0 && base == mutation.ref[0]) {
                alignment[alignment_index[pos]] = MatchingType::MATCH;

                if (mutation.alt.size()==1) {
                    remove_mutation(pos);
                }
            }
        }
    }
}

}   // Mutations

}   // CLONES
