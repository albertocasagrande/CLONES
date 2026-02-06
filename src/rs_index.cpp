/**
 * @file rs_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to compute the repeated substring index
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

#include "rs_index.hpp"

#include "utils.hpp"
#include "basic_IO.hpp"
#include "fasta_chr_reader.hpp"

namespace CLONES
{

namespace Mutations
{


RSIndex::RepetitionStorage::RepetitionStorage():
    total_number(0), stored_repetitions(0)
{}

const Repetition<RSIndex::RepetitionType>& RSIndex::RepetitionStorage::extract(const size_t& pos)
{
    if (pos >= stored_repetitions) {
        throw std::out_of_range("Position not available in the storage");
    }

    if (pos+1 != stored_repetitions) {
        std::swap(repetitions[pos], repetitions[stored_repetitions-1]);
    }

    --stored_repetitions;

    return repetitions[stored_repetitions];
}

void RSIndex::RepetitionStorage::push_back(Repetition<RSIndex::RepetitionType>&& repetition)
{
    repetitions.push_back(repetition);

    ++stored_repetitions;
    ++total_number;
}

size_t RSIndex::init_suffix_array(const char* s, std::vector<ChrPosition>& suffix_array,
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

void RSIndex::update_suffix_array(const size_t h,
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

bool RSIndex::add(RepetitionMap& r_map, const ChromosomeId& chr_id,
                  const ChrPosition& begin,
                  const uint8_t index,
                  const RepetitionType& num_of_repetitions,
                  const char*& unit, const uint8_t& unit_size,
                  const char& previous_base)
{
    auto r_it = r_map.find(index);

    if (r_it == r_map.end()) {
        r_it = r_map.emplace(index, RepetitionStorage()).first;
    }

    auto& r_storage = r_it->second;

    if (r_storage.total_number < max_stored_repetitions) {
        r_storage.push_back({chr_id, begin, num_of_repetitions, unit, unit_size,
                                previous_base});
        return true;
    }

    std::uniform_int_distribution<size_t> dist(0, ++r_storage.total_number);

    auto pos = dist(random_gen);

    if (pos<max_stored_repetitions) {
        r_storage[pos] = Repetition<RepetitionType>(chr_id, begin,
                                                    num_of_repetitions,
                                                    unit, unit_size,
                                                    previous_base);
        return true;
    }

    return false;
}

void RSIndex::add_repetition(const char* s, const ChromosomeId& chr_id,
                             const ChrPosition& begin, const size_t unit_size,
                             const ChrPosition& r_begin, const ChrPosition& r_end,
                             std::vector<bool>& covered)
{
    const auto rep_begin = r_begin+begin;
    if (rep_begin>1) {
        auto num_of_repetitions = 1+(r_end-r_begin)/unit_size;
        add_polymer(chr_id, rep_begin, static_cast<RepetitionType>(num_of_repetitions),
                    s+r_begin, unit_size, *(s+r_begin-1));

        fill(covered.begin()+r_begin, covered.begin()+r_end+unit_size, true);
    }
}

std::string RSIndex::build_k_mer(const uint8_t k)
{
    std::uniform_int_distribution<uint8_t> dist(0, 3);

    std::string k_mer(k, 'T');
    for (auto& nucleotide : k_mer) {
        switch(dist(random_gen)) {
            case 0:
                nucleotide = 'A';
                break;
            case 1:
                nucleotide = 'C';
                break;
            case 2:
                nucleotide = 'G';
                break;
            case 3:
            default:
                //nucleotide = 'T';
                break;
        }
    }

    return k_mer;
}

void RSIndex::add_null_heteropolymer(const char* s, const ChromosomeId& chr_id,
                                     const uint8_t unit_size,
                                     const ChrPosition& begin,
                                     const ChrPosition& r_begin)
{
    const auto rep_begin = r_begin+begin+1;
    auto& r_storage = (*hetero_map)[unit_size][0];

    if (r_storage.total_number < max_stored_repetitions) {
        std::string k_mer = build_k_mer(unit_size);
        r_storage.push_back(Repetition<RepetitionType>(chr_id, rep_begin,
                                                        0, k_mer.c_str(),
                                                        unit_size, *(s+r_begin)));
    } else {
        std::uniform_int_distribution<size_t> dist(0, ++r_storage.total_number);

        auto pos = dist(random_gen);

        if (pos<max_stored_repetitions) {
            std::string k_mer = build_k_mer(unit_size);
            r_storage[pos] = Repetition<RepetitionType>(chr_id, rep_begin,
                                                            0, k_mer.c_str(),
                                                            unit_size, *(s+r_begin));
        }
    }
}

void RSIndex::add_null_homopolymer(const size_t nucleotide_index, const char* s,
                                   const ChromosomeId& chr_id, const ChrPosition& begin,
                                   const ChrPosition& r_begin)
{
    const auto rep_begin = r_begin+begin+1;

    add_polymer(chr_id, rep_begin, static_cast<uint8_t>(0),
                s+nucleotide_index, 1, *(s+r_begin));
}

std::map<ChrPosition, std::map<size_t, ChrPosition>>
RSIndex::collect_candidates(const ChrPosition& begin, const size_t& h,
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

void RSIndex::collect_repetitions(const char* s, const ChromosomeId& chr_id,
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

                    add_repetition(s, chr_id, begin, unit_size, r_begin, r_end, covered);
                }
            } else {
                r_endings.insert({unit_size, r_end});

                add_repetition(s, chr_id, begin, unit_size, r_begin, r_end, covered);
            }
        }
    }
}

void RSIndex::collected_non_repeated_seq(const char* s, const ChromosomeId& chr_id,
                                        const ChrPosition& begin, std::vector<bool>& covered)
{
    ChrPosition begin_uncovered{0};
    std::vector<ChrPosition> last_char(1<<8, 0);
    for (ChrPosition i=0; i<covered.size(); ++i) {
        if (covered[i]) {
            if (begin_uncovered != i) {
                for (ChrPosition unit_size=2; unit_size<6; ++unit_size) {
                    for (ChrPosition j=begin_uncovered; j+unit_size<i; ++j) {
                        add_repetition(s, chr_id, begin, unit_size, j, j, covered);
                        add_null_heteropolymer(s, chr_id, unit_size, begin, j);
                    }
                }
            }
            begin_uncovered = i+1;
        } else {
            if (begin_uncovered == i) {
                last_char['A'] = last_char['C'] = last_char['G']
                                = last_char['T'] = i;
            }

            const char& curr_char = *(s+i);
            if (last_char[curr_char]+4<i) {
                for (ChrPosition j=last_char[curr_char]+2; j<i-2; ++j) {
                    add_null_homopolymer(i, s, chr_id, begin, j);
                }
            }
            last_char[curr_char] = i;

            add_repetition(s, chr_id, begin, 1, i, i, covered);
        }
    }
}

void RSIndex::collect_microhomologies(const char* subseq, const ChromosomeId& chr_id,
                                      const ChrPosition& begin, std::vector<bool>& covered)
{
    for (size_t i=1; i<covered.size()-2; ++i) {
        if (!covered[i]) {
            char const* head = subseq+i;
            const auto mh_begin = begin+i;

            size_t j=i+2;
            auto cover_it = covered.begin()+j;
            while (j<std::min(covered.size()-1, i+50) && !(*cover_it)) {
                char const* head_z = head;
                char const* tail_z = subseq+j;

                while (tail_z<subseq+std::min(covered.size()-1, i+50)
                        && !(*cover_it)
                        && *(head_z)==*(tail_z)
                        && head_z < subseq+j) {
                    ++head_z;
                    ++tail_z;
                    ++cover_it;
                }

                if (head < head_z && head_z < subseq+j) {
                    size_t homology_distance = j-i;
                    size_t homology_size = head_z-head;

                    add((*micro_map)[get_first_index(homology_distance)], chr_id, mh_begin,
                        (homology_size<5?homology_size:5), 1, head, homology_distance,
                        *(head-1));
                }
                ++j;
                cover_it = covered.begin()+j;
            }
        }
    }
}

std::vector<bool>
RSIndex::collect_repetitions(const char *subseq, const ChromosomeId& chr_id,
                             const ChrPosition begin, const size_t& length,
                             UI::ProgressBar* progress_bar)
{
    std::vector<bool> covered(length, false);

    std::vector<ChrPosition> suffix_array(length), classes(length),
                                tmp_a(length), tmp_b(length);

    size_t num_of_classes = init_suffix_array(subseq, suffix_array, classes);

    size_t h;
    size_t next_h;
    const auto h_max = std::min(ceil_div(max_unit_size, static_cast<size_t>(2)), length);
    for (h = 1; h < h_max; h=next_h) {
        next_h = (h>std::numeric_limits<size_t>::max()/2?
                    std::numeric_limits<size_t>::max():2*h);

        collect_repetitions(subseq, chr_id, begin, h, suffix_array, classes, covered);
        update_suffix_array(h, suffix_array, classes, num_of_classes, tmp_a, tmp_b);

        if (progress_bar != nullptr) {
            progress_bar->set_progress(progress_bar->get_progress());
        }
    }
    collect_repetitions(subseq, chr_id, begin, h, suffix_array, classes, covered);

    return covered;
}

void RSIndex::collect_data_from(const std::string& sequence, const ChromosomeId& chr_id,
                                const ChrPosition begin, size_t length,
                                UI::ProgressBar* progress_bar)
{
    if (length < 2) {
        return;
    }

    const char *subseq = sequence.c_str()+begin-1;
    length = std::min(sequence.size()-begin+1, length);
    auto covered = collect_repetitions(subseq, chr_id, begin, length, progress_bar);

    collect_microhomologies(subseq, chr_id, begin, covered);

    collected_non_repeated_seq(subseq, chr_id, begin, covered);
}

char RSIndex::canonize(const char& nucleotide)
{
    if (is_AT(nucleotide)) {
        return 'T';
    }

    return 'C';
}

void RSIndex::add_polymer(const ChromosomeId& chr_id, const ChrPosition& begin,
                          const RepetitionType& num_of_repetitions,
                          const char* unit, const size_t& unit_size,
                          const char& previous_base)
{
    if (unit_size==0) {
        throw std::domain_error("Only initialized repetitions can be added.");
    }

    RepetitionMap* r_map;
    if (unit_size == 1) {
        r_map = &((*homo_map)[canonize(unit[0])]);
    } else {
        r_map = &((*hetero_map)[get_first_index(unit_size)]);
    }

    const uint8_t sl_index = static_cast<uint8_t>(num_of_repetitions>6?6:num_of_repetitions);

    add(*r_map, chr_id, begin, sl_index, num_of_repetitions, unit,
        unit_size, previous_base);
}

RSIndex::RSIndex(const size_t max_unit_size, const size_t max_stored_repetitions,
                 const int seed):
    hetero_map(std::make_shared<std::map<uint8_t, RepetitionMap>>()),
    homo_map(std::make_shared<std::map<char, RepetitionMap>>()),
    micro_map(std::make_shared<std::map<uint8_t, RepetitionMap>>()),
    random_gen(seed), max_unit_size(max_unit_size),
    max_stored_repetitions(max_stored_repetitions)
{}

void RSIndex::collect_data_from(const ChromosomeId& chr_id, const std::string& chr_sequence,
                                UI::ProgressBar* progress_bar)
{
    ChrPosition begin=1;
    size_t length=0;
    for (ChrPosition i=0; i<static_cast<ChrPosition>(chr_sequence.size()); ++i) {
        if (chr_sequence[i] != 'N' && chr_sequence[i] != 'n') {
            if (length == 0) {
                begin = i+1;
            }
            ++length;
        } else {
            if (length > 0) {
                collect_data_from(chr_sequence, chr_id, begin, length, progress_bar);
                length = 0;
            }
        }
    }
    collect_data_from(chr_sequence, chr_id, begin, length, progress_bar);
}

size_t RSIndex::count_available_for(const IDType& id_type) const
{
    std::vector<RepetitionStorage const*> s_pointers;

    switch (id_type.ftype) {
        case IDType::FragmentType::HETEROPOLYMER:
            s_pointers = find_storages(static_cast<const std::map<uint8_t, RepetitionMap>&>(*hetero_map),
                                       get_first_index(id_type.fl_index),
                                       id_type.sl_index, id_type.insertion);
            break;
        case IDType::FragmentType::HOMOPOLYMER:
            s_pointers = find_storages(static_cast<const std::map<char, RepetitionMap>&>(*homo_map),
                                       canonize(id_type.fl_index),
                                       id_type.sl_index, id_type.insertion);
            break;
        case IDType::FragmentType::MICROHOMOLOGY:
            s_pointers = find_storages(static_cast<const std::map<uint8_t, RepetitionMap>&>(*micro_map),
                                       get_first_index(id_type.fl_index),
                                       id_type.sl_index, id_type.insertion);
            break;
        default:
            throw std::domain_error("Unsupported indel type.");
    }

    size_t total_size{0};
    for (const auto& s_pointer: s_pointers) {
        total_size += s_pointer->size();
    }

    return total_size;
}

const Repetition<RSIndex::RepetitionType>&
RSIndex::select(const IDType& id_type, const bool remove)
{
    std::pair<RepetitionStorage*, size_t> rep_reference;
    switch (id_type.ftype) {
        case IDType::FragmentType::HETEROPOLYMER:
            rep_reference = find_a_heteropolymer(id_type.fl_index, id_type.sl_index,
                                                 id_type.insertion);
            break;
        case IDType::FragmentType::HOMOPOLYMER:
            rep_reference = find_a_homopolymer(id_type.fl_index, id_type.sl_index,
                                               id_type.insertion);
            break;
        case IDType::FragmentType::MICROHOMOLOGY:
            rep_reference = find_a_microhomology(id_type.fl_index, id_type.sl_index);
            break;
        default:
            throw std::domain_error("Unsupported indel type.");
    }

    if (remove) {
        return rep_reference.first->extract(rep_reference.second);
    } else {
        return (*rep_reference.first)[rep_reference.second];
    }
}

void RSIndex::restore(const IDType& id_type)
{
    RepetitionMap* r_map;
    switch (id_type.ftype) {
        case IDType::FragmentType::HETEROPOLYMER:
        case IDType::FragmentType::MICROHOMOLOGY:
            r_map = &(hetero_map->at(get_first_index(id_type.fl_index)));
            break;
        case IDType::FragmentType::HOMOPOLYMER:
            r_map = &(homo_map->at(id_type.fl_index));
            break;
        default:
            throw std::domain_error("Microhomology is not supported yet.");
    }

    r_map->at(get_second_index(id_type.sl_index, id_type.insertion)).restore();
}

RSIndex RSIndex::clone() const
{
    RSIndex copy(max_unit_size, max_stored_repetitions);

    *(copy.hetero_map) = *hetero_map;
    *(copy.homo_map) = *homo_map;
    *(copy.micro_map) = *micro_map;

    return copy;
}

RSIndex RSIndex::build_index(const std::filesystem::path& genome_fasta,
                             const std::set<GenomicRegion>& regions_to_avoid,
                             const size_t max_unit_size,
                             const size_t max_stored_repetitions,
                             const int seed, UI::ProgressBar* progress_bar)
{
    if (regions_to_avoid.size()>0) {
        throw std::domain_error("The region to avoid parameter is not supported yet.");
    }

    using namespace IO::FASTA;

    RSIndex rs_index(max_unit_size, max_stored_repetitions, seed);

    Reader<ChromosomeData<Sequence>> chr_reader(genome_fasta);

    const auto streamsize = chr_reader.get_stream_size();

    ChromosomeData<Sequence> chr;
    while (chr_reader.read(chr, *progress_bar)) {
        if (progress_bar != nullptr) {
            progress_bar->set_progress(static_cast<uint8_t>(100*chr_reader.get_position()/streamsize),
                                        "Processing chr. " + GenomicPosition::chrtos(chr.chr_id));
        }
        rs_index.collect_data_from(chr.chr_id, chr.nucleotides, progress_bar);
    }

    if (progress_bar != nullptr) {
        progress_bar->set_progress(100, "RS index built");
    }

    return rs_index;
}

void RSIndex::restore()
{
    for (auto& [unit_nucleotipe, r_map]: *homo_map) {
        for (auto& [num_of_repeats, r_storage]: r_map) {
            r_storage.restore();
        }
    }

    for (auto& [unit_size, r_map]: *hetero_map) {
        for (auto& [num_of_repeats, r_storage]: r_map) {
            r_storage.restore();
        }
    }

    for (auto& [unit_size, r_map]: *micro_map) {
        for (auto& [num_of_repeats, r_storage]: r_map) {
            r_storage.restore();
        }
    }
}

}   // Mutations

}   // CLONES
