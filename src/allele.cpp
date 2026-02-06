/**
 * @file allele.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements allele representation
 * @version 1.7
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

#include "allele.hpp"

namespace CLONES
{

namespace Mutations
{

AlleleFragment::Data::Data()
{}

AlleleFragment::AlleleFragment():
    GenomicRegion(), _data(std::make_shared<AlleleFragment::Data>())
{}

AlleleFragment::AlleleFragment(const ChromosomeId& chromosome_id,
                               const ChrPosition& begin, const ChrPosition& end):
    GenomicRegion({chromosome_id, begin}, end-begin+1),
    _data(std::make_shared<AlleleFragment::Data>())
{}

AlleleFragment::AlleleFragment(const GenomicRegion& genomic_region):
    GenomicRegion(genomic_region),
    _data(std::make_shared<AlleleFragment::Data>())
{}

bool AlleleFragment::has_context_free(const SID& mutation) const
{
    if (!strictly_contains(mutation)) {
        return false;
    }

    GenomicPosition g_pos{mutation};
    if (g_pos.position>0) {
        g_pos.position -= 1;
    }

    auto end = g_pos.position + mutation.ref.size() + 3;
    auto it = _data->mutations.lower_bound(g_pos);
    if (it != _data->mutations.end() && it->first.position < end) {
        return false;
    }

    if (it != _data->mutations.begin()) {
        --it;

        return (it->first.position+(it->second)->ref.size() < g_pos.position+1);
    }

    return true;
}

bool AlleleFragment::apply(const SID& mutation)
{
    if (!contains(mutation)) {
        throw std::out_of_range("The allele fragment does not "
                                "contain mutation region");
    }

    if (has_context_free(mutation)) {
        make_data_exclusive();

        (_data->mutations)[mutation] = std::make_shared<SID>(mutation);

        return true;
    }

    return false;
}

bool AlleleFragment::remove_mutation(const GenomicPosition& genomic_position)
{
    if (!contains(genomic_position)) {
        throw std::out_of_range("The allele fragment does not "
                                "contain mutation region");
    }

    auto backup_data = make_data_exclusive();

    auto it = _data->mutations.find(genomic_position);
    if (it == _data->mutations.end()) {
        _data = backup_data;

        return false;
    }

    _data->mutations.extract(it);

    return true;
}

bool AlleleFragment::includes(const SID& mutation) const
{
    if (!contains(mutation)) {
        return false;
    }

    auto it = _data->mutations.find(mutation);
    if (it == _data->mutations.end()) {
        return false;
    }

    return ((it->second)->alt == mutation.alt)
            && ((it->second)->ref == mutation.ref);
}

AlleleFragment AlleleFragment::split(const GenomicPosition& split_point)
{
    if (!contains(split_point)) {
        throw std::out_of_range("The allele fragment does not "
                                "contain mutation region");
    }

    if (begin() == split_point.position) {
        throw std::domain_error("The split point cannot be the "
                                "initial position of the fragment");
    }

    make_data_exclusive();

    AlleleFragment new_fragment(GenomicRegion::split(split_point));

    auto it = _data->mutations.lower_bound(split_point);
    while (it != _data->mutations.end()) {
        new_fragment._data->mutations.insert(_data->mutations.extract(it++));
    }

    return new_fragment;
}

std::shared_ptr<AlleleFragment::Data> AlleleFragment::make_data_exclusive()
{
    // if the data are referenced by other AlleleFragment objects
    // copy them in an exclusive object to allow the copy
    if (_data.use_count()>1) {

        auto backup = _data;
        _data = std::make_shared<AlleleFragment::Data>(*_data);

        return backup;
    }

    return _data;
}

AlleleFragment AlleleFragment::copy(const GenomicRegion& genomic_region) const
{
    AlleleFragment new_fragment(genomic_region);

    auto final_position = genomic_region.end();
    for (auto it = _data->mutations.lower_bound(genomic_region.get_initial_position());
            it != _data->mutations.end() && it->first.position <= final_position; ++it) {
        new_fragment._data->mutations.insert(*it);
    }

    return new_fragment;
}

bool AlleleFragment::has_driver_mutations_in(const GenomicRegion& genomic_region) const
{
    auto final_position = genomic_region.end();
    for (auto it = _data->mutations.lower_bound(genomic_region.get_initial_position());
            it != _data->mutations.end() && it->first.position <= final_position; ++it) {
        if ((it->second)->nature == Mutation::DRIVER) {
            return true;
        }
    }

    return false;
}

Allele::Allele()
{}

Allele::Allele(const AlleleId& identifier, const std::list<AlleleId>& history):
    history(history)
{
    this->history.push_back(identifier);
}

Allele::Allele(const AlleleId& identifier, const ChromosomeId& chromosome_id,
               const ChrPosition& begin, const ChrPosition& end,
               const std::list<AlleleId>& history):
    fragments{{{chromosome_id, begin}, AlleleFragment(chromosome_id, begin, end)}}, history(history)
{
    this->history.push_back(identifier);
}

Allele::Allele(const AlleleId& identifier, const GenomicRegion& genomic_region,
               const std::list<AlleleId>& history):
    fragments{{genomic_region.get_initial_position(), AlleleFragment(genomic_region)}}, history(history)
{
    this->history.push_back(identifier);
}

template<typename KEY, typename VALUE>
typename std::map<KEY, VALUE>::const_iterator
find_not_after(const std::map<KEY, VALUE>& value_map, const KEY& key)
{
    auto it = value_map.upper_bound(key);

    if (it != value_map.begin()) {
        --it;
    }

    return it;
}

template<typename KEY, typename VALUE>
typename std::map<KEY, VALUE>::iterator
find_not_after(std::map<KEY, VALUE>& value_map, const KEY& key)
{
    auto it = value_map.upper_bound(key);

    if (it != value_map.begin()) {
        --it;
    }

    return it;
}

bool Allele::strictly_contains(const GenomicPosition& genomic_position) const
{
    auto it = find_not_after(fragments, genomic_position);

    return (it != fragments.end() && it->second.strictly_contains(genomic_position));
}

bool Allele::contains(const GenomicPosition& genomic_position) const
{
    auto it = find_not_after(fragments, genomic_position);

    return (it != fragments.end() && it->second.contains(genomic_position));
}

bool Allele::contains(const GenomicRegion& genomic_region) const
{
    auto it = find_not_after(fragments, genomic_region.get_initial_position());

    return (it != fragments.end() && it->second.contains(genomic_region));
}

bool Allele::includes(const SID& mutation) const
{
    auto it = find_not_after(fragments, static_cast<const GenomicPosition&>(mutation));

    return (it != fragments.end() && it->second.includes(mutation));
}

bool Allele::has_context_free(const SID& mutation) const
{
    GenomicPosition g_pos{mutation};
    g_pos.position -=1;

    auto it = find_not_after(fragments, g_pos);

    return it->second.has_context_free(mutation);
}

bool Allele::apply(const SID& mutation)
{
    if (!has_context_free(mutation)) {
        return false;
    }

    auto it = find_not_after(fragments, static_cast<const GenomicPosition&>(mutation));

    if (it != fragments.end() && it->second.contains(mutation)) {
        return it->second.apply(mutation);
    }

    return false;
}

bool Allele::remove_mutation(const GenomicPosition& genomic_position)
{
    auto it = find_not_after(fragments, genomic_position);

    if (it != fragments.end() && it->second.contains(genomic_position)) {
        return it->second.remove_mutation(genomic_position);
    }

    return false;
}

Allele Allele::copy(const AlleleId& new_allele_id, const GenomicRegion& genomic_region) const
{
    Allele new_sequence(new_allele_id, history);

    auto it = find_not_after(fragments, genomic_region.get_initial_position());

    if (it == fragments.end() || !it->second.contains(genomic_region)) {
        throw std::domain_error("The allele does not fully contain the genomic region.");
    }

    auto new_fragment = it->second.copy(genomic_region);

    new_sequence.fragments[new_fragment.get_initial_position()] = std::move(new_fragment);

    return new_sequence;
}

bool Allele::has_driver_mutations_in(const GenomicRegion& genomic_region) const
{
    auto it = find_not_after(fragments, genomic_region.get_initial_position());

    while (it != fragments.end() && !it->second.begins_after(genomic_region.get_final_position())) {
        if (it->second.has_driver_mutations_in(genomic_region)) {
            return true;
        }
        ++it;
    }

    return false;
}


bool Allele::remove(const GenomicRegion& genomic_region)
{
    auto it = find_not_after(fragments, genomic_region.get_initial_position());

    if (it == fragments.end() || !it->second.contains(genomic_region)) {
        return false;
    }

    if (genomic_region.begin() > it->first.position) {
        auto new_fragment = it->second.split(genomic_region.get_initial_position());

        if (new_fragment.end()>genomic_region.end()) {
            GenomicPosition g_pos(genomic_region.get_final_position());
            g_pos.position += 1;

            new_fragment = new_fragment.split(g_pos);
        }

        fragments[new_fragment.get_initial_position()] = std::move(new_fragment);
    } else {
        if (it->second.end()>genomic_region.end()) {
            GenomicPosition g_pos(genomic_region.get_final_position());
            g_pos.position += 1;

            auto new_fragment = it->second.split(g_pos);

            fragments.extract(it);

            fragments[new_fragment.get_initial_position()] = std::move(new_fragment);
        }
    }

    return true;
}

std::map<GenomicPosition, std::shared_ptr<SID>> Allele::get_mutations() const
{
    std::map<GenomicPosition, std::shared_ptr<SID>> mutations;

    for (const auto& [pos, fragment]: fragments) {
        for (const auto& [mutation_pos, mutation]: fragment.get_mutations()) {
            mutations.insert({mutation_pos, mutation});
        }
    }

    return mutations;
}


Allele::Length Allele::size() const
{
    Length total_size{0};

    for (const auto& [pos, fragment]: fragments) {
        total_size += fragment.size();
    }

    return total_size;
}

std::string Allele::format_id(const CLONES::Mutations::AlleleId& allele_id)
{
    if (allele_id == RANDOM_ALLELE) {
        return "NA";
    }

    return std::to_string(allele_id);
}

Allele Allele::copy_structure() const
{
    Allele copy;

    copy.history = history;

    for (const auto& [pos, fragment] : fragments) {
        AlleleFragment d_frag(static_cast<const GenomicRegion&>(fragment));
        copy.fragments.emplace(pos, std::move(d_frag));
    }

    return copy;
}

bool operator==(const CLONES::Mutations::AlleleFragment& lhs,
                const CLONES::Mutations::AlleleFragment& rhs)
{
    using namespace CLONES::Mutations;

    if (static_cast<const GenomicRegion&>(lhs) != static_cast<const GenomicRegion&>(rhs)) {
        return false;
    }

    if (lhs.get_mutations().size() != rhs.get_mutations().size()) {
        return false;
    }
    
    auto lhs_it = lhs.get_mutations().begin();
    for (const auto& [pos, rhs_sid_ptr] : rhs.get_mutations()) {
        if (*rhs_sid_ptr != *(lhs_it->second)) {
            return false;
        }

        ++lhs_it;
    }

    return true;
}

bool operator!=(const CLONES::Mutations::AlleleFragment& lhs,
                const CLONES::Mutations::AlleleFragment& rhs)
{
    using namespace CLONES::Mutations;

    if (static_cast<const GenomicRegion&>(lhs) != static_cast<const GenomicRegion&>(rhs)) {
        return true;
    }

    if (lhs.get_mutations().size() != rhs.get_mutations().size()) {
        return true;
    }
    
    auto lhs_it = lhs.get_mutations().begin();
    for (const auto& [pos, rhs_sid_ptr] : rhs.get_mutations()) {
        if (*rhs_sid_ptr != *(lhs_it->second)) {
            return true;
        }

        ++lhs_it;
    }

    return false;
}

}   // Mutations

}   // CLONES

namespace std
{

std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::AlleleFragment& allele_fragment)
{
    os << "["<< allele_fragment.begin() << "-"
       << allele_fragment.end() << "]{";

    std::string sep;
    for (const auto& [pos, mutation]: allele_fragment.get_mutations()) {
        os << sep << mutation;
        sep = ",";
    }

    os << "}";

    return os;
}

std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::Allele& allele)
{
    std::string sep;
    for (const auto& [pos, fragment]: allele.get_fragments()) {
        os << sep << fragment;

        sep = ",";
    }

    return os;
}

} // std
