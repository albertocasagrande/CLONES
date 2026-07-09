/**
 * @file mutant_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements the mutant properties
 * @version 1.8
 * @date 2026-07-09
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
#include <sstream>

#include "mutant_properties.hpp"
#include "cell_event.hpp"

#include "error.hpp"

namespace CLONES
{

namespace Mutants
{

std::map<std::string, MutantId> MutantProperties::mutant_ids{};
std::map<std::string, SpeciesId> SpeciesProperties::species_ids{};


MutantProperties::MutantProperties()
{}

MutantProperties::MutantProperties(const std::string& name):
    name(name)
{
    SpeciesName::assert_name(name);

    auto ids_id = mutant_ids.find(name);
    if (ids_id != mutant_ids.end()) {
        id = ids_id->second;
    } else {
        id = mutant_ids.size();

        mutant_ids.emplace(name, id);
    }
}

void SpeciesProperties::record_name(const std::string& name)
{
    auto ids_id = species_ids.find(name);
    if (ids_id != species_ids.end()) {
        id = ids_id->second;
    } else {
        id = species_ids.size();

        species_ids.emplace(name, id);
    }
}

SpeciesProperties::SpeciesProperties()
{}

SpeciesProperties::SpeciesProperties(const MutantProperties& mutant):
    mutant_properties{mutant}, epigenetic_name{""}
{
    const std::string name = SpeciesName(mutant.get_name());

    record_name(name);
}

SpeciesProperties::SpeciesProperties(const MutantProperties& mutant,
                                     const std::string& epistate_name):
    mutant_properties{mutant}, epigenetic_name{epistate_name}
{
    SpeciesName::assert_name(epistate_name);

    const std::string name = SpeciesName(mutant.get_name(), epistate_name);

    record_name(name);
}

std::string SpeciesProperties::get_name() const
{
    if (epigenetic_name=="") {
        return SpeciesName(mutant_properties.get_name());
    }

    return SpeciesName(mutant_properties.get_name(),
                        epigenetic_name);
}

double SpeciesProperties::get_rate(const CellEventType& event) const
{
    if (event != CellEventType::DEATH && event != CellEventType::DUPLICATION) {
        throw Error<std::domain_error>("SpeciesProperties::get_rate(const CellEventType&) const can "
                                       "only be called to get death and duplication rates.");
    }

    return get_rate(event, get_id());
}

double SpeciesProperties::get_rate(const CellEventType& event,
                                   const SpeciesId& dst_species) const
{
    EventRate::const_iterator event_rates_it = event_rates.find(event);

    if (event_rates_it != event_rates.end()) {
        auto rate_it = event_rates_it->second.find(dst_species);

        if (rate_it != event_rates_it->second.end()) {
            return rate_it->second;
        }
    }

    return 0;
}

SpeciesProperties& SpeciesProperties::set_rate(const CellEventType& event,
                                               const double rate)
{
    if (event != CellEventType::DEATH
            && event != CellEventType::DUPLICATION) {
        throw Error<std::domain_error>("Death or duplication event expected; got "
                                      + cell_event_names[event] + ".");
    }

    return set_rate(event, *this, rate);
}

SpeciesProperties& SpeciesProperties::set_rate(const CellEventType& event,
                                               const SpeciesProperties& dst_species,
                                               const double rate)
{
    if (dst_species.get_mutant_id() != get_mutant_id()) {
        std::ostringstream oss;

        oss << "Exclusively admitted epigenetic switches between species "
            << "belonging to the same mutant. Got " << get_name() << " and "
            << dst_species.get_name() << ".";

        throw Error<std::runtime_error>(oss.str());
    }

    switch(event) {
        case CellEventType::DEATH:
        case CellEventType::DUPLICATION:
        case CellEventType::DUP_AND_EPI_SWITCH:
            event_rates[event][dst_species.get_id()] = rate;

            return *this;
        default:
            throw Error<std::runtime_error>("Unsupported event \""
                                            + cell_event_names[event] + "\".");
    }
}

bool operator==(const SpeciesProperties& a, const SpeciesProperties& b)
{
    return (a.get_id()==b.get_id()
            && a.get_mutant_id()==b.get_mutant_id()
            && a.get_rates()==b.get_rates());
}


}   // Mutants

}   // CLONES

namespace std
{

std::ostream& operator<<(std::ostream& out, const CLONES::Mutants::SpeciesProperties& species)
{
    out << "{name: \""<< species.get_name() << "\", id: " << species.get_id()
        << ", event_rates: {";

    std::string sep="";
    for (const auto& [event, event_type_rates]: species.get_rates()) {
        for (const auto& [dst_event, rate]: event_type_rates) {
            out << sep
                << "(event: "<< CLONES::Mutants::cell_event_names[event]
                << ", dst: "<< dst_event << ", rate: " << rate << ")";
            sep = ", ";
        }
    }

    out << "}}";

    return out;
}

std::ostream& operator<<(std::ostream& out, const CLONES::Mutants::MutantProperties& mutant)
{
    out << "{name: \"" << mutant.get_name() << "\", id: "
        << mutant.get_id() << "}";

    return out;
}

}  // std