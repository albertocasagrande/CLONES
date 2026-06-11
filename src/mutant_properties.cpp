/**
 * @file mutant_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements the mutant properties
 * @version 1.4
 * @date 2026-06-11
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

unsigned int SpeciesProperties::counter = 0;
unsigned int MutantProperties::counter = 0;

SpeciesProperties::SpeciesProperties():
    id(0), mutant_id(0)
{}

SpeciesProperties::SpeciesProperties(const MutantProperties& mutant,
                                     const std::string& epistate_name):
    id(counter++), epistate_name{epistate_name}, mutant_id(mutant.get_id()),
    mutant_name{mutant.get_name()}
{}

std::string SpeciesProperties::get_name(const std::string& mutant_name,
                                        const std::string& epistate_name)
{
    if (epistate_name.size()==0) {
        return mutant_name;
    }

    return mutant_name + "[" + epistate_name + "]";
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
        case CellEventType::EPIGENETIC_SWITCH:
        case CellEventType::DUP_AND_EPI_SWITCH:
            event_rates[event][dst_species.get_id()] = rate;

            return *this;
        default:
            throw Error<std::runtime_error>("Unsupported event \""
                                            + cell_event_names[event] + "\".");
    }
}

MutantProperties::MutantProperties(const std::string& name, std::list<std::string>&& epi_states):
    MutantProperties(name, epi_states)
{}


MutantProperties::MutantProperties(const std::string& name, const std::list<std::string>& epi_states):
    id(counter++), name(name)
{
    for (const auto& epi_state : epi_states) {
        species.emplace(epi_state, SpeciesProperties{*this, epi_state});
    }
}

MutantProperties::MutantProperties(const std::string& name):
    MutantProperties(name, {""})
{}

const SpeciesProperties& MutantProperties::operator[](const std::string& epistate_name) const
{
    auto src_species_it = species.find(epistate_name);
    if (src_species_it == species.end()) {
        throw Error<std::out_of_range>("Unknown epistate \"" + epistate_name + "\"");
    }

    return src_species_it->second;
}

SpeciesProperties& MutantProperties::operator[](const std::string& epistate_name)
{
    auto src_species_it = species.find(epistate_name);
    if (src_species_it == species.end()) {
        throw Error<std::out_of_range>("Unknown epistate \"" + epistate_name + "\"");
    }

    return src_species_it->second;
}

bool MutantProperties::have_the_same_epistates(const MutantProperties& mutant_a,
                                               const MutantProperties& mutant_b)
{
    if (mutant_a.species.size() != mutant_b.species.size()) {
        return false;
    }

    auto species_a_it = mutant_a.species.begin();
    for (const auto& [epistate_b, species_b]: mutant_b.species) {
        if (species_a_it->first != epistate_b) {
            return false;
        }

        ++species_a_it;
    }

    return true;
}


bool operator==(const SpeciesProperties& a, const SpeciesProperties& b)
{
    return (a.get_id()==b.get_id()
            && a.get_mutant_id()==b.get_mutant_id()
            && a.get_mutant_name()==b.get_mutant_name()
            && a.get_epistate_name()==b.get_epistate_name()
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
        << mutant.get_id() << ", species=[";

    std::string sep = "";
    if (mutant.get_species().size()>1) {
        sep = "\n";
    }

    for (const auto& [epistate_name, species]: mutant.get_species()) {
        out << sep << species;
        sep = ",\n";
    }

    out << "]}";
    return out;
}

}  // std