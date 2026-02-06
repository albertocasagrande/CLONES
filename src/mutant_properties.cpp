/**
 * @file mutant_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements the mutant properties
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
#include <sstream>

#include "mutant_properties.hpp"
#include "cell_event.hpp"

namespace CLONES
{

namespace Mutants
{

unsigned int SpeciesProperties::counter = 0;
unsigned int MutantProperties::counter = 0;

EpigeneticRates::EpigeneticRates(const double methylation_rate, const double demethylation_rate):
        methylation(methylation_rate), demethylation(demethylation_rate)
{
    if (methylation_rate<0 || methylation_rate>1) {
        throw std::domain_error("the methylation rate does not belong to the interval [0,1]");
    }

    if (demethylation_rate<0 || demethylation_rate>1) {
        throw std::domain_error("the demethylation rate does not belong to the interval [0,1]");
    }
}

EpigeneticRates::EpigeneticRates(const double rate):
        EpigeneticRates(rate, rate)
{}

EpigeneticRates& EpigeneticRates::set_methylation_rate(const double& rate)
{
    if (rate<0 || rate>1) {
        throw std::domain_error("the rate does not belong to the interval [0,1]");
    }

    methylation = rate;

    return *this;
}

EpigeneticRates& EpigeneticRates::set_demethylation_rate(const double& rate)
{
    if (rate<0 || rate>1) {
        throw std::domain_error("the rate does not belong to the interval [0,1]");
    }

    demethylation = rate;

    return *this;
}

SpeciesProperties::SpeciesProperties():
    id(0), mutant_id(0)
{}

SpeciesProperties::SpeciesProperties(const MutantProperties& mutant,
                                     const size_t num_of_promoters):
    id(counter++), mutant_id(mutant.get_id())
{
    const size_t epigenetic_index = mutant.get_species().size();

    name = mutant.get_name();

    methylation_signature = MutantProperties::index_to_signature(epigenetic_index, num_of_promoters);
}

double SpeciesProperties::get_rate(const CellEventType& event) const
{
    if (event == CellEventType::EPIGENETIC_SWITCH) {
        throw std::domain_error("Use get_epigenetic_rate_to to get epigenetic rate.");
    }

    if (event_rates.count(event)==0) {
        return 0;
    }
    return event_rates.at(event);
}

const double& SpeciesProperties::set_rate(const CellEventType& event, const double rate)
{
    if (event == CellEventType::EPIGENETIC_SWITCH) {
        throw std::domain_error("Use set_epigenetic_rate_to to set epigenetic rate.");
    }

    return (event_rates[event] = rate);
}

double SpeciesProperties::get_epigenetic_rate_to(const SpeciesId& species_id) const
{
    const auto found = epigenetic_rates.find(species_id);

    if (found == epigenetic_rates.end()) {
        throw std::out_of_range("The species " + std::to_string(id) + " and "
                                + std::to_string(species_id) + "do not belong to the "
                                + "same mutant.");
    }

    return found->second;
}

const double& SpeciesProperties::set_epigenetic_rate_to(const SpeciesId& species_id, const double rate)
{
    auto found = epigenetic_rates.find(species_id);

    if (found == epigenetic_rates.end()) {
        throw std::out_of_range("The species " + std::to_string(id) + " and "
                                + std::to_string(species_id) + "do not belong to the "
                                + "same mutant.");
    }

    return (found->second = rate);
}

std::string SpeciesProperties::get_name() const
{
    return name + MutantProperties::signature_to_string(methylation_signature);
}

MutantProperties::MutantProperties(const std::string& name,
                                   const std::vector<EpigeneticRates>& epigenetic_event_rates):
    id(counter++), name(name)
{
    size_t epigenetic_switches = 1<<epigenetic_event_rates.size();

    for (size_t i=0; i<epigenetic_switches; ++i) {
        species.push_back(SpeciesProperties(*this, epigenetic_event_rates.size()));
    }

    for (auto& l_species: species) {
        auto e_signature = l_species.get_methylation_signature();
        auto& e_rates = l_species.epigenetic_rates;

        for (size_t i=0; i<epigenetic_event_rates.size(); ++i) {
            // get the signature of the mutant reachable by
            // methylating/demethylating the i-th promoter
            e_signature[i] = !e_signature[i];

            // get the index of the mutant reachable by
            // methylating/demethylating the i-th promoter
            size_t index = signature_to_index(e_signature);

            // get the identifier of the mutant reachable by
            // methylating/demethylating the i-th promoter
            SpeciesId dst_id = species[index].get_id();

            if (e_signature[i]) {  // if the promoter is methylated
                // a methylation event must occur
                e_rates[dst_id] = epigenetic_event_rates[i].get_methylation_rate();
            } else {
                // a demethylation event must occur
                e_rates[dst_id] = epigenetic_event_rates[i].get_demethylation_rate();
            }

            // revert to the original signature
            e_signature[i] = !e_signature[i];
        }
    }
}

MutantProperties::MutantProperties(const std::string& name):
    MutantProperties(name, std::vector<EpigeneticRates>())
{}

size_t MutantProperties::num_of_promoters() const
{
    if (species.size()==0) {
        return 0;
    }

    return species[0].get_methylation_signature().size();
}

SpeciesProperties& MutantProperties::operator[](const MethylationSignature& methylation_signature)
{
    validate_signature(methylation_signature);

    return species[MutantProperties::signature_to_index(methylation_signature)];
}

const SpeciesProperties& MutantProperties::operator[](const MethylationSignature& methylation_signature) const
{
    validate_signature(methylation_signature);

    return species[MutantProperties::signature_to_index(methylation_signature)];
}

SpeciesProperties& MutantProperties::operator[](const std::string& methylation_signature)
{
    validate_signature(methylation_signature);

    return species[MutantProperties::string_to_index(methylation_signature)];
}

const SpeciesProperties& MutantProperties::operator[](const std::string& methylation_signature) const
{
    validate_signature(methylation_signature);

    return species[MutantProperties::string_to_index(methylation_signature)];
}

size_t MutantProperties::string_to_index(const std::string& methylation_signature)
{
    size_t index=0, digit_value=(1 << methylation_signature.size());

    for (const char& c: methylation_signature) {
        digit_value >>= 1;
        if (c=='+') {
            index |= digit_value;
        }
    }

    return index;
}

std::string MutantProperties::index_to_string(const size_t& index, const size_t num_of_promoters)
{
    size_t digit_value=(1 << num_of_promoters);
    std::ostringstream oss;

    while (digit_value>1) {
        digit_value >>= 1;
        oss << (((digit_value&index) != 0) ? "+" : "-");
    }

    return oss.str();
}

size_t MutantProperties::signature_to_index(const MethylationSignature& methylation_signature)
{
    size_t index=0, digit_value=(1 << methylation_signature.size());

    for (const bool& b_value: methylation_signature) {
        digit_value >>= 1;
        if (b_value) {
            index |= digit_value;
        }
    }

    return index;
}

MethylationSignature MutantProperties::index_to_signature(const size_t& index, const size_t num_of_promoters)
{
    MethylationSignature signature;
    size_t digit_value=(1 << num_of_promoters);

    while (digit_value>1) {
        digit_value >>= 1;
        signature.push_back((digit_value&index) != 0);
    }

    return signature;
}

}   // Mutants

}   // CLONES

namespace std
{

std::ostream& operator<<(std::ostream& out, const CLONES::Mutants::EpigeneticRates& epigentic_rates)
{
    out << "{\"on\": " << epigentic_rates.get_methylation_rate()
            << ",\"off\": " << epigentic_rates.get_demethylation_rate() << "}";
    return out;
}

std::ostream& operator<<(std::ostream& out, const CLONES::Mutants::SpeciesProperties& species)
{
    out << "{name: \""<< species.get_name() << "\", id: " << species.get_id()
        << ", event_rates: {";

    std::string sep="";
    for (const auto& [event, rate]: species.get_rates()) {
        out << sep << CLONES::Mutants::cell_event_names[event] << ": " << rate;
        sep = ", ";
    }

    const auto& e_rates = species.get_epigenetic_switch_rates();
    if (e_rates.size()>0) {
        out << "}, epigenetic rates: {";
        sep="";
        for (const auto& [dst_id, rate]: e_rates) {
            out << sep << dst_id << ": " << rate;
            sep = ", ";
        }
        out << "}";
    }

    out << "}";

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

    for (const auto& species: mutant.get_species()) {
        out << sep << species;
        sep = ",\n";
    }

    out << "]}";
    return out;
}

}  // std