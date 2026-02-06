/**
 * @file mutational_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to represent the mutational properties
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

#include <sstream>

#include "mutational_properties.hpp"

namespace CLONES
{

namespace Mutations
{

PassengerRates::PassengerRates():
    indel{0}, snv{0}, cna{0}
{}

PassengerRates::PassengerRates(const double& indel_rate, const double& SNV_rate,
                               const double& CNA_rate):
    indel{indel_rate}, snv{SNV_rate}, cna{CNA_rate}
{}

MutationalProperties::MutationalProperties()
{}

MutationalProperties&
MutationalProperties::add_mutant(const std::string& mutant_name,
                                 const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                 const std::list<MutationSpec<SID>>& driver_SIDs,
                                 const std::list<CNA>& driver_CNAs,
                                 const bool& wg_doubling)
{
    auto application_order = DriverMutations::get_default_order(driver_SIDs, driver_CNAs, wg_doubling);

    return add_mutant(mutant_name, epistate_passenger_rates, driver_SIDs,
                      driver_CNAs, application_order);
}

MutationalProperties&
MutationalProperties::add_mutant(const std::string& mutant_name,
                                 const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                 const std::list<MutationSpec<SID>>& driver_SIDs,
                                 const std::list<CNA>& driver_CNAs,
                                 const std::list<DriverMutations::MutationType>& application_order)
{
    if (driver_mutations.count(mutant_name)>0) {
        throw std::domain_error("The mutational properties of mutant \""+mutant_name+
                                "\" has been already added.");
    }

    driver_mutations.insert({mutant_name, {mutant_name, driver_SIDs,
                                           driver_CNAs, application_order}});

    for (const auto& [epistate, passenger_rate] : epistate_passenger_rates) {
        auto species_name = mutant_name+epistate;

        passenger_rates[species_name] = passenger_rate;
    }

    return *this;
}

void MutationalProperties::set_mutant_drivers(const DriverMutations& mutant_drivers)
{
    auto found = this->driver_mutations.find(mutant_drivers.name);

    if (found == this->driver_mutations.end()) {
        throw std::runtime_error("Unknown mutant \"" + mutant_drivers.name + "\".");
    }

    found->second = mutant_drivers;
}

}   // Mutations

}   // CLONES
