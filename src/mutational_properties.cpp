/**
 * @file mutational_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to represent the mutational properties
 * @version 1.9
 * @date 2026-06-30
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

#include "error.hpp"

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
MutationalProperties::change_rates_from(const Time& time, const std::string& mutant_name,
                                        const std::map<std::string, PassengerRates>& epistate_passenger_rates)
{
    for (const auto& [epistate, p_rates] : epistate_passenger_rates) {
        auto& t_p_rates = get_passenger_rates(mutant_name, epistate);

        t_p_rates.from(time) = p_rates;
    }

    return *this;
}

MutationalProperties&
MutationalProperties::change_rates_from(const Time& time, const std::string& mutant_name,
                                        const PassengerRates& passenger_rates)
{
    const CLONES::Mutants::SpeciesName species_name{mutant_name};

    auto& t_p_rates = get_passenger_rates(species_name);

    t_p_rates.from(time) = passenger_rates;

    return *this;
}

void MutationalProperties::set_mutant_drivers(const DriverMutations& mutant_drivers)
{
    auto found = this->driver_mutations.find(mutant_drivers.name);

    if (found == this->driver_mutations.end()) {
        throw Error<std::runtime_error>("Unknown mutant \"" + mutant_drivers.name + "\".");
    }

    found->second = mutant_drivers;
}

}   // Mutations

}   // CLONES
