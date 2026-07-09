/**
 * @file species_name.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement species name representation and parsing
 * @version 1.5
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

#include <sstream>

#include "species_name.hpp"

#include "error.hpp"

namespace CLONES
{

namespace Mutants
{

SpeciesName::SpeciesName()
{}

SpeciesName::SpeciesName(const std::string mutant_name,
                         const std::string epistate_name):
    mutant_name{mutant_name},  epistate_name{epistate_name}
{
    assert_name(mutant_name);
    assert_name(epistate_name);
}

SpeciesName::SpeciesName(const std::string& species_name)
{
    if (species_name=="") {
        throw Error<std::domain_error>("A species name cannot be empty.");
    }

    const auto open_par_idx = species_name.find_first_of('[');

    if (open_par_idx == std::string::npos) {
        mutant_name = species_name;

        const auto closed_par_idx = species_name.find_first_of(']');

        if (closed_par_idx != std::string::npos) {
            throw Error<std::domain_error>("\"" + species_name
                                           + "\" is not the name of a species. "
                                           + "Unmatched ']'.");
        }

        return;
    }

    mutant_name = species_name.substr(0, open_par_idx);

    if (species_name.find_first_of('[', open_par_idx+1) != std::string::npos) {
        throw Error<std::domain_error>("\"" + species_name
                                       + "\" is not the name of a species. "
                                       + "Unexpected second '['.");
    }

    const auto closed_par_idx = species_name.find_first_of(']');

    if (open_par_idx+1==closed_par_idx) {
        throw Error<std::domain_error>("\"" + species_name
                                       + "\" is not the name of a species. "
                                       + "The epigenetic state name cannot be empty.");
    }

    if (open_par_idx>closed_par_idx) {
        throw Error<std::domain_error>("\"" + species_name
                                       + "\" is not the name of a species. "
                                       + "Unexpected ']' before '['.");
    }

    if (closed_par_idx == std::string::npos) {
        throw Error<std::domain_error>("\"" + species_name
                                       + "\" is not the name of a species. "
                                       + "Unmatched '['.");
    }

    for (size_t i=closed_par_idx+1; i<species_name.size(); ++i) {
        if (species_name[i] != ' ') {
            throw Error<std::domain_error>("\"" + species_name
                                           + "\" is not the name of a species. "
                                           + "Unexpected suffix \""
                                           + species_name.substr(i) + "\".");
        }
    }

    const auto epistate_len = closed_par_idx-open_par_idx-1;

    epistate_name = species_name.substr(open_par_idx+1, epistate_len);
}

SpeciesName::operator std::string() const
{
    if (epistate_name=="") {
        return mutant_name;
    }

    return mutant_name + "[" + epistate_name + "]";
}

void SpeciesName::assert_name(const std::string& s)
{
    if (s == "") {
         throw Error<std::domain_error>("Mutant and epigenetic state names must "
                                        "not be empty string.");
    }

    if (s.find_first_of('[') != std::string::npos) {
         throw Error<std::domain_error>("Mutant and epigenetic state names must "
                                        "cannot contain the symbol '['.");
    }

    if (s.find_first_of(']') != std::string::npos) {
         throw Error<std::domain_error>("Mutant and epigenetic state names must "
                                        "cannot contain the symbol ']'.");
    }

    if (s == "Wild-type") {
         throw Error<std::domain_error>("\"Wild-type\" is a reserved name.");
    }
}

bool SpeciesName::is_valid_name(const std::string& s)
{
    if (s == "") {
        return false;
    }

    if (s.find_first_of('[') != std::string::npos) {
        return false;
    }

    if (s.find_first_of(']') != std::string::npos) {
        return false;
    }

    if (s == "Wild-type") {
        return false;
    }

    return true;
}

}   // namespace Mutants

}   // namespace CLONES
