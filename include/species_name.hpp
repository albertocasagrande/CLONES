/**
 * @file species_name.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define species name representation and parsing
 * @version 1.3
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

#ifndef __CLONES_SPECIES_NAME__
#define __CLONES_SPECIES_NAME__

#include <string>

#include "archive.hpp"

namespace CLONES
{

namespace Mutants
{

/**
 * @brief A class to represent species names
 */
class SpeciesName
{
    std::string mutant_name;    //!< The mutant name
    std::string epistate_name;  //!< The epigenetic state name

public:
    /**
     * @brief An empty constructor
     */
    SpeciesName();

    /**
     * @brief A constructor
     *
     * @param mutant_name is the mutant name
     * @param epistate_name is the name of the epigenetic state
     */
    SpeciesName(const std::string mutant_name,
                const std::string epistate_name);

    /**
     * @brief Test whether a string is a valid mutant and epigenetic state name
     *
     * @param s is the string to be tested
     * @return `true` if and only if `s` is not an empty string, does not
     *      contain the symbols `[` and `]`, and it differs from `Wild-type`
     */
    static bool is_valid_name(const std::string& s);

    /**
     * @brief Assert that a string is possible mutant or epigenetic state name
     *
     * If the string is not a valid mutant or epigenetic state name, the method
     * throw a `Error<std::domain_error>` object.
     *
     * @param s is the string to be validate
     */
    static void assert_name(const std::string& s);

    /**
     * @brief A constructor
     *
     * This constructor parses a string to read the species name.
     *
     * @param species_name is the species name
     */
    explicit SpeciesName(const std::string& species_name);

    /**
     * @brief Cast a species name into a string
     *
     * @return The string corresponding to the species name
     */
    operator std::string() const;

    /**
     * @brief Get the mutant name
     *
     * @return A constant reference to the mutant name
     */
    inline const std::string& get_mutant_name() const
    {
        return mutant_name;
    }

    /**
     * @brief Get the epistate name
     *
     * @return A constant reference to the epistate name
     */
    inline const std::string& get_epistate_name() const
    {
        return epistate_name;
    }

    /**
     * @brief Save the species name in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & mutant_name
                & epistate_name;
    }

    /**
     * @brief Load a species from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded species
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SpeciesName load(ARCHIVE& archive)
    {
        SpeciesName species_name;

        archive & species_name.mutant_name
                & species_name.epistate_name;

        return species_name;
    }
};

}   // namespace Mutants

}   // namespace CLONES

#endif // __CLONES_SPECIES_NAME__