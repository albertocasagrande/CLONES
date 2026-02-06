/**
 * @file mutant_mutation.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines mutations
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

#ifndef __CLONES_MUTANT_MUTATION__
#define __CLONES_MUTANT_MUTATION__

#include "mutant_properties.hpp"
#include "simulation_event.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

/**
 * @brief A structure to represent genomic mutation
 */
struct Mutation : public SimulationEvent
{
    using Type = SimulationEvent::Type;

    MutantId initial_id;   //!< The initial mutant identifier
    MutantId final_id;     //!< The final mutant identifier

    /**
     * @brief A constructor
     *
     * @param initial_id is the identifier of the initial mutant
     * @param final_id is the identifier of the final mutant
     */
    Mutation(const MutantId& initial_id, const MutantId& final_id);

    /**
     * @brief A constructor
     *
     * @param initial_mutant is the initial mutant
     * @param final_mutant is the final mutant
     */
    Mutation(const MutantProperties& initial_mutant, const MutantProperties& final_mutant);

    /**
     * @brief Save a genomic mutation in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & initial_id
                & final_id;
    }

    /**
     * @brief Load a genomic mutation from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed genomic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static Mutation load(ARCHIVE& archive)
    {
        MutantId initial_id;
        MutantId final_id;

        archive & initial_id
                & final_id;

        return {initial_id, final_id};
    }

    inline Type type() const {
        return Type::MUTATION;
    }
};

}   // Evolutions

}   // Mutants

}   // CLONES

/**
 * @brief Test the equivalence between two mutations
 *
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two mutations represent
 *      the same event
 */
inline bool operator==(const CLONES::Mutants::Evolutions::Mutation& lhs,
                       const CLONES::Mutants::Evolutions::Mutation& rhs)
{
    return (lhs.initial_id == rhs.initial_id)
            && (lhs.final_id == rhs.final_id);
}

#endif // __CLONES_MUTANT_MUTATION__