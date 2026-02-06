/**
 * @file mutation_spec.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a wrapper class to specify mutation allele
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

#ifndef __CLONES_MUTATION_SPEC__
#define __CLONES_MUTATION_SPEC__

#include "allele.hpp"
#include "mutation.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief A wrapper for mutations to specify the mutation allele
 *
 * This class is meant to be used during mutation specification
 * (e.g., during mutant specification) to declare the allele in
 * which the mutation itself should be placed.
 */
template<typename MUTATION_TYPE,
         std::enable_if_t<std::is_base_of_v<Mutation, MUTATION_TYPE>, bool> = true>
struct MutationSpec : public MUTATION_TYPE
{
    AlleleId allele_id;      //!< The mutation allele

    /**
     * @brief The empty constructor
     */
    MutationSpec():
        MUTATION_TYPE(), allele_id(RANDOM_ALLELE)
    {}

    /**
     * @brief A constructor
     *
     * This constructor builds a mutation that must be placed
     * in the specified allele.
     *
     * @tparam Args are the types of the wrapped class constructor
     *      parameters
     * @param allele_id is the allele in which the mutation must
     *      be placed
     * @param args are the wrapped class constructor parameters
     */
    template<typename... Args>
    MutationSpec(const AlleleId allele_id, Args... args):
        MUTATION_TYPE(args...), allele_id(allele_id)
    {}

    /**
     * @brief Save a mutation with allele in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const MUTATION_TYPE&>(*this)
                & allele_id;
    }

    /**
     * @brief Load a mutation with allele from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded mutation with allele
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static MutationSpec<MUTATION_TYPE> load(ARCHIVE& archive)
    {
        MutationSpec<MUTATION_TYPE> mutation;

        archive & static_cast<MUTATION_TYPE&>(mutation)
                & mutation.allele_id;

        return mutation;
    }
};

}   // Mutations

}   // CLONES


namespace std
{

template<typename MUTATION_TYPE>
struct less<CLONES::Mutations::MutationSpec<MUTATION_TYPE>>
{
    inline bool operator()(const CLONES::Mutations::MutationSpec<MUTATION_TYPE> &lhs,
                           const CLONES::Mutations::MutationSpec<MUTATION_TYPE> &rhs) const
    {
        less<MUTATION_TYPE> cmp;

        if (cmp(lhs, rhs)) {
            return true;
        }

        if (cmp(rhs, lhs)) {
            return false;
        }

        return lhs.allele_id < rhs.allele_id;
    }
};

template<typename MUTATION_TYPE>
std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::MutationSpec<MUTATION_TYPE>& specification)
{
    os << static_cast<const MUTATION_TYPE&>(specification) << "(allele: ";

    if (specification.allele_id == RANDOM_ALLELE) {
        os << "random";
    } else {
        os << specification.allele_id;
    }

    os << ")";

    return os;
}

}   // std

#endif // __CLONES_MUTATION_SPEC__