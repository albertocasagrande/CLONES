/**
 * @file mutation.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a general class for mutations
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

#ifndef __CLONES_MUTATION__
#define __CLONES_MUTATION__

#include "genomic_position.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief The mutation type
 */
struct MutationType
{
    /**
     * @brief Type of the mutation
     */
    enum class Type
    {
        SBS = 0,
        INDEL = 1,
        MNP = 2,
        UNKNOWN = 3
    };

    /**
     * @brief Get the type of the mutation type
     *
     * @return the type of the mutation type
     */
    static inline constexpr Type type()
    {
        return Type::UNKNOWN;
    }

    /**
     * @brief Get the mutation type name
     *
     * @return the mutation type name
     */
    static inline const std::string name()
    {
        return "UNKNOWN";
    }
};

/**
 * @brief A general class to represent mutations
 */
struct Mutation : public GenomicPosition
{
    /**
     * @brief The nature of the mutations
     *
     * This enumeration defines the set of supported
     * mutation nature.
     * A mutation can either be DRIVER, PASSENGER,
     * PRE_NEOPLASTIC, GERMLINE, or UNDEFINED.
     */
    enum Nature {
        DRIVER,
        PASSENGER,
        PRE_NEOPLASTIC,
        GERMINAL,
        UNDEFINED
    };

    Nature nature;      //!< The mutation nature
    std::string cause;  //!< The cause of the mutation

    /**
     * @brief The empty constructor
     */
    Mutation();

    /**
     * @brief Construct a new Mutation object
     *
     * @param chr_id is the id of the SNV chromosome
     * @param chr_position is the SNV position in the chromosome
     * @param nature is the nature of the mutation
     * @param cause is the cause of the mutation
     */
    Mutation(const ChromosomeId& chr_id, const ChrPosition& chr_position,
             const Nature& nature=UNDEFINED, const std::string& cause="");

    /**
     * @brief Construct a new Mutation object
     *
     * @param position is the mutation position
     * @param nature is the nature of the mutation
     * @param cause is the cause of the mutation
     */
    Mutation(const GenomicPosition& position, const Nature& nature=UNDEFINED,
             const std::string& cause="");

    /**
     * @brief Get the mutation's inital position
     *
     * @return a constant reference to the mutation's initial position
     */
    inline const ChrPosition& begin() const
    {
        return position;
    }

    /**
     * @brief Get a string describing a mutation nature
     *
     * @param nature is the mutation nature whose description is required
     * @return a string describing the mutation
     */
    static std::string get_nature_description(const Nature& nature);

    /**
     * @brief Get a string describing the mutation nature
     *
     * @return a string describing the mutation nature
     */
    inline std::string get_nature_description() const
    {
        return get_nature_description(nature);
    }

    /**
     * @brief Save a mutation in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const GenomicPosition&>(*this)
                & cause
                & nature;
    }

    /**
     * @brief Load a mutation from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static Mutation load(ARCHIVE& archive)
    {
        Mutation mutation;

        archive & static_cast<GenomicPosition&>(mutation)
                & mutation.cause
                & mutation.nature;

        return mutation;
    }
};

}   // Mutations

}   // CLONES

#endif // __CLONES_MUTATION__