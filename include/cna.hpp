/**
 * @file cna.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class for copy number alterations
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

#ifndef __CLONES_CNA__
#define __CLONES_CNA__

#include "allele.hpp"
#include "mutation.hpp"
#include "genomic_region.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief Copy Number Alteration type
 */
struct CNA : public Mutation
{
    /**
     * @brief The CNA type
     */
    enum class Type {
        AMPLIFICATION,
        DELETION
    };

    using Length = GenomicRegion::Length;

    Length length;      //!< The CNA length
    AlleleId source;    //!< source allele id
    AlleleId dest;      //!< destination allele id
    Type type;          //!< amplification/deletion flag

    /**
     * @brief The empty constructor
     */
    CNA();

    /**
     * @brief A constructor
     *
     * @param initial_position is the CNA initial position
     * @param length is the CNA length
     * @param type is the CNA type
     * @param nature is the nature of the CNA
     */
    CNA(const GenomicPosition& initial_position, const Length& length,
        const Type& type, const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param initial_position is the CNA initial position
     * @param length is the CNA length
     * @param source is the allele from which the region is copied/removed
     * @param type is the CNA type
     * @param nature is the nature of the CNA
     */
    CNA(const GenomicPosition& initial_position, const Length& length,
        const Type& type, const AlleleId& source,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param initial_position is the CNA initial position
     * @param length is the CNA length
     * @param source is the allele from which the region is copied/removed
     * @param destination is the allele in which the region is copied/removed
     * @param type is the CNA type
     * @param nature is the nature of the CNA
     */
    CNA(const GenomicPosition& initial_position, const Length& length,
        const Type& type, const AlleleId& source,
        const AlleleId& destination,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Build a new amplification
     *
     * @param initial_position is the CNA initial position
     * @param length is the CNA length
     * @param source is the allele from which the region is copied
     * @param destination is the allele in which the region is copied
     * @param nature is the nature of the CNA
     */
    static inline CNA new_amplification(const GenomicPosition& initial_position,
                                        const Length& length,
                                        const AlleleId& source,
                                        const AlleleId& destination,
                                        const Mutation::Nature& nature=Mutation::UNDEFINED)
    {
        return CNA(initial_position, length, Type::AMPLIFICATION,
                   source, destination, nature);
    }

    /**
     * @brief Build a new amplification
     *
     * @param initial_position is the CNA initial position
     * @param length is the CNA length
     * @param source is the allele from which the region is copied
     * @param nature is the nature of the CNA
     */
    static inline CNA new_amplification(const GenomicPosition& initial_position,
                                        const Length& length,
                                        const AlleleId& source,
                                        const Mutation::Nature& nature=Mutation::UNDEFINED)
    {
        return CNA::new_amplification(initial_position, length,
                                      source, RANDOM_ALLELE, nature);
    }

    /**
     * @brief Build a new amplification
     *
     * @param initial_position is the CNA initial position
     * @param length is the CNA length
     * @param nature is the nature of the CNA
     */
    static inline CNA new_amplification(const GenomicPosition& initial_position,
                                        const Length& length,
                                        const Mutation::Nature& nature=Mutation::UNDEFINED)
    {
        return CNA::new_amplification(initial_position, length,
                                      RANDOM_ALLELE, RANDOM_ALLELE, nature);
    }

    /**
     * @brief Build a new amplification
     *
     * @param initial_position is the CNA initial position
     * @param length is the CNA length
     * @param allele is the allele from which the region is removed
     * @param nature is the nature of the CNA
     */
    static inline CNA new_deletion(const GenomicPosition& initial_position,
                                   const Length& length,
                                   const AlleleId& allele,
                                   const Mutation::Nature& nature=Mutation::UNDEFINED)
    {
        return CNA(initial_position, length, Type::DELETION, allele, allele, nature);
    }

    /**
     * @brief Build a new amplification
     *
     * @param initial_position is the CNA initial position
     * @param length is the CNA length
     * @param nature is the nature of the CNA
     */
    static inline CNA new_deletion(const GenomicPosition& initial_position,
                                   const Length& length,
                                   const Mutation::Nature& nature=Mutation::UNDEFINED)
    {
        return CNA::new_deletion(initial_position, length, RANDOM_ALLELE, nature);
    }

    /**
     * @brief Get the CNA region
     *
     * @return the CNA region
     */
    inline GenomicRegion get_region() const
    {
        return GenomicRegion(*this, length);
    }

    /**
     * @brief Get the CNA final position
     *
     * @return a constant reference to the CNA final position
     */
    inline ChrPosition end() const
    {
        return begin()+length-1;
    }

    /**
     * @brief Get the CNA's size in the reference
     *
     * @return the CNA's size in the reference
     */
    inline const Length& size_in_ref() const
    {
        return length;
    }

    /**
     * @brief Save a CNA in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const Mutation&>(*this)
                & type
                & length
                & dest;

        switch (type) {
            case Type::AMPLIFICATION:
                archive & source;
                return;
            case Type::DELETION:
                return;
            default:
                throw std::runtime_error("CNA::save: Unsupported CNA::type "+
                                         std::to_string(static_cast<int>(type)));
        }
    }

    /**
     * @brief Load a CNA from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load CNA
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static CNA load(ARCHIVE& archive)
    {
        CNA cna;

        archive & static_cast<Mutation&>(cna)
                & cna.type
                & cna.length
                & cna.dest;

        switch (cna.type) {
            case Type::AMPLIFICATION:
                archive & cna.source;
                break;
            case Type::DELETION:
                cna.source = cna.dest;
                break;
            default:
            {
                auto type = static_cast<int>(cna.type);
                throw std::runtime_error("CNA::save: Unsupported CNA::type "+
                                         std::to_string(type));
            }
        }

        return cna;
    }
};

}   // Mutations

}   // CLONES

namespace std
{

template<>
struct less<CLONES::Mutations::CNA>
{
    bool operator()(const CLONES::Mutations::CNA &lhs,
                    const CLONES::Mutations::CNA &rhs) const;
};

/**
 * @brief Write a CNA in a output stream
 *
 * @param out is the output stream
 * @param cna is the CNA to stream
 * @return a reference of the updated stream
 */
std::ostream& operator<<(std::ostream& out, const CLONES::Mutations::CNA& cna);

};

#endif // __CLONES_CNA__