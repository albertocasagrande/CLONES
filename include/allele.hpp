/**
 * @file allele.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines allele representation
 * @version 1.5
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

#ifndef __CLONES_ALLELE__
#define __CLONES_ALLELE__

#include <map>
#include <limits>

#include "sid.hpp"
#include "genomic_region.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief A identifier type for alleles
 */
typedef size_t AlleleId;

// define a random allele identifier
#define RANDOM_ALLELE std::numeric_limits<CLONES::Mutations::AlleleId>::max()

/**
 * @brief A class to represent a fragment of an allele
 *
 */
class AlleleFragment : public GenomicRegion
{
    /**
     * @brief A structure for `AlleleFragment` data
     *
     * The class `AlleleFragment` implements copy-of-write
     * by storing its members in a `Data` object pointed by
     * `AlleleFragment:_data`.
     */
    class Data
    {    
        std::map<GenomicPosition,
                 std::shared_ptr<SID>> mutations;  //!< the fragment SIDs
    public:
        /**
         * @brief The empty constructor
         */
        Data();

        /**
         * @brief Save an allele fragment data object in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & mutations;
        }

        /**
         * @brief Load an allele fragment data object from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the load allele fragment data object
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
        inline static Data load(ARCHIVE& archive)
        {
            Data af_data;

            archive & af_data.mutations;

            return af_data;
        }

        friend AlleleFragment;
    };

    std::shared_ptr<Data> _data;    //!< The allele fragment data

    /**
     * @brief Make data exclusive
     *
     * When the data pointer is not exclusive and is referenced by many different
     * alleles, the method copy of the original data member into a data
     * object exclusively pointed by the current `AlleleFragment` object.
     *
     * @warning In order to avoid unneccessary and time-consuming copies of 
     *     the allele fragment data, the allele fragment data pointer should
     *     be exclusively maintained by the `AlleleFragment` object during
     *     the call to this method. If needed, the returned pointer can be 
     *     used for backup purpouse.
     *
     * @return the shared pointer to the original data for backup
     */
    std::shared_ptr<Data> make_data_exclusive();
public:
    /**
     * @brief The allele fragment length
     */
    using Length = GenomicRegion::Length;

    /**
     * @brief The empty constructor
     */
    AlleleFragment();

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the chromosome identifier
     * @param begin is the initial position of the allele fragment
     * @param end is the final position of the allele fragment
     */
    AlleleFragment(const ChromosomeId& chromosome_id,
                   const ChrPosition& begin, const ChrPosition& end);

    /**
     * @brief A constructor
     *
     * @param genomic_region is the genomic region of the allele fragment
     */
    explicit AlleleFragment(const GenomicRegion& genomic_region);

    /**
     * @brief Get the allele fragment SID mutations
     *
     * @return a constant reference to the allele fragment SID mutations
     */
    inline const std::map<GenomicPosition, std::shared_ptr<SID>>&
    get_mutations() const
    {
        return _data->mutations;
    }

    /**
     * @brief Check whether a mutation context is free
     *
     * @param mutation is a mutation
     * @return `true` if and only if the allele fragment does not
     *      contains SID mutations in the one-base neighborhood of
     *      `mutation`
     */
    bool has_context_free(const SID& mutation) const;

    /**
     * @brief Applied a new SID mutation
     *
     * This method tries to apply a new SID mutation. It succeeds
     * if no other SID mutations are contained in the context.
     *
     * @param mutation is the SID mutation to be applied
     * @return `true` if and only if the mutation application
     *          has succeeded
     */
    bool apply(const SID& mutation);

    /**
     * @brief Remove a SID mutation
     *
     * This method tries to remove a SID mutation from a genomic
     * position. It succeeds if the allele fragment contains a SID
     * mutation in the specified position.
     *
     * @param genomic_position is the genomic position containing the
     *          SID mutation to be removed
     * @return `true` if and only if the removal has succeeded
     */
    bool remove_mutation(const GenomicPosition& genomic_position);

    /**
     * @brief Check whether a SID is included among the fragment allele mutations
     *
     * @param mutation is the SID mutation whose inclusion among the
     *          fragment allele mutations is tested
     * @return `true` if and only if `snv` is contained among the fragment allele
     *          mutations
     */
    bool includes(const SID& mutation) const;

    /**
     * @brief Copy part of an allele fragment
     *
     * @param genomic_region is the genomic region to copy
     * @return an allelic fragment corresponding to `genomic_region` that
     *      contains all the original allele fragment SID mutations laying
     *      in `genomic_region`
     */
    AlleleFragment copy(const GenomicRegion& genomic_region) const;

    /**
     * @brief Split an allele fragment
     *
     * This method cuts an allele fragment in a position,
     * generates a new allele fragment beginning at the specified
     * position, and updates the length of the original one so
     * that the two fragments are contiguous, i.e., one of the two
     * follows the other one.
     * The SID mutations of the original allele fragment are
     * distributed on the two allele fragments according to their
     * genomic position.
     *
     * @param split_point is the position of the new allele fragment
     * @return the allele fragment originated by the split
     * @throw std::domain_error the allele fragment does not contain
     *          `split_point` or `split_point` and the allele
     *          fragment initial point are the same
     */
    AlleleFragment split(const GenomicPosition& split_point);

    /**
     * @brief Split an allele fragment
     *
     * This method cuts an allele fragment in a position,
     * generates a new allele fragment beginning at the specified
     * position, and updates the length of the original one so
     * that the two fragments are contiguous, i.e., one of the two
     * follows the other one.
     * The SID mutations of the original allele fragment are
     * distributed on the two allele fragments according to their
     * genomic position.
     *
     * @param split_point is the position of the new allele fragment
     * @return the allele fragment originated by the split
     * @throw std::domain_error the allele fragment does not contain
     *          `split_point` or `split_point` and the allele
     *          fragment initial point are the same
     */
    inline AlleleFragment split(GenomicPosition&& split_point)
    {
        return split(split_point);
    }

    /**
     * @brief Check whether the fragment contains driver mutations in a genomic region
     *
     * @param genomic_region is the genomic region to check
     * @return `true` if and only if the allele fragment contains some driver
     *      mutations in `genomic_region`
     */
    bool has_driver_mutations_in(const GenomicRegion& genomic_region) const;

    /**
     * @brief Save an allele fragment in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "CLONES AlleleFragment", 0);

        archive & static_cast<const GenomicRegion&>(*this)
                & _data;
    }

    /**
     * @brief Load an allele fragment from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load allele fragment
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static AlleleFragment load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "CLONES AlleleFragment", 0);

        AlleleFragment a_fragment;

        archive & static_cast<GenomicRegion&>(a_fragment)
                & a_fragment._data;

        return a_fragment;
    }
};

/**
 * @brief A possible representation for an allele
 */
class Allele
{
    std::map<GenomicPosition, AlleleFragment> fragments;    //!< the sequence fragments

    std::list<AlleleId> history;    //!< the allele history

    /**
     * @brief The empty constructor
     */
    Allele();
public:
    /**
     * @brief The allele length
     */
    using Length = AlleleFragment::Length;

    /**
     * @brief A constructor
     *
     * @param identifier is the identifier of the allele
     * @param history is the history of the allele
     */
    explicit Allele(const AlleleId& identifier,
                    const std::list<AlleleId>& history={});

    /**
     * @brief A constructor
     *
     * @param identifier is the identifier of the allele
     * @param chromosome_id is the chromosome identifier
     * @param begin is the initial position of the allele
     * @param end is the final position of the allele
     * @param history is the history of the allele
     */
    explicit Allele(const AlleleId& identifier, const ChromosomeId& chromosome_id,
                    const ChrPosition& begin, const ChrPosition& end,
                    const std::list<AlleleId>& history={});

    /**
     * @brief A constructor
     *
     * @param identifier is the identifier of the allele
     * @param genomic_region is the genomic region of the allele
     * @param history is the history of the allele
     */
    explicit Allele(const AlleleId& identifier, const GenomicRegion& genomic_region,
                    const std::list<AlleleId>& history={});

    /**
     * @brief Get the allele identifier
     *
     * @return the allele identifier
     */
    inline const AlleleId& get_id() const
    {
        return history.back();
    }

    /**
     * @brief Get the allele history
     *
     * @return a constant reference to the allele history
     */
    inline const std::list<AlleleId>& get_history() const
    {
        return history;
    }

    /**
     * @brief Get the allele fragments
     *
     * @return a constant reference to the allele fragments
     */
    inline const std::map<GenomicPosition, AlleleFragment>& get_fragments() const
    {
        return fragments;
    }

    /**
     * @brief Test whether the allele strictly contains a genomic position
     *
     * @param genomic_position is a genomic position
     * @return `true` if and only if one of the allele fragments
     *          strictly contains `genomic_position`, i.e.,
     *          `genomic_position` belongs to the allele and it is
     *          not in one of its borders
     */
    bool strictly_contains(const GenomicPosition& genomic_position) const;

    /**
     * @brief Test whether the allele contains a genomic position
     *
     * @param genomic_position is a genomic position
     * @return `true` if and only if one of the allele fragments
     *          contains `genomic_position`
     */
    bool contains(const GenomicPosition& genomic_position) const;

    /**
     * @brief Test whether the allele contains a genomic region
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if one of the allele fragments
     *          contains `genomic_regions`
     */
    bool contains(const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether a SID is included among the allele mutations
     *
     * @param mutation is the SID whose inclusion among the allele mutations
     *      is tested
     * @return `true` if and only if `mutation` is contained among the
     *      allele mutations
     */
    bool includes(const SID& mutation) const;

    /**
     * @brief Check whether a mutation context is free
     *
     * @param mutation is a mutation
     * @return `true` if and only if the allele does not contains
     *      SID mutations in the one-base neighborhood of `mutation`
     */
    bool has_context_free(const SID& mutation) const;

    /**
     * @brief Apply a new SID mutation
     *
     * This method tries to apply a new SID mutation. It succeeds
     * if no other SIDs are contained in the context.
     *
     * @param mutation is the SID to applied
     * @return `true` if and only if the SID application
     *          has succeeded
     */
    bool apply(const SID& mutation);

    /**
     * @brief Remove a SID mutation
     *
     * This method tries to remove a SID mutation from a genomic position.
     * It succeeds if the allele fragment contains a SID mutation in the
     * specified position.
     *
     * @param genomic_position is the genomic position containing the
     *          SID mutation to be removed
     * @return `true` if and only if the removal has succeeded
     */
    bool remove_mutation(const GenomicPosition& genomic_position);

    /**
     * @brief Copy part of an allele
     *
     * @param new_allele_id is the identifier of the resulting allele
     * @param genomic_region is the genomic region to copy
     * @return an allelic fragment corresponding to `genomic_region`
     *      that contains all the original allele SID mutations
     *      laying in `genomic_region`
     */
    Allele copy(const AlleleId& new_allele_id,
                const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether the allele contains driver mutations in a genomic region
     *
     * @param genomic_region is the genomic region to check
     * @return `true` if and only if the allele contains some driver
     *      mutations in `genomic_region`
     */
    bool has_driver_mutations_in(const GenomicRegion& genomic_region) const;

    /**
     * @brief Remove part of an allele
     *
     * This method tries to remove the part of an allele
     * corresponding to a genomic region. It succeeds only if the
     * original allele contains the specified genomic
     * region.
     *
     * @param genomic_region is the genomic region to remove
     * @return `true` if and only if the allele contains
     *          `genomic_region`
     */
    bool remove(const GenomicRegion& genomic_region);

    /**
     * @brief Get the allele SID mutations
     *
     * @return the allele SID mutations
     */
    std::map<GenomicPosition, std::shared_ptr<SID>> get_mutations() const;

    /**
     * @brief Get the size of the allele
     *
     * @return the size of the allele
     */
    Length size() const;

    /**
     * @brief Get the string representation of an allele identifier
     *
     * @param allele_id is the allele identifier to be formatted
     * @return the string representation of `allele_id`
     */
    static std::string format_id(const CLONES::Mutations::AlleleId& allele_id);

    /**
     * @brief Copy genomic structure
     *
     * This method copies the genomic structure of the current objects.
     * It returns an `Allele` object that has the same fragments of the current
     * objects, but misses the original SID mutations and indels.
     *
     * @return an `Allele` object that has the same fragments of the current
     *      objects, but misses the original SID mutations
     */
    Allele copy_structure() const;

    /**
     * @brief Save an allele in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & fragments
                & history;
    }

    /**
     * @brief Load an allele from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load allele
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static Allele load(ARCHIVE& archive)
    {
        Allele allele(0);

        archive & allele.fragments
                & allele.history;

        return allele;
    }
};

}   // Mutations

}   // CLONES

namespace std
{

/**
 * @brief Write allele fragment data in a stream
 *
 * @param os is the output stream
 * @param allele_fragment is the allele fragment to be written
 * @return a reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::AlleleFragment& allele_fragment);


/**
 * @brief Write allele data in a stream
 *
 * @param os is the output stream
 * @param allele is the allele to be written
 * @return a reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::Allele& allele);

} // std

#endif // __CLONES_ALLELE__