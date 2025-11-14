/**
 * @file genomic_region.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines genomic regions
 * @version 1.4
 * @date 2025-11-14
 *
 * @copyright Copyright (c) 2023-2025
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

#ifndef __RACES_GENOMIC_REGION__
#define __RACES_GENOMIC_REGION__

#include <ostream>

#include "genomic_position.hpp"

namespace RACES
{

namespace Mutations
{

/**
 * @brief A class to represent SNVs, Insertions, and Deletions
 */
struct SID;

/**
 * @brief A genomic region is a interval of genomic positions
 *
 */
struct GenomicRegion
{
    /**
     * @brief The genomic region length
     */
    using Length = ChrPosition;
protected:
    GenomicPosition initial_pos;    //!< the genomic region initial position
    Length length;                  //!< the genomic region length
public:
    /**
     * @brief The empty constructor
     */
    GenomicRegion();

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the identifier of the chromosome on which the region lays
     * @param length is the length of the genomic region
     * @throw std::domain_error `length` is 0
     */
    GenomicRegion(const ChromosomeId chromosome_id, const Length length);

    /**
     * @brief A constructor
     *
     * @param initial_pos is the genomic region initial position
     * @param length is the genomic region length
     * @throw std::domain_error `length` is 0
     */
    GenomicRegion(const GenomicPosition initial_pos, const Length length);

    /**
     * @brief Get the identifier of the chromosome containing the genomic region
     *
     * @return a constant reference to the identifier of the chromosome
     *          containing the genomic region
     */
    inline const ChromosomeId& get_chromosome_id() const
    {
        return initial_pos.chr_id;
    }

    /**
     * @brief Get the genomic region initial position
     *
     * @return a constant reference to the genomic region initial position
     */
    inline const GenomicPosition& get_initial_position() const
    {
        return initial_pos;
    }

    /**
     * @brief Get the genomic region final position
     *
     * @return the genomic region final position
     */
    inline GenomicPosition get_final_position() const
    {
        return GenomicPosition(initial_pos.chr_id, initial_pos.position+length-1);
    }

    /**
     * @brief Get the genomic region initial position in the chromosome
     *
     * @return the genomic region initial position in the chromosome
     */
    inline const ChrPosition& begin() const
    {
        return initial_pos.position;
    }

    /**
     * @brief Get the genomic region final position in the chromosome
     *
     * @return the genomic region final position in the chromosome
     */
    inline ChrPosition end() const
    {
        return initial_pos.position+length-1;
    }

    /**
     * @brief Get the genomic region length
     *
     * @return a constant reference to the genomic region length
     */
    inline const Length& size() const
    {
        return length;
    }

    /**
     * @brief Test whether a region follows another region
     *
     * Let `A` and `B` two regions . "`A` follows `B`" when the initial position
     * of `A` immediately follows the final position of `B`.
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current genomic region begins just
     *          after the end of `genomic_region`
     */
    inline bool follows(const GenomicRegion& genomic_region) const
    {
        return ((initial_pos.chr_id == genomic_region.initial_pos.chr_id)
                && (genomic_region.initial_pos.position + genomic_region.length == initial_pos.position));
    }

    /**
     * @brief Test whether a region follows another region
     *
     * Let `A` and `B` two regions . "`A` follows `B`" when the initial position
     * of `A` immediately follows the final position of `B`.
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current genomic region begins just
     *          after the end of `genomic_region`
     */
    inline bool follows(GenomicRegion&& genomic_region) const
    {
        return follows(genomic_region);
    }

    /**
     * @brief Test whether a region precedes another region
     *
     * Let `A` and `B` two regions . "`A` precedes `B`" when the final position
     * of `A` immediately follows the initial position of `B`.
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current `genomic_region` ends just
     *          before the beginning of the current genomic region
     */
    inline bool precedes(const GenomicRegion& genomic_region) const
    {
        return genomic_region.follows(*this);
    }

    /**
     * @brief Test whether a region precedes another region
     *
     * Let `A` and `B` two regions . "`A` precedes `B`" when the final position
     * of `A` immediately follows the initial position of `B`.
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current `genomic_region` ends just
     *          before the beginning of the current genomic region
     */
    inline bool precedes(GenomicRegion&& genomic_region) const
    {
        return genomic_region.follows(*this);
    }

    /**
     * @brief Check whether two genomic region overlap
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region and
     *       `genomic_region`
     */
    bool overlaps(const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether a genomic region overlaps a region container
     *
     * @tparam Container is the container template
     * @param regions is the container of the regions
     * @return true iff any of the regions in `regions` overlaps the
     *      current object
     */
    template <template<typename, typename> class Container>
    bool overlaps(const Container<GenomicRegion, std::allocator<GenomicRegion>>& regions) const
    {
        for (const auto& region: regions) {
            if (overlaps(region)) {
                return true;
            }
        }

        return false;
    }

    /**
     * @brief Check whether two genomic region overlap
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region and
     *       `genomic_region`
     */
    inline bool overlaps(GenomicRegion&& genomic_region) const
    {
        return this->overlaps(genomic_region);
    }

    /**
     * @brief Check whether a region ends before another one begin
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region ends before
     *          `genomic_region`'s initial position
     */
    inline bool ends_before(const GenomicRegion& genomic_region) const
    {
        return ends_before(genomic_region.get_initial_position());
    }

    /**
     * @brief Check whether a region ends before another one begin
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region ends before
     *          `genomic_region`'s initial position
     */
    inline bool ends_before(GenomicRegion&& genomic_region) const
    {
        return ends_before(genomic_region);
    }

    /**
     * @brief Check whether a region begin after another one end
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region ends before
     *          `genomic_position`
     */
    inline bool begins_after(const GenomicRegion& genomic_region) const
    {
        return genomic_region.ends_before(*this);
    }

    /**
     * @brief Check whether a position is contained in the genomic region
     *
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `genomic_position` is contained in the
     *      genomic region
     */
    inline bool contains(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position <= genomic_position.position)
                && (genomic_position.position+1 <= initial_pos.position+length));
    }

    /**
     * @brief Check whether a position is contained in the genomic region
     *
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `genomic_position` is contained in the
     *      genomic region
     */
    inline bool contains(GenomicPosition&& genomic_position) const
    {
        return contains(genomic_position);
    }

    /**
     * @brief Check whether a mutation lays in the genomic region
     *
     * @param mutation is a SID mutation
     * @return `true` if and only if `mutation` completely lays in the current
     *      genomic region
     */
    bool contains(const SID& mutation) const;

    /**
     * @brief Check whether a genomic region is contained in another genomic region
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if `genomic_region` is contained in the
     *      genomic region
     */
    inline bool contains(const GenomicRegion& genomic_region) const
    {
        return ((genomic_region.get_chromosome_id()==initial_pos.chr_id)
                && (begin() <= genomic_region.begin())
                && (genomic_region.end() <= end()));
    }

    /**
     * @brief Check whether a position is strictly contained in the genomic region
     *
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `genomic_position` is strictly contained in the
     *      genomic region
     */
    inline bool strictly_contains(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position < genomic_position.position)
                && (genomic_position.position+1 < initial_pos.position+length));
    }

    /**
     * @brief Check whether a position is strictly contained in the genomic region
     *
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `genomic_position` is strictly contained in the
     *      genomic region
     */
    inline bool strictly_contains(GenomicPosition&& genomic_position) const
    {
        return strictly_contains(genomic_position);
    }

    /**
     * @brief Check whether a genomic region is strictly contained in the genomic region
     *
     * @param genomic_region is the genomic region whose inclusion is questioned
     * @return `true` if and only if `genomic_region` is strictly contained in the
     *      genomic region
     */
    inline bool strictly_contains(const GenomicRegion& genomic_region) const
    {
        return ((genomic_region.get_chromosome_id()==initial_pos.chr_id)
                && (begin()+1 <= genomic_region.begin())
                && (genomic_region.end()+1 <= end()));
    }

    /**
     * @brief Check whether a mutation is strictly contained in the genomic region
     *
     * @param mutation is a SID mutation
     * @return `true` if and only if `mutation` is strictly contained in the
     *      genomic region
     */
    bool strictly_contains(const SID& mutation) const;

    /**
     * @brief Check whether a region ends before a position
     *
     * @param genomic_position is a genomic position
     * @return `true` if and only if the current region ends before
     *          `genomic_position`
     */
    inline bool ends_before(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position+length < genomic_position.position+1));
    }

    /**
     * @brief Check whether a region ends before a position
     *
     * @param genomic_position is a genomic position
     * @return `true` if and only if the current region ends before
     *          `genomic_position`
     */
    inline bool ends_before(GenomicPosition&& genomic_position) const
    {
        return ends_before(genomic_position);
    }

    /**
     * @brief Check whether a region begins after a position
     *
     * @param genomic_position is a genomic position
     * @return `true` if and only if the current region begins after
     *          `genomic_position`
     */
    inline bool begins_after(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position > genomic_position.position));
    }

    /**
     * @brief Check whether a region begins after a position
     *
     * @param genomic_position is a genomic position
     * @return `true` if and only if the current region begins after
     *          `genomic_position`
     */
    inline bool begins_after(GenomicPosition&& genomic_position) const
    {
        return begins_after(genomic_position);
    }

    /**
     * @brief Split a genomic region
     *
     * This method cuts a genomic region in a position, generates a
     * new genomic region beginning at the specified position, and
     * updates the length of the considered genomic region so that
     * the two fragments are contiguous, i.e., one of the two
     * follows the other one.
     *
     * @param split_point is the position of the new genomic region
     * @return the fragment originated by the split
     * @throw std::domain_error the genomic region does not contain
     *          `split_point` or `split_point` and the genomic region
     *          initial point are the same
     */
    GenomicRegion split(const GenomicPosition& split_point);

    /**
     * @brief Split a genomic region
     *
     * This method cuts a genomic region in a position, generates a
     * new genomic region beginning at the specified position, and
     * updates the length of the considered genomic region so that
     * the two fragments are contiguous, i.e., one of the two
     * follows the other one.
     *
     * @param split_point is the position of the new genomic region
     * @return the fragment originated by the split
     * @throw std::domain_error the genomic region does not contain
     *          `split_point` or `split_point` and the genomic region
     *          initial point are the same
     */
    inline GenomicRegion split(GenomicPosition&& split_point)
    {
        return split(split_point);
    }

    /**
     * @brief Join two contiguous genomic regions
     *
     * Two genomic regions are contiguous when the inital position
     * of one of the two is the initial position of the other
     * one plus the latter's length plus one.
     * This method joins the current genomic region and a contiguous
     * one. It sets the initial position of the current object
     * to the first of the two initial positions, and the its
     * length to the sum of the two lengths.
     *
     * @param contiguous_region is a genomic region contiguous to
     *          the current one
     * @return a reference to the updated genomic region
     */
    GenomicRegion& join(GenomicRegion& contiguous_region);

    /**
     * @brief Save a genomic region in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & initial_pos
                & length;
    }

    /**
     * @brief Load a genomic region from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load genomic region
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static GenomicRegion load(ARCHIVE& archive)
    {
        GenomicRegion g_region;

        archive & g_region.initial_pos
                & g_region.length;

        return g_region;
    }
};


/**
 * @brief Test whether two `GenomicRegion` objects are the same
 *
 * @param lhs is the left-hand side of the relation
 * @param rhs is the right-hand side of the relation
 * @return `true` is and only if `lhs` and `rhs` represent the same
 *      allele fragment, i.e., have the same allele fragments and allele history
 */
inline bool operator==(const RACES::Mutations::GenomicRegion& lhs,
                       const RACES::Mutations::GenomicRegion& rhs)
{
    return (lhs.get_initial_position() == rhs.get_initial_position())
            && (lhs.size() == rhs.size());
}

/**
 * @brief Test whether two `GenomicRegion` objects differ
 *
 * @param lhs is the left-hand side of the relation
 * @param rhs is the right-hand side of the relation
 * @return `false` is and only if `lhs` and `rhs` represent the same
 *      allele fragment, i.e., have the same allele fragments and allele history
 */
inline bool operator!=(const RACES::Mutations::GenomicRegion& lhs,
                       const RACES::Mutations::GenomicRegion& rhs)
{
    return !(lhs == rhs);
}

/**
 * @brief Split a set of genomic regions by chromosome id
 *
 * @param[in] genomic_regions is the set of genomic region to be split
 * @return a map that associates a chromosome id to the the set of genomic regions
 *     laying in the corresponding chromosome
 */
std::map<ChromosomeId, std::set<GenomicRegion>>
split_by_chromosome_id(const std::set<GenomicRegion>& genomic_regions);

}   // Mutations

}   // RACES

namespace std
{

template<>
struct less<RACES::Mutations::GenomicRegion>
{
    inline bool operator()(const RACES::Mutations::GenomicRegion &lhs,
                           const RACES::Mutations::GenomicRegion &rhs) const
    {
        return ((lhs.get_chromosome_id()<rhs.get_chromosome_id())
                || ((lhs.get_chromosome_id()==rhs.get_chromosome_id())
                    && (lhs.begin()<rhs.begin()))
                || ((lhs.get_chromosome_id()==rhs.get_chromosome_id())
                    && (lhs.begin()==rhs.begin())
                    && (lhs.size()<rhs.size())));
    }
};

/**
 * @brief Write a genomic region in a output stream
 *
 * @param out is the output stream
 * @param genomic_region is the genomic region to stream
 * @return a reference of the updated stream
 */
std::ostream& operator<<(std::ostream& out, const RACES::Mutations::GenomicRegion& genomic_region);

}

#endif // __RACES_GENOMIC_REGION__
