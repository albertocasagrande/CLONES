/**
 * @file genome_fragment.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines genome DNA fragments
 * @version 1.2
 * @date 2026-07-07
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

#ifndef __CLONES_GENOME_FRAGMENTS__
#define __CLONES_GENOME_FRAGMENTS__

#include <string>
#include <vector>

#include "sid.hpp"
#include "genomic_position.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief Matching type
 *
 * This enum class represents the possible match/mismatch
 * types when aligning two sequences
 */
enum class MatchingType
{
    MATCH,      //!< A match
    MISMATCH,   //!< A mismatch
    INSERTION,  //!< An insertion
    DELETION    //!< A deletion
};

/**
 * @brief A class for representing genome DNA fragments
 */
class GenomeFragment
{
    /**
     * @brief Fragment direction
     */
    enum class Direction
    {
        FORWARD,    //!< Forward direction
        BACKWARD    //!< Backward direction
    };

    /**
     * @brief A class for indices of bases in mutations
     *
     * Any valid object of this class refers to a base in the
     * altered sequence of a mutation laying of a genome
     * fragment. The objects of this class may also be not
     * valid if it does not refer to any base.
     */
    struct MutationBaseIndex
    {
        size_t index;       //!< The mutation index
        size_t base_pos;    //!< The base position

        /**
         * @brief The empty constructor
         */
        inline MutationBaseIndex():
            index(std::numeric_limits<size_t>::max())
        {}

        /**
         * @brief A constructor
         *
         * @param index is the index of the mutation
         * @param base_pos
         */
        inline MutationBaseIndex(const size_t& index,
                                 const size_t& base_pos):
            index(index), base_pos(base_pos)
        {}

        /**
         * @brief Check if this mutation base index is valid
         *
         * @return `true` if and only if this object is referring
         *      to a base in the altered sequence of the
         *      `index`-th mutation
         */
        inline bool is_valid() const
        {
            return (index != std::numeric_limits<size_t>::max());
        }
    };

    GenomicPosition genomic_position;   //!< The genome fragment position

    std::vector<SID> mutations;         //!< The mutations in the fragment

    std::string nucleotides;            //!< The nucleotide sequence
    std::vector<MatchingType> alignment;    //!< The alignment to the reference
    std::vector<size_t> alignment_index;    //!< The alignment index of each sequence base
    std::vector<MutationBaseIndex> mutation_index;     //!< The mutation index of each sequence base

    /**
     * @brief An in-order forward iterator for mutations
     *
     * This class simulates an iterator over the merged somatic
     * and germline mutations.
     */
    class MutationIterator
    {
        std::map<GenomicPosition, std::shared_ptr<SID>> const* somatic;     //!< The somatic map
        std::map<GenomicPosition, std::shared_ptr<SID>> const* germline;      //!< The germline map

        std::map<GenomicPosition, std::shared_ptr<SID>>::const_iterator s_it;    //!< The somatic iterator
        std::map<GenomicPosition, std::shared_ptr<SID>>::const_iterator g_it;    //!< The germline iterator

        Direction direction;    //!< Iterator direction

        bool p_begin;   //!< A Boolean flag that holds when the somatic iterator has already reached the begin
        bool p_end;     //!< A Boolean flag that holds when the somatic iterator has already reached the end
        bool g_begin;   //!< A Boolean flag that holds when the germline iterator has already reached the begin
        bool g_end;     //!< A Boolean flag that holds when the germline iterator has already reached the end

        bool s_it_curr;   //!< A Boolean flag that holds when the referenced mutation is a somatic mutation

        /**
         * @brief Set the current mutation
         */
        void set_current_mutation();

        /**
         * @brief Construct a new Forward Mutation Iterator object
         *
         * @param germline is the germline SID map
         * @param somatic is the somatic SID map
         * @param germline_it is an germline SID map iterator
         * @param somatic_it is an somatic SID map iterator
         */
        MutationIterator(const std::map<GenomicPosition, std::shared_ptr<SID>>& germline,
                         const std::map<GenomicPosition, std::shared_ptr<SID>>& somatic,
                         const std::map<GenomicPosition, std::shared_ptr<SID>>::const_iterator& germline_it,
                         const std::map<GenomicPosition, std::shared_ptr<SID>>::const_iterator& somatic_it);

    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   std::map<GenomicPosition, std::shared_ptr<SID>>::value_type;
        using allocator_type	=   std::map<GenomicPosition, std::shared_ptr<SID>>::allocator_type;
        using pointer           =   std::allocator_traits<allocator_type>::pointer;
        using const_pointer     =   std::allocator_traits<allocator_type>::const_pointer;
        using reference         =   const value_type&;
        using iterator_category =   std::bidirectional_iterator_tag;

        /**
         * @brief The empty constructor
         */
        MutationIterator();

        /**
         * @brief
         *
         * @param germline
         * @param somatic
         * @param genomic_position
         * @return MutationIterator
         */
        static MutationIterator lower_bound(const std::map<GenomicPosition, std::shared_ptr<SID>>& germline,
                                            const std::map<GenomicPosition, std::shared_ptr<SID>>& somatic,
                                            const GenomicPosition& genomic_position);

        /**
         * @brief Get the pointed direction
         *
         * @return the pointed direction
         */
        inline reference operator*() const
        {
            return *(s_it_curr?s_it:g_it);
        }

        /**
         * @brief Get the pointer to the direction
         *
         * @return the pointer to the direction
         */
        inline const_pointer operator->() const
        {
            return (s_it_curr?s_it.operator->():g_it.operator->());
        }

        /**
         * @brief Prefix increment
         *
         * @return a reference to the updated
         *      constant iterator
         */
        MutationIterator& operator++();

        /**
         * @brief Postfix increment
         *
         * @return a copy of the original iterator
         */
        MutationIterator operator++(int);

        /**
         * @brief Prefix decrement
         *
         * @return a reference to the updated
         *      constant iterator
         */
        MutationIterator& operator--();

        /**
         * @brief Postfix decrement
         *
         * @return a copy of the original iterator
         */
        MutationIterator operator--(int);

        /**
         * @brief Test whether the begin of the iterator has been reached
         *
         * @return `true` if and only if the iterator reached the begin
         *      of the germline and somatic maps
         */
        inline bool is_begin() const
        {
            return p_begin && g_begin;
        }

        /**
         * @brief Test whether the end of the iterator has been reached
         *
         * @return `true` if and only if the iterator reached the end
         *      of the germline and somatic maps
         */
        inline bool is_end() const
        {
            return p_end && g_end;
        }

        /**
         * @brief Compare two iterators
         *
         * @param a is the iterator to compare
         * @return `true` if and only if the two mutation iterators
         *      refers to the same mutation
         */
        inline bool operator==(const MutationIterator& a) const
        {
            return s_it == a.s_it && g_it == a.g_it;
        }

        /**
         * @brief Compare two iterators
         *
         * @param a is the iterator to compare
         * @return `true` if and only if the two mutation iterators
         *      refers to two distinct mutations
         */
        inline bool operator!=(const MutationIterator& a) const
        {
            return s_it != a.s_it || g_it != a.g_it;
        }
    };

    /**
     * @brief Copy a fragment of the reference
     *
     * This method copies a reference fragment in the genome fragment.
     * It also update the positions of both the genome fragment and
     * the reference to the end so they correspond to the end of the
     * genome fragment and the corresponding position in the reference,
     * respectively.
     *
     * @param reference is the reference
     * @param up_to_pos is reference position up to which the reference
     *   must be copied
     * @param genome_pos is the first position of the genome fragment
     *   to be filled
     * @param reference_pos is the first position of the reference to be
     *   copied
     */
    void copy_reference(const std::string& reference, const size_t length,
                        size_t& genome_pos, size_t& reference_pos);

    /**
     * @brief Remove the mutation in a genome fragment position
     *
     * @param pos is the genome fragment position in which the
     *    to-be-removed mutation lays
     */
    void remove_mutation(const size_t& pos);
public:
    /**
     * @brief The empty constructor
     */
    GenomeFragment();

    /**
     * @brief A genome fragment constructor
     *
     * @param reference is the reference sequence
     * @param germline is the position-mutation map for germline mutations
     * @param somatic is the position-mutation map for somatic mutations
     * @param begin_pos is the aimed position of the first reference base
     *   of the genome fragment
     * @param size is the aimed genome fragment size
     */
    GenomeFragment(const std::string& reference,
                    const std::map<GenomicPosition, std::shared_ptr<SID>>& germline,
                    const std::map<GenomicPosition, std::shared_ptr<SID>>& somatic,
                    const GenomicPosition& begin_pos,
                    const size_t& size);

    /**
     * @brief A genome fragment constructor
     *
     * @param reference_fragment is a fragment of the reference sequence
     * @param fragment_offset is the offset of the reference fragment with
     *   respect to the whole reference sequence
     * @param germline is the position-mutation map for germline mutations
     * @param somatic is the position-mutation map for somatic mutations
     * @param begin_pos is the aimed position of the first reference base
     *   of the genome fragment
     * @param size is the aimed genome fragment size
     */
    GenomeFragment(const std::string& reference_fragment,
                    const size_t& fragment_offset,
                    const std::map<GenomicPosition, std::shared_ptr<SID>>& germline,
                    const std::map<GenomicPosition, std::shared_ptr<SID>>& somatic,
                    const GenomicPosition& begin_pos,
                    const size_t& size);

    /**
     * @brief Count the number of mismatched
     *
     * @return the number of mismatches in the genome
     *   fragment
     */
    size_t Hamming_distance() const;

    /**
     * @brief Get the genome fragment sequence
     *
     * @return a constant reference to the genome
     *   fragment sequence
     */
    inline const std::string& get_sequence() const
    {
        return nucleotides;
    }

    /**
     * @brief Get the n-th nucleotide in the genome fragment
     *
     * @param pos is the position of the aimed nucleotide
     * @return a constant reference to the aimed nucleotide
     */
    inline const char& operator[](const size_t& pos) const
    {
        return nucleotides[pos];
    }

    /**
     * @brief Get the covered reference region
     *
     * This method returns the region in the reference genome that the
     * genome fragment would cover if the latter contained no deletions.
     *
     * @see get_covered_reference_region()
     * @return The reference region covered by this genome fragment
     */
    GenomicRegion get_covered_reference_region() const;

    /**
     * @brief Get the reference regions in the genome fragment
     *
     * This method returns the list of regions in the reference genome with
     * bases in the current genome fragment. The returned regions are
     * contained in the region returned by the method
     * `get_covered_reference_region()`. However, there may be bases that
     * belong to the latter, but not to the output of this method, because
     * they are missing from the genome fragment due to deletions.
     *
     * @see get_covered_reference_region()
     * @return The list of the reference regions that are present
     *    in this this genome fragment
     */
    std::list<GenomicRegion> get_contained_reference_regions() const;

    /**
     * @brief Compute the CIGAR string
     *
     * The CIGAR string a code that represents the alignment of a sequence
     * over a reference in SAM files (see [1]). This method considers a set
     * of SIDs, a genomic position, and a size, and it produces the CIGAR
     * code corresponding to a genome fragment whose sequence correspond
     * to that of the reference genome from the specified position with the
     * exception of positions in which the given SIDs were applied.
     *
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format
     *     Specification Working Group, 22 August 2022,
     *     http://samtools.github.io/hts-specs/SAMv1.pdf
     *
     * @return the CIGAR code corresponding to the mismatch vector
     */
    std::string get_CIGAR() const;

    /**
     * @brief Alter one of the bases
     *
     * @param pos is the position of the base
     * @param base is the new base
     */
    void alter_base(const size_t pos, const char base);

    /**
     * @brief Get the mutation on this genome fragment
     *
     * @return a constant reference to the mutation laying on this genome
     *   fragment
     */
    inline const std::vector<SID>& get_mutations() const
    {
        return mutations;
    }

    /**
     * @brief Get the genomic position of this genome fragment
     *
     * @return a constant reference to the genomic position of the first
     *   base of this genome fragment
     */
    inline const GenomicPosition& get_genomic_position() const
    {
        return genomic_position;
    }

    /**
     * @brief Get the genome fragment size
     *
     * @return the genome fragment size
     */
    inline size_t size() const {
        return nucleotides.size();
    }
};

}   // Mutations

}   // CLONES


#endif // __CLONES_GENOME_FRAGMENTS__
