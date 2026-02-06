/**
 * @file read.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines sequencing reads
 * @version 1.6
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

#ifndef __CLONES_READ__
#define __CLONES_READ__

#include <string>
#include <vector>

#include "sid.hpp"
#include "genomic_position.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief The sequencing simulation namespace
 */
namespace SequencingSimulations
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
 * @brief A class for representing the simulated reads
 */
class Read
{
    /**
     * @brief Read build direction
     */
    enum class Direction
    {
        FORWARD,    //!< Forward direction
        BACKWARD    //!< Backward direction
    };

    /**
     * @brief A class for indices of bases in mutations
     *
     * Any valid object of this class refers to a base in
     * the altered sequence of a mutation laying a read.
     * The objects of this class may also be not valid if
     * it does not refer to any base.
     */
    struct MutationBaseIndex
    {
        size_t index;       //!< The read mutation index
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
         * @param index is the index of the mutation in the read
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
         *      `index`-th mutation laying in the read
         */
        inline bool is_valid() const
        {
            return (index != std::numeric_limits<size_t>::max());
        }
    };

    GenomicPosition genomic_position;   //!< The read position

    std::vector<SID> mutations;         //!< The read mutations

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
     * @brief Copy a fragment of the reference in the read
     *
     * This method copies a reference fragment in the read.
     * It also update the indices of both the read and the
     * reference.
     *
     * @param reference is the reference
     * @param up_to_index is reference index up to which
     *      the reference must be copied
     * @param read_idx is the first position of the read to
     *      be filled by the read
     * @param reference_idx is the first position of the
     *      reference to be copied
     */
    void copy_reference(const std::string& reference,
                        const size_t length,
                        size_t& read_idx, size_t& reference_idx);

    /**
     * @brief Remove the mutation in a read position
     *
     * @param pos is the read position in which the to-be-removed
     *    mutation lays
     */
    void remove_mutation(const size_t& pos);
public:
    /**
     * @brief The empty constructor
     */
    Read();

    /**
     * @brief A read constructor
     *
     * @param reference is the reference sequence
     * @param germline is the position-mutation map for germline mutations
     * @param somatic is the position-mutation map for somatic mutations
     * @param genomic_position is the aimed position of the first base in the read
     * @param read_size is the aimed read size
     */
    Read(const std::string& reference,
         const std::map<GenomicPosition, std::shared_ptr<SID>>& germline,
         const std::map<GenomicPosition, std::shared_ptr<SID>>& somatic,
         const GenomicPosition& genomic_position,
         const size_t& read_size);

    /**
     * @brief Count the number of mismatched
     *
     * @return the number of mismatches in the read
     */
    size_t Hamming_distance() const;

    /**
     * @brief Get the read sequence
     *
     * @return a constant reference to the read sequence
     */
    inline const std::string& get_sequence() const
    {
        return nucleotides;
    }

    /**
     * @brief Get the n-th nucleotide in the read
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
     * This method returns the region in the reference genome that the read
     * would cover if the latter contained no deletions.
     *
     * @see get_covered_reference_region()
     * @return The reference region covered by this read
     */
    GenomicRegion get_covered_reference_region() const;

    /**
     * @brief Get the reference regions in the read
     *
     * This method returns the list of regions in the reference genome with
     * bases in the current read. The returned regions are contained in the
     * region returned by the method `get_covered_reference_region()`.
     * However, there may be bases that belong to the latter, but not to the
     * output of this method, because they are missing from the read due to
     * deletion.
     *
     * @see get_reference_regions_in_read()
     * @return The list of the reference regions that are present
     *    in this this read
     */
    std::list<GenomicRegion> get_reference_regions_in_read() const;

    /**
     * @brief Compute the CIGAR string of a genomic region potentially affected by SIDs
     *
     * The CIGAR string a code that represents the alignment of a sequence over a reference
     * in SAM files (see [1]). This method considers a set of SIDs, a genomic position,
     * and a size, and it produces the CIGAR code corresponding to a simulated read whose
     * sequence correspond to that of the reference genome from the specified position with
     * the exception of positions in which the given SIDs were applied.
     *
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format Specification
     *     Working Group, 22 August 2022, http://samtools.github.io/hts-specs/SAMv1.pdf
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
     * @brief Get the mutation on this read
     *
     * @return a constant reference to the mutation laying on
     *      this read
     */
    inline const std::vector<SID>& get_mutations() const
    {
        return mutations;
    }

    /**
     * @brief Get the genomic position of this read
     *
     * @return a constant reference to the genomic position
     *      of the first base of this read
     */
    inline const GenomicPosition& get_genomic_position() const
    {
        return genomic_position;
    }

    /**
     * @brief Get the read size
     *
     * @return the read size
     */
    inline size_t size() const {
        return nucleotides.size();
    }
};

}   // SequencingSimulations

}   // Mutations

}   // CLONES


#endif // __CLONES_READ__
