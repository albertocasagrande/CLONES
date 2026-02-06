/**
 * @file genome_mutations.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines genome and chromosome data structures
 * @version 1.17
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

#ifndef __CLONES_GENOME_MUTATIONS__
#define __CLONES_GENOME_MUTATIONS__

#include <map>
#include <list>
#include <fstream>
#include <memory>   // for std::shared_ptr
#include <utility>  // for std::pair

#include "mutation.hpp"
#include "mutation_spec.hpp"
#include "mutation_list.hpp"

#include "allele.hpp"
#include "cna.hpp"

#include "cell.hpp"

#include "progress_bar.hpp"

namespace CLONES
{

/**
 * @brief A namespace for classes related to mutations
 */
namespace Mutations
{

/**
 * @brief A class to represent allelic type
 */
using AllelicType = std::vector<uint32_t>;

/**
 * @brief A class to represent the mutations in a chromosome
 */
class ChromosomeMutations;

class ChromosomeMutations
{
    /**
     * @brief The chromosome length type
     */
    using Length = GenomicRegion::Length;
private:

    /**
     * @brief A structure for `ChromosomeMutations` data
     *
     * The class `ChromosomeMutations` implements copy-of-write
     * by storing all its data in a `Data` object pointed by
     * `ChromosomeMutations:_data`.
     */
    class Data
    {
        ChromosomeId identifier;    //!< the chromosome identifier
        Length length;              //!< the chromosome length
        Length allelic_length;      //!< the sum of the lengths of all the alleles

        std::map<AlleleId, Allele> alleles;  //!< the chromosome alleles

        std::list<std::shared_ptr<CNA>> CNAs;        //!< the occurred CNAs

        AlleleId next_allele_id;    //!< the identifier of the next allele

    public:
        /**
         * @brief The empty constructor
         */
        Data();

        /**
         * @brief A constructor
         *
         * @param identifier is the chromosome identifier
         * @param size is the chromosome size
         * @param num_of_alleles is the initial number of alleles
         */
        Data(const ChromosomeId& identifier,
             const Length& size,
             const size_t& num_of_alleles);

        /**
         * @brief A constructor
         *
         * @param chromosome_region is the chromosome region
         * @param num_of_alleles is the initial number of alleles
         */
        Data(const GenomicRegion& chromosome_region,
             const size_t& num_of_alleles);

        /**
         * @brief Save chromosome mutation data in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & identifier
                    & length
                    & allelic_length
                    & alleles
                    & CNAs
                    & next_allele_id;
        }

        /**
         * @brief Load chromosome mutation data from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the loaded chromosome mutation data
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
        inline static Data load(ARCHIVE& archive)
        {
            Data data;

            archive & data.identifier
                    & data.length
                    & data.allelic_length
                    & data.alleles
                    & data.CNAs
                    & data.next_allele_id;

            return data;
        }

        friend ChromosomeMutations;
    };

    std::shared_ptr<Data> _data;   //!< The chromosome mutation data

    /**
     * @brief Apply a CNA
     *
     * This function applies a CNA. Whenever the `CNA::dest` member of the
     * parameter is set to `RANDOM_ALLELE`, the function change its value to
     * the destination allele identifier.
     *
     * @param[in,out] cna is the CNA to be applied
     * @return `true` if and only if the CNA has been successfully applied
     */
    bool apply_CNA(CNA& cna);

    /**
     * @brief Make data exclusive
     *
     * When the data pointer is not exclusive and is referenced by many different
     * chromosomes, the method copy of the original data member into a data
     * object exclusively pointed by the current `ChromosomeMutations` object.
     *
     * @warning In order to avoid unnecessary and time-consuming copies of the
     *     chromosome mutation data, the chromosome mutation data pointer should
     *     be exclusively maintained by the `ChromosomeMutations` object during
     *     the call to this method. If needed, the returned pointer can be used
     *     for backup purpose.
     *
     * @return the shared pointer to the original data for backup
     */
    std::shared_ptr<Data> make_data_exclusive();

    /**
     * @brief Get an allele private reference in the chromosome
     *
     * This method is meant to be used to retrieve an allele that is going
     * to be changed. It guarantees that the data member is exclusive.
     *
     * @warning In order to avoid unnecessary and time-consuming copies of the
     *     chromosome mutation data, the chromosome mutation data pointer should
     *     be exclusively maintained by the `ChromosomeMutations` object during
     *     the call to this method. For backup purpose, use the second component
     *     of the returned pair.
     *
     * @param allele_id is the identifier of the allele to find
     * @return a pair whose first component is a pointer to the allele and second
     *     component is the shared pointer to the original data for backup
     * @throw std::out_of_range `allele_id` is not a valid allele identifier for
     *     the chromosome
     */
    std::pair<Allele*, std::shared_ptr<Data>>
    get_modifiable_allele(const AlleleId& allele_id);
public:
    /**
     * @brief The empty constructor
     */
    ChromosomeMutations();

    /**
     * @brief A constructor
     *
     * @param identifier is the chromosome identifier
     * @param size is the chromosome size
     * @param num_of_alleles is the initial number of alleles
     */
    ChromosomeMutations(const ChromosomeId& identifier, const ChromosomeMutations::Length& size, const size_t& num_of_alleles);

    /**
     * @brief A constructor
     *
     * @param chromosome_region is the chromosome region
     * @param num_of_alleles is the initial number of alleles
     */
    ChromosomeMutations(const GenomicRegion& chromosome_region, const size_t& num_of_alleles);

    /**
     * @brief Get chromosome identifier
     *
     * @return the chromosome identifier
     */
    inline const ChromosomeId& id() const
    {
        return _data->identifier;
    }

    /**
     * @brief Get chromosome size
     *
     * @return the chromosome size
     */
    inline const ChromosomeMutations::Length& size() const
    {
        return _data->length;
    }

    /**
     * @brief Get chromosome allelic size
     *
     * @return the chromosome allelic size
     */
    inline const ChromosomeMutations::Length& allelic_size() const
    {
        return _data->allelic_length;
    }

    /**
     * @brief Get the chromosome alleles
     *
     * @return a constant reference to chromosome alleles
     */
    inline const std::map<AlleleId, Allele>& get_alleles() const
    {
        return _data->alleles;
    }

    /**
     * @brief Get an allele in the chromosome
     *
     * @param allele_id is the identifier of the allele to find
     * @return a constant reference to the allele
     * @throw std::out_of_range `allele_id` is not a valid allele identifier for the
     *          chromosome
     */
    const Allele& get_allele(const AlleleId& allele_id) const;

    /**
     * @brief Get the CNAs occurred in the chromosome
     *
     * @return the list of the CNAs occurred in the chromosome
     */
    inline const std::list<std::shared_ptr<CNA>>& get_CNAs() const
    {
        return _data->CNAs;
    }

    /**
     * @brief Get the identifiers of the alleles containing a genomic position
     *
     * @param genomic_position is a position of the chromosome
     * @return a list of the identifiers of the alleles containing `genomic_position`
     */
    std::list<AlleleId> get_alleles_containing(const GenomicPosition& genomic_position) const;

    /**
     * @brief Get the identifiers of the alleles that completely include a region
     *
     * @param genomic_region is a region of the chromosome
     * @return a list of the identifiers of the alleles that completely include `genomic_region`
     */
    std::list<AlleleId> get_alleles_containing(const GenomicRegion& genomic_region) const;

    /**
     * @brief Get the alleles with context free for a mutation
     *
     * @param mutation is a mutation
     * @return a list of the identifiers of the allele in which the context of `mutation`
     *      is free
     */
    std::list<AlleleId> get_alleles_with_context_free_for(const SID& mutation) const;

    /**
     * @brief Check whether an allele contains a chromosomic region
     *
     * @param allele_id is the identifier of an allele
     * @param genomic_region is the region on which the check is performed
     * @return `true` if and only if the chromosome have an allele having `allele_id`
     *          as allele identifier and the whole `genomic_region` is contained
     *          in this allele
     * @throw std::domain_error `genomic_region` does not lay in this chromosome
     * @throw std::out_of_range the chromosome has not the allele `allele_id`
     */
    bool allele_contains(const AlleleId& allele_id, const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether a genomic position lays in a chromosome
     *
     * @param genomic_position is a genomic position
     * @return `true` if and only if `genomic_position` lays in the current chromosome
     */
    inline bool contains(const GenomicPosition& genomic_position) const
    {
        return (genomic_position.chr_id==id()       // same chromosome
                && genomic_position.position>0      // lays in [1, chromosome size]
                && genomic_position.position<=size());
    }

    /**
     * @brief Check whether a whole mutation lays in a chromosome
     *
     * @param mutation is a SID mutation
     * @return `true` if and only if `mutation` completely lays in the current
     *      chromosome
     */
    bool contains(const SID& mutation) const;

    /**
     * @brief Check whether a whole genomic region is mapped in a chromosome
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if `genomic_region` lays in the current chromosome
     */
    inline bool contains(const GenomicRegion& genomic_region) const
    {
        return (genomic_region.get_chromosome_id()==id()     // same chromosome
                && genomic_region.begin()>0   // it lays in [1, chromosome size]
                && genomic_region.end()<=size());
    }

    /**
     * @brief Check whether a mutation context is free
     *
     * @param mutation is a mutation
     * @return `true` if and only if the allele does not contains
     *      SID mutations in the one-base neighborhood of `mutation`
     * @throw std::domain_error `mutation` does not lay in the chromosome
     */
    bool has_context_free(const SID& mutation) const;

    /**
     * @brief Amplify a genomic region
     *
     * This method increases the number of copies in a genomic region by selecting
     * the same allele in each fragment of the region. If not all the fragments touching
     * the genomic region have the specified allele, the computation ends and the
     * method returns `false`.
     *
     * @param[in] genomic_region is the genomic region whose copy number must be increased
     * @param[in] allele_id is the identifier of the allele that must be amplified
     * @param[in,out] new_allele_id is the identifier of the new allele if amplification
     *      succeeded
     * @param[in] nature is the nature of the amplification
     * @return `true` if and only if the all the fragments touching `genomic_region`
     *      have `allele_id` among their allele ids and the region has been amplified
     * @throw std::domain_error `genomic_region` does not lay in this chromosome
     */
    bool amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                        AlleleId& new_allele_id,
                        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Amplify a genomic region
     *
     * This method increases the number of copies in a genomic region by selecting
     * the same allele in each fragment of the region. If not all the fragments touching
     * the genomic region have the specified allele, the computation ends and the
     * method returns `false`.
     *
     * @param[in] genomic_region is the genomic region whose copy number must be increased
     * @param[in] allele_id is the identifier of the allele that must be amplified
     * @param[in] nature is the nature of the amplification
     * @return `true` if and only if the all the fragments touching `genomic_region`
     *      have `allele_id` among their allele ids and the region has been amplified
     * @throw std::domain_error `genomic_region` does not lay in this chromosome
     */
    inline bool amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                               const Mutation::Nature& nature=Mutation::UNDEFINED)
    {
        AlleleId new_allele_id{RANDOM_ALLELE};

        return amplify_region(genomic_region, allele_id, new_allele_id, nature);
    }

    /**
     * @brief Remove a genomic region
     *
     * This method decreases the number of copies in a genomic region by selecting
     * the same allele in each fragment of the region. If not all the fragments touching
     * the genomic region have the specified allele, the computation ends and the
     * method returns `false`.
     *
     * @param genomic_region is the genomic region whose copy number must be decrease
     * @param allele_id is the identifier of the allele that must be deleted
     * @param nature is the nature of the deletion
     * @return `true` if and only if the all the fragments touching `genomic_region`
     *      have `allele_id` among their allele ids and the allele has been removed
     * @throw std::domain_error `genomic_region` does not lay in this chromosome
     */
    bool remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                       const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Apply a CNA
     *
     * This function applies a CNA. Whenever the `CNA::dest` member of the
     * parameter is set to `RANDOM_ALLELE`, the function change its value to
     * the destination allele identifier.
     *
     * @param[in,out] cna is the CNA to be applied
     * @return `true` if and only if the CNA has been successfully applied
     */
    inline bool apply(CNA& cna)
    {
        return apply_CNA(cna);
    }

    /**
     * @brief Apply a CNA
     *
     * @param[in] cna is the CNA to be applied
     * @return `true` if and only if the CNA has been successfully applied
     */
    inline bool apply(const CNA& cna)
    {
        auto cna_copy{cna};

        return apply_CNA(cna_copy);
    }

    /**
     * @brief Apply a SID mutation in a context free position
     *
     * This method tries to apply a SID mutation. If the chromosome contains another
     * SID mutation occurring in the mutational context of the specified SID, then
     * the application fails.
     *
     * @param mutation is the SID mutation to be applied in the chromosome
     * @param allele_id is the identifier of the allele in which the SID mutation
     *      must be placed
     * @return `true` if and only if the chromosome do not contain any other
     *      SID mutations in `mutation`'s mutational context
     * @throw std::domain_error `mutation` does not lay in the chromosome
     * @throw std::out_of_range the chromosome has not the allele
     *      `allele_id` or `mutation` does not lay in the allele
     */
    bool apply(const SID& mutation, const AlleleId& allele_id);

    /**
     * @brief Apply a SID mutation in a context free position
     *
     * This method tries to apply a SID mutation. If the chromosome contains another
     * SID mutation occurring in the mutational context of the specified SID, then
     * the application fails.
     *
     * @param mutation_spec is the SID mutation to be applied in the chromosome
     * @return `true` if and only if the chromosome do not contain any other
     *      SID mutations in `mutation`'s mutational context
     * @throw std::domain_error `mutation` does not lay in the chromosome
     * @throw std::out_of_range the chromosome has not the allele
     *      `allele_id` or `mutation` does not lay in the allele
     */
    inline bool apply(const MutationSpec<SID>& mutation_spec)
    {
        return apply(static_cast<const SID&>(mutation_spec), mutation_spec.allele_id);
    }

    /**
     * @brief Apply a list of mutations
     *
     * This method applies the mutations that lays in this chromosome.
     * If a mutations in the list does not affect this chromosome,
     * it is not applied. If all mutations in the list are applied
     * the method returns `true`. Otherwise, it returns `false`.
     *
     * @param mutation_list is the list of mutations to be applied
     * @return `true` if and only if all mutations in the list are applied
     *      the method returns `true`
     * @throw std::out_of_range any mutation in `mutation_list` does not lay
     *      in the corresponding genome allele
     */
    bool apply_contained(const MutationList& mutation_list);

    /**
     * @brief Copy genomic structure
     *
     * This method copies the genomic structure of the current objects.
     * It returns a `ChromosomeMutations` object that has the same alleles of
     * the current objects, but misses the original SID mutations.
     *
     * @return a `ChromosomeMutations` object that has the same alleles of
     *      the current objects, but misses the original SID mutations
     */
    ChromosomeMutations copy_structure() const;

    /**
     * @brief Duplicate genome alleles
     *
     * This method duplicate all the alleles in the genome
     */
    void duplicate_alleles();

    /**
     * @brief Get the allelic map of some chromosome fragments
     *
     * @param break_points are the fragment break points
     * @param min_num_of_alleles is the minimum number of alleles
     *    to report
     * @return A map that associates each break point to the
     *    allelic type of the genomic fragment beginning at the
     *    break point and ending at the following break point
     */
    std::map<ChrPosition, AllelicType>
    get_allelic_map(const std::set<ChrPosition>& break_points,
                    AlleleCounter min_num_of_alleles=0) const;

    /**
     * @brief Remove a SID mutation
     *
     * This method tries to remove a SID mutation from the chromosome.
     *
     * @param genomic_position is the position of the SID mutation to be
     *      removed
     * @return `true` if and only if a SID occurred at `genomic_position`
     * @throw std::domain_error `genomic_position` does not lay in the
     *      chromosome
     * @throw std::out_of_range the chromosome has not the allele
     *      `allele_id` or `genomic_position` does not lay in the allele
     */
    bool remove_mutation(const GenomicPosition& genomic_position);

    /**
     * @brief Check whether a SID mutation is included among the chromosome mutations
     *
     * @param mutation is the SID mutation whose inclusion among the chromosome
     *      mutations is tested
     * @return `true` if and only if `mutation` is contained among the chromosome
     *      mutations
     */
    bool includes(const SID& mutation) const;

    /**
     * @brief Save chromosome mutations in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "CLONES Genome Mutations", 1);

        archive & _data;
    }

    /**
     * @brief Load chromosome mutations from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded chromosome mutations
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static ChromosomeMutations load(ARCHIVE& archive)
    {
        ChromosomeMutations chr_mutations;

        ARCHIVE::read_header(archive, "CLONES Genome Mutations", 1);

        archive & chr_mutations._data;

        return chr_mutations;
    }
};

/**
 * @brief Test whether two `ChromosomeMutations` objects are the same
 *
 * @param lhs is the left-hand side of the relation
 * @param rhs is the right-hand side of the relation
 * @return `true` is and only if `lhs` and `rhs` represent the same
 *      chromosome, i.e., have the same alleles, contain the same
 *      mutations, and have the same length and chromosome id
 */
inline bool operator==(const CLONES::Mutations::ChromosomeMutations& lhs,
                       const CLONES::Mutations::ChromosomeMutations& rhs)
{
    return (lhs.id() == rhs.id()) && (lhs.size() == rhs.size())
            && (lhs.allelic_size() == rhs.allelic_size())
            && (lhs.get_alleles() == rhs.get_alleles());
}

/**
 * @brief Test whether two `ChromosomeMutations` objects differ
 *
 * @param lhs is the left-hand side of the relation
 * @param rhs is the right-hand side of the relation
 * @return `false` is and only if `lhs` and `rhs` represent the same
 *      chromosome, i.e., have the same alleles, contain the same
 *      mutations, and have the same length and chromosome id
 */
inline bool operator!=(const CLONES::Mutations::ChromosomeMutations& lhs,
                       const CLONES::Mutations::ChromosomeMutations& rhs)
{
    return (lhs.id() != rhs.id()) || (lhs.size() != rhs.size())
            || (lhs.allelic_size() != rhs.allelic_size())
            || (lhs.get_alleles() != rhs.get_alleles());
}

/**
 * @brief A class to represent the mutations of a genome
 */
class GenomeMutations
{
    std::map<ChromosomeId, ChromosomeMutations> chromosomes;     //!< the chromosomes

    /**
     * @brief Apply a CNA
     *
     * This function applies a CNA. Whenever the `CNA::dest` member of the
     * parameter is set to `RANDOM_ALLELE`, the function change its value to
     * the destination allele identifier.
     *
     * @param[in,out] cna is the CNA to be applied
     * @return `true` if and only if the CNA has been successfully applied
     */
    bool apply_CNA(CNA& cna);
public:
    using AllelicMap = std::map<ChromosomeId, std::map<ChrPosition, AllelicType>>;

    /**
     * @brief The genome length type
     */
    using Length = size_t;

    /**
     * @brief The empty constructor
     */
    GenomeMutations();

    /**
     * @brief A constructor
     *
     * @param chromosomes is the vector of genome chromosomes
     */
    explicit GenomeMutations(const std::vector<ChromosomeMutations>& chromosomes);

    /**
     * @brief A constructor
     *
     * This constructor builds the genome mutations model corresponding
     * to a set of context positions.
     *
     * @param chromosome_regions is a list of the chromosome regions
     * @param num_of_alleles is the initial number of alleles
     */
    GenomeMutations(const std::list<GenomicRegion>& chromosome_regions,
                    const size_t& num_of_alleles);

    /**
     * @brief A constructor
     *
     * This constructor builds the genome mutations model corresponding
     * to a set of context positions.
     *
     * @param chromosome_regions is a list of the chromosome regions
     * @param alleles_per_chromosome is the initial number of alleles per chromosome
     */
    GenomeMutations(const std::list<GenomicRegion>& chromosome_regions,
                    const std::map<ChromosomeId, size_t>& alleles_per_chromosome);

    /**
     * @brief A constructor
     *
     * This constructor builds the genome mutations model corresponding
     * to a set of context positions.
     *
     * @param chromosome_regions is a list of the chromosome regions
     * @param num_of_alleles is the initial number of alleles
     */
    inline GenomeMutations(std::list<GenomicRegion>&& chromosome_regions,
                           const size_t& num_of_alleles):
        GenomeMutations(chromosome_regions, num_of_alleles)
    {}

    /**
     * @brief A constructor
     *
     * This constructor builds the genome mutations model corresponding
     * to a set of context positions.
     *
     * @param chromosome_regions is a list of the chromosome regions
     * @param alleles_per_chromosome is the initial number of alleles per chromosome
     */
    inline GenomeMutations(std::list<GenomicRegion>&& chromosome_regions,
                           const std::map<ChromosomeId, size_t>& alleles_per_chromosome):
        GenomeMutations(chromosome_regions, alleles_per_chromosome)
    {}

    /**
     * @brief Get genome size
     *
     * @return the genome size
     */
    GenomeMutations::Length size() const;

    /**
     * @brief Get genome allelic size
     *
     * @return the genome allelic size
     */
    GenomeMutations::Length allelic_size() const;

    /**
     * @brief Get the mutations of genome chromosomes
     *
     * @return the mutations of genome chromosomes
     */
    inline const std::map<ChromosomeId, ChromosomeMutations>& get_chromosomes() const
    {
        return chromosomes;
    }

    /**
     * @brief Get the absolute chromosome positions
     *
     * @return a map from chromosome identifier to the position of the first base
     *      of the chromosome
     */
    std::map<ChromosomeId, size_t> get_absolute_chromosome_positions() const;

    /**
     * @brief Get the mutations of a genome chromosome
     *
     * @param chromosome_id is the identifier of the aimed chromosome mutations
     * @return the mutations of a genome chromosome
     */
    const ChromosomeMutations& get_chromosome(const ChromosomeId& chromosome_id) const;

    /**
     * @brief Get the mutations of a genome chromosome
     *
     * @param chromosome_id is the identifier of the aimed chromosome mutations
     * @return the mutations of a genome chromosome
     */
    ChromosomeMutations& get_chromosome(const ChromosomeId& chromosome_id);

    /**
     * @brief Get the identifiers of the alleles containing a genomic position
     *
     * @param genomic_position is a position of the genome
     * @return a list of the identifiers of the alleles containing `genomic_position`
     */
    std::list<AlleleId> get_alleles_containing(const GenomicPosition& genomic_position) const;

    /**
     * @brief Get the identifiers of the alleles that completely include a region
     *
     * @param genomic_region is a region of the chromosome
     * @return a list of the identifiers of the allele that completely include `genomic_region`
     */
    std::list<AlleleId> get_alleles_containing(const GenomicRegion& genomic_region) const;

    /**
     * @brief Get the alleles with context free for a mutation
     *
     * @param mutation is a mutation
     * @return a list of the identifiers of the allele in which the context of `mutation`
     *      is free
     */
    std::list<AlleleId> get_alleles_with_context_free_for(const SID& mutation) const;

    /**
     * @brief Check whether an allele contains a chromosomic region
     *
     * @param allele_id is the identifier of an allele
     * @param genomic_region is the region on which the check is performed
     * @return `true` if and only if the whole `genomic_region` is contained
     *          in the allele `allele_id` of the corresponding chromosome
     * @throw std::domain_error `genomic_region` does not lay in this chromosome
     * @throw std::out_of_range the chromosome has not the allele `allele_id`
     */
    bool allele_contains(const AlleleId& allele_id, const GenomicRegion& genomic_region) const;

    /**
     * @brief Amplify a genomic region
     *
     * This method increases the number of copies in a genomic region by selecting
     * the same allele in each fragment of the region. If not all the fragments touching
     * the genomic region have the specified allele, the computation ends and the
     * method returns `false`.
     *
     * @param[in] genomic_region is the genomic region whose copy number must be increased
     * @param[in] allele_id is the identifier of the allele that must be amplified
     * @param[in,out] new_allele_id is the identifier of the new allele if amplification
     *      succeeded
     * @param[in] nature is the nature of the amplification
     * @return `true` if and only if the all the fragments touching `genomic_region`
     *      have `allele_id` among their allele ids and the region has been amplified
     * @throw std::domain_error `genomic_region` does not lay in this chromosome
     */
    bool amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                        AlleleId& new_allele_id,
                        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Amplify a genomic region
     *
     * This method increases the number of copies in a genomic region by selecting
     * the same allele in each fragment of the region. If not all the fragments touching
     * the genomic region have the specified allele, the computation ends and the
     * method returns `false`.
     *
     * @param[in] genomic_region is the genomic region whose copy number must be increased
     * @param[in] allele_id is the identifier of the allele that must be amplified
     * @param[in] nature is the nature of the amplification
     * @return `true` if and only if the all the fragments touching `genomic_region`
     *      have `allele_id` among their allele ids and the region has been amplified
     * @throw std::domain_error `genomic_region` does not lay in this chromosome
     */
    inline bool amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                               const Mutation::Nature& nature=Mutation::UNDEFINED)
    {
        AlleleId new_allele_id{RANDOM_ALLELE};

        return amplify_region(genomic_region, allele_id, new_allele_id, nature);
    }

    /**
     * @brief Remove a genomic region
     *
     * This method decreases the number of copies in a genomic region by selecting
     * the same allele in each fragment of the region. If not all the fragments touching
     * the genomic region have the specified allele, the computation ends and the
     * method returns `false`
     *
     * @param genomic_region is the genomic region whose copy number must be decrease
     * @param allele_id is the identifier of the allele from which the region must be removed
     * @param nature is the nature of the deletion
     * @return `true` if and only if deletion has been successful
     * @throw std::domain_error no allele can be removed
     */
    bool remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                       const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Apply a CNA
     *
     * This function applies a CNA. Whenever the `CNA::dest` member of the
     * parameter is set to `RANDOM_ALLELE`, the function change its value to
     * the destination allele identifier.
     *
     * @param[in,out] cna is the CNA to be applied
     * @return `true` if and only if the CNA has been successfully applied
     */
    inline bool apply(CNA& cna)
    {
        return apply_CNA(cna);
    }

    /**
     * @brief Apply a CNA
     *
     * @param[in] cna is the CNA to be applied
     * @return `true` if and only if the CNA has been successfully applied
     */
    inline bool apply(const CNA& cna)
    {
        auto cna_copy{cna};

        return apply_CNA(cna_copy);
    }

    /**
     * @brief Apply a SID mutation in a context free position
     *
     * This method tries to apply a SID mutation. If the genome contains another
     * SID mutation occurring in the mutational context of the specified SID, then
     * the application fails.
     *
     * @param mutation is the SID mutation to be applied to the genome
     * @param allele_id is the identifier of the allele in which the SID mutation
     *      must be placed
     * @return `true` if and only if the genome do not contain any other
     *      SID mutations in `mutation`'s mutational context
     * @throw std::domain_error `mutation` does not lay to the genome
     * @throw std::out_of_range the genome has not the allele
     *      `allele_id` or `mutation` does not lay in the allele
     */
    bool apply(const SID& mutation, const AlleleId& allele_id);

    /**
     * @brief Apply a SID mutation in a context free position
     *
     * This method tries to apply a SID mutation. If the genome contains another
     * SID mutation occurring in the mutational context of the specified SID, then
     * the application fails.
     *
     * @param mutation_spec is the SID mutation to be applied to the genome
     * @return `true` if and only if the genome do not contain any other
     *      SID mutations in `mutation`'s mutational context
     * @throw std::domain_error `mutation` does not lay into the genome
     * @throw std::out_of_range `mutation_spec` does not lay in the corresponding
     *      genome allele
     */
    bool apply(const MutationSpec<SID>& mutation_spec);

    /**
     * @brief Apply a list of mutations
     *
     * This method tries to apply a list of mutations. If all mutations in
     * the list are applied the method returns `true`. Otherwise, it returns
     * `false`.
     *
     * @param mutation_list is the list of mutations to be applied
     * @return `true` if and only if all mutations in the list are applied
     *      the method returns `true`
     * @throw std::domain_error any mutation in `mutation_list` does not lay
     *      into the genome
     * @throw std::out_of_range any mutation in `mutation_list` does not lay
     *      in the corresponding genome allele
     */
    bool apply(const MutationList& mutation_list);

    /**
     * @brief Apply a list of mutations
     *
     * This method tries to apply the mutations in a list that lays in this
     * genome. If all mutations in the list are applied the method returns
     * `true`. Otherwise, it returns `false`.
     *
     * @param mutation_list is the list of mutations to be applied
     * @return `true` if and only if all mutations in the list are applied
     *      the method returns `true`
     * @throw std::domain_error any mutation in `mutation_list` does not lay
     *      into the genome
     * @throw std::out_of_range any mutation in `mutation_list` does not lay
     *      in the corresponding genome allele
     */
    inline bool apply_contained(const MutationList& mutation_list)
    {
        return apply(mutation_list);
    }

    /**
     * @brief Remove a SID mutation
     *
     * This method tries to remove a SID mutation from the chromosome.
     *
     * @param genomic_position is the position of the SID to be removed
     * @return `true` if and only if a SID occurred at `genomic_position`
     */
    bool remove_mutation(const GenomicPosition& genomic_position);

    /**
     * @brief Check whether a mutation context is free
     *
     * @param mutation is a mutation
     * @return `true` if and only if the allele does not contains
     *      SID mutations in the one-base neighborhood of `mutation`
     * @throw std::domain_error `mutation` does not lay in the genome
     */
    bool has_context_free(const SID& mutation) const;

    /**
     * @brief Check whether a SID mutation is included among the genome mutations
     *
     * @param mutation is the SID mutation whose inclusion among the genome
     *      mutations is tested
     * @return `true` if and only if `mutation` is contained among the genome
     *      mutations
     */
    bool includes(const SID& mutation) const;

    /**
     * @brief Copy genomic structure
     *
     * This method copies the genomic structure of the current objects.
     * It returns a `GenomeMutations` object that has the same chromosomes and
     * alleles of the current objects, but misses the original SNVs and indels.
     *
     * @return a `GenomeMutations` object that has the same chromosomes and
     *      alleles of the current objects, but misses the original SNVs and
     *      indels
     */
    GenomeMutations copy_structure() const;

    /**
     * @brief Duplicate genome alleles
     *
     * This method duplicate all the alleles in the genome
     */
    void duplicate_alleles();

    /**
     * @brief Get the allelic map of the fragments
     *
     * @param min_num_of_alleles is the minimum number of alleles
     *    to report
     * @return A map that, for each chromosome identifier,
     *    associates each of the genome break points to the allelic
     *    type of the genomic fragment beginning at the break point
     *    and ending at the following break point
     */
    inline AllelicMap
    get_allelic_map(const AlleleCounter& min_num_of_alleles=0) const
    {
        const auto b_points = get_CNA_break_points();

        return get_allelic_map(b_points, min_num_of_alleles);
    }

    /**
     * @brief Get the allelic map of the genomic fragments
     *
     * @param break_points is a map from the set of chromosome
     *    identifiers to the fragment break points in the
     *    corresponding chromosome
     * @param min_allelic_size is the minimum number of alleles
     *    to report
     * @return A map that, for each chromosome identifier,
     *    associates each break point to the allelic type of the
     *    genomic fragment beginning at the break point and ending
     *    at the following break point
     */
    AllelicMap
    get_allelic_map(const std::map<ChromosomeId, std::set<ChrPosition>>& break_points,
                    const AlleleCounter& min_allelic_size=0) const;

    /**
     * @brief Get the CNA break points
     *
     * @return A map that associates a chromosome identifier
     *   to the set of CNA break points on the corresponding
     *   chromosome
     */
    std::map<ChromosomeId, std::set<ChrPosition>>
    get_CNA_break_points() const;

    /**
     * @brief Save genome mutations in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "CLONES GenomeMutations", 0);

        archive & chromosomes;
    }

    /**
     * @brief Load genome mutations from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded genome mutations
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static GenomeMutations load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "CLONES GenomeMutations", 0);

        GenomeMutations g_mutations;

        archive & g_mutations.chromosomes;

        return g_mutations;
    }
};

/**
 * @brief A class to represent the mutations of a specific cell
 */
struct CellGenomeMutations : public Mutants::Cell, public GenomeMutations
{
    /**
     * @brief The genome length type
     */
    using Length = size_t;

    /**
     * @brief The empty constructor
     */
    CellGenomeMutations();

    /**
     * @brief A constructor for wild-type cells
     *
     * @param germline_mutations are the germline mutations
     */
    explicit CellGenomeMutations(const GenomeMutations& germline_mutations);

    /**
     * @brief A constructor
     *
     * @param cell is the cell
     * @param genome_mutations are the genome mutations
     */
    explicit CellGenomeMutations(const Mutants::Cell& cell, const GenomeMutations& genome_mutations);

    /**
     * @brief Save cell genome mutations in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const Mutants::Cell&>(*this)
                & static_cast<const GenomeMutations&>(*this);
    }

    /**
     * @brief Load cell genome mutations from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded cell genome mutations
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static CellGenomeMutations load(ARCHIVE& archive)
    {
        CellGenomeMutations cg_mutations;

        archive & static_cast<Mutants::Cell&>(cg_mutations)
                & static_cast<GenomeMutations&>(cg_mutations);

        return cg_mutations;
    }
};

/**
 * @brief Test whether two `GenomeMutations` objects are the same
 *
 * @param lhs is the left-hand side of the relation
 * @param rhs is the right-hand side of the relation
 * @return `true` is and only if `lhs` and `rhs` represent the same
 *      genome, i.e., have the same chromosomes and contain the same
 *      mutations
 */
inline bool operator==(const CLONES::Mutations::GenomeMutations& lhs,
                       const CLONES::Mutations::GenomeMutations& rhs)
{
    return (lhs.get_chromosomes() == rhs.get_chromosomes());
}

/**
 * @brief Test whether two `GenomeMutations` objects differ
 *
 * @param lhs is the left-hand side of the relation
 * @param rhs is the right-hand side of the relation
 * @return `false` is and only if `lhs` and `rhs` represent the same
 *      genome, i.e., have the same chromosomes and contain the same
 *      mutations
 */
inline bool operator!=(const CLONES::Mutations::GenomeMutations& lhs,
                       const CLONES::Mutations::GenomeMutations& rhs)
{
    return (lhs.get_chromosomes() != rhs.get_chromosomes());
}

}   // Mutations

}   // CLONES


namespace std
{

/**
 * @brief Write chromosome mutations data in a stream
 *
 * @param os is the output stream
 * @param chromosome_mutations is a chromosome mutations
 * @return a reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::ChromosomeMutations& chromosome_mutations);

/**
 * @brief Write genome mutations data in a stream
 *
 * @param os is the output stream
 * @param genome_mutations is a genome mutations
 * @return a reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::GenomeMutations& genome_mutations);

}   // std

#endif // __CLONES_GENOME_MUTATIONS__
