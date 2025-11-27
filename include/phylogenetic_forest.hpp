/**
 * @file phylogenetic_forest.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for phylogenetic forests
 * @version 1.16
 * @date 2025-11-28
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

#ifndef __RACES_PHYLOGENETIC_FOREST__
#define __RACES_PHYLOGENETIC_FOREST__

#include <map>
#include <memory>

#include "mutation_list.hpp"
#include "mutational_properties.hpp"
#include "descendant_forest.hpp"

#include "genome_mutations.hpp"
#include "mutation_spec.hpp"

#include "label_tour.hpp"

namespace RACES
{

namespace Mutations
{

/**
 * @brief A class representing descendant forests
 */
class PhylogeneticForest : public Mutants::DescendantForest
{
public:

    /**
     * @brief Some statistics about samples
     *
     * This structure collects some statistics about one sample.
     * It maintains the quantity of DNA and the number of cells
     * in the sample.
     */
    struct SampleStatistics
    {
        size_t total_allelic_size;  //!< The quantity of DNA in the sample
        size_t number_of_cells;     //!< The number of cells in the sample

        /**
         * @brief The empty constructor
         */
        SampleStatistics();

        /**
         * @brief Save sample statistics in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param[out] archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & total_allelic_size
                    & number_of_cells;
        }

        /**
         * @brief Load sample statistics from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param[in,out] archive is the input archive
         * @return the load sample statistics
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
        inline static SampleStatistics load(ARCHIVE& archive)
        {
            SampleStatistics statistics;

            archive & statistics.total_allelic_size
                    & statistics.number_of_cells;

            return statistics;
        }
    };
private:
    using CellIdSet = std::set<Mutants::CellId>;

    using MutationsPerCell = std::map<Mutants::CellId, MutationList>;

    using TissueSampleId = Mutants::Evolutions::TissueSampleId;

    MutationsPerCell pre_neoplastic_mutations;   //!< The pre-neoplastic mutations per forest root
    MutationsPerCell arising_mutations;          //!< The non-pre-neoplastic mutations arising in the forest cells
    std::map<SID, CellIdSet> SID_first_cells;      //!< A map associating each SID to the first cells in which it occurred
    std::map<CNA, CellIdSet> CNA_first_cells;      //!< A map associating each CNA to the first cells in which it occurred

    std::map<TissueSampleId, SampleStatistics>  sample_statistics;  //!< The sample statistics

    std::shared_ptr<GenomeMutations> germline_mutations; //!< The germline mutations

    MutationalProperties mutational_properties; //!< The mutational properties used to build the forest
public:
    using AllelicCount = std::map<ChromosomeId, std::map<ChrPosition, std::map<AllelicType, size_t>>>;

    /**
     * @brief A constant node of the forest
     */
    class const_node : public Mutants::DescendantForest::_const_node<PhylogeneticForest>
    {
    public:
        /**
         * @brief A constructor for a constant node
         *
         * @param[in] forest is the forest of the node
         * @param[in] cell_id is the cell id of a cell in the forest
         */
        const_node(const PhylogeneticForest* forest, const Mutants::CellId cell_id);

        /**
         * @brief Get the genome mutations of the cell represented by the node
         *
         * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
         *      mutations
         * @param[in] with_germinal is a Boolean flag to add/avoid germinal mutations
         * @return a constant reference to the genome mutations of the cell
         *      represented by the node
         */
        inline CellGenomeMutations cell_mutations(const bool& with_pre_neoplastic=true,
                                                  const bool& with_germinal=false) const
        {
            return forest->get_cell_mutations(get_id(), with_pre_neoplastic, with_germinal);
        }

        /**
         * @brief Get the parent node
         *
         * @return the parent node of this node
         */
        inline const_node parent() const
        {
            return Mutants::DescendantForest::_const_node<PhylogeneticForest>::parent<const_node>();
        }

        /**
         * @brief Get the node direct descendants
         *
         * @return a vector containing the direct descendants of
         *         this node
         */
        inline std::vector<const_node> children() const
        {
            return Mutants::DescendantForest::_const_node<PhylogeneticForest>::children<const_node>();
        }

        /**
         * @brief Get the mutations arising in the cell represented by this node
         *
         * @return a constant reference to the list of mutations arising in the
         *      cell represented by this node
         */
        inline const MutationList& arising_mutations() const
        {
            return forest->arising_mutations.at(get_id());
        }

        /**
         * @brief Get the pre-neoplastic mutations of the cell represented by this node
         *
         * This method returns the pre-neoplastic mutations arising in the cell represented
         * by this node if this node is a root of the considered forest. Otherwise, it throws
         * an `std::runtime_error` exception.
         *
         * @return a constant reference to the list of pre-neoplastic mutations arising
         *      in the cell represented by this node
         * @throw `std::runtime_error` if the current node is not a forest root
         */
        const MutationList& pre_neoplastic_mutations() const;

        friend class PhylogeneticForest;
    };

    class node : public Mutants::DescendantForest::_node<PhylogeneticForest>
    {
        /**
         * @brief Add the current cell node as first cell occurrence for a SID
         *
         * @param[in] sid is the SID appearing for the first time in the current
         *      node cell
         */
        inline void add_first_cell(const MutationSpec<SID>& sid)
        {
            get_forest().SID_first_cells[sid].insert(cell_id);
        }

        /**
         * @brief Add the current cell node as first cell occurrence for a CNA
         *
         * @param[in] cna is the CNA appearing for the first time in the current
         *      node cell
         */
        inline void add_first_cell(const CNA& cna)
        {
            get_forest().CNA_first_cells[cna].insert(cell_id);
        }
    public:
        /**
         * @brief A constructor for a constant node
         *
         * @param[in,out] forest is the forest of the node
         * @param[in] cell_id is the cell id of a cell in the forest
         */
        node(PhylogeneticForest* forest, const Mutants::CellId cell_id);

        /**
         * @brief Get the parent node
         *
         * @return the parent node of this node
         */
        inline node parent()
        {
            return Mutants::DescendantForest::_node<PhylogeneticForest>::parent<node>();
        }

        /**
         * @brief Get the node direct descendants
         *
         * @return a vector containing the direct descendants of
         *         this node
         */
        inline std::vector<node> children()
        {
            return Mutants::DescendantForest::_node<PhylogeneticForest>::children<node>();
        }

        /**
         * @brief Add a newly introduced mutation
         *
         * @tparam MUTATION_TYPE is the type of mutation
         * @param[in] mutation is a mutation that has appeared for the first time in the cell
         *      corresponding to the current node
         */
        template<typename MUTATION_TYPE,
                 std::enable_if_t<std::is_same_v<MUTATION_TYPE,MutationSpec<SID>>
                                                 || std::is_same_v<MUTATION_TYPE,CNA>, bool> = true>
        void add_new_mutation(const MUTATION_TYPE& mutation)
        {
            if (mutation.nature == Mutation::PRE_NEOPLASTIC) {
                get_forest().pre_neoplastic_mutations[cell_id].insert(mutation);
            } else {
                get_forest().arising_mutations[cell_id].insert(mutation);
            }
            add_first_cell(mutation);
        }

        /**
         * @brief Add the whole genome doubling to the node mutations
         */
        inline void add_whole_genome_doubling()
        {
            get_forest().arising_mutations[cell_id].insert_WGD();
        }

        /**
         * @brief Get the mutations arising in the cell represented by this node
         *
         * @return a constant reference to the list of mutations arising in the
         *      cell represented by this node
         */
        inline const MutationList& arising_mutations() const
        {
            return forest->arising_mutations.at(get_id());
        }

        /**
         * @brief Get the mutations arising in the cell represented by this node
         *
         * @return a reference to the list of mutations arising in the
         *      cell represented by this node
         */
        inline MutationList& arising_mutations()
        {
            return forest->arising_mutations[get_id()];
        }

        friend class PhylogeneticForest;
    };

public:
    /**
     * @brief The empty constructor
     */
    PhylogeneticForest();

    /**
     * @brief Get a constant node with the given id
     *
     * @param[in] cell_id is the id of the aimed cell node
     * @return the corresponding constant cell node
     */
    const_node get_node(const Mutants::CellId& cell_id) const;

    /**
     * @brief Get a node with the given id
     *
     * @param[in] cell_id is the id of the aimed cell node
     * @return the corresponding cell node
     */
    node get_node(const Mutants::CellId& cell_id);

    /**
     * @brief Get the forest roots
     *
     * @return std::vector<const_node>
     */
    std::vector<const_node> get_roots() const;

    /**
     * @brief Get the forest roots
     *
     * @return std::vector<const_node>
     */
    std::vector<node> get_roots();

    /**
     * @brief Get the forest for a subset of the tissue samples
     *
     * @param[in] sample_names are the names of the samples to be considered
     * @return the descendant forest for the tissue samples whose name is
     *          in `sample_names`
     */
    PhylogeneticForest get_subforest_for(const std::vector<std::string>& sample_names) const;

    /**
     * @brief Get the mutations arising in the forest cells
     *
     * @return a constant reference to the mutations arising in the
     *      forest cells
     */
    inline const std::map<Mutants::CellId, MutationList>&
    get_arising_mutations() const
    {
        return arising_mutations;
    }

    /**
     * @brief Get the genome mutations of a cell represented in the forest
     *
     * @param[in] cell_id is a cell identifier of a cell represented in the forest
     * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
     *      mutations
     * @param[in] with_germinal is a Boolean flag to add/avoid germinal mutations
     * @return the genome mutations of the cell with identifier `cell_id`
     */
    CellGenomeMutations get_cell_mutations(const Mutants::CellId& cell_id,
                                           const bool& with_pre_neoplastic=true,
                                           const bool& with_germinal=false) const;

    /**
     * @brief Get map associating each SID to the first cell in which it occurs
     *
     * @return a constant reference to a map associating each SID in the phylogenetic
     *         forest to the identifier of the first cell in which the SID occurred
     */
    inline const std::map<SID, std::set<Mutants::CellId>>& get_mutation_first_cells() const
    {
        return SID_first_cells;
    }

    /**
     * @brief Get map associating each CNA to the first cell in which it occurs
     *
     * @return a constant reference to a map associating each CNA in the phylogenetic
     *         forest to the identifier of the first cell in which the CNA occurred
     */
    inline const std::map<CNA, std::set<Mutants::CellId>>& get_CNA_first_cells() const
    {
        return CNA_first_cells;
    }

    /**
     * @brief Get the germline mutations
     *
     * @return A constant reference to germline mutations
     */
    inline const GenomeMutations& get_germline_mutations() const
    {
        return *germline_mutations;
    }

    /**
     * @brief Get the forest mutational properties
     *
     * @return A constant reference to the forest mutational properties
     */
    inline const MutationalProperties& get_mutational_properties() const
    {
        return mutational_properties;
    }

    /**
     * @brief Get the pre-neoplastic mutations
     *
     * @return A constant reference to the map of pre-neoplastic mutations grouped by forest
     *      root cell identifier
     */
    inline const std::map<Mutants::CellId, MutationList>& get_pre_neoplastic_mutations() const
    {
        return pre_neoplastic_mutations;
    }

    /**
     * @brief Build the wild type genomes
     *
     * This method build the wild-type genome of the embryo cell and of every cell
     * represented by a root in the forest.
     *
     * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
     *      mutations
     * @param[in] with_germinal is a Boolean flag to include the germinal mutations
     * @return a map that associates to each forest root and to the embryo cell
     *      the corresponding wild type genome.
     */
    std::map<Mutants::CellId, CellGenomeMutations>
    get_wild_type_genomes(const bool with_pre_neoplastic=true,
                          const bool with_germinal=true) const;

    /**
     * @brief Get the CNA break points
     *
     * @return The list of the CNA break points grouped by chromosome
     *   identifier
     */
    std::map<ChromosomeId, std::set<ChrPosition>> get_CNA_break_points() const;

    /**
     * @brief Get the allelic count of all the leaves
     *å
     * @param[in] min_allelic_size is the minimum number of alleles to report
     * @return A map that, for each chromosome and for each break point, reports
     *   the number of leaves per allelic type
     */
    AllelicCount get_allelic_count(const AlleleCounter& min_allelic_size=0) const;

    /**
     * @brief Get the allelic count of a set of cells
     *
     * @param[in] cell_ids is a list of cell identifiers corresponding to leaves in
     *   the phylogenetic forest
     * @param[in] min_allelic_size is the minimum number of alleles to report
     * @return A map that, for each chromosome and for each break point, reports
     *   the number of cells among those in corresponding to `cell_ids` per
     *   allelic type
     */
    AllelicCount get_allelic_count(const std::list<Mutants::CellId>& cell_ids,
                                   const AlleleCounter& min_allelic_size=0) const;

    /**
     * @brief Get the allelic count of a sample
     *
     * @param[in] sample_name is the sample name
     * @param[in] min_allelic_size is the minimum number of alleles to report
     * @return A map that, for each chromosome and for each break point, reports
     *   the number of cells in the sample per allelic type
     */
    AllelicCount get_allelic_count(const std::string& sample_name,
                                   const size_t& min_allelic_size=0) const;

    /**
     * @brief Get the sample statistics
     *
     * @return a constant reference to the map associating the sample identifiers to
     *      the corresponding statistics
     */
    inline const std::map<TissueSampleId, SampleStatistics>& get_sample_statistics() const
    {
        return sample_statistics;
    }

    /**
     * @brief Clear the forest
     */
    void clear();

    /**
     * @brief Save a phylogenetic forest in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param[out] archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "RACES Phylogenetic Forest", 4);

        archive & static_cast<const Mutants::DescendantForest&>(*this)
                & pre_neoplastic_mutations
                & arising_mutations
                & SID_first_cells
                & CNA_first_cells
                & sample_statistics
                & *germline_mutations
                & mutational_properties;
    }

    /**
     * @brief Load a phylogenetic forest from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param[in,out] archive is the input archive
     * @return the load phylogenetic forest
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static PhylogeneticForest load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "RACES Phylogenetic Forest", 4);

        PhylogeneticForest forest;

        forest.germline_mutations = std::make_shared<GenomeMutations>();

        archive & static_cast<Mutants::DescendantForest&>(forest)
                & forest.pre_neoplastic_mutations
                & forest.arising_mutations
                & forest.SID_first_cells
                & forest.CNA_first_cells
                & forest.sample_statistics
                & *(forest.germline_mutations)
                & forest.mutational_properties;

        return forest;
    }

    template<typename GENOME_WIDE_POSITION, typename RANDOM_GENERATOR>
    friend class MutationEngine;
};

/**
 * @brief The mutation labelling functor
 *
 * @tparam MUTATION_CONTAINER is the type of container for the mutations. It can
 *      be either `GenomeMutations` or `ChromosomeMutations`
 */
template<typename MUTATION_CONTAINER,
            std::enable_if_t<std::is_same_v<MUTATION_CONTAINER, ChromosomeMutations>
                            || std::is_same_v<MUTATION_CONTAINER, GenomeMutations>, bool> = true>
class MutationLabellingFunctor
{
    bool with_pre_neoplastic;   //!< A Boolean flag to add/avoid pre-neoplastic mutations

public:
    /**
     * @brief The label type
     */
    using label_type = MUTATION_CONTAINER;

    /**
     * @brief The constructor
     *
     * @param with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic mutations
     */
    explicit MutationLabellingFunctor(const bool with_pre_neoplastic=true):
        with_pre_neoplastic(with_pre_neoplastic)
    {}

    /**
     * @brief The functor application
     *
     * @param parent_mutations is the mutation of the parent of current node
     * @param node is the current node
     * @return Either the `GenomeMutations` or `ChromosomeMutations` labelling the
     *      current node
     */
    MUTATION_CONTAINER operator()(const MUTATION_CONTAINER& parent_mutations,
                                    const PhylogeneticForest::const_node& node) const
    {
        MUTATION_CONTAINER mutations{parent_mutations};

        if (node.is_root() && with_pre_neoplastic) {
            mutations.apply_contained(node.pre_neoplastic_mutations());
        }

        mutations.apply_contained(node.arising_mutations());

        return mutations;
    }
};

/**
 * @brief The tour of the mutations of a phylogenetic forest
 *
 * @tparam MUTATION_TYPE is the type of container for the mutations. It can
 *      be either `GenomeMutations` or `ChromosomeMutations`
 */
template<typename MUTATION_TYPE>
using MutationTour = LabelTour<PhylogeneticForest, MutationLabellingFunctor<MUTATION_TYPE>>;

/**
 * @brief Get a tour over the genome mutations of the forest leaves
 *
 * @param[in] forest is a constant reference to the forest whose labels
 *          will be visited
 * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
 *      mutations
 * @param[in] with_germinal is a Boolean flag to add/avoid germinal mutations
 * @return a tour over the genome mutations of the forest leaves
 */
MutationTour<GenomeMutations>
get_leaf_mutation_tour(const PhylogeneticForest& forest,
                       const bool with_pre_neoplastic=true,
                       const bool with_germinal=false);

/**
 * @brief Get a tour over the chromosome mutations of the forest leaves
 *
 * @param[in] forest is a constant reference to the forest whose labels
 *          will be visited
 * @param[in] chromosome_id is the identifier of the chromosome of interest
 * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
 *      mutations
 * @param[in] with_germinal is a Boolean flag to add/avoid germinal mutations
 * @return a tour over the genome mutations of the forest leaves
 */
MutationTour<ChromosomeMutations>
get_leaf_mutation_tour(const PhylogeneticForest& forest,
                       const ChromosomeId& chromosome_id,
                       const bool with_pre_neoplastic=true,
                       const bool with_germinal=false);

}   // Mutations

}   // RACES

#endif // __RACES_PHYLOGENETIC_TREE__
