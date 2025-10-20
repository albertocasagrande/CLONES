/**
 * @file phylogenetic_forest.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for phylogenetic forests
 * @version 1.14
 * @date 2025-10-20
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

namespace RACES
{

namespace Mutations
{

/**
 * @brief The template dependent data for the template `PhylogeneticForest::MutationTour`
 *
 * @tparam MUTATION_CONTAINER is the structure that contains the mutations,
 *      i.e., `GenomeMutations` or `ChromosomeMutations`
 */
template<typename MUTATION_CONTAINER>
struct ExtraTourData
{};

/**
 * @brief The `PhylogeneticForest::MutationTour` data for `ChromosomeMutations`
 *
 * This is a specialization of `ExtraTourData` for `ChromosomeMutations`.
 */
template<>
struct ExtraTourData<ChromosomeMutations>
{
    ChromosomeId chr_id;    //!< The chromosome identificator
};

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
     * @brief Tours of the mutations of the forest nodes
     *
     * This class implements tours of the genome or chromosome mutations --
     * depending on the template parameter -- mutations of the forest nodes. It
     * implements a constant iterator that allows the online generation of the
     * genome or chromosome mutations of each node or leaf in the forest.
     *
     * @tparam MUTATION_CONTAINER is the structure that contains the mutations,
     *      i.e., `GenomeMutations` or `ChromosomeMutations`
     */
    template<typename MUTATION_CONTAINER,
             std::enable_if_t<std::is_same_v<MUTATION_CONTAINER, ChromosomeMutations>
                              || std::is_same_v<MUTATION_CONTAINER, GenomeMutations>, bool> = true>
    class MutationTour : ExtraTourData<MUTATION_CONTAINER>
    {
        PhylogeneticForest const* forest;   //!< A pointer to the forest
        bool only_leaves;       //!< A Boolean flag to enable/disable internal node visit
        bool with_pre_neoplastic;   //!< A Boolean flag to add/avoid pre-neoplastic mutations
        bool with_germinal;     //!< A Boolean flag to add/avoid germline mutations

        /**
         * @brief Construct a new mutation tour object
         *
         * @param[in] forest is a constant pointer to the forest
         * @param[in] only_leaves is a Boolean flag to enable/disable internal node visit
         * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
         *      mutations
         * @param[in] with_germinal is Boolean flag to add/avoid germline mutations in the
         *      produced mutations structure
         */
        MutationTour(const PhylogeneticForest* forest, const bool only_leaves,
                     const bool with_pre_neoplastic, const bool with_germinal):
            forest{forest}, only_leaves{only_leaves},
            with_pre_neoplastic{with_pre_neoplastic},
            with_germinal{with_germinal}
        {}
    public:
        /**
         * @brief A constant iterator for the mutation tour
         *
         * This class implements a constant iterator that visits all nodes or
         * leaves in the tour forest and generates the corresponding mutations
         * structure.
         * The asymptotic complexity of the successor operator is linear in the
         * number of forest nodes in the worst case. However, this asymptotic
         * complexity also holds for the full tour.
         * In the worst case, the memory required by each iterator is
         * logarithmic in the number of forest nodes.
         */
        class const_iterator
        {
        public:
            /**
             * @brief The dereferenced value type
             */
            using value_type = std::pair<Mutants::CellId, MUTATION_CONTAINER>;
        private:
            PhylogeneticForest const* forest;   //!< A pointer to the forest
            bool only_leaves;   //!< A Boolean flag to enable/disable internal node visit
            bool tour_end;      //!< A Boolean flag to mark the end of the tour

            std::stack<value_type> iterator_stack; //!< The recursion stack
            value_type node_mutations; //!< The genome mutation of the current node

            /**
             * @brief A constructor
             *
             * @param[in] forest is a constant pointer to the forest
             * @param[in] only_leaves is a Boolean flag to enable/disable internal node visit
             * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
             *      mutations
             * @param[in] with_germinal is Boolean flag to add/avoid germline mutations
             */
            const_iterator(const PhylogeneticForest* forest,
                           const MUTATION_CONTAINER& template_mutations,
                           const bool& only_leaves, const bool& with_pre_neoplastic):
                forest{forest}, only_leaves{only_leaves}, tour_end{false}
            {
                if (forest != nullptr) {
                    auto forest_roots = forest->get_roots();
                    for (auto root_it = forest_roots.rbegin();
                            root_it != forest_roots.rend(); ++root_it) {
                        MUTATION_CONTAINER mutations = template_mutations;
                        mutations.apply_contained(root_it->arising_mutations());

                        if (with_pre_neoplastic) {
                            mutations.apply_contained(root_it->pre_neoplastic_mutations());
                        }

                        iterator_stack.emplace(root_it->get_id(), std::move(mutations));
                    }

                    std::swap(node_mutations, iterator_stack.top());

                    iterator_stack.pop();

                    if (only_leaves) {
                        this->operator++();
                    }
                } else {
                    tour_end = true;
                }
            }

            const_iterator(const PhylogeneticForest* forest, const bool& only_leaves):
                forest{forest}, only_leaves{only_leaves}, tour_end{true}
            {}
        public:
            /**
             * @brief The empty constructor
             */
            const_iterator():
                const_iterator{nullptr, false}
            {}

            /**
             * @brief The successor operator
             *
             * This method moves the iterator to the next node of the tour
             * and, at the same time, builds the mutations structure of the
             * reached node.
             * The asymptotic complexity of this method is linear in the
             * number of the forest nodes. However, the full tour has the
             * same complexity.
             *
             * @return a reference to the updated object
             */
            const_iterator& operator++()
            {
                if (forest == nullptr) {
                    return *this;
                }

                const_node node{forest, node_mutations.first};

                if (node.is_leaf()) {
                    if (iterator_stack.empty()) {
                        tour_end = true;

                        return *this;
                    }

                    // take a new cell genome mutations object from the stack
                    std::swap(node_mutations, iterator_stack.top());
                    iterator_stack.pop();

                    if (!only_leaves) {
                        return *this;
                    }

                    // update the node
                    node = const_node{forest, node_mutations.first};
                }

                bool next_node_found = node.is_leaf();
                while (!next_node_found) {
                    auto children = node.children();

                    // place all children's cell genome mutations, but the first one, into the stack
                    for (auto child_it = children.rbegin();
                            child_it != children.rend()-1; ++child_it) {

                        MUTATION_CONTAINER child_mutations{node_mutations.second};
                        child_mutations.apply_contained(child_it->arising_mutations());

                        iterator_stack.emplace(child_it->get_id(), std::move(child_mutations));
                    }

                    // apply the first children mutations to the current cell genome mutations
                    std::swap(node, children.front());
                    node_mutations.second.apply_contained(node.arising_mutations());

                    next_node_found = node.is_leaf() || !only_leaves;
                }

                // update the node cell id
                node_mutations.first = node.get_id();

                return *this;
            }

            /**
             * @brief The dereference operator
             *
             * This method returns a constant reference to the pair
             * identifier/mutations structure of the current node in
             * the tour. It takes constant time.
             *
             * @return a constant reference to the pair identifier/mutations
             *      structure of the current node in the tour
             */
            inline const value_type& operator*() const
            {
                return node_mutations;
            }

            /**
             * @brief Check whether the tour end has been reached
             *
             * @return `true` if and only if the tour end has been
             *      reached
             */
            inline const bool& is_end() const
            {
                return tour_end;
            }

            /**
             * @brief Equality operator
             *
             * This method checks whether two `const_iterator` are
             * the same.
             *
             * @param[in] rhs is the right-hand side of the equality
             * @return `true` if and only if this object and `rhs`
             *      iterates over the same tour and refer to the
             *      same node
             */
            bool operator==(const const_iterator& rhs) const
            {
                if (this->forest != rhs.forest) {
                    return false;
                }

                if (this->iterator_stack.size() != rhs.iterator_stack.size()) {
                    return false;
                }

                if (this->iterator_stack.size()==0) {
                    if (this->tour_end != rhs.tour_end) {
                        return false;
                    }

                    if (this->tour_end) {
                        return true;
                    }
                }

                return this->iterator_stack.top() == rhs.iterator_stack.top();
            }

            /**
             * @brief inequality operator
             *
             * This method checks whether two `const_iterator` differ.
             *
             * @param[in] rhs is the right-hand side of the inequality
             * @return `true` if and only if this object and `rhs`
             *      iterates over different tours or refer to different
             *      nodes
             */
            inline bool operator!=(const const_iterator& rhs) const
            {
                return !(*this == rhs);
            }

            friend class MutationTour<MUTATION_CONTAINER>;
        };

        /**
         * @brief The empty constructor
         */
        MutationTour():
            forest{nullptr}, only_leaves{false}, with_pre_neoplastic{true},
            with_germinal{false}
        {}

        /**
         * @brief Construct a new mutation tour for `GermlineMutations`
         *
         * @param[in] forest is a constant reference to the forest
         * @param[in] only_leaves is a Boolean flag to enable/disable internal node visit
         * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
         *      mutations
         * @param[in] with_germinal is Boolean flag to add/avoid germline mutations in the
         *      produced genome mutations
         */

        template<typename U = MUTATION_CONTAINER,
                 std::enable_if_t<std::is_same_v<U, MUTATION_CONTAINER>
                                  && std::is_same_v<U, GenomeMutations>, bool> = true>
        MutationTour(const PhylogeneticForest& forest, const bool only_leaves,
                     const bool with_pre_neoplastic, const bool with_germinal):
            MutationTour<MUTATION_CONTAINER>(&forest, only_leaves, with_pre_neoplastic,
                                             with_germinal)
        {}

        /**
         * @brief Construct a new mutation tour for `ChromosomeMutations`
         *
         * @param[in] forest is a constant reference to the forest
         * @param[in] chr_id is the identifier of the chromosomes that will be
         *      produced
         * @param[in] only_leaves is a Boolean flag to enable/disable internal node visit
         * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
         *      mutations
         * @param[in] with_germinal is Boolean flag to add/avoid germline mutations in the
         *      produced chromosome mutations
         */
        template<typename U = MUTATION_CONTAINER,
                 std::enable_if_t<std::is_same_v<U, MUTATION_CONTAINER>
                                  && std::is_same_v<U, ChromosomeMutations>, bool> = true>
        MutationTour(const PhylogeneticForest& forest, const ChromosomeId& chr_id,
                     const bool only_leaves, const bool with_pre_neoplastic,
                     const bool with_germinal):
            MutationTour<MUTATION_CONTAINER>(&forest, only_leaves, with_pre_neoplastic,
                                             with_germinal)
        {
            this->chr_id = chr_id;
        }
private:
        /**
         * @brief Get the constant iterator to the tour begin
         *
         * @param template_mutations is the mutations structure that may or may not
         *      contains the germinal mutations
         * @return the constant iterator to the tour begin
         */
        const_iterator begin(const MUTATION_CONTAINER& germline_mutations) const
        {
            if (with_germinal) {
                return const_iterator{forest, germline_mutations,
                                      only_leaves, with_pre_neoplastic};
            }

            const auto germinal_structure = germline_mutations.copy_structure();
            return const_iterator{forest, germinal_structure,
                                  only_leaves, with_pre_neoplastic};
        }
public:
        /**
         * @brief Get the constant iterator to the tour begin
         *
         * @return the constant iterator to the tour begin
         */
        const_iterator begin() const
        {
            if (forest == nullptr) {
                return const_iterator{forest, only_leaves};
            }

            if (forest->num_of_nodes()==0) {
                return const_iterator{forest, only_leaves};
            }

            const auto& germline_mutations = forest->get_germline_mutations();
            if constexpr(std::is_same_v<MUTATION_CONTAINER, GenomeMutations>) {
                return begin(germline_mutations);
            }

            if constexpr(std::is_same_v<MUTATION_CONTAINER, ChromosomeMutations>) {
                return begin(germline_mutations.get_chromosome(this->chr_id));
            }

            static_assert(true, "MutationTour::begin(): unsupported template type.");
        }

        /**
         * @brief Get the constant iterator to the tour end
         *
         * @return the constant iterator to the tour end
         */
        inline const_iterator end() const
        {
            return const_iterator{forest, only_leaves};
        }

        /**
         * @brief Get the associated forest
         *
         * @return a constant reference to the forest to which
         *      the forest is associated
         */
        inline const PhylogeneticForest& get_forest() const
        {
            return *forest;
        }

        friend class PhylogeneticForest;
    };

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
     * @brief Get a tour over the genome mutations of the forest leaves
     *
     * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
     *      mutations
     * @param[in] with_germinal is a Boolean flag to add/avoid germinal mutations
     * @return a tour over the genome mutations of the forest leaves
     */
    inline MutationTour<GenomeMutations>
    get_leaf_mutation_tour(const bool with_pre_neoplastic=true,
                           const bool with_germinal=false) const
    {
        return MutationTour<GenomeMutations>{this, true, with_pre_neoplastic,
                                             with_germinal};
    }

    /**
     * @brief Get a tour over the chromosome mutations of the forest leaves
     *
     * @param[in] chromosome_id is the identifier of the chromosome of interest
     * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
     *      mutations
     * @param[in] with_germinal is a Boolean flag to add/avoid germinal mutations
     * @return a tour over the genome mutations of the forest leaves
     */
    inline MutationTour<ChromosomeMutations>
    get_leaf_mutation_tour(const ChromosomeId& chromosome_id,
                           const bool with_pre_neoplastic=true,
                           const bool with_germinal=false) const
    {
        return MutationTour<ChromosomeMutations>{*this, chromosome_id, true,
                                                 with_pre_neoplastic,
                                                 with_germinal};
    }

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

}   // Mutations

}   // RACES

#endif // __RACES_PHYLOGENETIC_TREE__
