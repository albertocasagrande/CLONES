/**
 * @file descendant_forest.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for descendant forests
 * @version 1.7
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

#ifndef __CLONES_DESCENDANT_FOREST__
#define __CLONES_DESCENDANT_FOREST__

#include <map>
#include <set>
#include <vector>
#include <queue>
#include <string>
#include <limits>
#include <ranges>
#include <functional>

#include "cell.hpp"
#include "binary_logger.hpp"
#include "tissue_sample.hpp"

namespace CLONES
{

namespace Mutants
{

/**
 * @brief The concept of cell id labelling function
 */
template <typename T>
concept CellIdLabellingFunction =
    std::is_constructible_v< std::function<std::string(const CellId&)>,
    std::decay_t<T> >;

/**
 * @brief A class representing descendant forests
 */
class DescendantForest
{
public:
    using SamplePosition = Evolutions::TissueSampleId;
private:
    /**
     * @brief Species data
     *
     * This structure privately stores species data.
     */
    struct SpeciesData
    {
        MutantId mutant_id;                //!< The species mutant id
        MethylationSignature signature;    //!< The species signature

        /**
         * @brief The constructor
         */
        SpeciesData(const MutantId& mutant_id, const MethylationSignature& signature);

        /**
         * @brief Save species data in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & mutant_id
                    & signature;
        }

        /**
         * @brief Load species data from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the load species data
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
        inline static SpeciesData load(ARCHIVE& archive)
        {
            SpeciesData data;

            archive & data.mutant_id
                    & data.signature;

            return data;
        }
private:

        /**
         * @brief Construct a new Species Data object
         *
         */
        SpeciesData();
    };

    using EdgeTail = std::set<CellId>;

    std::set<CellId> roots;                 //!< The cell ids of the forest roots
    std::map<CellId, Cell> cells;           //!< The forest cell id-cell map
    std::map<CellId, EdgeTail> branches;    //!< The descendant branches

    std::map<SpeciesId, SpeciesData> species_data;  //!< The species id to data map
    std::map<MutantId, std::string> mutant_names;   //!< The mutant id to mutant name map

    std::vector<Evolutions::TissueSample> samples;  //!< The vector of the samples that produced the forest
    std::map<CellId, SamplePosition> coming_from;   //!< The map associating each leaf to the sample which it comes from

    /**
     * @brief Grow a forest from the cells of a sample
     *
     * This function grows a forest from the cells in a sample.
     * The cells will be represented as the leaves of the forest and
     * their ancestors are loaded from a cell storage.
     *
     * @tparam CELL_STORAGE is the type of the cell storage
     * @param sample is the sample from which the forest must grow
     * @param cell_storage is the cell storage
     */
    template<typename CELL_STORAGE>
    void grow_from(const Evolutions::TissueSample& sample, CELL_STORAGE& cell_storage)
    {
        for (auto cell_id : sample.get_cell_ids()) {
            // the sampled cells are leaves and we can add them to branch
            branches.emplace(cell_id, EdgeTail{});

            bool new_branch{true};
            while (new_branch) {
                const auto cell = cell_storage.at(cell_id);
                cells.emplace(cell_id, cell);

                if (cell.has_parent()) {
                    const auto& parent_id = cell.get_parent_id();
                    auto found = branches.find(parent_id);

                    if (found == branches.end()) {
                        branches.emplace(parent_id, EdgeTail{cell_id});

                        cell_id = parent_id;
                    } else {
                        found->second.insert(cell_id);
                        new_branch = false;
                    }
                } else {
                    // it is a root
                    roots.insert(cell_id);

                    new_branch = false;
                }
            }
        }
    }

    /**
     * @brief Grow a forest from a list of tissue samples
     *
     * This method resets a forest and grows it from a list of tissue samples and
     * a cell storage. The cells in the samples will be the leaves of the forest
     * and their ancestors are loaded from a cell storage.
     *
     * @tparam CELL_STORAGE is the type of the cell storage
     * @param tissue_samples is a list of tissue samples coming from the simulation
     * @param cell_storage is the cell storage
     */
    template<typename CELL_STORAGE>
    void grow_from(const std::list<Evolutions::TissueSample>& tissue_samples,
                   CELL_STORAGE& cell_storage)
    {
        clear();

        samples = std::vector<Evolutions::TissueSample>(tissue_samples.begin(),
                                                        tissue_samples.end());

        SamplePosition i{0};
        for (const auto& sample: samples) {
            for (const auto& cell_id: sample.get_cell_ids()) {
                coming_from.emplace(cell_id, i);
            }
            ++i;
        }

        std::list<CellId> cell_ids;
        for (auto sample_it = tissue_samples.begin(); sample_it != tissue_samples.end(); ++sample_it) {
            grow_from(*sample_it, cell_storage);
        }
    }

    /**
     * @brief Collect sticks below a node
     *
     * See `DescendantForest::get_sticks() const` for the definition of stick.
     *
     * @param sticks is the list of sticks in the subtree rooted in `cell_id`
     * @param cell_id is the identifier of a cell represented in the descendant
     *              forest
     * @param birth_time_threshold is the maximum birth time for a
     *      cell associated to the returned sticks
     * @return the tail of a candidate stick that ends in the deepest crucial node
     *      in the subtree rooted in `cell_id`
     */
    std::list<CellId>
    collect_sticks_from(std::list<std::list<CellId>>& sticks, const CellId& cell_id,
                        const double& birth_time_threshold) const;

protected:

    /**
     * @brief Get the map associating each leaf to the sample which it comes from
     *
     * @return a constant reference to the map associating each leaf to the
     *      sample which it comes from
     */
    inline const std::map<CellId, SamplePosition>& get_coming_from() const
    {
        return coming_from;
    }

    /**
     * @brief A constant node of the forest
     *
     * @param FOREST_TYPE is the type of forest
     */

    template<typename FOREST_TYPE>
    class _const_node
    {
    public:
        FOREST_TYPE* forest;    //!< A pointer to the forest
        CellId cell_id;         //!< The cell id of a cell in the forest

        /**
         * @brief A constructor for a constant node
         *
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        _const_node(const FOREST_TYPE* forest, const CellId cell_id):
                forest(const_cast<FOREST_TYPE*>(forest)), cell_id(cell_id)
        {}

        /**
         * @brief Cast to Cell
         *
         * This method returns a constant reference to a cell.
         * Notice that the cells in the forest should be modified
         * exclusively by using `DescendantForest::node` methods
         *
         * @return a constant reference to a cell
         */
        inline operator const Cell&() const
        {
            return forest->cells.at(cell_id);
        }

        /**
         * @brief Get the parent node
         *
         * @tparam NODE_TYPE is the node type
         * @return the parent node of this node
         */
        template<typename NODE_TYPE = _const_node<FOREST_TYPE>>
        NODE_TYPE parent() const
        {
            if (forest==nullptr) {
                throw std::runtime_error("The forest node has not been initialized");
            }

            return NODE_TYPE(forest, forest->cells.at(cell_id).get_parent_id());
        }

        /**
         * @brief Get the node direct descendants
         *
         * @tparam NODE_TYPE is the node type
         * @return a vector containing the direct descendants of
         *         this node
         */
        template<typename NODE_TYPE = _const_node<FOREST_TYPE>>
        std::vector<NODE_TYPE> children() const
        {
            if (forest==nullptr) {
                throw std::runtime_error("The forest node has not been initialized");
            }

            std::vector<NODE_TYPE> nodes;

            for (const auto& child_id: forest->branches.at(cell_id)) {
                nodes.push_back(NODE_TYPE(forest, child_id));
            }

            return nodes;
        }

        /**
         * @brief Test whether this node is a leaf
         *
         * @return `true` if and only if this node is a leaf
         */
        inline bool is_leaf() const
        {
            return forest->branches.at(cell_id).size()==0;
        }

        /**
         * @brief Test whether this node is a root
         *
         * @return `true` if and only if this node is a root
         */
        inline bool is_root() const
        {
            const auto& cell = forest->cells.at(cell_id);

            return cell.get_id() == cell.get_parent_id();
        }

        /**
         * @brief Compute the node height
         *
         * @return the node height
         */
        size_t height() const
        {
            size_t curr_height{0};
            for (const auto& child: children()) {
                curr_height = std::max(curr_height, child.height()+1);
            }

            return curr_height;
        }

        /**
         * @brief Get the cell id of the node
         *
         * @return a constant reference to the cell id of the node
         */
        inline const CellId& get_id() const
        {
            return cell_id;
        }

        /**
         * @brief Get the sample that collected the cell
         *
         * @return the sample that collected the identifier
         *      of this node in the tissue
         * @throw `std::domain_error` when the node is not
         *      a leaf
         */
        const Evolutions::TissueSample& get_sample() const
        {
            if (forest==nullptr) {
                throw std::runtime_error("The forest node has not been initialized");
            }

            auto found = forest->coming_from.find(cell_id);

            if (found == forest->coming_from.end()) {
                throw std::domain_error("The node does not correspond to a sampled cell");
            }

            return forest->samples.at(found->second);
        }

        /**
         * @brief Get the node species id
         *
         * @return a constant reference to the node species id
         */
        inline const SpeciesId& get_species_id() const
        {
            return forest->cells.at(cell_id).get_species_id();
        }

        /**
         * @brief Get the node mutant id
         *
         * @return a constant reference to the node mutant id
         */
        inline const MutantId& get_mutant_id() const
        {
            return forest->species_data.at(get_species_id()).mutant_id;
        }

        /**
         * @brief Get the node species name
         *
         * @return the node species name
         */
        inline std::string get_species_name() const
        {
            return forest->get_species_name(get_species_id());
        }

        /**
         * @brief Get the node mutant name
         *
         * @return a constant reference to the node mutant name
         */
        inline const std::string& get_mutant_name() const
        {
            return forest->mutant_names.at(get_mutant_id());
        }

        /**
         * @brief Get the cell's birth time
         *
         * @return the cell's birth time
         */
        inline const Time& get_birth_time() const
        {
            return forest->cells.at(cell_id).get_birth_time();
        }

        /**
         * @brief Get the cell's death time
         *
         * @return the cell's death time
         */
        const Time& get_death_time() const
        {
            const auto& siblings = forest->branches.at(cell_id);
            if (siblings.empty()) {
                return get_sample().get_time();
            }

            const auto& a_child = *(siblings.begin());

            return forest->cells.at(a_child).get_birth_time();
        }

        /**
         * @brief Get the node methylation signature
         *
         * @return a constant reference to the node methylation signature
         */
        inline const MethylationSignature& get_methylation_signature() const
        {
            return forest->species_data.at(get_species_id()).signature;
        }

        /**
         * @brief Get the node forest
         *
         * @return a constant reference to the node forest
         */
        inline const FOREST_TYPE& get_forest() const
        {
            return *forest;
        }

        friend FOREST_TYPE;
    };

    /**
     * @brief A non-constant node of the forest
     *
     * @param FOREST_TYPE is the type of forest
     */
    template<typename FOREST_TYPE>
    class _node : public _const_node<FOREST_TYPE>
    {
    public:
        /**
         * @brief A constructor for a constant node
         *
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        _node(FOREST_TYPE* forest, const CellId cell_id):
            _const_node<FOREST_TYPE>(forest, cell_id)
        {}

        /**
         * @brief Cast to Cell
         *
         * @return a non-constant reference to a cell
         */
        inline operator Cell&()
        {
            return _const_node<FOREST_TYPE>::forest->cells.at(_const_node<FOREST_TYPE>::cell_id);
        }

        /**
         * @brief Get the parent node
         *
         * @tparam NODE_TYPE is the node type
         * @return the parent node of this node
         */
        template<typename NODE_TYPE = _node<FOREST_TYPE>>
        NODE_TYPE parent()
        {
            if (_const_node<FOREST_TYPE>::forest==nullptr) {
                throw std::runtime_error("The forest node has not been initialized");
            }

            return NODE_TYPE(_const_node<FOREST_TYPE>::forest,
                             _const_node<FOREST_TYPE>::forest->cells.at(_const_node<FOREST_TYPE>::cell_id).get_parent_id());
        }

        /**
         * @brief Get the node direct descendants
         *
         * @tparam NODE_TYPE is the node type
         * @return a vector containing the direct descendants of
         *         this node
         */
        template<typename NODE_TYPE = _node<FOREST_TYPE>>
        std::vector<NODE_TYPE> children()
        {
            if (_const_node<FOREST_TYPE>::forest==nullptr) {
                throw std::runtime_error("The forest node has not been initialized");
            }

            std::vector<NODE_TYPE> nodes;

            for (const auto& child_id: _const_node<FOREST_TYPE>::forest->branches.at(_const_node<FOREST_TYPE>::cell_id)) {
                nodes.push_back(NODE_TYPE(_const_node<FOREST_TYPE>::forest, child_id));
            }

            return nodes;
        }

        /**
         * @brief Get the node forest
         *
         * @return a constant reference to the node forest
         */
        inline FOREST_TYPE& get_forest()
        {
            return *(_const_node<FOREST_TYPE>::forest);
        }

        friend FOREST_TYPE;
    };

public:
    using const_node = _const_node<DescendantForest>;
    using node = _node<DescendantForest>;

    /**
     * @brief The empty constructor
     */
    DescendantForest();

    /**
     * @brief Construct a descendant forest by using a tissue sample of a simulation
     *
     * This method builds a descendant forest by using mutant simulation
     * pre-sampled cells as leaves.
     *
     * @param simulation is a simulation
     */
    DescendantForest(const Evolutions::Simulation& simulation);

    /**
     * @brief Construct a descendant forest by using a tissue sample of a simulation
     *
     * @param simulation is a simulation
     * @param tissue_samples is a list of tissue samples coming from the simulation
     */
    DescendantForest(const Evolutions::Simulation& simulation,
                      const std::list<Evolutions::TissueSample>& tissue_samples);

    /**
     * @brief Get a constant node with the given id
     *
     * @param cell_id is the id of the aimed cell node
     * @return the corresponding constant cell node
     */
    const_node get_node(const CellId& cell_id) const;

    /**
     * @brief Get a node with the given id
     *
     * @param cell_id is the id of the aimed cell node
     * @return the corresponding cell node
     */
    node get_node(const CellId& cell_id);

    /**
     * @brief Get the forest roots
     *
     * @return a vector containing constant nodes referring to
     *      forest roots
     */
    std::vector<const_node> get_roots() const;

    /**
     * @brief Get the number of forest leaves
     *
     * @return the number of forest leaves
     */
    inline size_t num_of_leaves() const
    {
        return coming_from.size();
    }

    /**
     * @brief Get the forest leaves
     *
     * @return a vector containing constant nodes referring to
     *      forest leaves
     */
    std::vector<const_node> get_leaves() const;

    /**
     * @brief Get the forest roots
     *
     * @return std::vector<const_node>
     */
    std::vector<node> get_roots();

    /**
     * @brief Get the number of forest nodes
     *
     * @return the number of forest nodes
     */
    inline size_t num_of_nodes() const
    {
        return cells.size();
    }

    /**
     * @brief Get the forest height
     *
     * @return the forest height
     */
    inline size_t height() const;

    /**
     * @brief Get the tissue samples that produced the forest
     *
     * @return a constant reference to a vector of tissue samples
     *      that produced the descendant forest
     */
    inline const std::vector<Evolutions::TissueSample>& get_samples() const
    {
        return samples;
    }

    /**
     * @brief Count how many leaves are progeny of each root
     *
     * @return A map that, for any sample and for any root, reports how many cells
     *      in the sample are represented by a leave that is progeny of the root.
     */
    std::map<Mutants::Evolutions::TissueSampleId, std::map<Mutants::CellId, size_t>>
    count_cells_per_root() const;

    /**
     * @brief Partition the forest samples
     *
     * This method partition the forest samples according to a labelling function that
     * associated a cell id to a string. The names of the resulting samples follow the
     * regular expression `<old sample name>_<label>`.
     *
     * @tparam LABELLING_FUNCTION is a cell id labelling function
     * @param labelling is a function labelling each forest cell identifier
     *      with a string
     */
    template<CellIdLabellingFunction LABELLING_FUNCTION>
    void partition_samples(const LABELLING_FUNCTION& labelling)
    {
        using TissueSample = Evolutions::TissueSample;
        std::map<std::string, TissueSample> new_samples;

        for (const auto& [cell_id, sample_pos] : coming_from) {
            const auto& old_sample = samples[sample_pos];

            std::string new_sample_name = old_sample.get_name()
                                            + "_" + labelling(cell_id);

            auto found = new_samples.find(new_sample_name);
            if (found == new_samples.end()) {
                TissueSample new_sample{new_sample_name,
                                        old_sample.get_time(),
                                        old_sample.get_bounding_box(),
                                        old_sample.get_tumour_cells_in_bbox()};

                new_sample.add_cell_id(cell_id);

                new_samples.emplace(new_sample_name, std::move(new_sample));
            } else {
                found->second.add_cell_id(cell_id);
            }
        }

        samples = std::vector<TissueSample>(new_samples.size());

        SamplePosition pos = 0;
        auto sample_it = samples.begin();
        for (auto& sample : std::views::values(new_samples)) {
            for (const auto& cell_id : sample.get_cell_ids()) {
                coming_from[cell_id] = pos;
            }

            std::swap(*sample_it,  sample);

            ++sample_it;
            ++pos;
        }
    }

    /**
     * @brief Get the most recent common ancestors of the forest leaves
     *
     * @return the most recent common ancestors of the forest leaves
     */
    std::vector<CellId> get_coalescent_cells() const;

    /**
     * @brief Get the most recent common ancestors
     *
     * @param cell_ids is the list of the cells whose most recent common
     *          ancestors are searched
     * @return the most recent common ancestors of the cells having
     *          the identifier among those in `cell_ids`
     */
    std::vector<CellId> get_coalescent_cells(const std::list<CellId>& cell_ids) const;

    /**
     * @brief Get the forest for a subset of the tissue samples
     *
     * @param sample_names are the names of the samples to be considered
     * @return the descendant forest for the tissue samples whose name is
     *          in `sample_names`
     */
    DescendantForest get_subforest_for(const std::vector<std::string>& sample_names) const;

    /**
     * @brief Check whether a cell is represented by a forest leaf
     *
     * @param cell_id is a cell identifier
     * @return `true` if and only if `cell_id` is the identifier of a cell represented
     *      by one of the forest leaves
     */
    bool is_leaf(const CellId& cell_id) const;

    /**
     * @brief Get the species data
     *
     * @return a constant reference to the species data.
     */
    inline const std::map<SpeciesId, SpeciesData>& get_species_data() const
    {
        return species_data;
    }

    /**
     * @brief Get the name of a mutant
     *
     * @param mutant_id is the identifier of the mutant whose name is aimed
     * @return a constant reference to the name of the mutant having
     *          `mutant_id` as identifier
     */
    inline const std::string& get_mutant_name(const MutantId& mutant_id) const
    {
        return mutant_names.at(mutant_id);
    }

    /**
     * @brief Get the name of a species
     *
     * @param species_id is the identifier of the species whose name is aimed
     * @return the name of the species having `species_id` as identifier
     */
    std::string get_species_name(const SpeciesId& species_id) const;

    /**
     * @brief Get the forest root cell identifiers
     *
     * @return a constant reference to the forest root cell identifiers
     */
    inline const std::set<CellId>& get_root_cell_ids() const
    {
        return roots;
    }

    /**
     * @brief Get the cell ids to cells maps
     *
     * @return a constant reference to the cell ids to cells maps
     */
    inline const std::map<CellId, Cell>& get_cells() const
    {
        return cells;
    }

    /**
     * @brief Get the forest sticks
     *
     * A _crucial node_ is a root of the forest, a node whose parent belongs
     * to a different species, or the most recent common ancestor of two
     * crucial nodes.
     * A _stick_ is a path of the forest in which the only crucial nodes are
     * the first and the last.
     *
     * @param birth_time_threshold is the maximum birth time for a cell
     *      associated to the returned sticks
     * @return a list of all the forest sticks whose associated cells have
     *      birth time smaller than or equal to `birth_time_threshold`
     */
    std::list<std::list<CellId>>
    get_sticks(const double birth_time_threshold) const;

    /**
     * @brief Get the forest sticks
     *
     * A _crucial node_ is a root of the forest, a node whose parent belongs
     * to a different species, or the most recent common ancestor of two
     * crucial nodes.
     * A _stick_ is a path of the forest in which the only crucial nodes are
     * the first and the last.
     *
     * @return a list of all the forest sticks
     */
    inline std::list<std::list<CellId>> get_sticks() const
    {
        return get_sticks(std::numeric_limits<double>::max());
    }

    /**
     * @brief Clear the forest
     */
    void clear();

    /**
     * @brief Save a descendant forest in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "CLONES Descendants Forest", 1);

        archive & roots
                & cells
                & branches
                & species_data
                & mutant_names
                & samples
                & coming_from;
    }

    /**
     * @brief Load a descendant forest from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load descendant forest
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static DescendantForest load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "CLONES Descendants Forest", 1);

        DescendantForest forest;

        archive & forest.roots
                & forest.cells
                & forest.branches
                & forest.species_data
                & forest.mutant_names
                & forest.samples
                & forest.coming_from;

        return forest;
    }
};

}   // Mutants

}   // CLONES

#endif // __CLONES_DESCENDANT_FOREST__
