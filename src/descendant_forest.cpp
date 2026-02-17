/**
 * @file descendant_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for descendant forests
 * @version 1.5
 * @date 2026-02-17
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

#include <set>
#include <algorithm>

#include "descendant_forest.hpp"

#include "simulation.hpp"

namespace CLONES
{

namespace Mutants
{

DescendantForest::SpeciesData::SpeciesData()
{}

DescendantForest::SpeciesData::SpeciesData(const MutantId& mutant_id, const MethylationSignature& signature):
    mutant_id(mutant_id), signature(signature)
{}

DescendantForest::const_node DescendantForest::get_node(const CellId& cell_id) const
{
    if (cells.count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return DescendantForest::const_node(this, cell_id);
}

DescendantForest::node DescendantForest::get_node(const CellId& cell_id)
{
    if (cells.count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return DescendantForest::node(this, cell_id);
}

DescendantForest::DescendantForest()
{}

DescendantForest::DescendantForest(const Evolutions::Simulation& simulation):
    DescendantForest(simulation, simulation.get_tissue_samples())
{}

DescendantForest::DescendantForest(const Evolutions::Simulation& simulation,
                                     const std::list<Evolutions::TissueSample>& tissue_samples)
{
    std::set<MutantId> mutant_ids;

    auto species_properties = simulation.tissue().get_species_properties();
    for (const auto& s_properties: species_properties) {
        const auto& mutant_id = s_properties.get_mutant_id();
        mutant_ids.insert(mutant_id);

        const auto& signature = s_properties.get_methylation_signature();
        species_data.insert({s_properties.get_id(), {mutant_id, signature}});
    }

    for (const auto& mutant_id : mutant_ids) {
        mutant_names[mutant_id] = simulation.find_mutant_name(mutant_id);
    }

    Evolutions::BinaryLogger::CellReader reader(simulation.get_logger().get_directory());

    grow_from(tissue_samples, reader);
}

std::vector<DescendantForest::const_node> DescendantForest::get_roots() const
{
    std::vector<DescendantForest::const_node> nodes;

    for (const auto& root_id : roots) {
        nodes.push_back(DescendantForest::const_node(this, root_id));
    }

    return nodes;
}

std::vector<DescendantForest::node> DescendantForest::get_roots()
{
    std::vector<DescendantForest::node> nodes;

    for (const auto& root_id : roots) {
        nodes.push_back(DescendantForest::node(this, root_id));
    }

    return nodes;
}

std::vector<DescendantForest::const_node> DescendantForest::get_leaves() const
{
    std::vector<DescendantForest::const_node> leaves;

    for (const auto& [leaf_id, sample_idx] : coming_from) {
        leaves.push_back(DescendantForest::const_node(this, leaf_id));
    }

    return leaves;
}

std::string DescendantForest::get_species_name(const SpeciesId& species_id) const
{
    const auto& species_it = species_data.find(species_id);
    if (species_it == species_data.end()) {
        throw std::runtime_error("Unknown species id "+std::to_string(species_id));
    }

    const auto& data = species_it->second;

    return get_mutant_name(data.mutant_id) +
                MutantProperties::signature_to_string(data.signature);
}

bool is_crucial(const DescendantForest::const_node& node, const std::list<std::list<CellId>>& sticks_from)
{
    if (sticks_from.size()>1) {
        // if sticks_from contains more than one stick, then this node is
        // the MRCA of two crucial nodes and it is a crucial node too.
        return true;
    }

    if (node.is_root()) {
            // if node is a root, then it is crucial.
        return true;
    }

    // if the node and its parent have different species id,
    // then it is crucial
    return node.get_species_id() != node.parent().get_species_id();
}

std::list<CellId>
DescendantForest::collect_sticks_from(std::list<std::list<CellId>>& sticks, const CellId& cell_id,
                                       const double& birth_time_threshold) const
{
    std::list<std::list<CellId>> sticks_from;

    if (birth_time_threshold < cells.at(cell_id).get_birth_time()) {
        return {};
    }

    for (const auto& child_id : branches.at(cell_id)) {
        auto stick_from = collect_sticks_from(sticks, child_id, birth_time_threshold);

        if (stick_from.size()>0) {
            stick_from.push_front(cell_id);

            sticks_from.push_back(std::move(stick_from));
        }
    }

    const auto node = get_node(cell_id);

    // if sticks_from contains more than one stick, then this node is the MRCA of two
    // crucial nodes and it is a crucial node too.
    // if node is the tree root
    if (is_crucial(node, sticks_from))  {
        std::move(std::begin(sticks_from), std::end(sticks_from),
                  std::back_inserter(sticks));

        return {cell_id};
    }

    if (sticks_from.size()==1) {
        return sticks_from.front();
    }

    return {};
}

std::list<std::list<CellId>> DescendantForest::get_sticks(const double birth_time_threshold) const
{
    std::list<std::list<CellId>> sticks;

    for (const auto& root_id : roots) {
        collect_sticks_from(sticks, root_id, birth_time_threshold);
    }

    return sticks;
}
 
std::map<CellId, size_t> DescendantForest::get_node_depths() const
{
    std::map<CellId, size_t> depths;
    std::list<CellId> next{roots.begin(), roots.end()};

    for (const auto& root: roots) {
        depths[root] = 0;
        next.push_back(root);
    }

    while (!next.empty()) {
        auto found = branches.find(next.front());

        if (found != branches.end()) {
            const auto front_depth = depths.at(next.front());
            for (const auto& child: found->second) {
                next.push_back(child);

                depths.emplace(child, front_depth+1);
            }
        }

        next.pop_front();
    }

    return depths;
}

std::map<CellId, size_t>
count_descendants_in_subtree(const DescendantForest& forest,
                             const std::list<CellId>& leaf_ids)
{
    std::map<CellId, size_t> counter;

    for (const auto& cell_id: leaf_ids) {
        auto found = counter.find(cell_id);
        if (found == counter.end()) {
            counter.insert({cell_id, 1});
        } else {
            ++found->second;
        }

        auto cell_node = forest.get_node(cell_id);
        while (!cell_node.is_root()) {
            cell_node = cell_node.parent();

            ++counter[cell_node.get_id()];
        }
    }

    return counter;
}

CellId find_tree_coalescent(const DescendantForest& forest, const CellId& root,
                            const std::map<CellId, size_t>& counter)
{
    const auto& num_of_leaves = counter.at(root);
    auto cell_node = forest.get_node(root);

    while (true) {
        auto children = cell_node.children();

        if (children.size()>0) {
            size_t i{0};
            auto counter_it = counter.find(children[i].get_id());
            while (counter_it==counter.end()) {
                counter_it = counter.find(children[++i].get_id());
            }

            if (counter_it->second != num_of_leaves) {
                return cell_node.get_id();
            }

            cell_node = children[i];
        } else {
            return cell_node.get_id();
        }
    }
}

std::vector<CellId>
DescendantForest::get_coalescent_cells(const std::list<CellId>& cell_ids) const
{
    const auto counter = count_descendants_in_subtree(*this, cell_ids);

    std::vector<CellId> coalescent_nodes;
    coalescent_nodes.reserve(roots.size());

    for (const auto& root: roots) {
        if (counter.find(root) != counter.end()) {
            coalescent_nodes.push_back(find_tree_coalescent(*this, root, counter));
        }
    }

    return coalescent_nodes;
}

size_t
DescendantForest::height() const
{
    size_t curr_height{0};

    for (const auto& root: get_roots()) {
        curr_height = std::max(curr_height, root.height());
    }

    return curr_height;
}


std::map<Mutants::Evolutions::TissueSampleId, std::map<Mutants::CellId, size_t>>
DescendantForest::count_cells_per_root() const
{
    std::map<Mutants::Evolutions::TissueSampleId, std::map<Mutants::CellId, size_t>> counter;

    for (const auto& sample : get_samples()) {
        std::map<Mutants::CellId, size_t> forest_counter;
        for (const auto& root_id : get_root_cell_ids()) {
            forest_counter.emplace(root_id, 0);
        }
        counter.emplace(sample.get_id(), std::move(forest_counter));
    }

    for (const auto& root_id : get_root_cell_ids()) {
        std::list<Mutants::CellId> queue{root_id};

        while (!queue.empty()) {
            const auto node = this->get_node(queue.front());
            queue.pop_front();

            if (node.is_leaf()) {
                auto& sample_counter = counter[node.get_sample().get_id()];

                ++(sample_counter[root_id]);
            } else {
                for (const auto& child_id : branches.at(node.get_id())) {
                    queue.push_back(child_id);
                }
            }
        }
    }

    return counter;
}

std::vector<CellId>
DescendantForest::get_coalescent_cells() const
{
    std::list<CellId> leaf_ids;

    for (const auto& sample: get_samples()) {
        for (const auto& cell_id: sample.get_cell_ids()) {
            leaf_ids.push_back(cell_id);
        }
    }

    return get_coalescent_cells(leaf_ids);
}

DescendantForest
DescendantForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
    std::set<std::string> names(sample_names.begin(), sample_names.end());
    std::set<std::string> found_names;

    std::list<Evolutions::TissueSample> tissue_samples;
    for (const auto& sample : samples) {
        if (names.count(sample.get_name())>0) {
            tissue_samples.push_back(sample);
            found_names.insert(sample.get_name());
        }
    }

    if (names.size() != found_names.size()) {
        std::set<std::string> missing_names;
        std::set_difference(names.begin(), names.end(),
                            found_names.begin(), found_names.end(),
                            std::inserter(missing_names, missing_names.end()));

        std::ostringstream oss;
        oss << "Unknown sample names:";
        std::string sep = " ";
        for (const auto& missing_name : missing_names) {
            oss << sep << missing_name;
            sep = ", ";
        }

        throw std::domain_error(oss.str());
    }

    DescendantForest subforest = *this;
    subforest.grow_from(tissue_samples, cells);

    return subforest;
}

bool DescendantForest::is_leaf(const CellId& cell_id) const
{
    auto branches_it = branches.find(cell_id);

    if (branches_it == branches.end()) {
        throw std::domain_error(std::to_string(cell_id)+" is not a forest cell");
    }
    return branches_it->second.size()==0;
}

void DescendantForest::clear()
{
    cells.clear();
    roots.clear();
    branches.clear();
    samples.clear();
    coming_from.clear();
}

}   // Mutants

}   // CLONES
