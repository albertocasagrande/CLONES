/**
 * @file phylogenetic_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for phylogenetic forests
 * @version 1.15
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

#include <set>
#include <ranges>
#include <algorithm>

#include "phylogenetic_forest.hpp"

#include "simulation.hpp"

namespace RACES
{

namespace Mutations
{

PhylogeneticForest::SampleStatistics::SampleStatistics():
    total_allelic_size{0}, number_of_cells{0}
{}

PhylogeneticForest::PhylogeneticForest():
    DescendantForest()
{}

PhylogeneticForest::const_node::const_node(const PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantForest::_const_node<PhylogeneticForest>(forest, cell_id)
{}

PhylogeneticForest::node::node(PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantForest::_node<PhylogeneticForest>(forest, cell_id)
{}

const MutationList& PhylogeneticForest::const_node::pre_neoplastic_mutations() const
{
    auto found = forest->get_pre_neoplastic_mutations().find(get_id());

    if (found != forest->get_pre_neoplastic_mutations().end()) {
        return found->second;
    }

    throw std::runtime_error("const_node::pre_neoplastic_mutations(): The "
                             "node is not a forest root.");
}

PhylogeneticForest::const_node PhylogeneticForest::get_node(const Mutants::CellId& cell_id) const
{
    if (get_cells().count(cell_id)==0) {
        throw std::domain_error("The cell " + std::to_string(cell_id)
                                + " is not in the forest.");
    }

    return PhylogeneticForest::const_node(this, cell_id);
}

PhylogeneticForest::node PhylogeneticForest::get_node(const Mutants::CellId& cell_id)
{
    if (get_cells().count(cell_id)==0) {
        throw std::domain_error("The cell " + std::to_string(cell_id)
                                + " is not in the forest.");
    }

    return PhylogeneticForest::node(this, cell_id);
}

std::vector<PhylogeneticForest::const_node> PhylogeneticForest::get_roots() const
{
    std::vector<PhylogeneticForest::const_node> nodes;

    for (const auto& root_id : get_root_cell_ids()) {
        nodes.push_back(PhylogeneticForest::const_node(this, root_id));
    }

    return nodes;
}

std::vector<PhylogeneticForest::node> PhylogeneticForest::get_roots()
{
    std::vector<PhylogeneticForest::node> nodes;

    for (const auto& root_id : get_root_cell_ids()) {
        nodes.push_back(PhylogeneticForest::node(this, root_id));
    }

    return nodes;
}

template<typename MUTATION>
std::map<MUTATION, std::set<Mutants::CellId>>
filter_by_cells_in(const std::map<MUTATION, std::set<Mutants::CellId>>& mutation_map,
                   const Mutants::DescendantForest& forest)
{
    std::map<MUTATION, std::set<Mutants::CellId>> filtered_map;

    for (const auto& [mutation, cell_ids]: mutation_map) {
        std::set<Mutants::CellId>* in_forest = nullptr;
        for (const auto& cell_id : cell_ids) {
            if (forest.get_cells().count(cell_id)>0) {
                if (in_forest == nullptr) {
                    in_forest = &(filtered_map[mutation]);
                }
                in_forest->insert(cell_id);
            }
        }
    }

    return filtered_map;
}

PhylogeneticForest
PhylogeneticForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
    PhylogeneticForest forest;

    static_cast<Mutants::DescendantForest&>(forest) = Mutants::DescendantForest::get_subforest_for(sample_names);

    for (const auto& root_id : forest.get_root_cell_ids()) {
        auto found = pre_neoplastic_mutations.find(root_id);

        if (found != pre_neoplastic_mutations.end()) {
            forest.pre_neoplastic_mutations.emplace(root_id, found->second);
        }
    }

    for (const auto& cell_id: std::views::keys(forest.get_cells())) {
        auto found = arising_mutations.find(cell_id);
        if (found != arising_mutations.end()) {
            forest.arising_mutations.emplace(cell_id, found->second);
        }
    }

    for (const auto& sample : forest.get_samples()) {
        const auto& sample_id = sample.get_id();

        auto found = sample_statistics.find(sample_id);
        if (found != sample_statistics.end()) {
            forest.sample_statistics.emplace(sample_id, found->second);
        }
    }

    forest.SID_first_cells = filter_by_cells_in(SID_first_cells, forest);
    forest.CNA_first_cells = filter_by_cells_in(CNA_first_cells, forest);

    forest.germline_mutations = germline_mutations;
    forest.mutational_properties = mutational_properties;

    return forest;
}

CellGenomeMutations PhylogeneticForest::get_cell_mutations(const Mutants::CellId& cell_id,
                                                           const bool& with_pre_neoplastic,
                                                           const bool& with_germinal) const
{
    if (arising_mutations.count(cell_id) == 0) {
        throw std::domain_error("The cell " + std::to_string(cell_id)
                                + " is not in the forest.");
    }

    auto node = get_node(cell_id);

    const Mutants::Cell cell = static_cast<const Mutants::Cell&>(node);

    // find the branch from root to the node
    std::stack<Mutants::CellId> cell_id_stack;
    do {
        cell_id_stack.push(node.get_id());

        node = node.parent();
    } while (!node.is_root());

    // initialize the node genome mutations to the wild-type with or
    // without germline
    CellGenomeMutations node_mutations;
    if (with_germinal) {
        node_mutations = CellGenomeMutations{cell, *germline_mutations};
    } else {
        node_mutations = CellGenomeMutations{cell, germline_mutations->copy_structure()};
    }

    if (with_pre_neoplastic) {
        node_mutations.apply(pre_neoplastic_mutations.at(node.get_id()));
    }
    node_mutations.apply(arising_mutations.at(node.get_id()));

    // browse the branch down to the node and apply the arising mutations
    while (!cell_id_stack.empty()) {
        node_mutations.apply(arising_mutations.at(cell_id_stack.top()));

        cell_id_stack.pop();
    }

    return node_mutations;
}

std::map<Mutants::CellId, CellGenomeMutations>
PhylogeneticForest::get_wild_type_genomes(const bool with_pre_neoplastic,
                                          const bool with_germinal) const
{
    GenomeMutations germinal;

    if (with_germinal) {
        germinal = get_germline_mutations();
    } else {
        germinal = get_germline_mutations().copy_structure();
    }

    std::map<Mutants::CellId, CellGenomeMutations> wild_type_genomes;
    wild_type_genomes.emplace(EMBRYO_CELL_ID,
                              CellGenomeMutations{germinal});

    for (const auto& root_id : get_root_cell_ids()) {
        PhylogeneticForest::const_node root{this, root_id};

        const Mutants::Cell& root_cell = static_cast<const Mutants::Cell&>(root);
        CellGenomeMutations root_mutations{root_cell, germinal};

        if (with_pre_neoplastic) {
            root_mutations.apply(root.pre_neoplastic_mutations());
        }

        wild_type_genomes.emplace(root_id, std::move(root_mutations));
    }

    return wild_type_genomes;
}

std::map<ChromosomeId, std::set<ChrPosition>>
PhylogeneticForest::get_CNA_break_points() const
{
    std::map<ChromosomeId, std::set<ChrPosition>> b_points;

    for (const auto& [leaf_id, leaf_mutations] : get_leaf_mutation_tour()) {
        for (const auto& [chr_id, cb_points] : leaf_mutations.get_CNA_break_points()) {
            auto& chr_b_points = b_points[chr_id];
            for (const auto& cb_point : cb_points) {
                chr_b_points.insert(cb_point);
            }
        }
    }

    return b_points;
}

void update_allelic_count_on(PhylogeneticForest::AllelicCount& allelic_count,
                             const GenomeMutations::AllelicMap& allelic_map)
{
    std::greater<uint32_t> cmp;

    for (const auto& [chr_id, chr_allelic_map] : allelic_map) {
        auto& chr_allelic_count = allelic_count[chr_id];
        for (auto [pos, atype] : chr_allelic_map) {
            std::sort(atype.begin(), atype.end(), cmp);
            ++(chr_allelic_count[pos][atype]);
        }
    }
}

PhylogeneticForest::AllelicCount
PhylogeneticForest::get_allelic_count(const AlleleCounter& min_allelic_size) const
{
    auto b_points = get_CNA_break_points();

    PhylogeneticForest::AllelicCount allelic_count;

    for (const auto& [leaf_id, leaf_mutations] : get_leaf_mutation_tour()) {
        const auto allelic_map = leaf_mutations.get_allelic_map(b_points, min_allelic_size);

        update_allelic_count_on(allelic_count, allelic_map);
    }

    return allelic_count;
}

PhylogeneticForest::AllelicCount
PhylogeneticForest::get_allelic_count(const std::list<Mutants::CellId>& cell_ids,
                                      const AlleleCounter& min_allelic_size) const
{
    auto b_points = get_CNA_break_points();

    PhylogeneticForest::AllelicCount allelic_count;

    for (const auto& cell_id : cell_ids) {
        const auto cell_mutations = get_cell_mutations(cell_id);

        auto allelic_map = cell_mutations.get_allelic_map(b_points, min_allelic_size);

        update_allelic_count_on(allelic_count, allelic_map);
    }

    return allelic_count;
}


PhylogeneticForest::AllelicCount
PhylogeneticForest::get_allelic_count(const std::string& sample_name,
                                      const size_t& min_allelic_size) const
{
    for (const auto& sample : get_samples()) {
        if (sample.get_name() == sample_name) {
            return get_allelic_count(sample.get_cell_ids(), min_allelic_size);
        }
    }

    throw std::domain_error("\"" + sample_name +"\" is not a sample of the "
                            + "phylogenetic forest");
}

void PhylogeneticForest::clear()
{
    arising_mutations.clear();
    SID_first_cells.clear();
    CNA_first_cells.clear();
    sample_statistics.clear();
    DescendantForest::clear();
}

}   // Mutants

}   // RACES
