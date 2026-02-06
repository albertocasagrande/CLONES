/**
 * @file phylogenetic_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for phylogenetic forests
 * @version 1.10
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

#include <set>
#include <algorithm>

#include "phylogenetic_forest.hpp"

#include "simulation.hpp"

namespace CLONES
{

namespace Mutations
{

PhylogeneticForest::PhylogeneticForest():
    DescendantsForest()
{}

PhylogeneticForest::const_node::const_node(const PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantsForest::_const_node<PhylogeneticForest>(forest, cell_id)
{}

PhylogeneticForest::node::node(PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantsForest::_node<PhylogeneticForest>(forest, cell_id)
{}

void PhylogeneticForest::node::add_new_mutation(const MutationSpec<SID>& mutation)
{
    get_forest().novel_mutations[cell_id].insert(mutation);
    get_forest().SID_first_cells[mutation].insert(cell_id);
}

void PhylogeneticForest::node::add_new_mutation(const CNA& cna)
{
    get_forest().novel_mutations[cell_id].insert(cna);
    get_forest().CNA_first_cells[cna].insert(cell_id);
}

PhylogeneticForest::const_node PhylogeneticForest::get_node(const Mutants::CellId& cell_id) const
{
    if (get_cells().count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return PhylogeneticForest::const_node(this, cell_id);
}

PhylogeneticForest::node PhylogeneticForest::get_node(const Mutants::CellId& cell_id)
{
    if (get_cells().count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return PhylogeneticForest::node(this, cell_id);
}

std::vector<PhylogeneticForest::const_node> PhylogeneticForest::get_roots() const
{
    std::vector<PhylogeneticForest::const_node> nodes;

    for (const auto& root_id : get_root_cells()) {
        nodes.push_back(PhylogeneticForest::const_node(this, root_id));
    }

    return nodes;
}

std::vector<PhylogeneticForest::node> PhylogeneticForest::get_roots()
{
    std::vector<PhylogeneticForest::node> nodes;

    for (const auto& root_id : get_root_cells()) {
        nodes.push_back(PhylogeneticForest::node(this, root_id));
    }

    return nodes;
}

template<typename MUTATION>
std::map<MUTATION, std::set<Mutants::CellId>>
filter_by_cells_in(const std::map<MUTATION, std::set<Mutants::CellId>>& mutation_map,
                   const Mutants::DescendantsForest& forest)
{
    std::map<MUTATION, std::set<Mutants::CellId>> filtered_map;

    for (const auto& [mutation, cell_ids]: mutation_map) {
        auto& in_forest = filtered_map[mutation];
        for (const auto& cell_id : cell_ids) {
            if (forest.get_cells().count(cell_id)>0) {
                in_forest.insert(cell_id);
            }
        }
    }

    return filtered_map;
}

PhylogeneticForest
PhylogeneticForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
    PhylogeneticForest forest;

    static_cast<Mutants::DescendantsForest&>(forest) = Mutants::DescendantsForest::get_subforest_for(sample_names);

    for (const auto& [cell_id, cell]: forest.get_cells()) {
        const auto novel_found = novel_mutations.find(cell_id);
        if (novel_found != novel_mutations.end()) {
            forest.novel_mutations.emplace(cell_id, novel_found->second);
        }

        if (forest.is_leaf(cell_id)) {
            const auto leaves_found = leaves_mutations.find(cell_id);
            if (leaves_found != leaves_mutations.end()) {
                forest.leaves_mutations.emplace(cell_id, leaves_found->second);
            }
        }
    }

    forest.SID_first_cells = filter_by_cells_in(SID_first_cells, forest);
    forest.CNA_first_cells = filter_by_cells_in(CNA_first_cells, forest);

    forest.germline_mutations = germline_mutations;
    forest.mutational_properties = mutational_properties;

    return forest;
}

const CellGenomeMutations&
PhylogeneticForest::get_leaf_mutations(const Mutants::CellId& cell_id) const
{
    if (is_leaf(cell_id)) {
        return *(leaves_mutations.at(cell_id));
    }

    throw std::domain_error(std::to_string(cell_id)+" is not a leaf of the forest");
}

std::list<SampleGenomeMutations> PhylogeneticForest::get_sample_mutations_list() const
{
    using namespace Mutants::Evolutions;

    std::list<SampleGenomeMutations> sample_mutations_list;
    std::map<TissueSampleId, SampleGenomeMutations*> sample_mutation_map;
    for (const auto& sample: get_samples()) {
        sample_mutations_list.emplace_back(sample.get_name(), germline_mutations);
        sample_mutation_map.insert({sample.get_id(), &(sample_mutations_list.back())});
    }

    for (const auto& [cell_id, mutations_ptr] : leaves_mutations) {
        const auto& sample = get_samples()[get_coming_from().at(cell_id)];
        auto& sample_genomic_mutations = *(sample_mutation_map.at(sample.get_id()));
        sample_genomic_mutations.mutations.push_back(mutations_ptr);
    }

    return sample_mutations_list;
}

size_t
find_sample_index(const std::vector<CLONES::Mutants::Evolutions::TissueSample>& sample_list,
                 const std::string& sample_name)
{
    size_t idx=0;
    for (const auto& sample: sample_list) {
        if (sample.get_name() == sample_name) {
            return idx;
        }
        ++idx;
    }

    throw std::runtime_error("Unknown sample \""+sample_name+"\".");
}

SampleGenomeMutations PhylogeneticForest::get_sample_mutations(const std::string& sample_name) const
{
    using namespace Mutants::Evolutions;

    SampleGenomeMutations sample_mutations(sample_name, germline_mutations);

    const size_t sample_idx = find_sample_index(get_samples(), sample_name);

    for (const auto& cell_id : get_samples()[sample_idx].get_cell_ids()) {
        sample_mutations.mutations.push_back(leaves_mutations.at(cell_id));
    }

    return sample_mutations;
}

std::map<Mutants::CellId, MutationList>
PhylogeneticForest::get_preneoplastic_mutations() const
{
    std::map<Mutants::CellId, MutationList> pnp_mutations;

    for (const auto& root: get_roots()) {
        auto it = novel_mutations.find(root.cell_id);

        if (it != novel_mutations.end()) {
            MutationList root_mutations;

            for (auto mut_it = it->second.begin(); mut_it != it->second.end();
                 ++mut_it) {

                switch(mut_it.get_type()) {
                    case MutationList::SID_TURN:
                    {
                        const auto& mutation = mut_it.get_last_SID();

                        if (mutation.nature == SID::PRENEOPLASTIC) {
                            root_mutations.insert(mutation);
                        }
                        break;
                    }
                    case MutationList::CNA_TURN:
                    {
                        const auto& mutation = mut_it.get_last_CNA();

                        if (mutation.nature == SID::PRENEOPLASTIC) {
                            root_mutations.insert(mutation);
                        }
                        break;
                    }
                    default:
                        break;
                }
            }

            pnp_mutations.insert({root.cell_id, std::move(root_mutations)});
        }
    }

    return pnp_mutations;
}

SampleGenomeMutations PhylogeneticForest::get_normal_sample(const std::string& name,
                                                            const bool& with_preneoplastic) const
{
    SampleGenomeMutations sample(name, germline_mutations);

    if (with_preneoplastic) {
        for (auto& [cell_id, cell_genome] : get_normal_genomes(with_preneoplastic)) {
            auto genome = std::make_shared<CellGenomeMutations>();

            std::swap(cell_genome, *genome);

            sample.mutations.push_back(genome);
        }
    } else {
        sample.mutations.push_back(std::make_shared<CellGenomeMutations>(germline_mutations->copy_structure()));
    }

    return sample;
}

std::map<Mutants::CellId, CellGenomeMutations>
PhylogeneticForest::get_normal_genomes(const bool& with_preneoplastic) const
{
    std::map<Mutants::CellId, CellGenomeMutations> normal_genomes;

    if (with_preneoplastic) {
        for (const auto& [cell_id, pnp_mutations] : get_preneoplastic_mutations()) {
            auto cell_genome = CellGenomeMutations(germline_mutations->copy_structure());
            for (auto mut_it = pnp_mutations.begin(); mut_it != pnp_mutations.end();
                ++mut_it) {

                switch(mut_it.get_type()) {
                    case MutationList::SID_TURN:
                    {
                        const auto& mutation = mut_it.get_last_SID();

                        cell_genome.apply(mutation);
                        break;
                    }
                    case MutationList::CNA_TURN:
                    {
                        const auto& mutation = mut_it.get_last_CNA();

                        cell_genome.apply(mutation);
                        break;
                    }
                    case MutationList::WGD_TURN:
                    {
                        cell_genome.duplicate_alleles();
                    }
                    default:
                        break;
                }
            }

            normal_genomes.insert({cell_id, std::move(cell_genome)});
        }
    } else {
        for (const auto& root : get_roots()) {
            normal_genomes.insert({root.get_id(),
                                   CellGenomeMutations(germline_mutations->copy_structure())});
        }
    }

    return normal_genomes;
}

std::map<ChromosomeId, std::set<ChrPosition>>
PhylogeneticForest::get_CNA_break_points() const
{
    std::map<ChromosomeId, std::set<ChrPosition>> b_points;

    for (const auto& [cell_id, cell_ptr] : get_leaves_mutations()) {
        for (const auto& [chr_id, cb_points] : cell_ptr->get_CNA_break_points()) {
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
PhylogeneticForest::get_allelic_count(const size_t& min_allelic_size) const
{
    auto b_points = get_CNA_break_points();

    PhylogeneticForest::AllelicCount allelic_count;

    for (const auto& [cell_id, mutations] : leaves_mutations) {
        auto allelic_map = mutations->get_allelic_map(b_points, min_allelic_size);

        update_allelic_count_on(allelic_count, allelic_map);
    }

    return allelic_count;
}

PhylogeneticForest::AllelicCount
PhylogeneticForest::get_allelic_count(const std::list<Mutants::CellId>& cell_ids,
                                      const size_t& min_allelic_size) const
{
    auto b_points = get_CNA_break_points();

    PhylogeneticForest::AllelicCount allelic_count;

    for (const auto& cell_id : cell_ids) {
        auto cell_it = leaves_mutations.find(cell_id);
        if (cell_it == leaves_mutations.end()) {
            throw std::domain_error("The cell "+std::to_string(cell_id)
                                    + " is not a leaf of the phylogenetic"
                                    + "forest.");
        }

        auto allelic_map = cell_it->second->get_allelic_map(b_points, min_allelic_size);

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
    DescendantsForest::clear();
    leaves_mutations.clear();
    novel_mutations.clear();
    SID_first_cells.clear();
    CNA_first_cells.clear();
}

}   // Mutants

}   // CLONES
