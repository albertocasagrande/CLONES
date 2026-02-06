/**
 * @file tissue.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define tissue class
 * @version 1.1
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

#include <random>

#include "tissue.hpp"
#include "mutant_properties.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

Tissue::SpeciesView::SpeciesView(const std::vector<Species>& species, const std::vector<size_t>& species_pos):
    species(species), species_pos(species_pos)
{}

Tissue::SpeciesView::const_iterator::const_iterator(const std::vector<Species>& species, const std::vector<size_t>::const_iterator it):
    species(&species), it(std::make_shared<std::vector<size_t>::const_iterator>(it))
{}

Tissue::SpeciesView::const_iterator Tissue::SpeciesView::const_iterator::operator++(int)
{
    Tissue::SpeciesView::const_iterator copy(*this);

    this->operator++();

    return copy;
}

Tissue::SpeciesView::const_iterator Tissue::SpeciesView::const_iterator::operator--(int)
{
    Tissue::SpeciesView::const_iterator copy(*this);

    this->operator--();

    return copy;
}

size_t Tissue::SpeciesView::num_of_cells() const
{
    size_t num_of_cells{0};
    for (auto species_it=begin(); species_it != end(); ++species_it) {
        num_of_cells += species_it->num_of_cells();
    }

    return num_of_cells;
}

Tissue::CellInTissueConstantProxy::CellInTissueConstantProxy(const Tissue &tissue, const PositionInTissue position):
    BaseCellInTissueProxy<const Tissue>(tissue, position)
{
}

Tissue::CellInTissueProxy::CellInTissueProxy(Tissue &tissue, const PositionInTissue position):
    BaseCellInTissueProxy<Tissue>(tissue, position)
{}

Tissue::CellInTissueProxy& Tissue::CellInTissueProxy::operator=(const Cell& cell)
{
    // erase the cell in the current position
    erase();

    Species& species = tissue->get_species(cell.get_species_id());

    // if the new cell is already in its species
    if (species.contains(cell.get_id())) {
        auto cell_ptr = &species(cell.get_id());

        // reset the corresponding tissue position
        tissue->cell_pointer(*cell_ptr) = nullptr;

        // update the cell position in the species
        *cell_ptr = position;
        tissue->cell_pointer(position) = cell_ptr;
    } else { // if, otherwise, the new cell is not in its species
        // add the cell to the species
        tissue->cell_pointer(position) = species.add(CellInTissue(cell, position));
    }

    return *this;
}

void Tissue::CellInTissueProxy::erase()
{
    // if the position already contains a cell with mutations
    if (!is_wild_type()) {
        CellInTissue*& space_ptr = tissue->cell_pointer(position);

        // remove the cell from its species
        Species& former_species = tissue->get_species(space_ptr->get_species_id());
        former_species.erase(space_ptr->get_id());

        space_ptr = nullptr;
    }
}

CellInTissue Tissue::CellInTissueProxy::copy_and_erase()
{
    CellInTissue copy = *this;

    erase();

    return copy;
}

void Tissue::CellInTissueProxy::switch_duplication(const bool duplication_on)
{
    // if the position already contains a cell
    if (!is_wild_type()) {
        CellInTissue*& space_ptr = tissue->cell_pointer(position);

        // switch duplication behaviour of the cell
        Species& species = tissue->get_species(space_ptr->get_species_id());
        species.switch_duplication_for(space_ptr->get_id(), duplication_on);
    }
}

Tissue::CellInTissueProxy::operator CellInTissue&()
{
    const auto ptr = tissue->cell_pointer(position);

    if (ptr!=nullptr) {
        return *ptr;
    }

    throw std::runtime_error("Wild-type cell");
}

Tissue::Tissue(const std::string& name, const std::vector<AxisSize>& sizes):
    name(name), dimensions(sizes.size())
{
    AxisSize z_size{1};
    if (dimensions==3) {
        z_size = sizes[2];
    } else {
        if (dimensions!=2) {
            throw std::domain_error("The tissue must be a either a 3D or 2D shape");
        }
    }

    std::vector<CellInTissue *> z_vector(z_size, nullptr);
    std::vector<std::vector<CellInTissue *>> y_vector(sizes[1], z_vector);
    space = TissueSpace(sizes[0], y_vector);
}

Tissue::Tissue(const std::string& name, const AxisSize x_size,const AxisSize y_size, const AxisSize z_size):
    Tissue(name, {x_size, y_size, z_size})
{
}

Tissue::Tissue(const std::string& name, const AxisSize x_size, const AxisSize y_size):
    Tissue(name, {x_size, y_size})
{
}

Tissue::Tissue(const std::string& name, const std::vector<MutantProperties>& mutants,
               const AxisSize  x_size, const AxisSize y_size, const AxisSize z_size):
    Tissue(name, {x_size, y_size, z_size})
{
    for (const auto& mutant: mutants) {
        add_mutant_species(mutant);
    }

    register_species_cells();
}

Tissue::Tissue(const std::string& name, const std::vector<MutantProperties>& mutants,
               const AxisSize x_size, const AxisSize y_size):
    Tissue(name, {x_size, y_size})
{
    for (const auto& mutant: mutants) {
        add_mutant_species(mutant);
    }

    register_species_cells();
}

Tissue::Tissue(const std::vector<AxisSize>& sizes):
    Tissue("", sizes)
{
}

Tissue::Tissue(const AxisSize x_size, const AxisSize y_size, const AxisSize z_size):
    Tissue("", {x_size, y_size, z_size})
{
}

Tissue::Tissue(const std::vector<MutantProperties>& mutants,
               const AxisSize x_size, const AxisSize y_size, const AxisSize z_size):
    Tissue("", mutants, x_size, y_size, z_size)
{
}

Tissue::Tissue(const std::vector<MutantProperties>& mutants,
               const AxisSize  x_size, const AxisSize  y_size):
    Tissue("", mutants, x_size, y_size)
{
}

size_t Tissue::num_of_mutated_cells() const
{
    size_t mutated{0};

    for (const auto& S: species) {
        mutated += S.num_of_cells();
    }

    return mutated;
}

std::vector<SpeciesProperties> Tissue::get_species_properties() const
{
    std::vector<SpeciesProperties> species_vector;

    for (const auto& species: *this) {
        species_vector.push_back(species);
    }

    return species_vector;
}

const Species& Tissue::get_species(const SpeciesId& species_id) const
{
    const auto pos_it = id_pos.find(species_id);
    if (pos_it == id_pos.end()) {
        throw std::out_of_range("Species identifier \""+
                                std::to_string(static_cast<int>(species_id))+"\" is unknown");
    }

    return species[pos_it->second];
}

Species& Tissue::get_species(const SpeciesId& species_id)
{
    const auto pos_it = id_pos.find(species_id);
    if (pos_it == id_pos.end()) {
        throw std::out_of_range("Species identifier \""+
                                std::to_string(static_cast<int>(species_id))+"\" is unknown");
    }

    return species[pos_it->second];
}


const Species& Tissue::get_species(const std::string& species_name) const
{
    const auto pos_it = name_pos.find(species_name);
    if (pos_it == name_pos.end()) {
        throw std::out_of_range("Species \""+species_name+"\" is unknown");
    }

    return species[pos_it->second];
}

Species& Tissue::get_species(const std::string& species_name)
{
    const auto pos_it = name_pos.find(species_name);
    if (pos_it == name_pos.end()) {
        throw std::out_of_range("Species \""+species_name+"\" is unknown");
    }

    return species[pos_it->second];
}

const CellInTissue& Tissue::place_cell(const CellId& id, const SpeciesId& species_id,
                                       const PositionInTissue position)
{
    if (!is_valid(position)) {
        throw std::runtime_error("The position is not in the tissue");
    }

    auto*& cell_ptr = space[position.x][position.y][position.z];

    if (cell_ptr!=nullptr) {
        throw std::runtime_error("The position is not free");
    }

    Species& species = get_species(species_id);

    cell_ptr = species.add(CellInTissue(id, species_id, position));

    return *cell_ptr;
}

void Tissue::register_species_cells()
{
    for (const auto& single_species : species) {
        for (auto& cell : single_species) {
            cell_pointer(cell) = const_cast<CellInTissue*>(&cell);
        }
    }
}

Tissue& Tissue::add_mutant(const MutantProperties& mutant)
{
    add_mutant_species(mutant);

    register_species_cells();

    return *this;
}

Tissue& Tissue::add_mutant_species(const MutantProperties& mutant)
{
    // check whether the mutant is already in the tissue
    if (mutant_pos.count(mutant.get_id())>0) {
        throw std::runtime_error("Clone already in the tissue");
    }

    // check whether any of the species is already in the tissue
    for (const auto& species: mutant.get_species()) {
        if (id_pos.count(species.get_id())>0) {
            throw std::runtime_error("Species id "
                                     + std::to_string(static_cast<int>(species.get_id()))
                                     + " already in the tissue");
        }
        if (name_pos.count(species.get_name())>0) {
            throw std::runtime_error("Species \""
                                     + species.get_name()
                                     + "\" already in the tissue");
        }
    }

    // insert the mutants in the tissue
    auto& pos = mutant_pos[mutant.get_id()];
    for (const auto& in_species: mutant.get_species()) {
        // place the new species at the end of the species vector
        pos.push_back(species.size());

        id_pos[in_species.get_id()] = species.size();

        name_pos[in_species.get_name()] = species.size();
        species.push_back(Species(in_species));
    }

    return *this;
}

Tissue::CellInTissueProxy Tissue::operator()(const PositionInTissue& position)
{
    return CellInTissueProxy(*this, position);
}

const Tissue::CellInTissueConstantProxy Tissue::operator()(const PositionInTissue& position) const
{
    return CellInTissueConstantProxy(*this, position);
}

size_t Tissue::count_mutated_cells_from(PositionInTissue position,
                                        const Direction& direction) const
{
    size_t counter{0};

    PositionDelta delta(direction);
    while (is_valid(position)&&cell_pointer(position)!=nullptr) {
        position += delta;
        ++counter;
    }

    return counter;
}

size_t Tissue::count_neighbors_in(const PositionInTissue position, const SpeciesId& species_id) const
{
    int16_t min_x = std::max(0, position.x-1);
    int16_t max_x = std::min(static_cast<int>(space.size()),
                             position.x+1);
    int16_t min_y = std::max(0, position.y-1);
    int16_t max_y = std::min(static_cast<int>(space[0].size()),
                             position.y+1);
    int16_t min_z = 0;
    int16_t max_z = 0;

    size_t counter{0};

    if (dimensions>2) {
        min_z = std::max(0, position.z-1);
        max_y = std::min(static_cast<int>(space[0][0].size()),
                         position.z+1);
    }

    for (int16_t x = min_x; x < max_x; ++x) {
        for (int16_t y = min_y; y < max_y; ++y) {
            for (int16_t z = min_z; z < max_z; ++z) {
                const auto& cell_ptr = space[x][y][z];
                if (cell_ptr != nullptr && cell_ptr->get_species_id() == species_id) {
                    ++counter;
                }
            }
        }
    }

    return counter;
}

std::list<Cell> Tissue::push_cells(const PositionInTissue from_position, const Direction& direction)
{
    PositionDelta delta(direction);
    PositionInTissue to_position(from_position+delta);
    CellInTissue* to_be_moved = cell_pointer(from_position);
    while (to_be_moved!=nullptr && is_valid(to_position)) {
        CellInTissue*& dest_ptr = cell_pointer(to_position);

        std::swap(dest_ptr, to_be_moved);

        *dest_ptr = to_position;

        to_position += delta;
    }

    std::list<Cell> lost_cell;

    if (to_be_moved!=nullptr) {
        lost_cell.push_back(*to_be_moved);

        Species &species = get_species(to_be_moved->get_species_id());

        species.erase(to_be_moved->get_id());
    }

    cell_pointer(from_position) = nullptr;

    return lost_cell;
}

std::vector<AxisSize> Tissue::size() const
{
    std::vector<AxisSize> sizes(dimensions);

    sizes[0] = static_cast<AxisSize>(space.size());
    sizes[1] = static_cast<AxisSize>(space[0].size());

    if (dimensions==3) {
        sizes[2] = static_cast<AxisSize>(space[0][0].size());
    }
    return sizes;
}

}   // Evolutions

}   // Mutants

}   // CLONES
