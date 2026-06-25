/**
 * @file tissue.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define tissue class
 * @version 1.10
 * @date 2026-06-25
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

#include "error.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

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

    throw Error<std::runtime_error>("Wild-type cell");
}

Tissue::Tissue(const std::string& name, const std::vector<AxisSize>& sizes):
    name(name), dimensions(sizes.size())
{
    AxisSize z_size{1};
    if (dimensions==3) {
        z_size = sizes[2];
    } else {
        if (dimensions!=2) {
            throw Error<std::domain_error>("The tissue must be a either a 3D or 2D shape. "
                                           + std::to_string(sizes.size()) + " dimensions "
                                           + "have been specified.");
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

Tissue::Tissue(const std::vector<AxisSize>& sizes):
    Tissue("", sizes)
{
}

Tissue::Tissue(const AxisSize x_size, const AxisSize y_size, const AxisSize z_size):
    Tissue("", {x_size, y_size, z_size})
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
        throw Error<std::out_of_range>("Species identifier \""
                                       + std::to_string(static_cast<uint>(species_id))
                                       + "\" is unknown.");
    }

    return species[pos_it->second];
}

Species& Tissue::get_species(const SpeciesId& species_id)
{
    const auto pos_it = id_pos.find(species_id);
    if (pos_it == id_pos.end()) {
        throw Error<std::out_of_range>("Species identifier \""
                                       + std::to_string(static_cast<uint>(species_id))
                                       + "\" is unknown.");
    }

    return species[pos_it->second];
}

const Species& Tissue::get_species(const std::string& species_name) const
{
    const auto pos_it = name_pos.find(species_name);
    if (pos_it == name_pos.end()) {
        throw Error<std::out_of_range>("Species \"" + species_name
                                       + "\" is unknown.");
    }

    return species[pos_it->second];
}

Species& Tissue::get_species(const std::string& species_name)
{
    const auto pos_it = name_pos.find(species_name);
    if (pos_it == name_pos.end()) {
        throw Error<std::out_of_range>("Species \"" + species_name
                                       + "\" is unknown.");
    }

    return species[pos_it->second];
}

const Species& Tissue::operator[](const std::string& species_name) const
{
    const auto pos_it = name_pos.find(species_name);
    if (pos_it == name_pos.end()) {
        throw Error<std::out_of_range>("Species \"" + species_name
                                       + "\" is unknown.");
    }

    return species[pos_it->second];
}

Species& Tissue::operator[](const std::string& species_name)
{
    const auto pos_it = name_pos.find(species_name);
    if (pos_it == name_pos.end()) {
        throw Error<std::out_of_range>("Species \"" + species_name
                                       + "\" is unknown.");
    }

    return species[pos_it->second];
}

const CellInTissue& Tissue::place_cell(const CellId& id, const SpeciesId& species_id,
                                       const PositionInTissue position)
{
    if (!is_valid(position)) {
        std::ostringstream oss;

        const auto size_sizes = size();

        oss << "Position " << position << " does not belong to the "
            << static_cast<uint>(size_sizes[0]) << "x"
            << static_cast<uint>(size_sizes[1]);

        if (size_sizes.size()>2) {
            oss << "x" << static_cast<uint>(size_sizes[2]);
        }

        oss << "-tissue.";
        throw Error<std::out_of_range>(oss.str());
    }

    auto*& cell_ptr = space[position.x][position.y][position.z];

    if (cell_ptr!=nullptr) {
        std::ostringstream oss;

        oss << "Position " << position << " is already taken by "
            << "a cancer cell.";
        throw Error<std::runtime_error>(oss.str());
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

Tissue& Tissue::add_species(const SpeciesProperties& species_properties)
{
    const auto& mutant_properties = species_properties.get_mutant_properties();

    auto mutant_it = mutant_pos.find(mutant_properties.get_id());
    if (mutant_it == mutant_pos.end()) {
        throw Error<std::domain_error>("Unknown mutant \""
                                       + mutant_properties.get_name()
                                       + "\".");
    }

    auto& species_pos = mutant_it->second;

    species_pos.push_back(species.size());

    Species new_species{species_properties};
    id_pos.emplace(new_species.get_id(), species_pos.back());
    name_pos.emplace(new_species.get_name(), species_pos.back());
    new_species.disable_death();

    species.push_back(std::move(new_species));

    return *this;
}

Tissue& Tissue::add_mutant(const MutantProperties& mutant)
{
    SpeciesName::validate_name(mutant.get_name());

    if (mutant_pos.find(mutant.get_id()) != mutant_pos.end()) {
        throw Error<std::domain_error>("The mutant \"" + mutant.get_name()
                                       + "\" has been already added.");
    }

    mutant_pos.emplace(mutant.get_id(), std::vector<size_t>{});

    if (epistate_names.size()==0) {
        SpeciesProperties properties{mutant};

        add_species(properties);
    } else {
        for (const auto& epistate_name: epistate_names) {

            SpeciesProperties properties{mutant, epistate_name};

            add_species(properties);
        }
    }

    return *this;
}

bool Tissue::knowns_species(const std::string& species_name) const
{
    SpeciesName name{species_name};

    MutantProperties mutant_prop{name.get_mutant_name()};

    if (mutant_pos.find(mutant_prop.get_id()) == mutant_pos.end()) {
        return false;
    }

    if (name.get_epistate_name()=="") {
        return true;
    }

    return (name_pos.find(species_name) != name_pos.end());
}

bool Tissue::knowns_mutant(const std::string& mutant_name) const
{
    for (const auto& [mutant_id, species_pos]: mutant_pos) {
        if (mutant_name == species[species_pos[0]].get_mutant_name()) {
            return true;
        }
    }

    return false;
}

bool Tissue::knowns_epigenetic_state(const std::string& epigenetic_state) const
{
    return (epistate_names.find(epigenetic_state) != epistate_names.end());
}

Tissue& Tissue::add_epigenetic_state(const std::string& epistate_name)
{
    SpeciesName::validate_name(epistate_name);

    if (epistate_names.count(epistate_name)>0) {
        throw Error<std::runtime_error>("\"" + epistate_name +
                                        "\" has been already added.");
    }

    std::list<MutantProperties> mutant_list;

    for (auto& [mutant_id, species_pos] : mutant_pos) {
        auto& base_species = species[species_pos.front()];

        mutant_list.push_back(base_species.get_mutant_properties());
    }

    if (epistate_names.size()==0) {
        for (const auto& curr_species : species) {
            if (curr_species.num_of_cells()>0) {
                throw Error<std::runtime_error>("Epigenetic states can be added to the simulation "
                                                "exclusively before placing any cell.");
            }
        }

        // delete all species
        for (const auto& [mutant_id, species_pos] : mutant_pos) {
            if (species_pos.size()!=1) {
                throw Error<std::runtime_error>("The species position should contains 1 position.");
            }
        }

        species.clear();
        id_pos.clear();
        name_pos.clear();

        for (auto& [mutant_id, species_pos] : mutant_pos) {
            species_pos = std::vector<size_t>();
        }
    }

    for (const auto& mutant_properties: mutant_list) {
        SpeciesProperties properties{mutant_properties, epistate_name};

        add_species(properties);
    }

    epistate_names.insert(epistate_name);

    return *this;
}

Tissue::const_mutant_view Tissue::get_mutant_view(const std::string& mutant_name) const
{
    MutantProperties properties(mutant_name);
    try {
        return get_mutant_view(properties.get_id());
    } catch (const Error<std::domain_error>& ex) {
        throw Error<std::domain_error>("Unknown mutant name \""
                                       + mutant_name + "\".");
    }
}

Tissue::const_mutant_view Tissue::get_mutant_view(const MutantId mutant_id) const
{
    auto mutant_it = mutant_pos.find(mutant_id);

    if (mutant_it == mutant_pos.end()) {
        throw Error<std::domain_error>("Unknown mutant id "
                                       + std::to_string(mutant_id) + ".");
    }

    return const_mutant_view(&species, mutant_it->second);
}

Tissue::mutant_view Tissue::get_mutant_view(const std::string& mutant_name)
{
    MutantProperties properties(mutant_name);
    try {
        return get_mutant_view(properties.get_id());
    } catch (const Error<std::domain_error>& ex) {
        throw Error<std::domain_error>("Unknown mutant name \""
                                       + mutant_name + "\".");
    }
}

Tissue::mutant_view Tissue::get_mutant_view(const MutantId mutant_id)
{
    auto mutant_it = mutant_pos.find(mutant_id);

    if (mutant_it == mutant_pos.end()) {
        throw Error<std::domain_error>("Unknown mutant id "
                                       + std::to_string(mutant_id) + ".");
    }

    return mutant_view(&species, mutant_it->second);
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

std::list<PositionInTissue>
Tissue::get_neighborhood_positions(const PositionInTissue& position) const
{
    std::list<PositionInTissue> positions;

    for (int16_t x=(position.x==0?0:position.x-1); x<position.x+2; ++x) {
        for (int16_t y=(position.y==0?0:position.y-1); y<position.y+2; ++y) {
            for (int16_t z=(position.z==0?0:position.z-1); z<position.z+2; ++z) {
                PositionInTissue pos(x, y, z);
                if (is_valid(pos)) {
                    positions.push_back(pos);
                }
            }
        }
    }

    return positions;
}

std::list<Cell> Tissue::push_cells(const PositionInTissue from_position, const Direction& direction,
                                   const bool duplicate_internal_cells)
{
    PositionDelta delta(direction);
    PositionInTissue to_position(from_position+delta);
    CellInTissue* to_be_moved = cell_pointer(from_position);
    bool old_on_border;
    if (duplicate_internal_cells) {
        old_on_border = operator()(from_position).is_on_border();
    }
    while (to_be_moved!=nullptr && is_valid(to_position)) {
        CellInTissue*& dest_ptr = cell_pointer(to_position);

        std::swap(dest_ptr, to_be_moved);

        *dest_ptr = to_position;

        if (duplicate_internal_cells) {
            Tissue::CellInTissueProxy cell_in_tissue = operator()(to_position);
            const bool on_border = cell_in_tissue.is_on_border();

            // switch cell status only when on border status changed
            if (old_on_border != on_border) {
                cell_in_tissue.switch_duplication(cell_in_tissue.is_on_border());

                old_on_border = on_border;
            }
        }

        to_position += delta;
    }

    if (duplicate_internal_cells && is_valid(to_position)) {
        to_position -= delta;
        Tissue::CellInTissueProxy cell_in_tissue = operator()(to_position);
        for (auto pos : cell_in_tissue.get_neighborhood_positions()) {
            Tissue::CellInTissueProxy neighbor = operator()(pos);
            if (!neighbor.is_wild_type()) {
                neighbor.switch_duplication(neighbor.is_on_border());
            }
        }
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
