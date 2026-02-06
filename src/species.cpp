/**
 * @file species.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements species representation methods
 * @version 1.2
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

#include "species.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

Species::Species():
    SpeciesProperties()
{}

void Species::reset()
{
    for (auto& [cell_id, cell_ptr]: cells) {
        delete cell_ptr;
    }

    cells.clear();
    duplication_enabled.clear();
}

Species::Species(const SpeciesProperties& species_properties):
    SpeciesProperties(species_properties), last_insertion_time(0),
    simulated_cells(0)
{}

Species::Species(const Species& orig):
    SpeciesProperties(orig), last_insertion_time(orig.last_insertion_time),
    simulated_cells(orig.simulated_cells)
{
    for (const auto& [cell_id, cell_ptr]: orig.cells) {
        this->add(*cell_ptr);
    }

    for (const auto& [cell_id, cell_ptr]: orig.duplication_enabled) {
        this->enable_duplication_for(cell_id);
    }
}

Species& Species::operator=(const Species& orig)
{
    reset();

    static_cast<SpeciesProperties&>(*this) = static_cast<SpeciesProperties>(orig);

    last_insertion_time = orig.last_insertion_time;
    simulated_cells = orig.simulated_cells;

    for (const auto& [cell_id, cell_ptr]: orig.cells) {
        this->add(*cell_ptr);
    }

    for (const auto& [cell_id, cell_ptr]: orig.duplication_enabled) {
        this->enable_duplication_for(cell_id);
    }

    return *this;
}

Species& Species::operator=(Species&& orig)
{
    reset();

    std::swap(*this, orig);

    return *this;
}

size_t Species::num_of_cells_available_for(const CellEventType& event_type) const
{
    switch(event_type) {
        case CellEventType::DEATH:
            return cells.size();
        case CellEventType::DUPLICATION:
        case CellEventType::EPIGENETIC_SWITCH:
        case CellEventType::MUTATION:
            return duplication_enabled.size();
        case CellEventType::ANY:
            return cells.size();
        default:
            throw std::domain_error("Unsupported event type");
    }
}

void Species::erase(const CellId& cell_id)
{
    // find the cell to be removed
    auto cell_it = cells.find(cell_id);

    if (cell_it == cells.end()) {
        throw std::domain_error("Unknown cell "+std::to_string(cell_id));
    }

    // delete it from the duplication enabled map
    duplication_enabled.erase(cell_id);

    // delete the cell
    delete cell_it->second;

    // remove the cell pointer
    cells.erase(cell_it);
}

CellInTissue* Species::add(CellInTissue* cell)
{
    // update the species id
    cell->species_id = get_id();

    // update `last_insertion_time`
    if (cell->get_birth_time()>last_insertion_time) {
        last_insertion_time = cell->get_birth_time();
    }

    cells.insert({cell->get_id(), cell});
    duplication_enabled.insert({cell->get_id(), cell});

    ++simulated_cells;

    return cell;
}

CellInTissue* Species::add(CellInTissue& cell)
{
    cell.species_id = get_id();

    return add(new CellInTissue(cell));
}

void Species::switch_duplication_for(const CellId& cell_id,
                                     const bool duplication_on)
{
    if (duplication_on) {
        enable_duplication_for(cell_id);
    } else {
        disable_duplication_for(cell_id);
    }
}

void Species::enable_duplication_for(const CellId& cell_id)
{
    auto cell_it = cells.find(cell_id);
    if (cell_it == cells.end()) {
        throw std::domain_error("Unknown cell "+std::to_string(cell_id));
    }

    duplication_enabled.insert({cell_id, cell_it->second});
}

void Species::disable_duplication_for(const CellId& cell_id)
{
    auto cell_it = cells.find(cell_id);
    if (cell_it == cells.end()) {
        throw std::domain_error("Unknown cell "+std::to_string(cell_id));
    }

    duplication_enabled.erase(cell_id);
}

CellInTissue& Species::operator()(const CellId& cell_id)
{
    auto cell_it = cells.find(cell_id);

    if (cell_it == cells.end()) {
        throw std::domain_error("Unknown cell "+std::to_string(cell_id));
    }

    return *(cell_it->second);
}

const CellInTissue& Species::operator()(const CellId& cell_id) const
{
    auto cell_it = cells.find(cell_id);

    if (cell_it == cells.end()) {
        throw std::domain_error("Unknown cell "+std::to_string(cell_id));
    }

    return *(cell_it->second);
}

Species::~Species()
{
    reset();
}

Species::const_iterator Species::begin() const
{
    return const_iterator(std::begin(cells));
}

Species::const_iterator Species::end() const
{
    return const_iterator(std::end(cells));
}

Species::const_iterator::const_iterator(const Species::CellIdToCell::const_iterator it): it(it)
{}

Species::const_iterator::const_iterator(): it()
{}

Species::const_iterator Species::const_iterator::operator++(int)
{
    const_iterator copy = *this;
    ++(*this);
    return copy;
}

Species::const_iterator Species::const_iterator::operator--(int)
{
    const_iterator copy = *this;
    --(*this);
    return copy;
}


void swap(Species& a, Species& b)
{
    std::swap(static_cast<SpeciesProperties&>(a),
              static_cast<SpeciesProperties&>(b));
    std::swap(a.cells,b.cells);
    std::swap(a.duplication_enabled,b.duplication_enabled);
    std::swap(a.last_insertion_time,b.last_insertion_time);
    std::swap(a.simulated_cells,b.simulated_cells);
}

}   // Evolutions

}   // Mutants

}   // CLONES

namespace std
{

std::ostream& operator<<(std::ostream& out, const CLONES::Mutants::Evolutions::Species& species)
{
    out << "{species_properties: " << static_cast<CLONES::Mutants::SpeciesProperties>(species)
        << ", cells: {";
    std::string sep{""};
    for (const auto& cell: species) {
        out << sep << cell;
        sep = ",";
    }
    out <<"}}";

    return out;
}

}   // std
