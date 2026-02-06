/**
 * @file cell.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements cell representation
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

#include "cell.hpp"
#include "species.hpp"
#include "tissue.hpp"
#include "mutant_properties.hpp"

namespace CLONES
{

namespace Mutants
{

Cell::Cell():
    Cell(0, WILD_TYPE_SPECIES, 0, 0)
{}

Cell::Cell(const CellId& id, const SpeciesId species_id):
    Cell(id, species_id, id, 0)
{}

Cell::Cell(const CellId& id, const SpeciesId species_id, const CellId parent_id):
    Cell(id, species_id, parent_id, 0)
{}

Cell::Cell(const CellId& id, const SpeciesId species_id, const CellId parent_id, const Time birth_time):
    id(id), parent(parent_id), birth_time(birth_time), species_id(species_id)
{}

void swap(Cell& a, Cell &b)
{
    std::swap(a.id,b.id);
    std::swap(a.parent, b.parent);
    std::swap(a.birth_time, b.birth_time);
    std::swap(a.species_id, b.species_id);
}

namespace Evolutions
{


CellInTissue::CellInTissue():
    Cell(), PositionInTissue()
{}

CellInTissue::CellInTissue(const CellId& id, const SpeciesId& species_id, const PositionInTissue& position):
    Cell(id, species_id), PositionInTissue(position)
{
}

CellInTissue::CellInTissue(const PositionInTissue& position):
    Cell(), PositionInTissue(position)
{
}

CellInTissue::CellInTissue(const AxisPosition& x, const AxisPosition& y, const AxisPosition& z):
    Cell(), PositionInTissue(x, y, z)
{
}

CellInTissue::CellInTissue(const Cell& cell, const PositionInTissue& position):
    Cell(cell), PositionInTissue(position)
{
}

CellInTissue& CellInTissue::operator=(const Cell& cell)
{
    Cell::operator=(cell);

    return *this;
}


CellInTissue& CellInTissue::operator=(const PositionInTissue& position)
{
    x = position.x;
    y = position.y;
    z = position.z;

    return *this;
}

}   // Evolutions

}   // Mutants

}   // CLONES

namespace std
{

std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Cell& cell)
{
    os << "Cell{id: "<< cell.get_id()
       << ", parent_id: "<< cell.get_parent_id()
       << ", birth_time: "<< cell.get_birth_time()
       << ", species_id: " << cell.get_species_id()
       << "}";

    return os;
}

std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::CellInTissue& cell)
{
    os << "Cell{id: "<< cell.get_id()
       << ", parent_id: "<< cell.get_parent_id()
       << ", birth_time: "<< cell.get_birth_time()
       << ", species_id: " << cell.get_species_id()
       << ", position: " << static_cast<CLONES::Mutants::Evolutions::PositionInTissue>(cell)
       << "}";

    return os;
}

}  // std
