/**
 * @file position.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a position class in a tissue
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

#include "position.hpp"
#include "tissue.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

Direction::Direction():
    Direction::Direction(Direction::X_DOWN)
{}

Direction::Direction(const uint8_t value):
    Direction::Direction(value, 0)
{}

Direction::Direction(const uint8_t value, const uint8_t& dimension):
    bit_vector(value)
{
    if (!bit_vector) {
        throw std::domain_error("Unsupported direction");
    }

    for (size_t i=0; i<3; ++i) {
        const uint8_t axis = (value >> 2*i)&0x03;

        if (axis == 0x03) {
            throw std::domain_error("Unsupported direction");
        }
    }

    switch (dimension) {
        case 2:
            if (value>>4) {
                throw std::domain_error(("Direction specification contains "
                                         "z-axis component, but the "
                                         "specified dimension is ")
                                        + std::to_string(dimension));
            }
            break;
        case 0:
            if (value>>4) {
                bit_vector = bit_vector | (0x01 << 6);
            }
            break;
        case 3:
            bit_vector = bit_vector | (0x01 << 6);
            break;
        default:
            throw std::domain_error("Supported dimension are 2 and 3. "
                                    "Specified "
                                    + std::to_string(dimension) + ".");
    }
}

uint8_t Direction::num_of_dimensions() const
{
    if ((bit_vector>>6)&0x03) {
        return 3;
    }
    return 2;
}

int Direction::get_delta(const size_t& axis_index) const
{
    const uint8_t axis = get_component(axis_index)>>(2*axis_index);

    if (axis == 0x00) {
        return 0;
    }

    if (axis == 0x01) {
        return -1;
    }

    return 1;
}

//! @private
void collect_positive_dot_products(std::list<Direction>& directions, const uint8_t bitvector,
                                   const std::vector<std::list<uint8_t>>& components,
                                   const uint8_t axis_index,
                                   const std::vector<uint8_t>& admitted_changes,
                                   const uint8_t total_changes)
{
    if (axis_index<components.size()) {
        const uint8_t mask = (0x03 << 2*axis_index);

        if (admitted_changes[axis_index]
                 && static_cast<size_t>(total_changes+1)<components.size()) {
            for (const auto& component : components[axis_index]) {
                collect_positive_dot_products(directions, (bitvector & ~mask) | component,
                                              components, axis_index+1, admitted_changes,
                                              total_changes+1);
            }
        }
        collect_positive_dot_products(directions, bitvector,
                                      components, axis_index+1, admitted_changes,
                                      total_changes);
    } else {
        if (bitvector != 0) {
            directions.emplace_back(bitvector, components.size());
        }
    }
}

//! @private
std::list<Direction> collect_positive_dot_products(const std::vector<std::list<uint8_t>>& components,
                                                   const uint8_t& bitvector,
                                                   const std::vector<uint8_t>& admitted_changes)
{
    std::list<Direction> directions;

    collect_positive_dot_products(directions, bitvector, components, 0, admitted_changes, 0);

    return directions;
}

std::list<Direction> Direction::get_near_by_directions() const
{
    std::list<Direction> directions;

    const auto num_of_dims = num_of_dimensions();
    std::vector<std::list<uint8_t>> components(num_of_dims);
    std::vector<uint8_t> admitted_changes(num_of_dims, false);
    for (uint8_t axis_index=0; axis_index < num_of_dims; ++axis_index) {
        const auto shift = 2*axis_index;
        const auto component = get_component(axis_index) >> shift;
        auto& component_list = components[axis_index];

        switch (component) {
            case 0x00:
                component_list.push_back(0x01 << shift);
                component_list.push_back(0x02 << shift);
                break;
            case 0x01:
            case 0x02:
                component_list.push_back(0x00 << shift);

                if (num_of_dims==3) {
                    admitted_changes[(axis_index+1)%3] = true;
                    admitted_changes[(axis_index+2)%3] = true;
                } else {
                    admitted_changes[(axis_index+1)%2] = true;
                }
                break;
            default:
                throw std::out_of_range("Unsupported direction");
        }
    }

    return collect_positive_dot_products(components, bit_vector & ~(0x01 << 6),
                                         admitted_changes);
}

bool operator==(const Direction& a, const Direction& b)
{
    const uint8_t dim = a.num_of_dimensions();

    if (dim != b.num_of_dimensions()) {
        return false;
    }

    for (uint8_t i=0; i<dim; ++i) {
        if (a.get_component(i) != b.get_component(i)) {
            return false;
        }
    }

    return true;
}

DirectionGenerator::DirectionGenerator():
    initial_dir{Direction::X_DOWN}
{}

DirectionGenerator::DirectionGenerator(const uint8_t& direction_value, const uint8_t& direction):
    initial_dir{direction_value, direction}
{}

DirectionGenerator::DirectionGenerator(const Direction& initial_direction):
    initial_dir{initial_direction}
{}

DirectionGenerator::DirectionGenerator(Direction&& initial_direction):
    initial_dir{std::move(initial_direction)}
{}

DirectionGenerator::const_iterator::const_iterator(const Direction& initial_direction,
                                                   const uint32_t& counter):
    direction{initial_direction}, counter(counter)
{}

Direction DirectionGenerator::next(const Direction& direction)
{
    const uint8_t dim = direction.num_of_dimensions();

    uint8_t reminder{1};
    uint8_t bitvector{0};

    for (uint8_t i=0; i<dim; ++i) {
        uint8_t value = direction.get_component(i)>>(2*i);

        value += reminder;

        if (value == 0x03) {
            value = 0x00;

            reminder = 1;
        } else {
            reminder = 0;
        }

        bitvector = bitvector | (value << 2*i);
    }

    if (reminder>0) {
        bitvector = 0x01;
    }

    return Direction(bitvector, dim);
}

Direction DirectionGenerator::prev(const Direction& direction)
{
    const uint8_t dim = direction.num_of_dimensions();

    uint8_t reminder{1};
    uint8_t bitvector{0};

    for (uint8_t i=0; i<dim; ++i) {
        uint8_t value = direction.get_component(i)>>(2*i);

        if (reminder>0) {
            if (value == 0x00) {
                value = 0x02;

                reminder = 1;
            } else {
                value -= reminder;

                reminder = 0;
            }
        }

        bitvector = bitvector | (value << 2*i);
    }

    if (bitvector==0) {
        bitvector = 0x02 | (0x02 << 2);
        if (dim==3) {
            bitvector = bitvector | (0x02 << 4);
        }
    }
    return Direction(bitvector, dim);
}

DirectionGenerator::const_iterator& DirectionGenerator::const_iterator::operator++()
{
    const uint8_t dim = direction.num_of_dimensions();

    if (!((dim==2 && counter==8) || (dim==3 && counter==26))) {
        ++counter;

        direction = next(direction);
    }

    return *this;
}

DirectionGenerator::const_iterator DirectionGenerator::const_iterator::operator++(int)
{
    const_iterator prev_it(*this);

    this->operator++();

    return prev_it;
}

DirectionGenerator::const_iterator& DirectionGenerator::const_iterator::operator--()
{
    if (counter!=0) {
        --counter;

        direction = prev(direction);
    }

    return *this;
}

DirectionGenerator::const_iterator DirectionGenerator::const_iterator::operator--(int)
{
    const_iterator prev_it(*this);

    this->operator--();

    return prev_it;
}

DirectionGenerator::const_iterator DirectionGenerator::end() const
{
    switch (initial_dir.num_of_dimensions()) {
        case 2:
            return const_iterator(initial_dir, 8);
        case 3:
            return const_iterator(initial_dir, 26);
        default:
            throw std::runtime_error("Unsupported dimension "
                                     + std::to_string(initial_dir.num_of_dimensions()));
    }
}

bool operator==(const DirectionGenerator::const_iterator& a,
                const DirectionGenerator::const_iterator& b)
{
    const auto num_of_dimensions = a.direction.num_of_dimensions();
    if ((num_of_dimensions != b.direction.num_of_dimensions())
            || (a.counter!=b.counter)) {
        return false;
    }

    return a.direction==b.direction;
}

PositionDelta::PositionDelta(const int x, const int y, const int z):
    x(x), y(y), z(z)
{
}

PositionDelta::PositionDelta(const Direction& direction):
    PositionDelta(direction.get_delta_x(),direction.get_delta_y(),direction.get_delta_z())
{
}

PositionDelta PositionDelta::operator-() const
{
    return PositionDelta{-this->x,-this->y,-this->z};
}

PositionInTissue::PositionInTissue():
    x(0), y(0), z(0)
{}

PositionInTissue::PositionInTissue(const AxisPosition x, const AxisPosition y, const AxisPosition z):
    x(x), y(y), z(z)
{}

PositionInTissue& PositionInTissue::operator+=(const PositionDelta delta)
{
    x += delta.x;
    y += delta.y;
    z += delta.z;

    return *this;
}

PositionInTissue& PositionInTissue::operator-=(const PositionDelta delta)
{
    x -= delta.x;
    y -= delta.y;
    z -= delta.z;

    return *this;
}


PositionInTissue operator+(const PositionInTissue& position, const PositionDelta delta)
{
    PositionInTissue res{position};

    res += delta;

    return res;
}


PositionInTissue operator-(const PositionInTissue& position, const PositionDelta delta)
{
    PositionInTissue res{position};

    res -= delta;

    return res;
}

bool operator==(const PositionInTissue& p1, const PositionInTissue& p2)
{
    return p1.x==p2.x && p1.y==p2.y && p1.z==p2.z;
}


bool operator!=(const PositionInTissue& p1, const PositionInTissue& p2)
{
    return !(p1==p2);
}

size_t Manhattan_distance(const PositionInTissue& p1, const PositionInTissue& p2)
{
    size_t distance{0};

    distance += (p1.x>p2.x ? p1.x-p2.x : p2.x-p1.x);
    distance += (p1.y>p2.y ? p1.y-p2.y : p2.y-p1.y);
    distance += (p1.z>p2.z ? p1.z-p2.z : p2.z-p1.z);

    return distance;
}

Position::Position():
    PositionInTissue(), tissue(nullptr)
{}

Position::Position(Tissue& tissue, const AxisPosition& x, const AxisPosition& y, const AxisPosition& z):
    PositionInTissue(x, y, z), tissue(&tissue)
{
    if (!tissue.is_valid(*this)) {
        throw std::domain_error("Invalid position");
    }
}

Position::Position(Tissue& tissue, const PositionInTissue& pos):
    PositionInTissue(pos), tissue(&tissue)
{
    if (!tissue.is_valid(pos)) {
        throw std::domain_error("Invalid position");
    }
}

}   // Evolutions

}   // Mutants

}   // CLONES

namespace std
{

std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::Direction& direction)
{
    for (size_t i=0; i<3; ++i) {
        switch(direction.get_delta(i)) {
            case 1:
                os << "U";
                break;
            case -1:
                os << "D";
                break;
            case 0:
                os << "N";
                break;
            default:
                throw std::runtime_error("Unknown direction");
        }
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::PositionDelta& delta)
{
    os << "(" << delta.x << "," << delta.y << "," << delta.z << ")";

    return os;
}

std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::PositionInTissue& position)
{
    os << "(" << position.x <<","<< position.y <<","<< position.z <<")";

    return os;
}

std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::Position& position)
{
    os << "\""<< position.tissue->get_name() << "\""
       << static_cast<CLONES::Mutants::Evolutions::PositionInTissue>(position);

    return os;
}

}  // std
