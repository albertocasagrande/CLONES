/**
 * @file position_set.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes to represent tissue position set
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

#include "position_set.hpp"

namespace CLONES
{

namespace Mutants
{

BasicPositionSet::const_iterator BasicPositionSet::begin() const
{
    throw std::runtime_error("The BasicPositionSet class must be inherited");
}

BasicPositionSet::const_iterator BasicPositionSet::end() const
{
    throw std::runtime_error("The BasicPositionSet class must be inherited");
}

BasicPositionSet::~BasicPositionSet()
{}

RectangleSet::const_iterator::const_iterator(const RectangleSet* rectangle,
                                             const Evolutions::PositionInTissue& position):
    rectangle(rectangle), pos(position)
{
}

RectangleSet::const_iterator::const_iterator():
    rectangle(nullptr)
{
}

template<typename T>
bool increase_and_update(T& value, const T& min, const T& max)
{
    if (value >= max) {
        value = min;

        return true;
    }

    ++value;
    return false;
}

void set_to_invalid_position(Evolutions::PositionInTissue& pos, const RectangleSet& rectangle)
{
    pos = rectangle.upper_corner;
    ++pos.x;
}

inline bool is_invalid_in(Evolutions::PositionInTissue& pos, const RectangleSet& rectangle)
{
    return ((pos.x<rectangle.lower_corner.x || pos.x>rectangle.upper_corner.x
            || pos.y<rectangle.lower_corner.y || pos.y>rectangle.upper_corner.y
            || pos.z<rectangle.lower_corner.z || pos.z>rectangle.upper_corner.z));
}

RectangleSet::const_iterator& RectangleSet::const_iterator::operator++()
{
    if (rectangle==nullptr) {
        return *this;
    }

    if (is_invalid_in(pos, *rectangle)) {
        throw std::runtime_error("The iterator has already reached the end of the set");
    }

    if (increase_and_update(pos.z, rectangle->lower_corner.z, rectangle->upper_corner.z)) {
        if (increase_and_update(pos.y, rectangle->lower_corner.y, rectangle->upper_corner.y)) {
            if (increase_and_update(pos.x, rectangle->lower_corner.x, rectangle->upper_corner.x)) {
                set_to_invalid_position(pos, *rectangle);

                return *this;
            }
        }
    }

    return *this;
}

template<typename T>
bool decrease_and_update(T& value, const T& min, const T& max)
{
    if (value <= min) {
        value = max;

        return true;
    }

    --value;
    return false;
}

RectangleSet::const_iterator& RectangleSet::const_iterator::operator--()
{
    if (rectangle==nullptr) {
        return *this;
    }

    if (is_invalid_in(pos, *rectangle)) {
        throw std::runtime_error("The iterator has already reached the end of the set");
    }

    if (decrease_and_update(pos.z, rectangle->lower_corner.z, rectangle->upper_corner.z)) {
        if (decrease_and_update(pos.y, rectangle->lower_corner.y, rectangle->upper_corner.y)) {
            if (decrease_and_update(pos.x, rectangle->lower_corner.x, rectangle->upper_corner.x)) {
                set_to_invalid_position(pos, *rectangle);

                return *this;
            }
        }
    }

    return *this;
}

RectangleSet::RectangleSet():
    RectangleSet({0,0,0})
{}

RectangleSet::RectangleSet(const Evolutions::PositionInTissue& position):
    RectangleSet(position, position)
{}

RectangleSet::RectangleSet(const Evolutions::PositionInTissue& lower_corner,
                           const Evolutions::PositionInTissue& upper_corner):
    lower_corner(lower_corner), upper_corner(upper_corner)
{
}

RectangleSet::RectangleSet(const Evolutions::PositionInTissue& lower_corner,
                            const Evolutions::AxisSize& x_size,
                            const Evolutions::AxisSize& y_size,
                            const Evolutions::AxisSize& z_size):
    lower_corner(lower_corner),
    upper_corner(lower_corner.x+x_size-1,lower_corner.y+y_size-1,lower_corner.z+z_size-1)
{
}

RectangleSet::RectangleSet(const Evolutions::PositionInTissue& lower_corner,
                           const Evolutions::AxisSize& x_size,
                           const Evolutions::AxisSize& y_size):
    lower_corner(lower_corner),
    upper_corner(lower_corner.x+x_size-1,lower_corner.y+y_size-1)
{
    if (lower_corner.z!=0) {
        throw std::domain_error("The lower corner must be a 2D position.");
    }
}

RectangleSet::const_iterator RectangleSet::begin() const
{
    return const_iterator(this, lower_corner);
}

RectangleSet::const_iterator RectangleSet::end() const
{
    const_iterator it(this, upper_corner);

    return ++it;
}

}   // Mutants

}   // CLONES

namespace std
{

std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::RectangleSet& rectangle)
{
    os << "RectangleSet{" << rectangle.lower_corner
       << "," << rectangle.upper_corner << "}";

    return os;
}

}   // std
