/**
 * @file position.hpp
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

#ifndef __CLONES_POSITION__
#define __CLONES_POSITION__

#include <iostream>
#include <type_traits>

#include "archive.hpp"



namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

struct Direction;

}

}

}

namespace std
{

/**
 * @brief Write a direction in an output stream
 *
 * @param os is the output stream
 * @param direction is the direction to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::Direction& direction);

};

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

class Tissue;

/**
 * @brief A class to represent 2D/3D directions
 */
struct Direction {
    /**
     * @brief Canonical directions
     */
    enum Value {
        X_NULL = 0x00,
        X_DOWN = 0x01,
        X_UP = 0x02,
        Y_NULL = 0x00 << 2,
        Y_DOWN = 0x01 << 2,
        Y_UP = 0x02 << 2,
        Z_NULL = 0x00 << 4,
        Z_DOWN = 0x01 << 4,
        Z_UP = 0x02 << 4
    };

    /**
     * @brief The empty constructor
     */
    Direction();

    /**
     * @brief A constructor
     *
     * @param value is a direction obtained by composing `Value`'s
     */
    Direction(const uint8_t value);

    /**
     * @brief A constructor
     *
     * @param value is a direction obtained by composing `Value`'s
     * @param dimension is the space dimension
     */
    Direction(const uint8_t value, const uint8_t& dimension);

    /**
     * @brief Get the direction space dimension
     *
     * @return the direction space dimension
     */
    uint8_t num_of_dimensions() const;

    /**
     * @brief Get the x-axis delta
     *
     * @return the x-axis delta
     */
    inline int get_delta_x() const
    {
        return get_delta(0);
    }

    /**
     * @brief Get the y-axis delta
     *
     * @return the y-axis delta
     */
    inline int get_delta_y() const
    {
        return get_delta(1);
    }

    /**
     * @brief Get the z-axis delta
     *
     * @return the z-axis delta
     */
    inline int get_delta_z() const
    {
        return get_delta(2);
    }

    /**
     * @brief Get the x-axis component of the direction
     *
     * @return the x-axis component of the direction
     */
    inline Value get_component_x() const
    {
        return get_component(0);
    }

    /**
     * @brief Get the y-axis component of the direction
     *
     * @return the y-axis component of the direction
     */
    inline Value get_component_y() const
    {
        return get_component(1);
    }

    /**
     * @brief Get the z-axis component of the direction
     *
     * @return the z-axis component of the direction
     */
    inline Value get_component_z() const
    {
        return get_component(2);
    }

    /**
     * @brief Get the near-by directions
     *
     * Two directions are near-by if their dot-produce is
     * positive.
     *
     * @return a list of the near-by directions of the
     *      this direction
     */
    std::list<Direction> get_near_by_directions() const;

    /**
     * @brief Get the direction component on one axis
     *
     * @param axis_index is the index of an axis
     * @return the direction component on one axis
     */
    inline Value get_component(const size_t& axis_index) const
    {
        return static_cast<Value>(bit_vector&(0x03<<(2*axis_index)));
    }

    friend std::ostream& std::operator<<(std::ostream& os, const Direction& direction);
private:
    uint8_t bit_vector;     //!< the bit vector representation of the direction

    /**
     * @brief Get the direction delta on one axis
     *
     * @param axis_index is the index of an axis
     * @return the direction delta on one axis
     */
    int get_delta(const size_t& axis_index) const;
};

/**
 * @brief Test whether two directions are the same
 *
 * @param a is a direction
 * @param b is a direction
 * @return `true` if and only if `a` and `b` are the same
 */
bool operator==(const Direction& a, const Direction& b);

/**
 * @brief Test whether two directions differ
 *
 * @param a is a direction
 * @param b is a direction
 * @return `true` if and only if `a` and `b` differ
 */
inline bool operator!=(const Direction& a, const Direction& b)
{
    return !(a==b);
}

/**
 * @brief Generate all the directions
 *
 * This class defines a direction order and supports iterations
 * over it.
 */
class DirectionGenerator
{
    Direction initial_dir;  //!< the initial direction of the generator
public:

    /**
     * @brief The constant iterator over all the directions
     */
    class const_iterator
    {
        Direction direction;    //!< the current direction

        uint32_t counter;      //!< the nu

        /**
         * @brief Construct a new const iterator
         *
         * @param direction is the direction
         * @param counter is the number of already visited directions
         */
        const_iterator(const Direction& direction, const uint32_t& counter);
    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   Direction;
        using pointer           =   const Direction*;
        using reference         =   const Direction&;
        using iterator_category =   std::bidirectional_iterator_tag;

        /**
         * @brief Get the pointed direction
         *
         * @return the pointed direction
         */
        inline reference operator*() const
        {
            return direction;
        }

        /**
         * @brief Get the pointer to the direction
         *
         * @return the pointer to the direction
         */
        inline pointer operator->() const
        {
            return &direction;
        }

        /**
         * @brief Prefix increment
         *
         * @return an updated reference to the updated
         *      constant iterator
         */
        const_iterator& operator++();

        /**
         * @brief Postfix increment
         *
         * @return a copy of the original iterator
         */
        const_iterator operator++(int);

        /**
         * @brief Prefix decrement
         *
         * @return an updated reference to the updated
         *      constant iterator
         */
        const_iterator& operator--();

        /**
         * @brief Postfix decrement
         *
         * @return a copy of the original iterator
         */
        const_iterator operator--(int);

        friend bool operator==(const const_iterator& a, const const_iterator& b);

        friend class DirectionGenerator;
    };

    /**
     * @brief Build a direction generator
     */
    DirectionGenerator();

    /**
     * @brief Build a direction generator
     *
     * @param direction_value is a disjunction of `Direction::Value`s
     * @param dimension is the direction space dimension
     */
    DirectionGenerator(const uint8_t& direction_value, const uint8_t& dimension);

    /**
     * @brief Build a direction generator
     *
     * @param initial_direction is the initial direction of the generator
     */
    DirectionGenerator(const Direction& initial_direction);

    /**
     * @brief Build a direction generator
     *
     * @param initial_direction is the initial direction of the generator
     */
    DirectionGenerator(Direction&& initial_direction);

    /**
     * @brief Get the successive direction in the direction order
     *
     * @param direction is a direction
     * @return the direction that comes just after `direction` in
     *      the direction order
     */
    static Direction next(const Direction& direction);

    /**
     * @brief Get the previous direction in the direction order
     *
     * @param direction is a direction
     * @return the direction that comes just before `direction` in
     *      the direction order
     */
    static Direction prev(const Direction& direction);

    /**
     * @brief Get the initial iterator of the generator
     *
     * @return the initial iterator of the generator
     */
    inline const_iterator begin() const
    {
        return const_iterator(initial_dir, 0);
    }

    /**
     * @brief Get the final iterator of the generator
     *
     * @return the final iterator of the generator
     */
    const_iterator end() const;
};

/**
 * @brief Test whether two iterators of a direction generator are equal
 *
 * @param a is a direction generator constant iterator
 * @param b is a direction generator constant iterator
 * @return `true` if and only if `a` and `b` are the same
 */
bool operator==(const DirectionGenerator::const_iterator& a,
                const DirectionGenerator::const_iterator& b);

/**
 * @brief Test whether two iterators of a direction generator differ
 *
 * @param a is a direction generator constant iterator
 * @param b is a direction generator constant iterator
 * @return `true` if and only if `a` and `b` differ
 */
inline bool operator!=(const DirectionGenerator::const_iterator& a,
                       const DirectionGenerator::const_iterator& b)
{
    return !(a==b);
}

/**
 * @brief A class representing the difference between two positions
 */
struct PositionDelta {
    int x;  //!< x axis
    int y;  //!< y axis
    int z;  //!< z axis

    /**
     * @brief Create a new position delta
     *
     * @param x is the x-axis delta
     * @param y is the y-axis delta
     * @param z is the z-axis delta
     */
    PositionDelta(const int x, const int y, const int z);

    /**
     * @brief Create a new position delta towards a direction
     *
     * @param direction
     */
    explicit PositionDelta(const Direction& direction);

    /**
     * @brief Negate position delta;
     *
     */
    PositionDelta operator-() const;
};

/**
 * @brief The type of a position on an axis
 */
using AxisPosition = uint16_t;

/**
 * @brief A 3D position in a tissues
 */
struct PositionInTissue {
    int16_t x;  //!< x axis
    int16_t y;  //!< y axis
    int16_t z;  //!< z axis

    /**
     * @brief A constructor
     */
    PositionInTissue();

    /**
     * @brief A new constructor
     *
     * @param x is the x-axis position
     * @param y is the y-axis position
     * @param z is the z-axis position
     */
    PositionInTissue(const AxisPosition x, const AxisPosition y, const AxisPosition z=0);

    /**
     * @brief Add a delta to the position
     *
     * @param delta is the delta to be added
     * @return a reference to the updated position
     */
    PositionInTissue& operator+=(const PositionDelta delta);

    /**
     * @brief Subtract a delta from a position
     *
     * @param delta is the delta to be subtracted to
     * @return a reference to the updated position
     */
    PositionInTissue& operator-=(const PositionDelta delta);

    /**
     * @brief Save a position in tissue in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & x & y & z;
    }

    /**
     * @brief Load a timed genomic mutation from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed genomic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static PositionInTissue load(ARCHIVE& archive)
    {
        PositionInTissue pos;

        archive & pos.x & pos.y & pos.z;

        return pos;
    }
};

/**
 * @brief Add a delta to a position
 *
 * @param position is the position to which a delta should be added to
 * @param delta is the delta to be added to
 * @return the position reachable from `position` adding `delta`
 */
PositionInTissue operator+(const PositionInTissue& position, const PositionDelta delta);

/**
 * @brief Subtract a delta from a position
 *
 * @param position is the position to which a delta should be subtracted to
 * @param delta is the delta to be subtracted to
 * @return the position reachable from `position` subtracting `delta`
 */
PositionInTissue operator-(const PositionInTissue& position, const PositionDelta delta);

/**
 * @brief Check whether two positions are the same
 *
 * @param p1 is a position
 * @param p2 is a position
 * @return true if and only if the two positions have the same coordinates
 */
bool operator==(const PositionInTissue& p1, const PositionInTissue& p2);

/**
 * @brief Check whether two positions differ
 *
 * @param p1 is a position
 * @param p2 is a position
 * @return true if and only if the two positions have different coordinates
 */
bool operator!=(const PositionInTissue& p1, const PositionInTissue& p2);

/**
 * @brief Compute the Manhattan distance between two positions
 *
 * @param p1 is a position
 * @param p2 is a position
 * @return the Manhattan distance
 */
size_t Manhattan_distance(const PositionInTissue& p1, const PositionInTissue& p2);

/**
 * @brief The class representing cell position
 */
struct Position : public PositionInTissue
{
    Tissue* tissue;         //!< tissue

    /**
     * @brief An empty constructor
     */
    Position();

    /**
     * @brief Construct a new Position object
     *
     * @param tissue is the tissue referred by the position
     * @param x is the x axis position in the tissue
     * @param y is the y axis position in the tissue
     * @param z is the z axis position in the tissue
     */
    Position(Tissue& tissue, const AxisPosition& x, const AxisPosition& y, const AxisPosition& z);

    /**
     * @brief A constructor
     *
     * @param tissue is the tissue referred by the position
     * @param position is the position in the tissue
     */
    Position(Tissue& tissue, const PositionInTissue& position);
};


}   // Evolutions

}   // Mutants

}   // CLONES


namespace std
{

/**
 * @brief Write a position delta in an output stream
 *
 * @param os is the output stream
 * @param delta is the position delta to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::PositionDelta& delta);

/**
 * @brief Write a position in a tissue in an output stream
 *
 * @param os is the output stream
 * @param position is the position in a tissue to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::PositionInTissue& position);

/**
 * @brief Write a position in an output stream
 *
 * @param os is the output stream
 * @param position is the position to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::Position& position);

}   // std

#endif // __CLONES_POSITION__
