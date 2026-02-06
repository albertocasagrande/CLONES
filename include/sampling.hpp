/**
 * @file sampling.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines tissue samplings
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

#ifndef __CLONES_SAMPLING__
#define __CLONES_SAMPLING__

#include "position_set.hpp"
#include "simulation_event.hpp"
#include "sample_specification.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

/**
 * @brief A structure to represent tissue sampling event
 */
struct Sampling : public SampleSpecification, public SimulationEvent
{
    using Type = SimulationEvent::Type;

    /**
     * @brief A constructor
     *
     * @param orig is a template sample specification for the new sampling event
     */
    explicit Sampling(const SampleSpecification& orig);

    /**
     * @brief A constructor
     *
     * @param bounding_box is the sample bounding box
     * @param num_of_cells is the maximum number of cells to sample
     */
    Sampling(const RectangleSet& bounding_box, const size_t& num_of_cells);

    /**
     * @brief A constructor
     *
     * @param bounding_box is the sample bounding box
     */
    explicit Sampling(const RectangleSet& bounding_box);

    /**
     * @brief A constructor
     *
     * @param name is the specification name
     * @param bounding_box is the sample bounding box
     * @param num_of_cells is the maximum number of cells to sample
     */
    Sampling(const std::string& name, const RectangleSet& bounding_box,
             const size_t& num_of_cells);

    /**
     * @brief A constructor
     *
     * @param name is the specification name
     * @param bounding_box is the sample bounding box
     */
    Sampling(const std::string& name, const RectangleSet& bounding_box);

    inline Type type() const {
        return Type::SAMPLING;
    }
};

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_RATE_UPDATE__
