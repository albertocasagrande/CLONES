/**
 * @file sample_specification.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements sample specification
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

#include <sstream>

#include "sample_specification.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

std::string SampleSpecification::get_default_name(const RectangleSet& region,
                                                  const size_t num_of_cells)
{
    std::ostringstream oss;

    oss << "S_" << region.lower_corner << "-" << region.upper_corner << "_" << num_of_cells;

    return oss.str();
}

SampleSpecification::SampleSpecification(const RectangleSet& bounding_box):
    name(), bbox(bounding_box), num_of_cells(bounding_box.size())
{
    name = get_default_name(bounding_box, num_of_cells);
}

SampleSpecification::SampleSpecification(const std::string& name, const RectangleSet& bounding_box):
    name(name), bbox(bounding_box), num_of_cells(bounding_box.size())
{}


SampleSpecification::SampleSpecification(const RectangleSet& bounding_box, const size_t num_of_cells):
    name(), bbox(bounding_box), num_of_cells(num_of_cells)
{
    name = get_default_name(bounding_box, num_of_cells);
}

SampleSpecification::SampleSpecification(const std::string& name, const RectangleSet& bounding_box,
                                         const size_t num_of_cells):
    name(name), bbox(bounding_box), num_of_cells(num_of_cells)
{}

}   // Evolutions

}   // Mutants

}   // CLONES
