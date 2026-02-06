/**
 * @file sampling.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements samplings
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

#include "sampling.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

Sampling::Sampling(const SampleSpecification& orig):
    SampleSpecification(orig)
{}

Sampling::Sampling(const RectangleSet& bounding_box):
    SampleSpecification(bounding_box)
{}

Sampling::Sampling(const RectangleSet& bounding_box,
                   const size_t& num_of_cells):
    SampleSpecification(bounding_box, num_of_cells)
{}

Sampling::Sampling(const std::string& name, const RectangleSet& bounding_box):
    SampleSpecification(name, bounding_box)
{}

Sampling::Sampling(const std::string& name, const RectangleSet& bounding_box,
                   const size_t& num_of_cells):
    SampleSpecification(name, bounding_box, num_of_cells)
{}

}   // Evolutions

}   // Mutants

}   // CLONES
