/**
 * @file mutant_id.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines mutant ans species identifiers
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

#ifndef __CLONES_MUTANT_ID__
#define __CLONES_MUTANT_ID__

#include <limits>
#include <cstdint>

namespace CLONES
{

/**
 * @brief The namespace of classes handing mutants
 */
namespace Mutants
{

/**
 * @brief The mutant identifier type
 */
using MutantId = uint16_t;

/**
 * @brief The species identifier type
 */
using SpeciesId = uint16_t;


}   // Mutants

}   // CLONES

/**
 * @brief An identifier for wild-type species
 *
 * This macro represents the species identifier of the wild-type cells.
 */
#define WILD_TYPE_SPECIES std::numeric_limits<CLONES::Mutants::SpeciesId>::max()

#endif // __CLONES_MUTANT_ID__