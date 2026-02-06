/**
 * @file mutant_mutation.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements mutations
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

#include "mutant_mutation.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

Mutation::Mutation(const MutantId& initial_id, const MutantId& final_id):
    initial_id(initial_id), final_id(final_id)
{}

Mutation::Mutation(const MutantProperties& initial_mutant,
                   const MutantProperties& final_mutant):
    Mutation(initial_mutant.get_id(), final_mutant.get_id())
{}


}   // Evolutions

}   // Mutants

}   // CLONES
