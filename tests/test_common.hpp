/**
 * @file test_common.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Header of the functions common to many tests
 * @version 1.0
 * @date 2026-06-19
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

#ifndef __CLONES_TEST_COMMON__
#define __CLONES_TEST_COMMON__

#include <map>
#include <vector>

#include "tissue_simulation.hpp"

template<typename T>
bool operator==(const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size()!=b.size()) {
        return false;
    }

    auto a_it=a.begin(), b_it=b.begin();
    for (; a_it!=a.end(); ++a_it, ++b_it) {
        if (*a_it!=*b_it) {
            return false;
        }
    }

    return true;
}

template<typename KEY, typename T>
bool operator==(const std::map<KEY,T>& a, const std::map<KEY,T>& b)
{
    if (a.size()!=b.size()) {
        return false;
    }

    auto a_it=a.begin(), b_it=b.begin();
    for (; a_it!=a.end(); ++a_it, ++b_it) {
        if (a_it->first!=b_it->first || a_it->second!=b_it->second) {
            return false;
        }
    }

    return true;
}

bool operator==(const CLONES::Mutants::Cell& a, const CLONES::Mutants::Cell& b);

bool operator==(const CLONES::Mutants::Evolutions::CellInTissue& a,
                const CLONES::Mutants::Evolutions::CellInTissue& b);

bool operator==(const CLONES::Mutants::Evolutions::Species& a,
                const CLONES::Mutants::Evolutions::Species& b);

bool operator==(const CLONES::Mutants::Evolutions::Tissue& a,
                const CLONES::Mutants::Evolutions::Tissue& b);

bool operator==(const CLONES::Mutants::Evolutions::TissueSimulation& a,
                const CLONES::Mutants::Evolutions::TissueSimulation& b);

#endif  // __CLONES_TEST_COMMON__