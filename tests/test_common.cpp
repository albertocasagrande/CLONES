/**
 * @file test_common.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Functions common to many tests
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

#include "test_common.hpp"

bool operator==(const CLONES::Mutants::Cell& a, const CLONES::Mutants::Cell& b)
{
    return (a.get_id()==b.get_id() &&
            a.get_parent_id()==b.get_parent_id() &&
            a.get_species_id()==b.get_species_id());
}

bool operator==(const CLONES::Mutants::Evolutions::CellInTissue& a,
                const CLONES::Mutants::Evolutions::CellInTissue& b)
{
    using namespace  CLONES::Mutants;

    return (static_cast<const Cell&>(a)==static_cast<const Cell&>(b) &&
            static_cast<const Evolutions::PositionInTissue&>(a)
                == static_cast<const Evolutions::PositionInTissue&>(b));
}

bool operator==(const CLONES::Mutants::Evolutions::Species& a,
                const CLONES::Mutants::Evolutions::Species& b)
{
    using namespace std;
    using namespace CLONES::Mutants;

    if (static_cast<const SpeciesProperties&>(a)
            != static_cast<const SpeciesProperties&>(b)) {
        return false;
    }

    auto cell_a_it = a.begin(), cell_b_it = b.begin();

    for (; cell_a_it != a.end(); ++cell_a_it, ++cell_b_it) {
        if (*cell_a_it != *cell_b_it) {
            return false;
        }
    }

    return true;
}

bool operator==(const CLONES::Mutants::Evolutions::Tissue& a,
                const CLONES::Mutants::Evolutions::Tissue& b)
{
    using namespace CLONES::Mutants;

    if (a.get_name()!=b.get_name()) {
        return false;
    }

    std::set<MutantId> mutant_ids;

    auto a_it = a.begin(), b_it = b.begin();

    for (; a_it != a.end(); ++a_it, ++b_it) {

        if (*a_it != *b_it) {
            return false;
        }
        mutant_ids.insert(a_it->get_mutant_id());
    }

    for (const auto& mutant_id : mutant_ids) {
        auto a_view = a.get_mutant_view(mutant_id);
        auto b_view = b.get_mutant_view(mutant_id);

        auto a_it = a_view.begin(), b_it = b_view.begin();

        for (; a_it != a_view.end(); ++a_it, ++b_it) {

            if (a_it->get_id() != b_it->get_id()) {
                return false;
            }
        }
    }

    return true;
}

bool operator==(const CLONES::Mutants::Evolutions::TissueSimulation& a,
                const CLONES::Mutants::Evolutions::TissueSimulation& b)
{
    return a.get_time()==b.get_time() &&
           a.tissue()==b.tissue() &&
           a.death_activation_level==b.death_activation_level;
}
