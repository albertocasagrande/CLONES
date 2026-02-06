/**
 * @file lineage_graph.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements lineage graphs
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

#include "lineage_graph.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

LineageEdge::LineageEdge():
    ancestor(WILD_TYPE_SPECIES), progeny(WILD_TYPE_SPECIES)
{}

LineageEdge::LineageEdge(const SpeciesId& ancestor, const SpeciesId& progeny):
    ancestor(ancestor), progeny(progeny)
{}

LineageGraph::LineageGraph()
{}

LineageGraph& LineageGraph::add_edge(const LineageEdge& edge, const Time& time)
{
    if (has_edge(edge)) {
        throw std::domain_error("The edge is already present in the lineage graph.");
    }

    first_occurrence.insert({edge, time});

    return *this;
}

}   // Evolutions

}   // Mutants

}   // CLONES
