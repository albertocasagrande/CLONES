/**
 * @file lineage_graph.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines lineage graphs
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

#ifndef __CLONES_LINEAGE_GRAPH__
#define __CLONES_LINEAGE_GRAPH__

#include <map>

#include "time.hpp"
#include "mutant_id.hpp"

#include "archive.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

/**
 * @brief Lineage of a species
 *
 * This class represents the descendance of a species
 * from another species.
 */
class LineageEdge
{
    SpeciesId ancestor;   //!< Ancestor species
    SpeciesId progeny;  //!< Progeny species
public:
    /**
     * @brief The empty constructor
     */
    LineageEdge();

    /**
     * @brief Construct a new lineage edge
     *
     * @param ancestor is the identifier of the ancestor species
     * @param progeny is the identifier of the progeny species
     */
    LineageEdge(const SpeciesId& ancestor, const SpeciesId& progeny);

    /**
     * @brief Get the ancestor of the lineage edge
     *
     * @return a constant reference to the ancestor identifier
     */
    inline const SpeciesId& get_ancestor() const
    {
        return ancestor;
    }

    /**
     * @brief Get the progeny of the lineage edge
     *
     * @return a constant reference to the progeny identifier
     */
    inline const SpeciesId& get_progeny() const
    {
        return progeny;
    }

    /**
     * @brief Save a lineage edge in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & ancestor & progeny;
    }

    /**
     * @brief Load a lineage edge from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded lineage edge
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static LineageEdge load(ARCHIVE& archive)
    {
        LineageEdge lineage_edge;

        archive & lineage_edge.ancestor
                & lineage_edge.progeny;

        return lineage_edge;
    }
};


/**
 * @brief A class to represent lineage between species
 *
 * This class represents lineage
 */
class LineageGraph
{
    /**
     * @brief A comparator for `LineageEdge`
     */
    struct EdgeCmp
    {
        inline
        bool operator()(const LineageEdge& a, const LineageEdge& b) const
        {
            return (a.get_ancestor()<b.get_ancestor()
                    || (a.get_ancestor()==b.get_ancestor()
                        && a.get_progeny()<b.get_progeny()));
        }
    };

    std::map<LineageEdge, Time, EdgeCmp> first_occurrence; //!< The time of the edge first occurrence

public:
    using const_iterator = std::map<LineageEdge, Time, EdgeCmp>::const_iterator;

    /**
     * @brief The empty constructor
     */
    LineageGraph();

    /**
     * @brief Test whether the graph contains an edge
     *
     * @param ancestor is the lineage edge ancestor
     * @param progeny is the lineage edge progeny
     * @return `true` if and only if the graph contains the edge
     *          `{ancestor, progeny}`.
     */
    inline bool has_edge(const SpeciesId& ancestor, const SpeciesId& progeny) const
    {
        return has_edge({ancestor, progeny});
    }

    /**
     * @brief Test whether the graph contains an edge
     *
     * @param edge is the lineage edge
     * @return `true` if and only if the graph contains the edge
     *          `edge`.
     */
    inline bool has_edge(const LineageEdge& edge) const
    {
        return first_occurrence.count(edge)>0;
    }

    /**
     * @brief Add a lineage edge to the graph and label it
     *
     * @param ancestor is the lineage edge ancestor
     * @param progeny is the lineage edge progeny
     * @param time is the the time of first occurrence for the edge
     * @return a reference to the updated lineage graph
     */
    inline LineageGraph& add_edge(const SpeciesId& ancestor, const SpeciesId& progeny,
                                  const Time& time)
    {
        return add_edge({ancestor, progeny}, time);
    }

    /**
     * @brief Add a lineage edge to the graph and label it
     *
     * @param edge is the lineage edge
     * @param time is the the time of first occurrence for the edge
     * @return a reference to the updated lineage graph
     */
    LineageGraph& add_edge(const LineageEdge& edge, const Time& time);

    /**
     * @brief Get the number of graph edges
     *
     * @return the number of graph edges
     */
    inline size_t num_of_edges() const
    {
        return first_occurrence.size();
    }

    /**
     * @brief Get the initial constant iterator over the edges
     *
     * @return the initial constant iterator
     */
    inline const_iterator begin() const
    {
        return first_occurrence.begin();
    }

    /**
     * @brief Get the final constant iterator over the edges
     *
     * @return the final constant iterator
     */
    inline const_iterator end() const
    {
        return first_occurrence.end();
    }

    /**
     * @brief Save a lineage graph in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & first_occurrence;
    }

    /**
     * @brief Load a lineage graph from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded lineage graph
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static LineageGraph load(ARCHIVE& archive)
    {
        LineageGraph lineage_graph;

        archive & lineage_graph.first_occurrence;

        return lineage_graph;
    }
};

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_LINEAGE_GRAPH__
