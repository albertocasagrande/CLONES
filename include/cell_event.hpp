/**
 * @file cell_event.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines cell events
 * @version 1.2
 * @date 2026-06-10
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

#ifndef __CLONES_CELL_EVENT__
#define __CLONES_CELL_EVENT__

#include <string>
#include <map>

#include "mutant_id.hpp"
#include "position.hpp"
#include "time.hpp"

namespace CLONES
{

namespace Mutants
{

/**
 * @brief Supported event type
 *
 */
enum class CellEventType {
    DEATH,                  //!< death of a cell
    DUPLICATION,            //!< duplication of a cell
    MUTATION,               //!< duplication and mutation of a cell
    EPIGENETIC_SWITCH,      //!< epigenetic switch for a cell
    DUP_AND_EPI_SWITCH,     //!< duplication and epigenetic switch for a cell
    ANY                     //!< unspecified
};

/**
 * @brief A map associating a name to each cell event
 */
extern std::map<CellEventType, std::string> cell_event_names;

namespace Evolutions
{

/**
 * @brief A description of a simulation event
 */
struct CellEvent
{
    CellEventType type;     //!< event type
    Position position;      //!< event position
    SpeciesId src_species;  //!< species of the cell on which event occurs
    SpeciesId dst_species;  //!< destination species for the events
    Time delay;             //!< event delay

    /**
     * @brief The empty constructor
     */
    CellEvent();

    /**
     * @brief Get the event identifier by name
     *
     * @param name is the event name
     * @return the identifier of the event having `name` as the name
     */
    static CellEventType get_event_id(const std::string& name);
};

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_CELL_EVENT__
