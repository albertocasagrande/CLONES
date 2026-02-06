/**
 * @file simulation_event.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a basic simulation class
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

#ifndef __CLONES_SIMULATION_EVENT__
#define __CLONES_SIMULATION_EVENT__

#include <map>

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

/**
 * @brief A structure to represent basic simulation events
 */
struct SimulationEvent
{
    /**
     * @brief The simulation event types
     */
    enum class Type {
        MUTATION,
        SAMPLING,
        LIVENESS_RATE_UPDATE
    };

    virtual Type type() const = 0;

    virtual ~SimulationEvent();
};

extern const std::map<SimulationEvent::Type, const char*> simulation_event_names;

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_SIMULATION_EVENT__
