/**
 * @file cell_event.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements cell events
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

#include <string>
#include <map>

#include "cell_event.hpp"

namespace CLONES
{

namespace Mutants
{

std::map<CellEventType, std::string> cell_event_names = {
    {CellEventType::DEATH, "deaths"},
    {CellEventType::DUPLICATION, "duplications"},
    {CellEventType::EPIGENETIC_SWITCH, "switches"},
    {CellEventType::MUTATION, "mutations"}
};

namespace Evolutions
{

CellEvent::CellEvent():
    type(), position(), initial_species(), final_species(), delay()
{}

CellEventType CellEvent::get_event_id(const std::string& event_name)
{
  for (const auto& [e_id, e_name]: CLONES::Mutants::cell_event_names) {
    if (e_name == event_name) {
      return e_id;
    }
  }

  throw std::domain_error("Unknown event name \"" + event_name + "\".");
}

}

}   // Mutants

}   // CLONES
