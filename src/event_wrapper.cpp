/**
 * @file event_wrapper.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement a simulation event wrapper
 * @version 1.3
 * @date 2026-06-11
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

#include "event_wrapper.hpp"

#include "error.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

TissueSimulationEventWrapper::TissueSimulationEventWrapper(const Mutation& mutation):
    type(TissueSimulationEvent::Type::MUTATION),
    event(std::make_shared<Mutation>(mutation))
{}

TissueSimulationEventWrapper::TissueSimulationEventWrapper(const RateUpdate& rate_update):
    type(TissueSimulationEvent::Type::LIVENESS_RATE_UPDATE),
    event(std::make_shared<RateUpdate>(rate_update))
{}

TissueSimulationEventWrapper::TissueSimulationEventWrapper(const Sampling& sampling):
    type(TissueSimulationEvent::Type::SAMPLING),
    event(std::make_shared<Sampling>(sampling))
{}

}   // Evolutions

}   // Mutants

}   // CLONES


bool operator==(const CLONES::Mutants::Evolutions::TissueSimulationEventWrapper& lhs,
                const CLONES::Mutants::Evolutions::TissueSimulationEventWrapper& rhs)
{
    using namespace CLONES::Mutants::Evolutions;

    if (lhs.type != rhs.type) {
        return false;
    }

    switch(lhs.type) {
        case TissueSimulationEvent::Type::MUTATION:
            {
                const auto& m_lhs = lhs.get_event<Mutation>();
                const auto& m_rhs = rhs.get_event<Mutation>();

                return (m_lhs == m_rhs);
            }
        case TissueSimulationEvent::Type::LIVENESS_RATE_UPDATE:
            {
                const auto& ru_lhs = lhs.get_event<RateUpdate>();
                const auto& ru_rhs = rhs.get_event<RateUpdate>();

                return (ru_lhs == ru_rhs);
            }
        default:
            throw CLONES::Error<std::runtime_error>("Unsupported event type "
                                                   + std::to_string(static_cast<uint>(lhs.type))
                                                   + ".");
    }
}
