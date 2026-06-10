/**
 * @file rate_update.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines liveness rate updates
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

#ifndef __CLONES_RATE_UPDATE__
#define __CLONES_RATE_UPDATE__

#include "mutant_properties.hpp"
#include "cell_event.hpp"
#include "simulation_event.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

/**
 * @brief A structure to represent liveness rate update
 */
struct RateUpdate : public TissueSimulationEvent
{
    using Type = TissueSimulationEvent::Type;

    SpeciesId src_id;           //!< The event source species id
    SpeciesId dst_id;           //!< The event destination species id
    CellEventType event_type;   //!< The event type
    double new_rate;            //!< The new rate

    /**
     * @brief A constructor
     *
     * @param species_id is the identifier of the species in which the event occurs
     * @param event_type is the event type whose rate is changed. The event type
     *      must be either `CellEventType::DEATH` or `CellEventType::DUPLICATION`
     * @param new_rate is the new rate for `event_type` in `species_id`
     */
    RateUpdate(const SpeciesId& species_id, const CellEventType& event_type,
               const double& new_rate);

    /**
     * @brief A constructor
     *
     * @param species is the species in which the event occurs
     * @param event_type is the event type whose rate is changed by the event. The
     *      event type must be either `CellEventType::DEATH` or
     *      `CellEventType::DUPLICATION`
     * @param new_rate is the new rate for `event_type` in `species`
     */
    RateUpdate(const SpeciesProperties& species,
               const CellEventType& event_type, const double& new_rate);

    /**
     * @brief A constructor
     *
     * @param src_id is the identifier of the species from which the event occurs
     * @param dst_id is the identifier of the species to which the event occurs
     * @param event_type is the event type whose rate is changed
     * @param new_rate is the new rate for `event_type` from `src_id` to `dst_id`
     */
    RateUpdate(const SpeciesId& src_id, const SpeciesId& dst_id,
               const CellEventType& event_type, const double& new_rate);

    /**
     * @brief A constructor
     *
     * @param src is the species from which the event occurs.
     * @param dst is the species to which the event occurs.
     * @param event_type is the event type whose rate is changed by the event.
     * @param new_rate is the new rate for `event_type`  from `src` to `dst`.
     */
    RateUpdate(const SpeciesProperties& src, const SpeciesProperties& dst,
               const CellEventType& event_type, const double& new_rate);

    /**
     * @brief Save a timed genomic mutation in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & src_id
                & dst_id
                & event_type
                & new_rate;
    }

    /**
     * @brief Load a timed genomic mutation from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed genomic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static RateUpdate load(ARCHIVE& archive)
    {
        SpeciesId src_id, dst_id;
        CellEventType event_type;
        double new_rate;

        archive & src_id
                & dst_id
                & event_type
                & new_rate;

        return {src_id, dst_id, event_type, new_rate};
    }

    inline Type type() const {
        return Type::LIVENESS_RATE_UPDATE;
    }
};

}   // Evolutions

}   // Mutants

}   // CLONES


/**
 * @brief Test the equivalence between two liveness rate updates
 *
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two liveness rate updates represent
 *      the same event
 */
inline bool operator==(const CLONES::Mutants::Evolutions::RateUpdate& lhs,
                       const CLONES::Mutants::Evolutions::RateUpdate& rhs)
{
    return (lhs.src_id == rhs.src_id)
            && (lhs.dst_id == rhs.dst_id)
            && (lhs.event_type == rhs.event_type)
            && (lhs.new_rate == rhs.new_rate);
}

#endif // __CLONES_RATE_UPDATE__
