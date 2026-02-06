/**
 * @file ending_conditions.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines simulation ending conditions
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

#ifndef __CLONES_ENDING_CONDITIONS__
#define __CLONES_ENDING_CONDITIONS__

#include "simulation.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

/**
 * @brief Simulation time test
 *
 * The objects of this class have a method testing the time
 * simulated by a simulation. If it is below a threshold,
 * then the test returns `false`; otherwise, it returns `true`.
 */
struct TimeTest : public Simulation::BasicTest
{
    const Time threshold; //!< The time threshold

    /**
     * @brief A constructor
     *
     * @param threshold is the threshold for simulation time test
     */
    explicit TimeTest(const Time& threshold);

    /**
     * @brief Test whether the simulated time is not below the test threshold
     *
     * @param simulation is the considered simulation
     * @return `false` if and only if the time simulated by `simulation`
     *          is lower than the test threshold
     */
    inline bool operator()(const Simulation& simulation) override
    {
        return threshold < simulation.get_time();
    }

    /**
     * @brief Return the percentage of the completed simulation
     *
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    inline uint8_t percentage(const Simulation& simulation) override
    {
        return static_cast<uint8_t>((100*simulation.get_time()/threshold));
    }
};

/**
 * @brief Species count test
 *
 * The objects of this class have a method testing the number of
 * cells in a species. If this number is below a threshold,
 * then the test returns `false`; otherwise, it returns `true`.
 */
struct SpeciesCountTest : public Simulation::BasicTest
{
    const SpeciesId species_id; //!< The species id
    const size_t threshold;  //!< The size threshold

    /**
     * @brief A constructor
     *
     * @param species_id is the identifier of the species whose
     *          number of cells is counted
     * @param threshold is the threshold for the count test
     */
    SpeciesCountTest(const SpeciesId& species_id, const size_t& threshold);

    /**
     * @brief Test whether the number of cells is below the threshold
     *
     * @param simulation is the considered simulation
     * @return `false` if and only if the number of cells of the
     *          considered species is below the test threshold
     */
    bool operator()(const Simulation& simulation) override;

    /**
     * @brief Return the percentage of the completed simulation
     *
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) override;
};

/**
 * @brief Clone count test
 *
 * The objects of this class have a method testing the number of
 * cells of a mutant. If this number is below a threshold,
 * then the test returns `false`; otherwise, it returns `true`.
 */
struct CloneCountTest : public Simulation::BasicTest
{
    const MutantId mutant_id;  //!< The mutant id
    const size_t threshold;  //!< The size threshold

    /**
     * @brief A constructor
     *
     * @param mutant_id is the identifier of the mutant whose
     *          number of cells is counted
     * @param threshold is the threshold for the count test
     */
    CloneCountTest(const MutantId& mutant_id, const size_t& threshold);

    /**
     * @brief Test whether the number of cells is below the threshold
     *
     * @param simulation is the considered simulation
     * @return `false` if and only if the number of cells of the
     *          considered mutant is below the test threshold
     */
    bool operator()(const Simulation& simulation) override;

    /**
     * @brief Return the percentage of the completed simulation
     *
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) override;
};

/**
 * @brief Event count test
 *
 * The objects of this class have a method testing the number of
 * cells stochastic events occurred to a species. If this number
 * is below a threshold, then the test returns `false`; otherwise,
 * it returns `true`.
 */
class EventCountTest : public Simulation::BasicTest
{
    /**
     * @brief Get the number of events tested by this object
     *
     * @param simulation is the considered simulation
     * @return the number of events tested by this object
     */
    size_t get_event_number(const Simulation& simulation) const;

    const CellEventType event_type; //!< The cell event
    const SpeciesId species_id;     //!< The id of the species in which the events occurred
    const SpeciesId dst_id;         //!< The id of the destination species in case of epigenetic events
    const size_t threshold;         //!< The threshold on the number of events

public:

    /**
     * @brief A constructor
     *
     * Whenever `event_type` is `CellEventType::EPIGENETIC_SWITCH`
     * the built object compare the threshold and the sum of all
     * the epigenetic events.
     *
     * @param event_type is the aimed event type
     * @param species_id is the identifier of the species whose
     *          number of events is counted
     * @param threshold is the threshold for the number of events
     */
    EventCountTest(const CellEventType& event_type, const SpeciesId& species_id,
                   const size_t& threshold);

    /**
     * @brief A constructor
     *
     * This constructor builds an object specifically designed to
     * test whether the number of epigenetic event from `src_id`
     * to `dst_id` has reached a threshold.
     *
     * @param src_id is the identifier of the source species of
     *          counted epigenetic events
     * @param dst_id is the identifier of the destination species of
     *          counted epigenetic events
     * @param threshold is the threshold for the number of events
     */
    EventCountTest(const SpeciesId& src_id, const SpeciesId& dst_id,
                   const size_t& threshold);

    /**
     * @brief Test whether the number of cells is below the threshold
     *
     * @param simulation is the considered simulation
     * @return `false` if and only if the number of cells of the
     *          considered species is below the test threshold
     */
    bool operator()(const Simulation& simulation) override;

    /**
     * @brief Return the percentage of the completed simulation
     *
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) override;
};

class FormulaTest : public Simulation::BasicTest
{
    Logics::Formula formula;
    size_t formula_distance;

public:

    /**
     * @brief The constructor
     *
     * @param formula is the formula to be satisfied
     */
    FormulaTest(const Logics::Formula& formula);

    /**
     * @brief Test whether the formula has been satisfied
     *
     * @param simulation is the considered simulation
     * @return `false` if and only if the simulation does not
     *      satisfies the formula of the test
     */
    bool operator()(const Simulation& simulation) override;

    /**
     * @brief Return the percentage of the completed simulation
     *
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) override;
};

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_ENDING_CONDITIONS__
