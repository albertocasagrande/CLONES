/**
 * @file tissue_simulation.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a tumour evolution simulation
 * @version 1.14
 * @date 2026-07-23
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

#ifndef __CLONES_SIMULATOR__
#define __CLONES_SIMULATOR__

#include <random>
#include <list>
#include <queue>
#include <vector>
#include <chrono>

#include "tissue.hpp"
#include "binary_logger.hpp"
#include "plot_2D.hpp"
#include "tissue_plotter.hpp"
#include "statistics.hpp"
#include "timed_event.hpp"
#include "tissue_sample.hpp"

#include "lineage_graph.hpp"

#include "descendant_forest.hpp"

#include "progress_bar.hpp"

namespace CLONES
{

namespace Mutants
{

/**
 * @brief The namespace of mutants evolutions
 */
namespace Evolutions
{

/**
 * @brief A tumour evolution simulation
 */
class TissueSimulation
{
public:

    /**
     * @brief A class to represent simulation status at snapshot
     */
    struct StatusAtSnapshot
    {
        using Clock = std::chrono::system_clock;
        using TimePoint = std::chrono::time_point<Clock>;
        using Duration = TimePoint::duration;

        TimePoint time;         //!< The time at snapshot
        Time clock;             //!< The simulated time at snapshot
        uint64_t num_of_cells;  //!< The number of cells at snapshot

        /**
         * @brief The empty constructor
         */
        StatusAtSnapshot();

        /**
         * @brief Construct a new Status At Snapshot object
         *
         * @param simulation is the tissue simulation whose snapshot
         *   has been taken
         */
        explicit StatusAtSnapshot(const TissueSimulation& simulation);

        /**
         * @brief Get the simulated time at the snapshot
         *
         * @return the simulated time at the snapshot
         */
        inline const Time& get_clock() const
        {
            return clock;
        }

        /**
         * @brief Set the time
         *
         * @return a constant reference to the updated object
         */
        inline const StatusAtSnapshot& set_time()
        {
            time = Clock::now();

            return *this;
        }

        /**
         * @brief Get the snapshot time
         *
         * @return the snapshot time
         */
        inline const TimePoint& get_time() const
        {
            return time;
        }

        /**
         * @brief Get the number of cells at snapshot
         *
         * @return the number of cells at snapshot
         */
        inline const uint64_t& get_num_of_cells() const
        {
            return num_of_cells;
        }

        /**
         * @brief Save an object of the class `StatusAtSnapshot` in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE>
          requires std::is_base_of_v<Archive::Basic::Out, ARCHIVE>
        inline void save(ARCHIVE& archive) const
        {
            archive & time.time_since_epoch().count()
                    & clock
                    & num_of_cells;
        }

        /**
         * @brief Load an object of the class `StatusAtSnapshot` from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the loaded object
         */
        template<typename ARCHIVE>
          requires std::is_base_of_v<Archive::Basic::In, ARCHIVE>
        inline static StatusAtSnapshot load(ARCHIVE& archive)
        {
            StatusAtSnapshot status;

            TimePoint::duration::rep ticks;

            archive & ticks
                    & status.clock
                    & status.num_of_cells;

            status.time = TimePoint(TimePoint::duration(ticks));

            return status;
        }
    };

    /**
     * @brief A snapshot trigger type
     *
     * The object of this class establish whether a snapshot must
     * be taken.
     */
    struct SnapshotTrigger
    {
        using Duration = StatusAtSnapshot::Duration;

        Duration time_delta;            //!< The maximum time difference between two successive snapshots
        Time clock_delta;               //!< The maximum simulated time difference between two successive snapshots
        uint64_t cardinality_delta;     //!< The maximum difference in the number of cells between two successive snapshots

        /**
         * @brief The empty constructor
         */
        SnapshotTrigger();

        /**
         * @brief Set the time trigger
         *
         * @param time_delta is the new maximum time interval between two
         *   consecutive snapshots.
         * @return A reference to the updated object.
         */
        inline SnapshotTrigger& set_time_trigger(const Duration& time_delta)
        {
            this->time_delta = time_delta;

            return *this;
        }

        /**
         * @brief Get the time trigger
         *
         * @return the maximum time interval between two
         *   consecutive snapshots.
         */
        inline const Duration& get_time_trigger() const
        {
            return this->time_delta;
        }

        /**
         * @brief Set the clock trigger
         *
         * @param clock_delta is the new maximum simulated time interval
         *   between two consecutive snapshots.
         * @return A reference to the updated object.
         */
        inline SnapshotTrigger& set_clock_trigger(const Time& clock_delta)
        {
            this->clock_delta = clock_delta;

            return *this;
        }

        /**
         * @brief Get the clock trigger
         *
         * @return the maximum simulated time interval
         *   between two consecutive snapshots.
         */
        inline const Time& get_clock_trigger() const
        {
            return clock_delta;
        }

        /**
         * @brief Set the cardinality trigger
         *
         * @param cardinality_delta is the new maximum difference in the number of cells
         *   between two consecutive snapshots.
         * @return A reference to the updated object.
         */
        inline SnapshotTrigger& set_cardinality_trigger(const uint64_t& cardinality_delta)
        {
            this->cardinality_delta = cardinality_delta;

            return *this;
        }

        /**
         * @brief Get the cardinality trigger
         *
         * @return The maximum difference in the number of cells
         *   between two consecutive snapshots.
         */
        inline const uint64_t& get_cardinality_trigger() const
        {
            return cardinality_delta;
        }

        /**
         * @brief Test whether the simulation triggers a new snapshot
         *
         * @param simulation is the tissue simulation to be tested.
         * @return `true` if and only if the current simulation triggers
         *   a new snapshot.
         */
        bool is_triggered_by(const TissueSimulation& simulation) const;

        /**
         * @brief Save an object of the class `SnapshotTrigger` in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE>
          requires std::is_base_of_v<Archive::Basic::Out, ARCHIVE>
        inline void save(ARCHIVE& archive) const
        {
            archive & time_delta.count()
                    & clock_delta
                    & cardinality_delta;
        }

        /**
         * @brief Load an object of the class `SnapshotTrigger` from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the loaded object
         */
        template<typename ARCHIVE>
          requires std::is_base_of_v<Archive::Basic::In, ARCHIVE>
        inline static SnapshotTrigger load(ARCHIVE& archive)
        {
            SnapshotTrigger trigger;

            Duration::rep ticks;

            archive & ticks
                    & trigger.clock_delta
                    & trigger.cardinality_delta;

            trigger.time_delta = Duration(ticks);

            return trigger;
        }
    };

    /**
     * @brief Information about the cells manually added to the simulation
     */
    struct AddedCell : public PositionInTissue
    {
        SpeciesId species_id;       //!< The species of the added cell
        Time time;                  //!< The adding time

        /**
         * @brief The empty constructor
         */
        AddedCell();

        /**
         * @brief Construct a new `AddedCell` object
         *
         * @param species_id is the identifier of the added cell
         * @param position is the position of the added cell
         * @param time is the adding time
         */
        AddedCell(const SpeciesId& species_id, const PositionInTissue& position, const Time& time);

        /**
         * @brief Save an object of the class `AddedCell` in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & *(static_cast<const PositionInTissue*>(this))
                    & species_id
                    & time;
        }

        /**
         * @brief Load an object of the class `AddedCell` from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the loaded object
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
        inline static AddedCell load(ARCHIVE& archive)
        {
            AddedCell added_cell;

            archive & static_cast<PositionInTissue&>(added_cell)
                    & added_cell.species_id
                    & added_cell.time;

            return added_cell;
        }
    };

protected:
    struct EventAffectedCells
    {
        std::list<CellInTissue> new_cells;
        std::list<Cell> lost_cells;
    };

    using TimedEventQueue = std::priority_queue<TimedEvent,
                                                std::vector<TimedEvent>,
                                                std::greater<TimedEvent>>;
    using system_clock = std::chrono::system_clock;

    std::vector<Tissue> tissues;     //!< Simulated tissues
    std::map<std::string, MutantId> mutant_name2id; //!< A map associating the mutant name to its id

    LineageGraph lineage_graph;     //!< The lineage graph of the simulation

    BinaryLogger logger;             //!< Event logger

    std::vector<Direction> valid_directions;   //!< valid simulation tissue directions

    StatusAtSnapshot status_at_snapshot;    //!< The simulation status at the last snapshot
    SnapshotTrigger snapshot_trigger;       //!< The object testing whether a snapshot is needed

    TissueStatistics statistics;     //!< The tissue simulation statistics

    Time time;                       //!< TissueSimulation time
    std::mt19937_64 random_gen;      //!< Pseudo-random generator

    TimedEventQueue timed_event_queue;   //!< The timed event queue

    std::list<AddedCell> added_cells;   //!< The list of the manually added cells

    std::list<TissueSample> samples;    //!< The list of the tissue samples

    using SampleConstIt = std::list<TissueSample>::const_iterator;

    std::map<std::string, SampleConstIt> name2sample; //!< A map associating each sample name to the sample

    CellId next_cell_id;    //!< The next cell id

    /**
     * @brief Simulate a cell duplication
     *
     * This method simulates the duplication of a cell in the
     * tissue. The duplicating cell is identified by its position.
     * One of the two sibling cells is placed in the parent cell
     * position; the other one is randomly placed near to the
     * former parent. This is done by selecting one direction among
     * the 6 possible directions (i.e., left or right on x-axis,
     * y-axis, or z-axis) and by pushing all the cells on that
     * direction one position await from the parent cell.
     * If the cell in the provided position has wild-type
     * genotype, then nothing is done.
     *
     * @param position is the position of the cell to be duplicate
     * @return list of the affected cells. Whenever the specified position
     *      is not in the tissue or does correspond to a wild-type cell,
     *      the returned list is empty. If one of the two children is
     *      pushed outside the tissue borders, the returned list contains
     *      only the new cell in position `position`. In the remaining
     *      case, the returned list contains two cells
     */
    EventAffectedCells simulate_duplication(const Position& position);

    /**
     * @brief Simulate the death of a cell
     *
     * This method simulates the death of a cell in tissue.
     * If the cell in the provided position has wild-type
     * genotype, then nothing is done.
     *
     * @param position is the position of the cell to be killed
     * @return list of the affected cells, i.e., a list containing
     *      the former status of the killed cell at most
     */
    EventAffectedCells simulate_death(const Position& position);

    /**
     * @brief Simulate a cell mutation
     *
     * This method simulates a mutation on a cell in tissue.
     * If the cell in the provided position has wild-type genotype,
     * then nothing is done.
     *
     * @param position is the position of the cell that will mutate
     * @param final_id is the resulting species identifier of the cell
     * @return list of the affected cells, i.e., a list containing
     *      the status of the mutated cell at most
     */
    EventAffectedCells simulate_mutation(const Position& position, const SpeciesId& final_id);

    /**
     * @brief Set the species of a tumour cell
     *
     * @param position is the position of a tumour cell
     * @param species_id is the new species identifier
     * @return `true` if and only if the cell in `position` is a tumour cell
     *      and its species is different from `species_id`
     */
    bool set_cell_species(const Position& position, const SpeciesId& species_id);

    /**
     * @brief Simulate a duplication and a mutation event
     *
     * This method simulates the duplication of a cell in the
     * tissue and change the species of the child in the provided position.
     *
     * @param position is the position of the cell to be duplicate
     * @param final_id is the resulting species identifier of the cell
     * @return The list of the affected cells. Whenever the specified position
     *      is not in the tissue or does correspond to a wild-type cell,
     *      the returned list is empty. If one of the two children is
     *      pushed outside the tissue borders, the returned list contains
     *      only the new cell in position `position`. In the remaining
     *      case, the returned list contains two cells
     */
    EventAffectedCells simulate_duplication_and_set_one_child_species(const Position& position,
                                                                      const SpeciesId& final_id);

    /**
     * @brief Record tissue initial cells in log
     */
    void log_initial_cells();

    /**
     * @brief Initialize the valid simulation tissue direction vector
     *
     * This method must be called every time the simulation tissue
     * changes.
     */
    void init_valid_directions();

    /**
     * @brief Handle a mutation
     *
     * This method tries to handle a mutation during the cell event selection and
     * it succeeds if and only if there exists at least one cell in the origin.
     *
     * @param timed_mutation is the timed mutation to be applied
     * @param candidate_event is the current candidate cell event
     * @return `true` if and only if there exists at least one cell in the origin and
     *          the candidate cell event has been updated
     */
    bool handle_timed_mutation(const TimedEvent& timed_mutation, CellEvent& candidate_event);

    /**
     * @brief Apply a rate update
     *
     * @param timed_rate_update is the timed rate update to be applied
     */
    void handle_timed_rate_update(const TimedEvent& timed_rate_update);

    /**
     * @brief Sample the tissue
     *
     * This method also update next candidate event.
     *
     * @param time_sampling is the timed sampling to be applied
     * @param candidate_event is the current candidate cell event
     */
    void handle_timed_sampling(const TimedEvent& time_sampling, CellEvent& candidate_event);

    /**
     * @brief Handle time event queue during cell event selection
     *
     * This method handles the time event queue during the cell event selection
     * by extracting the timed events that occurs before the candidate cell
     * event and applying them. If necessary (for instance when the next time
     * event is a mutation), the method also updates the candidate event.
     *
     * @param candidate_event is a candidate cell event
     */
    void handle_timed_event_queue(CellEvent& candidate_event);

    /**
     * @brief Select a cell event among those due to cell liveness
     *
     * @return a cell event among those due to cell liveness
     */
    CellEvent select_next_cell_event();

    /**
     * @brief Collect proxies of cells in a rectangle
     *
     * @param rectangle is the rectangle of tissue in which the proxies
     *      should be collected
     * @return a vector of constant proxies referencing the non-wild cells
     *      in the tissue region bounded by the rectangle
     */
    std::vector<Tissue::CellInTissueConstantProxy>
    collect_cell_proxies_in(const CLONES::Mutants::RectangleSet& rectangle) const;

    /**
     * @brief Collect proxies of cells in a rectangle
     *
     * @param rectangle is the rectangle of tissue in which the proxies
     *      should be collected
     * @return a vector of proxies referencing the non-wild cells in the
     *      tissue region bounded by the rectangle
     */
    std::vector<Tissue::CellInTissueProxy>
    collect_cell_proxies_in(const CLONES::Mutants::RectangleSet& rectangle);
public:
    /**
     * @brief TissueSimulation test
     */
    struct BasicTest
    {
        /**
         * @brief Test a simulation
         *
         * @param simulation is the considered simulation
         * @return a Boolean value
         */
        virtual bool operator()(const TissueSimulation& simulation) = 0;

        /**
         * @brief Return the percentage of the completed simulation
         *
         * @param simulation is the considered simulation
         * @return the percentage of the completed simulation
         */
        virtual uint8_t percentage(const TissueSimulation& simulation) = 0;

        virtual ~BasicTest()
        {}
    };

    size_t death_activation_level;  //!< The minimum number of cells required to activate death

    bool duplicate_internal_cells; //!< A flag to enable/disable duplication in internal cells
    bool storage_enabled;          //!< A flag to enable/disable storage

    /**
     * @brief The basic simulation constructor
     *
     * @param random_seed is the simulation random seed
     */
    explicit TissueSimulation(int random_seed=0);

    /**
     * @brief A simulation constructor
     *
     * @param log_directory is the simulation log directory
     * @param random_seed is the simulation random seed
     */
    explicit TissueSimulation(const std::filesystem::path& log_directory, int random_seed=0);

    /**
     * @brief A swap constructor
     *
     * @param orig is the original simulation
     */
    TissueSimulation(TissueSimulation&& orig);

    /**
     * @brief A copy operator
     *
     * @param orig is the original simulation
     * @return A reference of the updated object
     */
    TissueSimulation& operator=(TissueSimulation&& orig);

    /**
     * @brief Schedule a timed mutation
     *
     * @param src is the source mutant
     * @param dst is the destination mutant
     * @param time is the mutation timing
     * @return a reference to the updated simulation
     */
    inline TissueSimulation& schedule_mutation(const Tissue::const_mutant_view& src,
                                               const Tissue::const_mutant_view& dst,
                                               const Time time)
    {
        return schedule_mutation(src.get_mutant(), dst.get_mutant(), time);
    }

    /**
     * @brief Schedule a timed mutation
     *
     * @param src is the source mutant
     * @param dst is the destination mutant
     * @param time is the mutation timing
     * @return a reference to the updated simulation
     */
    TissueSimulation& schedule_mutation(const MutantProperties& src,
                                        const MutantProperties& dst, const Time time);

    /**
     * @brief Schedule a timed mutation
     *
     * @param src is the source mutant name
     * @param dst is the destination mutant name
     * @param time is the mutation timing
     * @return a reference to the updated simulation
     */
    TissueSimulation& schedule_mutation(const std::string& src, const std::string& dst, const Time time);

    /**
     * @brief Schedule a timed event
     *
     * @param timed_event is the timed event to schedule
     * @return a reference to the updated simulation
     */
    TissueSimulation& schedule_timed_event(const TimedEvent& timed_event);

    /**
     * @brief Select the next event
     *
     * This method select the next event by using the Gillespie's
     * first reaction method as detailed at page 42 of:
     *   Gillespie DT. Stochastic simulation of chemical kinetics.
     *   Annu Rev Phys Chem. 2007;58:35-55.
     *   doi: 10.1146/annurev.physchem.58.032806.104637.PMID: 17037977.
     *
     * @return a cell event
     */
    CellEvent select_next_event();

    /**
     * @brief  Simulate a mutation on a position
     *
     * This method simulates both the duplication of the cell in the
     * specified position and the birth of a cells of a different
     * mutant preserving the epigenetic status of the original cell.
     *
     * @param position is the position in which the
     * @param dst_mutant_name is the name of the mutated cell mutant
     * @return a reference to the updated simulation
     */
    TissueSimulation& simulate_mutation(const PositionInTissue& position,
                                        const std::string& dst_mutant_name);

    /**
     * @brief  Simulate a mutation on a position
     *
     * This method simulates both the duplication of the cell in the
     * specified position and the birth of a cells of a different
     * mutant preserving the epigenetic status of the original cell.
     *
     * @param position is the position in which the
     * @param dst_mutant_id is the identifier of the mutated cell mutant
     * @return a reference to the updated simulation
     */
    TissueSimulation& simulate_mutation(const PositionInTissue& position,
                                        const MutantId& dst_mutant_id);

    /**
     * @brief Simulate an event
     *
     * This method simulates an event on the tissue. If `plotter` is not
     * `nullptr`, then the simulation is also plotted in a graphical window.
     *
     * @tparam PLOT_WINDOW is the plotting window type
     * @param event is the event to be simulated
     * @param plotter is a tissue plotter pointer
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    TissueSimulation& simulate(const CellEvent& event, UI::TissuePlotter<PLOT_WINDOW>* plotter);

    /**
     * @brief Simulate an event
     *
     * This method simulates an event on the tissue. If `plotter` is not
     * `nullptr`, then the simulation is also plotted in a graphical window.
     *
     * @tparam PLOT_WINDOW is the plotting window type
     * @param event is the event to be simulated
     * @param plotter is a tissue plotter pointer
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    inline TissueSimulation& simulate(const CellEvent& event, UI::TissuePlotter<PLOT_WINDOW>& plotter)
    {
        return simulate(event, &plotter);
    }

    /**
     * @brief Simulate an event
     *
     * @param event is the event to be simulated
     * @return a reference to the updated simulation
     */
    inline TissueSimulation& simulate(const CellEvent& event)
    {
        return simulate<UI::Plot2DWindow>(event, nullptr);
    }

    /**
     * @brief Simulate up to the next event
     *
     * This method simulates a tissue up to the next event. If `plotter` is not
     * `nullptr`, then the simulation is also plotted in a graphical window.
     *
     * @tparam PLOT_WINDOW is the plotting window type
     * @param plotter is a tissue plotter pointer
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    TissueSimulation& run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>* plotter);

    /**
     * @brief Simulate up to the next event
     *
     * This method simulates a tissue up to the next event. The simulation is also
     * plotted in a graphical window.
     *
     * @tparam PLOT_WINDOW is the plotting window type
     * @param plotter is a tissue plotter
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    inline TissueSimulation& run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>& plotter)
    {
        return run_up_to_next_event(&plotter);
    }

    /**
     * @brief Simulate up to the next event
     *
     * This method simulates a tissue up to the next
     * event.
     *
     * @return a reference to the updated simulation
     */
    inline TissueSimulation& run_up_to_next_event()
    {
        return run_up_to_next_event<UI::Plot2DWindow>(nullptr);
    }

    /**
     * @brief Simulate a tissue
     *
     * This method simulates a tissue until a test on the simulation
     * object holds. If the user provide a pointer to a plotter, then
     * the simulation is also plotted in a graphical window.
     *
     * @tparam SIMULATION_TEST is a class whose objects can test whether
     *              the simulation is concluded
     * @tparam PLOT_WINDOW is the plotting window type
     * @tparam INDICATOR_TYPE is the progress indicator type
     * @param done is a test verifying whether the simulation has come to an end
     * @param plotter is a tissue plotter pointer
     * @param indicator is the progress indicator pointer
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST, typename PLOT_WINDOW, typename INDICATOR_TYPE,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST> &&
                              std::is_base_of_v<UI::Plot2DWindow, PLOT_WINDOW>, bool> = true>
    TissueSimulation& run(SIMULATION_TEST& done, UI::TissuePlotter<PLOT_WINDOW>* plotter,
                    INDICATOR_TYPE* indicator)
    {
        // the tissue() call checks whether a tissue has been
        // associated to the simulation and, if this is not the
        // case, it throws an std::runtime_error
        (void)tissue();

        while ((plotter == nullptr || !plotter->closed()) && !done(*this)
            && tissue().num_of_mutated_cells()>0) {

            run_up_to_next_event(plotter);

            if (storage_enabled) {
                if (snapshot_trigger.is_triggered_by(*this)) {
                    make_snapshot(indicator);
                }
            }

            if (indicator != nullptr) {
                indicator->set_progress(done.percentage(*this),
                                        "Cells: " + std::to_string(tissue().num_of_mutated_cells()));
            }
        }

        statistics.store_current_in_history(time);
        make_snapshot(indicator);

        if (indicator != nullptr) {
            indicator->set_message("Cells: " + std::to_string(tissue().num_of_mutated_cells()));
        }

        return *this;
    }

    /**
     * @brief Simulate a tissue
     *
     * This method simulates a tissue until a test on the simulation
     * object holds. The simulation is also plotted in a graphical window.
     *
     * @tparam SIMULATION_TEST is a class whose objects can test whether
     *              the simulation is concluded
     * @tparam PLOT_WINDOW is the plotting window type
     * @tparam INDICATOR_TYPE is the progress indicator type
     * @param done is a test verifying whether the simulation has come to an end
     * @param plotter is a tissue plotter
     * @param indicator is the progress indicator pointer
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST, typename PLOT_WINDOW, typename INDICATOR_TYPE,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST> &&
                              std::is_base_of_v<UI::Plot2DWindow, PLOT_WINDOW>, bool> = true>
    inline TissueSimulation& run(SIMULATION_TEST& done, UI::TissuePlotter<PLOT_WINDOW>& plotter,
                                 INDICATOR_TYPE& indicator)
    {
        return run(done, &plotter, &indicator);
    }

    /**
     * @brief Simulate a tissue
     *
     * This method simulates a tissue until a test on the simulation holds.
     *
     * @tparam SIMULATION_TEST is a class whose objects can test whether
     *              the simulation is concluded
     * @tparam INDICATOR_TYPE is the progress indicator type
     * @param done is a test verifying whether the simulation has come to an end
     * @param indicator is the progress indicator
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST, typename INDICATOR_TYPE,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST>, bool> = true>
    inline TissueSimulation& run(SIMULATION_TEST& done, INDICATOR_TYPE& indicator)
    {
        return run<SIMULATION_TEST, UI::Plot2DWindow, INDICATOR_TYPE>(done, nullptr, &indicator);
    }

    /**
     * @brief Simulate a tissue
     *
     * This method simulates a tissue until a test on the simulation
     * object holds. The simulation is also plotted in a graphical window.
     *
     * @tparam SIMULATION_TEST is a class whose objects can test whether
     *              the simulation is concluded
     * @tparam PLOT_WINDOW is the plotting window type
     * @param done is a test verifying whether the simulation has come to an end
     * @param plotter is a tissue plotter
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST, typename PLOT_WINDOW,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST> &&
                              std::is_base_of_v<UI::Plot2DWindow, PLOT_WINDOW>, bool> = true>
    inline TissueSimulation& run(SIMULATION_TEST& done, UI::TissuePlotter<PLOT_WINDOW>& plotter)
    {
        return run<SIMULATION_TEST, PLOT_WINDOW, UI::ProgressBar>(done, &plotter, nullptr);
    }

    /**
     * @brief Simulate a tissue
     *
     * This method simulates a tissue until a test on the simulation
     * object holds.
     *
     * @tparam SIMULATION_TEST is a class whose objects can test whether
     *              the simulation is concluded
     * @param done is a test verifying whether the simulation has come to an end
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST>, bool> = true>
    inline TissueSimulation& run(SIMULATION_TEST& done)
    {
        return run<SIMULATION_TEST, UI::Plot2DWindow, UI::ProgressBar>(done, nullptr, nullptr);
    }

    /**
     * @brief Get the last snapshot status
     *
     * @return the last snapshot status
     */
    inline const StatusAtSnapshot& get_last_snapshot_status() const
    {
        return status_at_snapshot;
    }

    /**
     * @brief Get the snapshot trigger
     *
     * @return A reference to the snapshot trigger
     */
    inline SnapshotTrigger& get_snapshot_trigger()
    {
        return snapshot_trigger;
    }

    /**
     * @brief Get the snapshot trigger
     *
     * @return A constant reference to the snapshot trigger
     */
    inline const SnapshotTrigger& get_snapshot_trigger() const
    {
        return snapshot_trigger;
    }

    /**
     * @brief Get the species identifiers from a name
     *
     * This method returns all the species identifiers corresponding to a
     * mutant or a species name.
     *
     * @param name is a mutant or a species name
     * @return a list of species identifiers corresponding to species whose name
     *      or mutant name is `name`
     */
    std::list<SpeciesId> get_species_ids_from(const std::string& name) const;

    /**
     * @brief Get the species identifiers from a list of names
     *
     * This method returns all the species identifiers corresponding to a
     * list of mutant or a species names.
     *
     * @param names is a list of mutant or a species names
     * @return a list of species identifiers corresponding to species whose name
     *      or mutant name is among those contained in `names`
     */
    std::list<SpeciesId> get_species_ids_from(const std::list<std::string>& names) const;

    /**
     * @brief Randomly select a cell
     *
     * This method randomly select a cell available for an event among those
     * belonging to specified species, and being available for an event.
     *
     * @param species_ids is the list of the species identifiers among those
     *      the cell must be chosen
     * @param event_type is the event type for which the choosen cell must be
     *          available
     * @return whenever the set of cells whose species identifier belongs
     *      to `species_ids` and available for `event_type` is not empty, a
     *      randomly selected cell in it. Otherwise, if the set is empty, a
     *      domain error is thrown.
     */
    const CellInTissue&
    choose_cell_in(const std::list<SpeciesId>& species_ids,
                   const CellEventType& event_type=CellEventType::ANY);

    /**
     * @brief Randomly select a cell in a tissue rectangle
     *
     * This method randomly select a cell available for an event among those
     * belonging to specified species, laying in a tissue rectangle, and being
     * available for an event.
     *
     * @param species_ids is the list of the species identifiers among those
     *      the cell must be chosen
     * @param rectangle is the tissue rectangle in which the cell must be selected
     * @param event_type is the event type for which the choosen cell must be
     *          available
     * @return whenever the set of cells whose species identifier belongs
     *      to `species_ids`, laying in `rectangle`, and available for
     *      `event_type` is not empty, a randomly selected cell in it.
     *      Otherwise, if the set is empty, a domain error is thrown.
     */
    const CellInTissue&
    choose_cell_in(const std::list<SpeciesId>& species_ids,
                   const RectangleSet& rectangle,
                   const CellEventType& event_type=CellEventType::ANY);

    /**
     * @brief Randomly select a cell on the border of the non-wild-type mass
     *
     * This method randomly select a cell on the external border of the
     * non-wild-type mass among those belonging to some species.
     *
     * @param species_ids is the list of the species identifiers among those
     *      the cell must be chosen
     * @return whenever the set of cells whose species identifier belongs
     *      to `species_ids` is not empty, a randomly selected cell in it.
     *      Otherwise, if the set is empty, a runtime error is thrown.
     */
    const CellInTissue&
    choose_border_cell_in(const std::list<SpeciesId>& species_ids);

    /**
     * @brief Randomly select a cell on the border of the non-wild-type mass
     *
     * This method randomly select a cell on the external border of the
     * non-wild-type mass among those belonging the some species.
     *
     * @param species_ids is the list of the species identifiers among those
     *      the cell must be chosen
     * @param rectangle is the tissue rectangle in which the cell must be
     *      selected
     * @return whenever the set of cells whose species identifier belongs
     *      to `species_ids` and laying in `rectangle` is not empty, a
     *      randomly selected cell in it.
     *      Otherwise, if the set is empty, a runtime error is thrown.
     */
    const CellInTissue&
    choose_border_cell_in(const std::list<SpeciesId>& species_ids,
                          const RectangleSet& rectangle);

    /**
     * @brief Randomly select a cell
     *
     * This method randomly select a cell available for an event among those
     * having a specified mutant or species name.
     *
     * @tparam STRINGS is the type of the first parameter
     * @tparam ARGS is the type of the parameter pack of the parameters from the
     *      second to the last one
     * @param names is the list of the mutant or species names that must contain
     *      the selected cell
     * @param args is the parameter pack of the remaining arguments of the
     *      function
     * @return whenever the set of cells whose mutant or species name is
     *      among those in `names` is not empty, a randomly selected cell in it.
     *      If one of the names does not correspond to either a mutant or a
     *      species or if any cell is available to be selected, a domain error
     *      is thrown.
     */
    template<typename STRINGS, typename... ARGS>
      requires (std::is_same_v<STRINGS, std::string>
                || std::is_same_v<STRINGS, std::list<std::string>>)
    const CellInTissue&
    choose_cell_in(const STRINGS& names, ARGS... args)
    {
        const std::list<SpeciesId> species_ids = get_species_ids_from(names);

        return choose_cell_in(species_ids, args...);
    }

    /**
     * @brief Randomly select a cell on the border of the non-wild-type mass
     *
     * This method randomly select a cell on the external border of the
     * non-wild-type mass among those belonging to one of the specified mutants
     * and species in a specified rectangle.
     *
     * @tparam STRINGS is the type of the first parameter
     * @tparam ARGS is the type of the parameter pack of the parameters from the
     *      second to the last one
     * @param names is the list of the mutant or species names that must contain
     *      the selected cell
     * @param args is the parameter pack of the remaining arguments of the
     *      function
     * @return whenever the set of cells whose mutant or species name is
     *      among those in `names` is not empty, a randomly selected cell in it.
     *      If one of the names does not correspond to either a mutant or a
     *      species or if any cell is available to be selected, a domain error
     *      is thrown.
     */
    template<typename STRINGS, typename... ARGS>
      requires (std::is_same_v<STRINGS, std::string>
                || std::is_same_v<STRINGS, std::list<std::string>>)
    const CellInTissue&
    choose_border_cell_in(const STRINGS& names, ARGS... args)
    {
        const std::list<SpeciesId> species_ids = get_species_ids_from(names);

        return choose_border_cell_in(species_ids, args...);
    }

    /**
     * @brief Get the current simulation time
     *
     * @return a constant reference to the simulation time
     */
    inline const Time& get_time() const
    {
        return time;
    }

    /**
     * @brief Get the simulation statistics
     *
     * @return the simulation statistics
     */
    inline const TissueStatistics& get_statistics() const
    {
        return statistics;
    }

    /**
     * @brief Get the simulation statistics
     *
     * @return the simulation statistics
     */
    inline TissueStatistics& get_statistics()
    {
        return statistics;
    }

    /**
     * @brief Get the cells manually added to the simulation
     *
     * @return a constant reference to the list of cells
     *          manually added to the simulation
     */
    inline const std::list<AddedCell>& get_added_cells() const
    {
        return added_cells;
    }

    /**
     * @brief Get the simulation lineage graph
     *
     * @return a constant reference to the simulation lineage graph
     */
    inline const LineageGraph& get_lineage_graph() const
    {
        return lineage_graph;
    }

    /**
     * @brief Add a mutant to the tissue
     *
     * @param mutant_properties is the mutant properties of the mutant
     * @return a reference to the updated object
     */
    TissueSimulation& add_mutant(const MutantProperties& mutant_properties);

    /**
     * @brief Add a mutant to the tissue
     *
     * @param mutant_name is the name of the mutant to add
     * @return a reference to the updated object
     */
    TissueSimulation& add_mutant(const std::string& mutant_name);

    /**
     * @brief Add a list of mutants to the tissue
     *
     * @param mutant_names is the list of the names of the mutants to add
     * @return a reference to the updated object
     */
    TissueSimulation& add_mutants(const std::list<std::string>& mutant_names);

    /**
     * @brief Add an epigenetic state to the tissue
     *
     * @param epistate_name is the name of the epigenetic state
     * @return a reference to the updated object
     */
    TissueSimulation& add_epigenetic_state(const std::string& epistate_name);

    /**
     * @brief Add a list of epigenetic states to the tissue
     *
     * @param epistate_names is the list of the epigenetic states to add
     * @return a reference to the updated object
     */
    TissueSimulation& add_epigenetic_states(const std::list<std::string>& epistate_names);

    /**
     * @brief Test whether the tissue knowns a species
     *
     * @param species_name is the species name
     * @return `true` if and only if the species has been added to the
     *      tissue
     */
    inline bool knowns_species(const std::string& species_name) const
    {
        return tissue().knowns_species(species_name);
    }

    /**
     * @brief Test whether the tissue knowns a mutant
     *
     * @param mutant_name is the mutant name
     * @return `true` if and only if the mutant has been added to the
     *      tissue
     */
    inline bool knowns_mutant(const std::string& mutant_name) const
    {
        return (mutant_name2id.find(mutant_name) != mutant_name2id.end());
    }

    /**
     * @brief Test whether the tissue knowns an epigenetic state
     *
     * @param epistate_name is the epigenetic state name
     * @return `true` if and only if the epigenetic state has been
     *      added to the tissue
     */
    inline bool knowns_epigenetic(const std::string& species_name) const
    {
        return tissue().knowns_species(species_name);
    }

    /**
     * @brief Get the mutant names
     *
     * @return A list of the mutant names
     */
    std::list<std::string> get_mutant_names() const;

    /**
     * @brief Get the epigenetic state names
     *
     * @return A list of the epigenetic  names
     */
    std::list<std::string> get_epigenetic_state_names() const
    {
        const auto& epistate_names = tissue().get_epigenetic_state_names();

        return {epistate_names.begin(), epistate_names.end()};
    }

    /**
     * @brief Place a cell in the simulated tissue
     *
     * @param species_name is the species name of the new cell
     * @param position is the cell position in the tissue
     * @return a reference to the updated object
     */
    TissueSimulation& place_cell(const std::string& species_name,
                                 const PositionInTissue& position);

    /**
     * @brief Place a cell in the simulated tissue
     *
     * @param species_properties is the species properties of the new cell
     * @param position is the cell position in the tissue
     * @return a reference to the updated object
     */
    inline TissueSimulation& place_cell(const SpeciesProperties& species_properties,
                                        const PositionInTissue& position)
    {
        return place_cell(species_properties.get_id(), position);
    }

    /**
     * @brief Place a cell in the simulated tissue
     *
     * @param species_id is the species identifier of the new cell
     * @param position is the cell position in the tissue
     * @return a reference to the updated object
     */
    TissueSimulation& place_cell(const SpeciesId& species_id, const PositionInTissue& position);

    /**
     * @brief Set a new simulation tissue
     *
     * This method resets the simulation and sets a
     * new simulation tissue.
     *
     * @param name is the tissue name
     * @param sizes are the sizes of the tissue
     * @return a reference to the updated object
     */
    TissueSimulation& set_tissue(const std::string& name, const std::vector<AxisSize>& sizes);

    /**
     * @brief Get the simulation tissue
     *
     * This method returns a reference to the simulation tissue.
     * A `std::runtime_error` object is throws if no tissue has been
     * associated to the simulation yet.
     *
     * @return A constant reference to the associated tissue
     */
    const Tissue& tissue() const;

    /**
     * @brief Get the simulated tissue
     *
     * This method returns a reference to the simulated tissue.
     * A `std::runtime_error` object is throws if no tissue has been
     * associated to the simulation yet.
     *
     * @return A constant reference to the associated tissue
     */
    Tissue& tissue();

    /**
     * @brief Reset a simulation
     */
    void reset();

    /**
     * @brief Update the log directory
     *
     * @param log_directory is the new simulation log directory
     */
    inline void rename_log_directory(const std::filesystem::path& log_directory)
    {
        logger.rename_directory(log_directory);
    }

    /**
     * @brief Inject random generator seed
     *
     * @param random_seed is the simulation random seed
     */
    inline TissueSimulation& random_generator_seed(int random_seed)
    {
        random_gen.seed(random_seed);

        return *this;
    }

    /**
     * @brief Performs a simulation snapshot
     *
     * This method performs a simulation snapshot.
     *
     * @tparam INDICATOR is the type of progress bar
     * @param indicator is the progress bar
     */
    template<typename INDICATOR>
    void make_snapshot(INDICATOR *indicator);

    /**
     * @brief Simulate tissue sampling
     *
     * This method simulates a tissue sampling collecting cells
     * data, but avoiding cell removal. The resulting sample is
     * not stored among those collected during the simulation.
     *
     * @param specification is the sample specification
     * @return the sample of the tissue in `rectangle`
     */
    TissueSample simulate_sampling(const SampleSpecification& specification);

    /**
     * @brief Sample the simulation tissue
     *
     * This method performs a tissue sampling collecting cells
     * data and removing them from them the tissue. The resulting
     * sample is stored among those collected during the
     * simulation.
     *
     * @param specification is the sample specification
     * @return the sample of the tissue in `rectangle`
     */
    TissueSample sample_tissue(const SampleSpecification& specification);

    /**
     * @brief Return all the simulation samples
     *
     * This method returns all the simulation samples collected during
     * the computation.
     *
     * @return the list of all the tissue samples collected during the
     *      computation sorted by collection time
     */
    inline
    const std::list<TissueSample>& get_tissue_samples() const
    {
        return samples;
    }

    /**
     * @brief Get the simulation logger
     *
     * @return a constant reference to the simulation logger
     */
    inline const BinaryLogger& get_logger() const
    {
        return logger;
    }

    /**
     * @brief Find a mutant id by name
     *
     * @param mutant_name is the name of the mutant whose id is aimed
     * @return a constant reference to the identifier of the mutant having
     *         `mutant_name` as name
     */
    const MutantId& find_mutant_id(const std::string& mutant_name) const;

    /**
     * @brief Find a mutant name by id
     *
     * @param mutant_id is the identifier of the mutant whose name is aimed
     * @return a constant reference to the name of the mutant having
     *         `mutant_id` as identifier
     */
    const std::string& find_mutant_name(const MutantId& mutant_id) const;

    /**
     * @brief Get the logic variable of a species cardinality
     *
     * @param species_name is the name of a species
     * @return the logic variable associated to the cardinality of the species
     *      whose name is `species_name`
     */
    inline Logics::Variable get_cardinality_variable(const std::string& species_name) const
    {
        return tissue().get_cardinality_variable(species_name);
    }

    /**
     * @brief Get the logic variable representing the number of a species event
     *
     * @param species_name is the name of a species
     * @param event_type is the variable event
     * @return the logic variable associated to the number event of type `event_type`
     *      occurring in the species whose name is `species_name`
     */
    inline Logics::Variable get_event_variable(const std::string& species_name,
                                               const CellEventType& event_type) const
    {
        return tissue().get_event_variable(species_name, event_type);
    }

    /**
     * @brief Get the logic variable of a species cardinality
     *
     * @param species_id is the identifier of a species
     * @return the logic variable associated to the cardinality of the species
     *      whose identifier is `species_id`
     */
    inline Logics::Variable get_cardinality_variable(const SpeciesId& species_id) const
    {
        return tissue().get_cardinality_variable(species_id);
    }

    /**
     * @brief Get the logic variable representing the number of a species event
     *
     * @param species_id is the identifier of a species
     * @param event_type is the variable event
     * @return the logic variable associated to the number event of type `event_type`
     *      occurring in the species whose name is `species_name`
     */
    inline Logics::Variable get_event_variable(const SpeciesId& species_id,
                                               const CellEventType& event_type) const
    {
        return tissue().get_event_variable(species_id, event_type);
    }

    /**
     * @brief Get the logic variable representing the evolution time
     *
     * @return the logic variable representing the evolution time
     */
    inline Logics::Variable get_time_variable() const
    {
        return Logics::Variable();
    }

    /**
     * @brief Evaluate a variable
     *
     * @param variable is the variable to be valuated
     * @return the value of the variable
     */
    double evaluate(const Logics::Variable& variable) const;

    /**
     * @brief Save a simulation in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "CLONES Tissue Simulation", 0);

        archive & tissues
                & lineage_graph
                & mutant_name2id
                & logger
                & status_at_snapshot
                & snapshot_trigger
                & statistics
                & time
                & timed_event_queue
                & death_activation_level
                & duplicate_internal_cells
                & storage_enabled
                & samples
                & next_cell_id;
    }

    /**
     * @brief Load a simulation from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded simulation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static TissueSimulation load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "CLONES Tissue Simulation", 1);

        TissueSimulation simulation;

        archive & simulation.tissues
                & simulation.lineage_graph
                & simulation.mutant_name2id
                & simulation.logger
                & simulation.status_at_snapshot
                & simulation.snapshot_trigger
                & simulation.statistics
                & simulation.time
                & simulation.timed_event_queue
                & simulation.death_activation_level
                & simulation.duplicate_internal_cells
                & simulation.storage_enabled
                & simulation.samples
                & simulation.next_cell_id;

        for (auto sample_it=simulation.samples.begin();
                sample_it != simulation.samples.end(); ++sample_it) {
            simulation.name2sample[sample_it->get_name()] = sample_it;
        }

        simulation.init_valid_directions();

        simulation.status_at_snapshot.set_time();

        return simulation;
    }

    ~TissueSimulation();
};

/* Template implementations */

template<typename INDICATOR>
void TissueSimulation::make_snapshot(INDICATOR *indicator)
{
    if (storage_enabled) {
        if (indicator != nullptr) {
            indicator->set_message("Saving snapshot");
        }

        logger.snapshot(*this);
        logger.flush_archives();
    }

    status_at_snapshot = StatusAtSnapshot(*this);
}

template<typename PLOT_WINDOW>
inline
TissueSimulation& TissueSimulation::run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>* plotter)
{
    CellEvent event = select_next_event();

    return simulate(event, plotter);
}

template<typename PLOT_WINDOW>
TissueSimulation& TissueSimulation::simulate(const CellEvent& event, UI::TissuePlotter<PLOT_WINDOW>* plotter)
{
    time += event.delay;

    EventAffectedCells affected;

    switch(event.type) {
        case CellEventType::DEATH:
            affected = simulate_death(event.position);
            break;
        case CellEventType::DUPLICATION:
            affected = simulate_duplication(event.position);
            break;
        case CellEventType::DUP_AND_EPI_SWITCH:
        case CellEventType::MUTATION:
            affected = simulate_duplication_and_set_one_child_species(event.position,
                                                                      event.dst_species);
            break;
        default:
            throw Error<std::runtime_error>("Unhandled event type \""
                                            + std::string(cell_event_names[event.type])
                                            + "\".");
    }

    for (const auto& cell : affected.new_cells) {
        if (storage_enabled) {
            logger.record(event.type, cell, time);
        }

        // if death has not been enabled yet
        const auto species_id = cell.get_species_id();
        Species& species = tissue().get_species(species_id);

        if ((!species.is_death_enabled())
                && (species.num_of_cells()>=death_activation_level)) {
            species.enable_death();
        }
    }

    statistics.record_event(event, time);
    for (const auto& cell: affected.lost_cells) {
        statistics.record_lost(cell.get_species_id(), time);
    }

    if (plotter != nullptr && !plotter->closed()) {
        plotter->plot(statistics);
    }

    return *this;
}

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_SIMULATOR__
