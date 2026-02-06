/**
 * @file statistics.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines simulation statistics
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

#ifndef __CLONES_STATISTICS__
#define __CLONES_STATISTICS__

#include <list>

#include "species.hpp"
#include "cell_event.hpp"


namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

class TissueStatistics;

/**
 * @brief A class for accounting simulation species statistics
 */
struct SpeciesStatistics
{
    Time rise_time;         //!< the time of the first appearance
    Time extinction_time;   //!< the time of the last appearance

    size_t total_cells;     //!< the number of cells that have been in the species
    size_t curr_cells;      //!< the number of cells currently in the species
    size_t killed_cells;    //!< the number of species cells that were killed
    size_t lost_cells;      //!< the number of cells that have overcome the tissue border
    size_t num_duplications;    //!< the number of duplications in the species
    std::map<SpeciesId, size_t> epigenetic_events; //!< the number of epigenetic event per species destination

    /**
     * @brief The empty constructor
     */
    SpeciesStatistics();

    /**
     * @brief A constructor
     *
     * @param num_of_cells is the number of cells in the species
     */
    explicit SpeciesStatistics(const size_t& num_of_cells);

    /**
     * @brief Compute the number of epigenetic event from the species
     *
     * @return the number of epigenetic event from the species
     */
    size_t num_of_epigenetic_events() const;

    /**
     * @brief Compute the number of epigenetic event leading to a species
     *
     * @param dst_id is the identifier of the species reached by the epigenetic
     *      events that must be counted
     * @return the number of epigenetic event leading to the species having
     *      `dst_id` as identifier
     */
    size_t num_of_epigenetic_events(const SpeciesId& dst_id) const;

    /**
     * @brief Save species statistics in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & rise_time
                & extinction_time
                & total_cells
                & curr_cells
                & killed_cells
                & lost_cells
                & num_duplications
                & epigenetic_events;
    }

    /**
     * @brief Load species statistics from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded species statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SpeciesStatistics load(ARCHIVE& archive)
    {
        SpeciesStatistics stats;

        archive & stats.rise_time
                & stats.extinction_time
                & stats.total_cells
                & stats.curr_cells
                & stats.killed_cells
                & stats.lost_cells
                & stats.num_duplications
                & stats.epigenetic_events;

        return stats;
    }

    friend class TissueStatistics;
};

/**
 * @brief A class to collect simulation statistics about a tissue
 */
class TissueStatistics
{
    using time_point = std::chrono::steady_clock::time_point;

    /**
     * @brief A structure that associate a species identifier to the species statistics
     */
    using SpeciesMap = std::map<SpeciesId, SpeciesStatistics>;

    std::map<Time, SpeciesMap> history;     //!< The statistics history
    Time history_delta;                     //!< The history sampling delta

    SpeciesMap s_statistics;                //!< The current species statistics

    std::list<Time> sim_times;              //!< The simulated times of the last recorded events
    std::list<time_point> real_times;       //!< The recording times of the last events

    size_t max_stored_times;                //!< The maximum number of times to store

    size_t total_events;                    //!< The number of recorded events
    time_point first_event_time;            //!< The first recoded event time

    /**
     * @brief Record a cell duplication and avoid history saving
     *
     * @param species_id is the species id of the duplicating cell
     */
    void record_duplication_no_save(const SpeciesId& species_id);

    /**
     * @brief Record a species change
     *
     * @param src_species is the source species identifier
     * @param dst_species is the destination species identifier
     * @param time is the epigenetic event time
     */
    void record_species_change(const SpeciesId& src_species,
                               const SpeciesId& dst_species, const Time &time);

    /**
     * @brief Record a duplication and a species change
     *
     * @param src_species is the source species identifier
     * @param dst_species is the destination species identifier
     * @param time is the epigenetic event time
     */
    inline
    void record_duplication_and_species_change(const SpeciesId& src_species,
                                               const SpeciesId& dst_species,
                                               const Time &time)
    {
        record_duplication(src_species, time);
        record_species_change(src_species, dst_species, time);
    }

    /**
     * @brief Save the current statistics if it is needed
     *
     * @param time is the time of the to-be-recorded event
     */
    void save_in_history_if_needed(const Time &time);

public:

    /**
     * @brief An empty constructor
     */
    TissueStatistics();

    /**
     * @brief Construct a new Tissue Statistics
     *
     * @param delta is the history sampling time
     */
    TissueStatistics(const Time& delta);

    /**
     * @brief Get the statistics of a species
     *
     * @param species is the species whose statistics are aimed
     * @return a non-constant reference to the statistics of `species`
     */
    inline SpeciesStatistics& operator[](const Species& species)
    {
        return s_statistics.at(species.get_id());
    }

    /**
     * @brief Get the statistics of a species
     *
     * @param species_id is the identifier of the species whose statistics are aimed
     * @return a non-constant reference to the statistics of species having
     *      `species_id` as identifier
     */
    inline SpeciesStatistics& operator[](const SpeciesId& species_id)
    {
        return s_statistics.at(species_id);
    }

    /**
     * @brief Get the statistics of a species
     *
     * @param species is the species whose statistics are aimed
     * @return a non-constant reference to the statistics of `species`
     */
    inline const SpeciesStatistics& at(const Species& species) const
    {
        return s_statistics.at(species.get_id());
    }

    /**
     * @brief Get the statistics of a species
     *
     * @param species_id is the identifier of the species whose statistics are aimed
     * @return a constant reference to the statistics of species having
     *      `species_id` as identifier
     */
    inline const SpeciesStatistics& at(const SpeciesId& species_id) const
    {
        return s_statistics.at(species_id);
    }

    /**
     * @brief Test whether the object contains statistics for a species
     *
     * @param species_id is the identifier of the species whose statistics are aimed
     * @return `true` if and only if the object contains statistics for the
     *          specified species
     */
    inline bool contains_data_for(const SpeciesId& species_id) const
    {
        return s_statistics.count(species_id)==1;
    }

    /**
     * @brief Test whether the object contains statistics for a species
     *
     * @param species is the species whose statistics are aimed
     * @return `true` if and only if the object contains statistics for the
     *          specified species
     */
    inline bool contains_data_for(const Species& species) const
    {
        return s_statistics.count(species.get_id())==1;
    }

    /**
     * @brief Record a death
     *
     * @param species_id is the species id of the dying cell
     * @param time is the death time
     */
    void record_death(const SpeciesId& species_id, const Time &time);

    /**
     * @brief Record a lost cell
     *
     * A cell is lost whenever it is pushed outside the tissue border.
     * This method record in the statistics that a cell is lost.
     *
     * @param species_id is the species id of the lost cell
     * @param time is the time in which the lost cell has been pushed
     *       outside the tissue border
     */
    void record_lost(const SpeciesId& species_id, const Time &time);

    /**
     * @brief Record a cell duplication
     *
     * @param species_id is the species id of the duplicating cell
     * @param time is the duplication time
     */
    inline
    void record_duplication(const SpeciesId& species_id, const Time &time)
    {
        record_duplication_no_save(species_id);

        save_in_history_if_needed(time);
    }

    /**
     * @brief Record a mutation
     *
     * @param src_species is the source species identifier
     * @param dst_species is the destination species identifier
     * @param time is the epigenetic switch time
     */
    inline
    void record_mutation(const SpeciesId& src_species,
                         const SpeciesId& dst_species,
                         const Time &time)
    {
        record_duplication_and_species_change(src_species, dst_species, time);

        save_in_history_if_needed(time);
    }

    /**
     * @brief Record an epigenetic switch
     *
     * @param src_species is the source species identifier
     * @param dst_species is the destination species identifier
     * @param time is the epigenetic switch time
     */
    void record_epigenetic_switch(const SpeciesId& src_species, const SpeciesId& dst_species,
                                  const Time &time);

    /**
     * @brief Record placed cell
     *
     * This method records a placed cell
     *
     * @param species_id is the species of the placed cell
     * @param time is the placing d time
     */
    void record_placed_cell(const SpeciesId& species_id, const Time &time);

    /**
     * @brief Record the last event
     *
     * This method records a new event and, whenever the
     * queues of the recorded times overlaps the maximum size,
     * deletes the oldest recorded times.
     *
     * @param event is the event to record
     * @param time is the event time
     */
    void record_event(const CellEvent& event, const Time &time);

    /**
     * @brief Get the number of recorded times
     *
     * @return the number of recorded times
     */
    const size_t& get_recorded_time_number() const;

    /**
     * @brief Get the elapsed time
     *
     * @return the elapsed time
     */
    inline std::chrono::steady_clock::duration get_elapsed_time() const
    {
        return std::chrono::steady_clock::now()-first_event_time;
    }

    /**
     * @brief Get the simulated time
     *
     * @return the simulated time
     */
    inline Time get_simulated_time() const
    {
        return sim_times.back();
    }

    /**
     * @brief Get the number of the last recorded events over time
     *
     * @return the number of the last recorded events over time
     */
    template<class ToDuration>
    double get_last_recorded_events_over_time() const
    {
        using namespace std::chrono;

        auto time = duration_cast<ToDuration>(real_times.back()-real_times.front()).count();

        return static_cast<double>(real_times.size())/time;
    }

    /**
     * @brief Store current statistical data in history
     *
     * @param time is the current simulation time
     */
    inline void store_current_in_history(const Time &time)
    {
        history[time] = s_statistics;
    }

    /**
     * @brief Reset the history
     */
    inline void clear_history()
    {
        history.clear();
    }

    /**
     * @brief Get the statistical history
     *
     * @return return a constant reference to the statistics history
     */
    inline const std::map<Time, std::map<SpeciesId, SpeciesStatistics>>&
    get_history() const
    {
        return history;
    }

    /**
     * @brief Get the last sample time
     *
     * @return the last sample time
     */
    Time get_last_time_in_history() const;

    /**
     * @brief Get the history time delta
     *
     * @return a constant reference to the the history time delta
     */
    inline const Time& get_history_delta() const
    {
        return history_delta;
    }

    /**
     * @brief Set the history time delta
     *
     * @param delta is the non-negative value to be set as history time delta
     */
    inline void set_history_delta(const Time& delta)
    {
        if (delta<0) {
            throw std::domain_error("The history time delta must be non-negative");
        }

        history_delta = delta;
    }

    /**
     * @brief Count fired event
     *
     * @param species_id is the identifier of the species in which the
     *   to-be-counted events occurred
     * @param event_type is the type of the event to-be-counted
     * @return the number of `event_type` events that occurred in the
     *   species whose identifier is `species_id`
     */
    size_t count_fired_events(const SpeciesId& species_id,
                              const CellEventType& event_type) const;

    /**
     * @brief Save tissue statistics in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & history
                & history_delta
                & s_statistics
                & sim_times
                & real_times
                & max_stored_times
                & total_events
                & first_event_time;
    }

    /**
     * @brief Load tissue statistics from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded tissue statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static TissueStatistics load(ARCHIVE& archive)
    {
        TissueStatistics stats;

        archive & stats.history
                & stats.history_delta
                & stats.s_statistics
                & stats.sim_times
                & stats.real_times
                & stats.max_stored_times
                & stats.total_events
                & stats.first_event_time;

        return stats;
    }
};

}   // Evolutions

}   // Mutants

}   // CLONES


#endif // __CLONES_STATISTICS__
