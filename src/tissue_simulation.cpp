/**
 * @file simulation.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define a tumour evolution simulation
 * @version 1.18
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

#include <set>
#include <list>
#include <string>
#include <limits>

#include "tissue_simulation.hpp"

#include "error.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

TissueSimulation::StatusAtSnapshot::StatusAtSnapshot():
    time{StatusAtSnapshot::Clock::now()},
    clock{static_cast<Time>(0)}, num_of_cells{0}
{}

TissueSimulation::StatusAtSnapshot::StatusAtSnapshot(const TissueSimulation& simulation):
    time{StatusAtSnapshot::Clock::now()},
    clock{simulation.get_time()},
    num_of_cells{simulation.tissue().num_of_mutated_cells()}
{}

TissueSimulation::SnapshotTrigger::SnapshotTrigger():
    time_delta{std::numeric_limits<Duration::rep>::max()},
    clock_delta{std::numeric_limits<Duration::rep>::max()},
    cardinality_delta{std::numeric_limits<uint64_t>::max()}
{}

bool TissueSimulation::SnapshotTrigger::is_triggered_by(const TissueSimulation& simulation) const
{
    const auto& last_status = simulation.get_last_snapshot_status();

    auto now = StatusAtSnapshot::Clock::now();

    return ((last_status.get_num_of_cells() + cardinality_delta <=
                simulation.tissue().num_of_cells())
            || (last_status.get_clock() + clock_delta <=
                simulation.get_time())
            || (now - last_status.get_time() > time_delta));
}

TissueSimulation::AddedCell::AddedCell():
    species_id(WILD_TYPE_SPECIES)
{}

TissueSimulation::AddedCell::AddedCell(const SpeciesId& species, const PositionInTissue& position,
                                 const Time& time):
    PositionInTissue(position), species_id{species}, time{time}
{}

TissueSimulation::TissueSimulation(int random_seed):
    logger{}, status_at_snapshot{}, snapshot_trigger{},
    time{static_cast<Time>(0)},
    next_cell_id{Cell::first_tumour_cell_id()}, death_activation_level{1},
    duplicate_internal_cells{false}, storage_enabled{true}
{
    random_gen.seed(random_seed);

    // Create a default tissue
    set_tissue("A tissue", {1000, 1000});
}

TissueSimulation::TissueSimulation(const std::filesystem::path& log_directory, int random_seed):
    logger{log_directory}, status_at_snapshot{}, snapshot_trigger{},
    time{static_cast<Time>(0)},
    next_cell_id{Cell::first_tumour_cell_id()}, death_activation_level{1},
    duplicate_internal_cells{false}, storage_enabled{true}
{
    random_gen.seed(random_seed);

    // Create a default tissue
    set_tissue("A tissue", {1000, 1000});
}

TissueSimulation::TissueSimulation(TissueSimulation&& orig):
    TissueSimulation()
{
    *this = std::move(orig);
}

TissueSimulation& TissueSimulation::operator=(TissueSimulation&& orig)
{
    std::swap(tissues, orig.tissues);
    std::swap(lineage_graph, orig.lineage_graph);
    std::swap(mutant_name2id, orig.mutant_name2id);
    std::swap(logger, orig.logger);
    std::swap(status_at_snapshot, orig.status_at_snapshot);
    std::swap(snapshot_trigger, orig.snapshot_trigger);
    std::swap(statistics, orig.statistics);
    std::swap(time, orig.time);
    std::swap(timed_event_queue, orig.timed_event_queue);
    std::swap(death_activation_level, orig.death_activation_level);
    std::swap(duplicate_internal_cells, orig.duplicate_internal_cells);
    std::swap(storage_enabled, orig.storage_enabled);
    std::swap(samples, orig.samples);

    std::swap(valid_directions, orig.valid_directions);
    std::swap(random_gen, orig.random_gen);
    std::swap(next_cell_id, orig.next_cell_id);

    return *this;
}

const Tissue& TissueSimulation::tissue() const
{
    if (tissues.size()==0) {
        throw Error<std::runtime_error>("No tissue has been associated "
                                        "to the simulation yet.");
    }
    return tissues[0];
}

Tissue& TissueSimulation::tissue()
{
    if (tissues.size()==0) {
        throw Error<std::runtime_error>("No tissue has been associated "
                                        "to the simulation yet.");
    }
    return tissues[0];
}


template<typename GENERATOR_TYPE>
void select_next_event_in_species(CellEvent& event, Tissue& tissue,
                                  const Species& species,
                                  std::uniform_real_distribution<double>& uni_dist,
                                  GENERATOR_TYPE& random_gen)
{
    if (species.num_of_cells()==0) {
        return;
    }

    bool selected_new_event = false;
    for (const auto& [event_type, event_type_rates]: species.get_rates()) {
        for (const auto& [dst_event, rate]: event_type_rates) {
            if (rate > 0) {
                const auto available_cells = species.num_of_cells_available_for(event_type);

                if (available_cells > 0) {
                    const auto r_value = uni_dist(random_gen);

                    // sample an exponential distribution with lambda=available_cells*rate
                    const auto candidate_delay = -log(r_value) / (available_cells * rate);

                    if (event.delay > candidate_delay) {
                        event.delay = candidate_delay;
                        event.dst_species = dst_event;
                        event.type = event_type;

                        // delaying cell selection because is expensive

                        selected_new_event = true;
                    }
                }
            }
        }
    }

    if (selected_new_event) {
        // at the end select cell selection

        event.position = Position(tissue, species.choose_a_cell(random_gen, event.type));
        event.src_species = species.get_id();
    }
}

std::list<SpeciesId> TissueSimulation::get_species_ids_from(const std::string& name) const
{
    auto mutant_it = mutant_name2id.find(name);
    if (mutant_it != mutant_name2id.end()) {
        std::list<SpeciesId> species_ids;

        for (const auto& species : tissue().get_mutant_view(mutant_it->second)) {
            species_ids.push_back(species.get_id());
        }

        return species_ids;
    }

    if (tissue().knowns_species(name)) {
        const Species& species = tissue().get_species(name);
        return {species.get_id()};
    }

    throw Error<std::runtime_error>("\"" + name +"\" is neither a mutant nor a species name.");
}

std::list<SpeciesId> TissueSimulation::get_species_ids_from(const std::list<std::string>& names) const
{
    std::list<SpeciesId> species_ids;

    for (const auto& name : names) {
        species_ids.splice(species_ids.end(), get_species_ids_from(name));
    }

    return species_ids;
}

const CellInTissue&
TissueSimulation::choose_cell_in(const std::list<SpeciesId>& species_ids, const CellEventType& event_type)
{
    if (species_ids.size()==0) {
        throw Error<std::domain_error>("The parameter \"species_ids\" must be non-empty.");
    }

    size_t total = 0;
    std::vector<size_t> num_of_cells;
    std::vector<SpeciesId> species_vec{species_ids.begin(), species_ids.end()};

    for (const SpeciesId& species_id : species_ids) {
        const Species& species = tissue().get_species(species_id);

        num_of_cells.push_back(species.num_of_cells_available_for(event_type));
        total += num_of_cells.back();
    }

    if (total == 0) {
        throw Error<std::runtime_error>("No cell satisfies the request.");
    }

    std::discrete_distribution<size_t> distribution(num_of_cells.begin(),
                                                    num_of_cells.end());
    const SpeciesId& species_id = species_vec[distribution(random_gen)];

    return tissue().get_species(species_id).choose_a_cell(random_gen, event_type);
}

const CellInTissue&
TissueSimulation::choose_cell_in(const std::list<SpeciesId>& species_ids,
                                 const RectangleSet& rectangle,
                                 const CellEventType& event_type)
{
    if (species_ids.size()==0) {
        throw Error<std::domain_error>("The parameter \"species_ids\" must be non-empty.");
    }

    std::set<SpeciesId> species_set{species_ids.begin(), species_ids.end()};

    std::list<PositionInTissue> position_vector;
    for (const auto& position : rectangle) {
        if (tissue()(position).is_available_for(event_type)) {
            const auto cell_proxy = tissue()(position);
            const auto species_id = static_cast<const CellInTissue&>(cell_proxy).get_species_id();

            if (species_set.count(species_id)>0) {
                position_vector.push_back(position);
            }
        }
    }

    std::uniform_int_distribution<size_t> distribution(1,position_vector.size());

    auto selected = distribution(random_gen);

    size_t i{0};
    for (const auto& position: position_vector) {
        if (++i==selected) {
            return tissue()(position);
        }
    }

    throw Error<std::runtime_error>("No cell satisfies the request.");
}

bool border_visible_from(PositionInTissue pos, const Tissue& tissue, const PositionDelta& delta)
{
    pos -= delta;
    while (tissue.is_valid(pos)) {
        const auto cell_proxy = tissue(pos);

        if (cell_proxy.is_wild_type()) {
            pos -= delta;
        } else {
            return false;
        }
    }

    return true;
}

bool search_tumour_cell(PositionInTissue& pos, const RectangleSet& rectangle,
                        const Tissue& tissue, const PositionDelta& delta,
                        const std::set<SpeciesId>& species_ids)
{
    while (tissue.is_valid(pos) && rectangle.contains(pos)) {
        const auto cell_proxy = tissue(pos);

        if (cell_proxy.is_wild_type()) {
            pos += delta;
        } else {
            const CellInTissue &cell = cell_proxy;

            return (species_ids.count(cell.get_species_id())>0);
        }
    }

    return false;
}

bool choose_border_cell_in(PositionInTissue& pos, const Tissue& tissue, const Direction& dir,
                           const std::set<SpeciesId>& species_ids, const RectangleSet& rectangle,
                           std::mt19937_64& random_gen)
{
    PositionDelta delta(dir);

    uint16_t random_init_z = static_cast<int16_t>(random_gen() % rectangle.depth());
    if (delta.x != 0) {
        uint16_t random_init_y = static_cast<int16_t>(random_gen() % rectangle.height());
        for (uint16_t y=0; y < rectangle.height(); ++y) {
            pos.y = static_cast<int16_t>((y+random_init_y) % rectangle.height())
                    + rectangle.lower_corner.y;
            for (uint16_t z=0; z < rectangle.depth(); ++z) {
                pos.z = static_cast<int16_t>((z+random_init_z) % rectangle.depth())
                        + rectangle.lower_corner.z;
                pos.x = (delta.x > 0?rectangle.lower_corner.x:rectangle.upper_corner.x);

                PositionInTissue new_pos{pos};
                if (search_tumour_cell(new_pos, rectangle, tissue, delta, species_ids)
                        && border_visible_from(pos, tissue, delta)) {
                    pos = new_pos;
                    return true;
                }
            }
        }
    }
    if (delta.y != 0) {
        uint16_t random_init_x = static_cast<int16_t>(random_gen() % rectangle.width());
        for (uint16_t x=0; x < rectangle.width(); ++x) {
            pos.x = static_cast<int16_t>((x+random_init_x) % rectangle.width())
                    + rectangle.lower_corner.x;

            for (uint16_t z=0; z < rectangle.depth(); ++z) {
                pos.z = static_cast<int16_t>((z+random_init_z) % rectangle.depth())
                        + rectangle.lower_corner.z;
                pos.y = (delta.y > 0?rectangle.lower_corner.y:rectangle.upper_corner.y);

                PositionInTissue new_pos{pos};
                if (search_tumour_cell(new_pos, rectangle, tissue, delta, species_ids)
                        && border_visible_from(pos, tissue, delta)) {
                    pos = new_pos;
                    return true;
                }
            }
        }
    }

    return false;
}

const CellInTissue&
TissueSimulation::choose_border_cell_in(const std::list<SpeciesId>& species_ids,
                                        const RectangleSet& rectangle)
{
    if (species_ids.size()==0) {
        throw Error<std::domain_error>("The parameter \"species_ids\" must be non-empty.");
    }

    std::set<SpeciesId> species_set{species_ids.begin(), species_ids.end()};

    PositionInTissue pos;
    size_t dir_offset = random_gen()%(valid_directions.size());
    for (size_t dir_idx=0; dir_idx < valid_directions.size(); ++dir_idx) {
        const auto& dir = valid_directions[(dir_offset+dir_idx)%valid_directions.size()];

        if (Evolutions::choose_border_cell_in(pos, tissue(), dir, species_set, rectangle, random_gen)) {
            return tissue()(pos);
        }
    }

    throw Error<std::runtime_error>("No border cells satisfies the request.");
}

const CellInTissue& TissueSimulation::choose_border_cell_in(const std::list<SpeciesId>& species_ids)
{
    RectangleSet rectangle;

    auto sizes = tissue().size();
    rectangle.lower_corner.x = 0;
    rectangle.upper_corner.x = sizes[0]-1;
    rectangle.lower_corner.y = 0;
    rectangle.upper_corner.y = sizes[1]-1;
    rectangle.lower_corner.z = 0;
    if (sizes.size()==3) {
        rectangle.upper_corner.y = sizes[2]-1;
    } else {
        rectangle.lower_corner.z = 0;
    }

    return choose_border_cell_in(species_ids, rectangle);
}

/**
 * @brief Create a mutation event
 *
 * @param tissue is the tissue in which the event will occurs
 * @param position is the position of the parent cell that will give birth to the mutated cell
 * @param final_id is the mutant identifier of the mutated cell
 * @param delay is the delay of the mutation with respect to the current simulation clock
 * @return the created cell event
 */
CellEvent create_mutation_event(Tissue& tissue, const PositionInTissue& position,
                                const MutantId& final_id, const Time& delay)
{
    CellEvent event;

    event.delay = delay;
    event.type = CellEventType::MUTATION;
    event.position = Position(tissue, position);
    event.src_species = static_cast<const CellInTissue&>(tissue(position)).get_species_id();

    const Species& src_species = tissue.get_species(event.src_species);
    const auto src_epistate = src_species.get_epistate_name();

    const auto& m_view = tissue.get_mutant_view(final_id);

    event.dst_species = m_view.get_species_by_epistate(src_epistate).get_id();

    return event;
}

bool TissueSimulation::handle_timed_mutation(const TimedEvent& timed_mutation, CellEvent& candidate_event)
{
    const auto& mutation = timed_mutation.get_event<Mutation>();

    try {
        std::list<SpeciesId> species_ids;
        for (const auto& species : tissue().get_mutant_view(mutation.initial_id)) {
            species_ids.push_back(species.get_id());
        }

        const auto& cell = choose_cell_in(species_ids);

        auto delay = (timed_mutation.time >= time ?
                        timed_mutation.time-time : 0);
        candidate_event = create_mutation_event(tissue(), cell, mutation.final_id, delay);

        return true;
    } catch (const std::runtime_error&) {
        return false;
    }
}

TissueSimulation& TissueSimulation::simulate_mutation(const PositionInTissue& position,
                                                      const std::string& dst_mutant_name)
{
    auto dst_mutant_id = find_mutant_id(dst_mutant_name);

    return simulate_mutation(position, dst_mutant_id);
}

TissueSimulation& TissueSimulation::simulate_mutation(const PositionInTissue& position,
                                                      const MutantId& dst_mutant_id)
{
    auto mutation_event = create_mutation_event(tissue(), position, dst_mutant_id, 0);

    simulate(mutation_event);

    if (storage_enabled) {
        logger.snapshot(*this);
        logger.flush_archives();
    }

    return *this;
}

void TissueSimulation::handle_timed_rate_update(const TimedEvent& timed_rate_update)
{
    const auto& rate_update = timed_rate_update.get_event<RateUpdate>();
    Species& species = tissue().get_species(rate_update.src_id);
    const Species& dst_species = tissue().get_species(rate_update.dst_id);

    species.set_rate(rate_update.event_type, dst_species,
                     rate_update.new_rate);
}

void TissueSimulation::handle_timed_sampling(const TimedEvent& timed_sampling, CellEvent& candidate_event)
{
    const auto& sampling = timed_sampling.get_event<Sampling>();

    const Time previous_time = time;

    time = timed_sampling.time;

    // sample tissue
    sample_tissue(sampling);

    // update candidate event to avoid removed regions
    candidate_event = select_next_cell_event();

    if (candidate_event.delay+previous_time < time) {
        candidate_event.delay = 0;
    } else {
        candidate_event.delay -= (time - previous_time);
    }
}

void TissueSimulation::handle_timed_event_queue(CellEvent& candidate_event)
{
    // if the timed event queue is not empty and the next time event occurs before the candidate cell event
    while (!timed_event_queue.empty() && (timed_event_queue.top().time <= candidate_event.delay + time)) {
        TimedEvent timed_event = timed_event_queue.top();

        timed_event_queue.pop();

        if (storage_enabled) {
            if (timed_event_queue.top().time != timed_event.time
                || timed_event_queue.top().type != timed_event.type) {
                logger.snapshot(*this);
            }
        }

        switch(timed_event.type) {
            case TimedEvent::Type::MUTATION:
                {
                    auto candidate_updated = handle_timed_mutation(timed_event, candidate_event);

                    if (candidate_updated) {
                        return;
                    }
                }
                break;
            case TimedEvent::Type::LIVENESS_RATE_UPDATE:
                {
                    handle_timed_rate_update(timed_event);
                }
                break;
            case TimedEvent::Type::SAMPLING:
                {
                    handle_timed_sampling(timed_event, candidate_event);
                }
                break;
            default:
                throw Error<std::runtime_error>("Unsupported timed event "
                                                + std::to_string(static_cast<uint>(timed_event.type))
                                                + ".");
        }
    }
}

CellEvent TissueSimulation::select_next_cell_event()
{
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

    CellEvent event;
    event.delay = std::numeric_limits<Time>::max();

    for (const Species& species: tissue()) {
        select_next_event_in_species(event, tissue(), species, uni_dist, random_gen);
    }

    if (event.delay == std::numeric_limits<Time>::max()) {
        throw Error<std::runtime_error>("No event available for the selection.");
    }

    return event;
}

CellEvent TissueSimulation::select_next_event()
{
    // the tissue() call checks whether a tissue has been
    // associated to the simulation and, if this is not the
    // case, it throws an std::runtime_error
    (void)tissue();

    CellEvent event = select_next_cell_event();

    handle_timed_event_queue(event);

    return event;
}

template<typename T>
inline
T get_min(const T v1, const size_t& v2)
{
    return std::min(v1, static_cast<T>(v2));
}

void enable_duplication_on_neighborhood_externals(Tissue& tissue, const PositionInTissue& position)
{
    auto sizes = tissue.size();
    PositionInTissue pos;
    pos.x = (position.x>0?position.x-1:0);
    auto max_x = get_min(position.x+2, sizes[0]);
    for (; pos.x < max_x; ++pos.x) {
        pos.y = (position.y>0?position.y-1:0);
        auto max_y = get_min(position.y+2, sizes[1]);
        for (; pos.y < max_y; ++pos.y) {
            pos.z = (position.z>0?position.z-1:0);
            auto max_z = (sizes.size()==2?1:get_min(position.z+2, sizes[2]));
            for (; pos.z < max_z; ++pos.z) {
                Tissue::CellInTissueProxy cell_in_tissue = tissue(pos);
                if (!cell_in_tissue.is_wild_type() && cell_in_tissue.is_on_border()) {
                    cell_in_tissue.enable_duplication();
                }
            }
        }
    }
}

typename TissueSimulation::EventAffectedCells
TissueSimulation::simulate_death(const Position& position)
{
    auto cell = (*(position.tissue))(position);

    if (cell.is_wild_type()) {
        return {{},{}};
    }

    TissueSimulation::EventAffectedCells affected = {{cell.copy_and_erase()},{}};

    enable_duplication_on_neighborhood_externals(*(position.tissue), position);

    return affected;
}

template<typename GENERATOR>
inline const Direction& select_push_direction(GENERATOR& random_gen,
                                              const Tissue& tissue,
                                              const PositionInTissue& position,
                                              const std::vector<Direction>& directions)
{
    std::vector<double> cells_to_push;
    double total{0.0};
    for (const auto& direction : directions) {
        size_t num_of_cells = tissue.count_mutated_cells_from(position, direction);
        cells_to_push.push_back(static_cast<double>(1.0)/num_of_cells);
        total += cells_to_push.back();
    }

    std::uniform_real_distribution<double> distribution(0,total);

    auto selected = distribution(random_gen);

    total = 0.0;
    auto it = cells_to_push.begin();
    for (const auto& direction : directions) {
        total += *it;
        ++it;

        if (total >= selected) {
            return direction;
        }
    }

    return directions.back();
}

template<typename GENERATOR>
inline const Direction& select_inverse_min_direction(GENERATOR& random_gen,
                                                     const Tissue& tissue,
                                                     const PositionInTissue& position,
                                                     const std::vector<Direction>& directions)
{
    std::list<size_t> cells_to_push;
    size_t total = 0;
    for (const auto& direction : directions) {
        cells_to_push.push_back(tissue.count_mutated_cells_from(position, direction));
        total += cells_to_push.back();
    }

    std::vector<double> ratios;
    ratios.reserve(directions.size());

    double sum{0};
    for (const auto& count : cells_to_push) {
        ratios.push_back(static_cast<double>(total)/count);
        sum += ratios.back();
    }
    for (auto& ratio : ratios) {
        ratio /= sum;
    }

    std::uniform_real_distribution<double> distribution(0,1);

    auto selected = distribution(random_gen);

    auto dir_it = directions.begin();
    sum = 0;
    for (const auto& ratio : ratios) {
        sum += ratio;
        if (sum >= selected) {
            return *dir_it;
        }

        ++dir_it;
    }

    return directions.back();
}

template<typename GENERATOR>
inline const Direction& select_min_push_direction(GENERATOR& random_gen,
                                                  const Tissue& tissue,
                                                  const PositionInTissue& position,
                                                  const std::vector<Direction>& directions)
{
    std::list<size_t> cells_to_push;
    size_t min_cells_to_push = std::numeric_limits<size_t>::max();
    uint8_t num_of_mins = 0;
    for (const auto& direction : directions) {
        cells_to_push.push_back(tissue.count_mutated_cells_from(position, direction));
        if (min_cells_to_push>cells_to_push.back()) {
            num_of_mins = 1;
            min_cells_to_push = cells_to_push.back();
        } else {
            if (min_cells_to_push==cells_to_push.back()) {
                ++num_of_mins;
            }
        }
    }

    std::uniform_int_distribution<uint8_t> distribution(1,num_of_mins);

    auto selected = distribution(random_gen);

    auto it = cells_to_push.begin();
    uint8_t counter = 0;
    for (const auto& direction : directions) {
        if (*it==min_cells_to_push) {
            if (++counter==selected) {
                return direction;
            }
        }
        ++it;
    }

    throw Error<std::runtime_error>("No direction has not been selected.");
}

template<typename GENERATOR, typename T>
inline const T& select_random_value(GENERATOR& random_gen, const std::vector<T>& values)
{
    std::uniform_int_distribution<size_t> distribution(0,values.size()-1);

    return values[distribution(random_gen)];
}

typename TissueSimulation::EventAffectedCells
TissueSimulation::simulate_mutation(const Position& position, const SpeciesId& final_id)
{
    Tissue& tissue = *(position.tissue);

    Cell parent_cell = tissue(position);

    parent_cell.species_id = final_id;

    tissue(position) = parent_cell.generate_descendent(next_cell_id, time);
    ++next_cell_id;

    return {{tissue(position)},{}};
}


void switch_duplication_on_neighborhood_internals(Tissue& tissue, const PositionInTissue& position)
{
    for (const auto pos : tissue.get_neighborhood_positions(position)) {
        Tissue::CellInTissueProxy cell_in_tissue = tissue(pos);
        if (!cell_in_tissue.is_wild_type()) {
            cell_in_tissue.switch_duplication(cell_in_tissue.is_on_border());
         }
     }
}

typename TissueSimulation::EventAffectedCells
TissueSimulation::simulate_duplication(const Position& position)
{
    Tissue& tissue = *(position.tissue);
    EventAffectedCells affected;

    if (tissue(position).is_wild_type()) {
        return affected;
    }

    Cell parent_cell = tissue(position);

    // push the cell in position towards a random direction
    const Direction& push_dir = select_inverse_min_direction(random_gen, tissue, position,
                                                             valid_directions);

    affected.lost_cells = tissue.push_cells(position, push_dir,
                                            !duplicate_internal_cells);

    Tissue::CellInTissueProxy cell_in_tissue = tissue(position);

    cell_in_tissue = parent_cell.generate_descendent(next_cell_id, time);
    ++next_cell_id;

    affected.new_cells.push_back(cell_in_tissue);

    PositionInTissue new_cell_position{position + PositionDelta(push_dir)};
    if (tissue.is_valid(new_cell_position)) {
        Tissue::CellInTissueProxy cell_in_tissue = tissue(new_cell_position);

        cell_in_tissue = parent_cell.generate_descendent(next_cell_id, time);
        ++next_cell_id;

        if (!duplicate_internal_cells) {
            switch_duplication_on_neighborhood_internals(tissue, new_cell_position);
        }

        affected.new_cells.push_back(cell_in_tissue);
    }

    return affected;
}

bool TissueSimulation::set_cell_species(const Position& position, const SpeciesId& species_id)
{
    Tissue& tissue = *(position.tissue);

    if (tissue(position).is_wild_type()) {
        return false;
    }

    Cell cell = tissue(position);

    if (cell.get_species_id() != species_id) {
        if (!lineage_graph.has_edge(cell.get_species_id(), species_id)) {
            lineage_graph.add_edge(cell.get_species_id(), species_id, time);
        }

        cell.species_id = species_id;

        tissue(position) = cell;

        return true;
    }

    return false;
}

typename TissueSimulation::EventAffectedCells
TissueSimulation::simulate_duplication_and_set_one_child_species(const Position& position, const SpeciesId& final_id)
{
    auto affected = simulate_duplication(position);

    if (affected.new_cells.size()==0) {
        return affected;
    }

    if (set_cell_species(position, final_id)) {
        Tissue& tissue = *(position.tissue);

        for (auto& cell : affected.new_cells) {
            cell = static_cast<Cell>(tissue(cell));
        }
    }

    return affected;
}

TissueSimulation&
TissueSimulation::schedule_mutation(const MutantProperties& src, const MutantProperties& dst, const Time time)
{
    // the tissue() call checks whether a tissue has been
    // associated to the simulation and, if this is not the
    // case, it throws an std::runtime_error
    (void)tissue();

    Mutation mutation(src, dst);

    timed_event_queue.emplace(time, mutation);

    return *this;
}

const MutantId& TissueSimulation::find_mutant_id(const std::string& mutant_name) const
{
    auto found_mutant = mutant_name2id.find(mutant_name);
    if (found_mutant == mutant_name2id.end()) {
        throw Error<std::out_of_range>("Unknown mutant name \"" + mutant_name + "\".");
    }

    return found_mutant->second;
}

const std::string& TissueSimulation::find_mutant_name(const MutantId& mutant_id) const
{
    for (const auto& [name, id] : mutant_name2id) {
        if (id == mutant_id) {
            return name;
        }
    }

    throw Error<std::out_of_range>("Unknown mutant id "+std::to_string(mutant_id) + ".");
}

TissueSimulation&
TissueSimulation::schedule_mutation(const std::string& src, const std::string& dst, const Time time)
{
    // the tissue() call checks whether a tissue has been
    // associated to the simulation and, if this is not the
    // case, it throws an std::runtime_error
    (void)tissue();

    auto src_id = find_mutant_id(src);
    auto dst_id = find_mutant_id(dst);

    Mutation mutation(src_id, dst_id);

    timed_event_queue.emplace(time, mutation);

    return *this;
}

TissueSimulation&
TissueSimulation::schedule_timed_event(const TimedEvent& timed_event)
{
    // the tissue() call checks whether a tissue has been
    // associated to the simulation and, if this is not the
    // case, it throws an std::runtime_error
    (void)tissue();

    timed_event_queue.push(timed_event);

    return *this;
}

void TissueSimulation::init_valid_directions()
{
    if (tissues.size()==0) {
        return;
    }

    valid_directions.clear();
    DirectionGenerator directions(Direction::X_UP, tissue().num_of_dimensions());
    for (const auto dir : directions) {
        valid_directions.push_back(dir);
    }
}

void TissueSimulation::reset()
{
    tissues.clear();
    added_cells.clear();

    statistics = TissueStatistics();
    time = 0;

    next_cell_id = Cell::first_tumour_cell_id();

    logger = BinaryLogger(logger.get_directory());

    //death_activation_level = 1;
}

TissueSimulation& TissueSimulation::add_mutant(const std::string& mutant_name)
{
    MutantProperties properties{mutant_name};

    return add_mutant(properties);
}

TissueSimulation& TissueSimulation::add_mutant(const MutantProperties& mutant)
{
    if (mutant_name2id.count(mutant.get_name())>0) {
        throw Error<std::runtime_error>("Mutant \"" + mutant.get_name()
                                        + "\" already in the simulation");
    }

    tissue().add_mutant(mutant);

    mutant_name2id[mutant.get_name()] = mutant.get_id();

    return *this;
}

TissueSimulation& TissueSimulation::add_mutants(const std::list<std::string>& mutant_names)
{
    for (const auto& name: mutant_names) {
        add_mutant(name);
    }

    return *this;
}

TissueSimulation& TissueSimulation::add_epigenetic_state(const std::string& epistate_name)
{
    tissue().add_epigenetic_state(epistate_name);

    return *this;
}

TissueSimulation& TissueSimulation::add_epigenetic_states(const std::list<std::string>& epistate_names)
{
    for (const auto& name: epistate_names) {
        add_epigenetic_state(name);
    }

    return *this;
}

std::list<std::string> TissueSimulation::get_mutant_names() const
{
    std::list<std::string> mutant_names;

    for (const auto& [name, id]: mutant_name2id) {
        mutant_names.push_back(name);
    }

    return mutant_names;
}

TissueSimulation& TissueSimulation::place_cell(const std::string& species_name,
                                               const PositionInTissue& position)
{
    if (!knowns_species(species_name)) {
        throw Error<std::runtime_error>("\"" + species_name
                                        + "\" has not been added to the simulation.");
    }

    SpeciesName name(species_name);

    SpeciesProperties s_prop(name.get_mutant_name(), name.get_epistate_name());

    return place_cell(s_prop.get_id(), position);
}

TissueSimulation& TissueSimulation::place_cell(const SpeciesId& species_id, const PositionInTissue& position)
{
    const auto& cell = tissue().place_cell(next_cell_id, species_id, position);

    if (storage_enabled) {
        logger.record_initial_cell(cell);
    }

    ++next_cell_id;

    added_cells.emplace_back(species_id, position, time);

    LineageEdge edge{WILD_TYPE_SPECIES, species_id};
    if (!lineage_graph.has_edge(edge)) {
        lineage_graph.add_edge(edge, time);
    }

    statistics.record_placed_cell(species_id, time);

    return *this;
}

TissueSimulation& TissueSimulation::set_tissue(const std::string& name, const std::vector<AxisSize>& sizes)
{
    if (time>0) {
        std::runtime_error("The tissue properties can be exclusively changed before"
                           "the simulation beginning");
    }

    if (tissues.size()>0) {
        reset();
    }

    tissues.emplace_back(name, sizes);

    init_valid_directions();

    return *this;
}

void TissueSimulation::log_initial_cells()
{
    if (storage_enabled) {
        for (const auto& species: tissue()) {
            for (const auto& cell: species) {
                logger.record_initial_cell(cell);
            }
        }
    }
}

std::vector<Tissue::CellInTissueConstantProxy>
TissueSimulation::collect_cell_proxies_in(const CLONES::Mutants::RectangleSet& rectangle) const
{
    std::vector<Tissue::CellInTissueConstantProxy> proxies;

    for (const auto& position: rectangle) {
        if (tissue().is_valid(position)) {
            auto proxy = tissue()(position);
            if (!proxy.is_wild_type()) {
                proxies.push_back(proxy);
            }
        }
    }

    return proxies;
}

std::vector<Tissue::CellInTissueProxy>
TissueSimulation::collect_cell_proxies_in(const CLONES::Mutants::RectangleSet& rectangle)
{
    std::vector<Tissue::CellInTissueProxy> proxies;

    for (const auto& position: rectangle) {
        if (tissue().is_valid(position)) {
            auto proxy = tissue()(position);
            if (!proxy.is_wild_type()) {
                proxies.push_back(proxy);
            }
        }
    }

    return proxies;
}

TissueSample
TissueSimulation::simulate_sampling(const SampleSpecification& specification)
{
    auto proxies = static_cast<const TissueSimulation*>(this)->collect_cell_proxies_in(specification.get_bounding_box());

    size_t num_of_proxies = proxies.size();

    TissueSample sample(specification.get_name(), time, specification.get_bounding_box(),
                        num_of_proxies);
    size_t to_collect = std::min(specification.get_num_of_cells(), num_of_proxies);
    while (to_collect>0) {
        std::uniform_int_distribution<size_t> dist(0, --num_of_proxies);

        size_t pos = dist(random_gen);

        auto& cell_proxy = proxies[pos];
        auto cell_in_tissue = static_cast<const CellInTissue&>(cell_proxy);
        sample.add_cell_id(cell_in_tissue.get_id());

        cell_proxy = proxies[num_of_proxies];
        --to_collect;
    }

    return sample;
}

TissueSample
TissueSimulation::sample_tissue(const SampleSpecification& specification)
{
    if (name2sample.count(specification.get_name())>0) {
        throw Error<std::runtime_error>("Sample name \"" + specification.get_name()
                                        + "\" has been already used.");
    }

    auto proxies = collect_cell_proxies_in(specification.get_bounding_box());
    size_t num_of_proxies = proxies.size();

    TissueSample sample(specification.get_name(), time, specification.get_bounding_box(),
                        num_of_proxies);
    size_t to_collect = std::min(specification.get_num_of_cells(), num_of_proxies);
    while (to_collect>0) {
        std::uniform_int_distribution<size_t> dist(0, --num_of_proxies);

        size_t pos = dist(random_gen);

        auto& cell_proxy = proxies[pos];
        auto cell_in_tissue = static_cast<const CellInTissue&>(cell_proxy);
        sample.add_cell_id(cell_in_tissue.get_id());
        ++(statistics[cell_in_tissue.get_species_id()].lost_cells);
        cell_proxy.erase();

        cell_proxy = proxies[num_of_proxies];
        --to_collect;
    }

    samples.push_back(sample);
    name2sample[specification.get_name()]=(samples.end()--);

    make_snapshot<CLONES::UI::ProgressBar>(nullptr);

    return sample;
}

double TissueSimulation::evaluate(const Logics::Variable& variable) const
{
    const auto& species_id = variable.get_species_id();
    switch(variable.get_type()) {
        case Logics::Variable::Type::CARDINALITY:
            return tissue().get_species(species_id).num_of_cells();
        case Logics::Variable::Type::EVENT:
            {
                const auto& t_stats = get_statistics();

                return t_stats.count_fired_events(species_id,
                                                  variable.event_type);
            }
        case Logics::Variable::Type::TIME:
            return get_time();
        default:
            throw Error<std::runtime_error>("Unsupported variable type "
                                            + std::to_string(static_cast<size_t>(variable.get_type()))
                                            + ".");
    }
}

TissueSimulation::~TissueSimulation()
{
}

}   // Evolutions

}   // Mutants

}   // CLONES
