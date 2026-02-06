/**
 * @file species.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines species representation
 * @version 1.2
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

#ifndef __CLONES_SPECIES__
#define __CLONES_SPECIES__

#include <iostream>

#include <vector>
#include <random>
#include <cstddef> // std::ptrdiff_t

#include "imap.hpp"
#include "archive.hpp"
#include "time.hpp"
#include "cell.hpp"
#include "mutant_properties.hpp"


namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

class Tissue;

/**
 * @brief Cell species
 *
 * This class represents the set of cells in a species.
 */
class Species: public SpeciesProperties
{
    /**
     * @brief A map from cell ids to the corresponding pointers
     */
    using CellIdToCell = imap<CellId, CellInTissue*>;

    CellIdToCell cells;                 //!< species cells sorted by birth time/id
    CellIdToCell duplication_enabled;   //!< species cells that are duplication enabled

    Time last_insertion_time;     //!< the last insertion time

    size_t simulated_cells;          //!< Total number of cells along the computation

    /**
     * @brief A constant iterator for the species cells
     */
    class const_iterator {
        CellIdToCell::const_iterator it; //!< the cell pointer vector iterator

        const_iterator(const CellIdToCell::const_iterator it);
    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   CellInTissue;
        using pointer           =   const CellInTissue*;
        using reference         =   const CellInTissue&;
        using iterator_category =   std::bidirectional_iterator_tag;

        const_iterator();

        inline reference operator*() const
        {
            return *(it->second);
        }

        inline pointer operator->() const
        {
            return it->second;
        }

        // Prefix increment
        inline const_iterator& operator++()
        {
            ++it;
            return *this;
        }

        // Postfix increment
        const_iterator operator++(int);

        // Prefix decrement
        inline const_iterator& operator--()
        {
            --it;
            return *this;
        }

        // Postfix decrement
        const_iterator operator--(int);

        friend inline bool operator==(const const_iterator& a, const const_iterator& b)
        {
            return a.it == b.it;
        }

        friend inline bool operator!=(const const_iterator& a, const const_iterator& b)
        {
            return a.it != b.it;
        }

        friend class Species;
    };

    /**
     * @brief Compute the cell age from the last addition to the species
     *
     * @param cell
     * @return Time
     */
    inline Time age(const CellInTissue *cell) const
    {
        return last_insertion_time-cell->get_birth_time();
    }

    /**
     * @brief Add a cell to the species
     *
     * This method add the cell pointer by the parameter to the
     * species. The species representation uses the memory pointed
     * by the parameter.
     *
     * @param cell is a pointer to the cell to be added to the species
     * @return pointer to the added cell
     */
    CellInTissue* add(CellInTissue* cell);

    /**
     * @brief Randomly select a cell among those in a id-cell pointer map
     *
     * This method select a cell by using an exponential over
     * cell ages.
     *
     * @tparam GENERATOR is a random number generator type
     * @param generator is the random number generator used to select
     *      the cell
     * @param id2cell_ptrs is the map id-cell pointers among which the
     *      cell must be selected
     * @return a non-constant reference to a randomly selected
     *       cell in the species
     */
    template<typename GENERATOR>
    const CellInTissue& choose_a_cell(GENERATOR& generator,
                                      const CellIdToCell& id2cell_ptrs) const
    {
        if (id2cell_ptrs.size()==0) {
            throw std::domain_error("No cells can be selected");
        }

        std::uniform_int_distribution<size_t> dist(0, id2cell_ptrs.size()-1);

        const auto pos = dist(generator);
        return *(id2cell_ptrs.get(pos)->second);
    }

    /**
     * @brief An empty constructor
     */
    Species();

    /**
     * @brief Reset the species
     *
     * This method resets the species and removes and deletes all
     * its cells.
     */
    void reset();
public:
    /**
     * @brief A constructor
     *
     * @param species_properties is the species of the species
     */
    explicit Species(const SpeciesProperties& species_properties);

    /**
     * @brief A copy constructor
     *
     * @param orig is the template object
     */
    Species(const Species& orig);

    /**
     * @brief The copy operator
     *
     * @param orig is the original object
     * @return a reference to the updated object
     */
    Species& operator=(const Species& orig);

    /**
     * @brief The replace operator
     *
     * @param orig is the original object
     * @return a reference to the updated object
     */
    Species& operator=(Species&& orig);

    /**
     * @brief Get the initial iterator for the species cells
     *
     * @return the initial iterator for the species cells
     */
    const_iterator begin() const;

    /**
     * @brief Get the final iterator for the species cells
     *
     * @return the final iterator for the species cells
     */
    const_iterator end() const;

    /**
     * @brief Get the number of cells in the species
     *
     * @return the number of cells in the species
     */
    inline size_t num_of_cells() const
    {
        return cells.size();
    }

    /**
     * @brief Get the total number of cells in the species along the whole computation
     *
     * @return the total number of cells in the species along the whole computation
     */
    inline const size_t& num_of_simulated_cells() const
    {
        return simulated_cells;
    }

    /**
     * @brief Test whether a cell in the species is avaiable for an event
     *
     * @param cell_id is the id of the cell to be tested
     * @param event_type is the event type for which the availability is tested
     * @return `true` if and only if the cell having id `cell_id` belong to the
     *          species and it is available for the the event `event_type`
     */
    inline bool cell_is_available_for(const CellId& cell_id,
                                      const CellEventType event_type) const
    {
        if (cells.count(cell_id)==0) {
            throw std::domain_error("The cell does not belong to the species");
        }

        switch (event_type) {
            case CellEventType::DEATH:
                return true;
            case CellEventType::DUPLICATION:
            case CellEventType::MUTATION:
            case CellEventType::EPIGENETIC_SWITCH:
                return (duplication_enabled.count(cell_id)>0);
            case CellEventType::ANY:
                return true;
            default:
                throw std::domain_error("Unsupported event type");
        }
    }

    /**
     * @brief Get the number of cells available for an event
     *
     * @param event_type is the event for which the cells are needed
     * @return the number of cells available for `event_type`
     */
    size_t num_of_cells_available_for(const CellEventType& event_type) const;

    /**
     * @brief Randomly select a cell for a specific event type
     *
     * This method select a cell by using an exponential over
     * cell ages.
     *
     * @tparam GENERATOR is a random number generator type
     * @param generator is the random number generator used to select
     *      the cell
     * @param event_type is the type of the event that will occurs
     *      to the selected cell
     * @return a non-constant reference to a randomly selected
     *       cell in the species
     */
    template<typename GENERATOR>
    const CellInTissue& choose_a_cell(GENERATOR& generator,
                                      const CellEventType event_type) const
    {
        switch (event_type) {
            case CellEventType::DEATH:
                return choose_a_cell(generator, cells);
            case CellEventType::DUPLICATION:
            case CellEventType::MUTATION:
            case CellEventType::EPIGENETIC_SWITCH:
                return choose_a_cell(generator, duplication_enabled);
            case CellEventType::ANY:
                return choose_a_cell(generator, cells);
            default:
                throw std::domain_error("Unsupported event type");
        }
    }

    /**
     * @brief Randomly select a cell in the species
     *
     * This method select a cell by using an exponential over
     * cell ages.
     *
     * @tparam GENERATOR is a random number generator type
     * @param generator is the random number generator used to select
     *      the cell
     * @return a non-constant reference to a randomly selected
     *       cell in the species
     */
    template<typename GENERATOR>
    inline const CellInTissue& choose_a_cell(GENERATOR& generator) const
    {
        return choose_a_cell(generator, cells);
    }

    /**
     * @brief Remove a cell from the species and delete it
     *
     * @param cell_id is the id of the cell to be deleted
     */
    void erase(const CellId& cell_id);

    /**
     * @brief Add a cell to the species
     *
     * By default the duplication of the cell is switched on.
     *
     * @param cell is the cell to be added to the species
     * @return pointer to the added cell
     */
    inline CellInTissue* add(CellInTissue&& cell)
    {
        return add(cell);
    }

    /**
     * @brief Add a cell to the species
     *
     * @param cell is the cell to be added to the species
     * @return pointer to the added cell
     */
    CellInTissue* add(CellInTissue& cell);

    /**
     * @brief Switch on/off the duplication of a cell
     *
     * @param cell_id is the identifier of the cell on which
     *          the duplication will be enabled/disabled
     * @param duplication_on is a Boolean flag to enable
     *          (True)/disable(False) duplication on `cell`
     */
    void switch_duplication_for(const CellId& cell_id,
                                const bool duplication_on);

    /**
     * @brief Enable the duplication of a cell
     *
     * @param cell_id is the identifier of the cell to be enabled
     *          for the duplication
     */
    void enable_duplication_for(const CellId& cell_id);

    /**
     * @brief Disable the duplication of a cell
     *
     * @param cell_id is the identifier of the cell to be disabled
     *          for the duplication
     */
    void disable_duplication_for(const CellId& cell_id);

    /**
     * @brief Get cell by identifier
     *
     * @param cell_id is the identifier of the aimed cell
     * @return a non-constant reference to the aimed cell
     */
    CellInTissue& operator()(const CellId& cell_id);

    /**
     * @brief Get cell by identifier
     *
     * @param cell_id is the identifier of the aimed cell
     * @return a constant reference to the aimed cell
     */
    const CellInTissue& operator()(const CellId& cell_id) const;

    /**
     * @brief Test whether a cell is already contained in a species
     *
     * @param cell_id is the identifier of the cell to be tested
     * @return `true` is an only if the cell is already in the species
     */
    inline bool contains(const CellId& cell_id) const
    {
        return cells.find(cell_id) != cells.end();
    }

    /**
     * @brief The destroyer
     */
    ~Species();

    /**
     * @brief Save a species in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "CLONES Species", 0);

        archive & static_cast<const SpeciesProperties &>(*this);

        archive & last_insertion_time
                & simulated_cells;

        // save species cells
        archive & cells.size();
        for (const auto& [cell_id, cell_ptr]: cells) {
            archive & *cell_ptr;
        }

        // save duplicated ids
        archive & duplication_enabled.size();
        for (const auto& [cell_id, cell_ptr]: duplication_enabled) {
            archive & cell_id;
        }
    }

    /**
     * @brief Load a species from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded species
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static Species load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "CLONES Species", 0);

        auto species_properties = SpeciesProperties::load(archive);
        Species species(species_properties);

        archive & species.last_insertion_time
                & species.simulated_cells;

        size_t num_of_cells;

        // load species cells
        archive & num_of_cells;
        for (size_t i=0; i<num_of_cells; ++i) {
            CellInTissue *cell = new CellInTissue();

            archive & *cell;

            species.cells.insert({cell->get_id(), cell});
        }

        // load duplicated enabled ids
        archive & num_of_cells;
        for (size_t i=0; i<num_of_cells; ++i) {
            CellId cell_id;

            archive & cell_id;

            species.enable_duplication_for(cell_id);
        }

        return species;
    }

    friend class Simulation;

    friend void swap(Species& a, Species& b);
};


/**
 * @brief Swap two species
 *
 * @param a is a species
 * @param b is a species
 */
void swap(Species& a, Species& b);

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_SPECIES__
