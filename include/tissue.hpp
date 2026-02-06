/**
 * @file tissue.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines tissue class
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

#ifndef __CLONES_TISSUE__
#define __CLONES_TISSUE__

#include <vector>
#include <list>
#include <map>
#include <set>
#include <string>
#include <memory> // SpeciesView::const_iterator

#include "archive.hpp"
#include "time.hpp"
#include "species.hpp"
#include "cell.hpp"
#include "logics.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

/**
 * @brief The type of tissue axis sizes
 */
using AxisSize = uint16_t;

class Simulation;

/**
 * @brief A class to represent tissues
 */
class Tissue {
    using ClonePosition = std::map<MutantId, std::vector<size_t>>;

    std::string name;                       //!< The tissue name
    std::vector<Species> species;           //!< The species in the tissue
    std::map<SpeciesId, size_t> id_pos;     //!< The identifier to position map
    std::map<std::string, size_t> name_pos; //!< The name to position map
    ClonePosition mutant_pos;     //!< The positions of the species associated to the same mutant

    uint8_t dimensions;   //!< The number of space dimension for the tissue

    using TissueSpace = std::vector<std::vector<std::vector<CellInTissue*>>>;

    TissueSpace space;     //!< Space in the tissue

    /**
     * @brief Get the pointer to the cell in a position
     *
     * This method returns a pointer to a cell in a given position.
     * The returned value is undefined when the position is not
     * a valid position for the tissue.
     *
     * @param position is the position of the aimed cell pointer
     * @return a non-constant reference to a pointer to the cell
     *          in `position`
     */
    inline CellInTissue*& cell_pointer(const PositionInTissue& position)
    {
        return this->space[position.x][position.y][position.z];
    }

    /**
     * @brief Get the pointer to the cell in a position
     *
     * This method returns a pointer to a cell in a given position.
     * The returned value is undefined when the position is not
     * a valid position for the tissue.
     *
     * @param position is the position of the aimed cell pointer
     * @return a constant reference to a pointer to the cell
     *          in `position`
     */
    inline const CellInTissue* cell_pointer(const PositionInTissue& position) const
    {
        return this->space[position.x][position.y][position.z];
    }

    /**
     * @brief Count the number of cells in a neighborwood in a species
     *
     * @param position is the central position of the neighborwood
     * @param species_id is the identifier of the searched species
     * @return the number of cells in a neighborwood of `position` having
     *      `species_id` as the identifier of their species
     */
    size_t count_neighbors_in(const PositionInTissue position, const SpeciesId& species_id) const;

    /**
     * @brief Register the species cells in the tissue space
     *
     * This method records the species cells in the tissue space by
     * copying each cell pointer in the cell position in the space.
     * This method must be called when either:
     * 1. we want to add the cells of newly added species to the tissue
     * 2. the species vector has been resized and we want to
     *    re-register the cell memory locations
     */
    void register_species_cells();

    /**
     * @brief Add mutant species to the tissue species
     *
     * This method add mutant species to the tissue species, but do not
     * register their cells into tissue space
     *
     * @param mutant_properties is the mutant properties of the mutant
     * @return a reference to the updated object
     */
    Tissue& add_mutant_species(const MutantProperties& mutant_properties);

    /**
     * @brief Place a cell in the tissue
     *
     * @param id is the cell identifier
     * @param species_id is the species identifier of the new cell
     * @param position is the cell position in the tissue
     * @return a constant reference to the new cell in the tissue
     */
    const CellInTissue& place_cell(const CellId& id, const SpeciesId& species_id,
                                   const PositionInTissue position);
public:

    /**
     * @brief A view class for species
     *
     * This class allows to have a partial view of the
     * tissue's species. For instance, it allows to
     * list all the species in the same mutant.
     */
    class SpeciesView {
        const std::vector<Species>& species;       //!< The species vector
        const std::vector<size_t>& species_pos;    //!< The species position vector

        /**
         * @brief A constructor
         *
         * @param species is a reference to the tissue's species vector
         * @param species_pos is a vector of valid positions for `vector`
         */
        SpeciesView(const std::vector<Species>& species, const std::vector<size_t>& species_pos);
    public:

        /**
         * @brief A constant iterator over species views
         */
        class const_iterator {
            const std::vector<Species>* species;    //!< The tissue's species vector
            std::shared_ptr<std::vector<size_t>::const_iterator> it; //!< A constant iterator over a vector of position in `species`

            /**
             * @brief A private constructor
             *
             * @param species is the vector of the tissue's species
             * @param it is a constant iterator over a vector of the positions in `species`
             */
            const_iterator(const std::vector<Species>& species, const std::vector<size_t>::const_iterator it);
        public:
            using difference_type   =   std::ptrdiff_t;
            using value_type        =   Species;
            using pointer           =   const Species*;
            using reference         =   const Species&;
            using iterator_category =   std::random_access_iterator_tag;

            /**
             * @brief An empty construtor
             */
            const_iterator();

            /**
             * @brief Reference operator
             *
             * @return a reference to the species pointer by the iterator
             */
            inline reference operator*() const
            {
                return (*species)[**it];
            }

            /**
             * @brief Pointer operator
             *
             * @return a pointer to the species pointer by the iterator
             */
            inline pointer operator->()
            {
                return &((*species)[**it]);
            }

            /**
             * @brief The prefix increment
             *
             * @return a reference to the updated object
             */
            inline const_iterator& operator++()
            {
                ++(*it);
                return *this;
            }

            /**
             * @brief The postfix increment
             *
             * @return a copy of the original object
             */
            const_iterator operator++(int);

            /**
             * @brief The prefix decrement
             *
             * @return a reference to the updated object
             */
            inline const_iterator& operator--()
            {
                --(*it);
                return *this;
            }

            /**
             * @brief The postfix decrement
             *
             * @return a copy of the original object
             */
            const_iterator operator--(int);

            /**
             * @brief Add operator
             *
             * @param delta is the value to add
             * @return a new iterator that points `delta` position ahead
             *      with respect to the original object
             */
            inline const_iterator operator+(const int delta)
            {
                return const_iterator(*species, *it + delta);
            }

            /**
             * @brief Subtract operator
             *
             * @param delta is the value to subtract
             * @return a new iterator that points `delta` position backwards
             *      with respect to the original object
             */
            inline const_iterator operator-(const int delta)
            {
                return const_iterator(*species, *it - delta);
            }

            /**
             * @brief Inplace add operator
             *
             * @param delta is the value to add
             * @return a reference to the update object
             */
            inline const_iterator& operator+=(const int& delta) {
                (*it) += delta;

                return *this;
            }

            /**
             * @brief Inplace subtract operator
             *
             * @param delta is the value to subtract
             * @return a reference to the update object
             */
            inline const_iterator& operator-=(const int& delta) {
                (*it) -= delta;

                return *this;
            }

            /**
             * @brief Index operator
             *
             * @param index is the index to access
             * @return a reference to the `index`-th objects after
             *      that pointer by the current iterator
             */
            inline reference operator[](const int& index) const
            {
                return (*species)[(*it)[index]];
            }

            /**
             * @brief Test whether two iterators are the same
             *
             * @param a is the first iterator to compare
             * @param b is the second iterator to compare
             * @return `true` if and only if the two iterators
             *      refer to the same object
             */
            friend inline bool operator==(const const_iterator& a, const const_iterator& b)
            {
                return (*a.it == *b.it) && (a.species == b.species);
            }

            friend class Tissue::SpeciesView;
        };

        /**
         * @brief Index operator
         *
         * @param index is the index to access
         * @return a constant reference to the `index`-th species
         *      in the view
         */
        inline const Species& operator[](const size_t& index) const
        {
            return species[species_pos[index]];
        }

        /**
         * @brief Get the view size
         *
         * @return the view size
         */
        inline size_t size() const
        {
            return species_pos.size();
        }

        /**
         * @brief Get the view begin
         *
         * @return a constant iterator to the first element
         *      in the view
         */
        inline const_iterator begin() const
        {
            return const_iterator(species, species_pos.begin());
        }

        /**
         * @brief Get the view end
         *
         * @return a constant iterator to the view end
         */
        inline const_iterator end() const
        {
            return const_iterator(species, species_pos.end());
        }

        /**
         * @brief Get the numer of cells in the species view
         *
         * @return the number of cells in the species view
         */
        size_t num_of_cells() const;

        friend class Tissue;
    };

    /**
     * @brief This class wraps pointer to constant cells in tissue space
     */
    template<typename TISSUE_TYPE>
    class BaseCellInTissueProxy
    {
    protected:
        TISSUE_TYPE *tissue;          //!< tissue
        PositionInTissue position;    //!< position of the cell

        BaseCellInTissueProxy(TISSUE_TYPE &tissue, const PositionInTissue position):
            tissue(&tissue), position(position)
        {
            if (!tissue.is_valid(position)) {
                std::ostringstream oss;

                oss << "The position (" << position.x << "," << position.y;

                if (tissue.num_of_dimensions()==3) {
                    oss << "," << position.z;
                }

                oss << ") does not belong to the tissue";

                throw std::out_of_range(oss.str());
            }
        }
    public:

        /**
         * @brief Test whether the referenced cell is of wild-type cell
         *
         * @return `true` if and only if the referenced cell is of wild-type cell
         */
        inline bool is_wild_type() const
        {
            return tissue->cell_pointer(position)==nullptr;
        }

        /**
         * @brief Test whether the referenced cell is available for an event
         *
         * @param event_type is the event type for which the avaiability is tested
         * @return `true` if and only if the referenced cell is available for an event
         */
        bool is_available_for(const CellEventType& event_type) const
        {
            if (is_wild_type()) {
                return false;
            }

            const CellInTissue& cell = *(tissue->cell_pointer(position));
            const Species& species = tissue->species[tissue->id_pos.at(cell.get_species_id())];

            return species.cell_is_available_for(cell.get_id(), event_type);
        }

        /**
         * @brief Check whether a cell is on the cancer edge
         *
         * @param species_border a Boolean flag
         * @return `true` if and only if either one of the
         *      eight neighbor cells is a wild-type cell and
         *      `species_border` is `false` or one of the
         *      eight neighbor cells and the pointed cell
         *      belong to different species and
         *      `species_border` is `true`.
         */
        bool is_on_border(const bool species_border=false) const
        {
            if (is_wild_type()) {
                return false;
            }

            auto sizes = tissue->size();

            const SpeciesId species_id = tissue->cell_pointer(position)->get_species_id();
            PositionInTissue pos;
            pos.x = (position.x>0?position.x-1:0);
            for (; (pos.x < position.x+2 &&  pos.x < sizes[0]); ++pos.x) {
                pos.y = (position.y>0?position.y-1:0);
                for (; (pos.y < position.y+2 &&  pos.y < sizes[1]); ++pos.y) {
                    pos.z = (position.z>0?position.z-1:0);
                    for (; ((sizes.size()==3 && pos.z < position.z+2 &&  pos.z < sizes[2])
                            || (sizes.size()==2 && pos.z==0)); ++pos.z) {
                        const auto* cell_ptr = tissue->cell_pointer(pos);
                        if (cell_ptr==nullptr) {
                            return true;
                        }
                        if (species_border && cell_ptr->get_species_id()!=species_id) {
                            return true;
                        }
                    }
                }
            }

            return false;
        }

        /**
         * @brief Get a constant reference to the referenced cell
         *
         * This method returns a constant reference of the referenced cell. When
         * the referenced cell is a wild-type cell, the method throws
         * a `std::runtime_error`.
         *
         * @return a constant reference of the referenced cell
         * @throws `std::runtime_error` if the referenced cell is a wild-type cell
         */
        operator const CellInTissue&() const
        {
            const auto ptr = tissue->cell_pointer(position);

            if (ptr!=nullptr) {
                return *ptr;
            }

            throw std::runtime_error("Wild-type cell");
        }

        friend class Tissue;
    };

    /**
     * @brief This class wraps pointer to constant cells in tissue space
     */
    class CellInTissueConstantProxy : public BaseCellInTissueProxy<const Tissue>
    {
        CellInTissueConstantProxy(const Tissue &tissue, const PositionInTissue position);

        friend class Tissue;
    };

    /**
     * @brief This class wraps pointer to cells in tissue space
     */
    class CellInTissueProxy : public BaseCellInTissueProxy<Tissue>
    {
        CellInTissueProxy(Tissue &tissue, const PositionInTissue position);
    public:
        CellInTissueProxy& operator=(const Cell& cell);

        operator CellInTissue&();

        void erase();

        CellInTissue copy_and_erase();

        void switch_duplication(const bool duplication_on);

        inline void enable_duplication()
        {
            switch_duplication(true);
        }

        void disable_duplication()
        {
            switch_duplication(false);
        }

        friend class Tissue;
    };

    /**
     * @brief A constructor
     *
     * @param sizes are the sizes of the tissue
     */
    explicit Tissue(const std::vector<AxisSize>& sizes);

    /**
     * @brief A constructor for a 3D tissue
     *
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const AxisSize x_size, const AxisSize y_size, const AxisSize z_size);

    /**
     * @brief A constructor for a 2D tissue
     *
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const AxisSize x_size, const AxisSize y_size);

    /**
     * @brief A constructor for a 3D tissue
     *
     * @param mutants is the vector of mutants
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::vector<MutantProperties>& mutants, const AxisSize x_size, const AxisSize y_size, const AxisSize z_size);

    /**
     * @brief A constructor for a 2D tissue
     *
     * @param mutants is the vector of mutants
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const std::vector<MutantProperties>& mutants, const AxisSize x_size, const AxisSize y_size);

    /**
     * @brief A constructor
     *
     * @param name is the tissue name
     * @param sizes are the sizes of the tissue
     */
    Tissue(const std::string& name, const std::vector<AxisSize>& sizes);

    /**
     * @brief A constructor for a 3D tissue
     *
     * @param name is the tissue name
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::string& name, const AxisSize x_size, const AxisSize y_size, const AxisSize z_size);

    /**
     * @brief A constructor for a 2D tissue
     *
     * @param name is the tissue name
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const std::string& name, const AxisSize x_size, const AxisSize y_size);

    /**
     * @brief A constructor for a 3D tissue
     *
     * @param name is the tissue name
     * @param mutants is the vector of mutants
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::string& name, const std::vector<MutantProperties>& mutants, const AxisSize x_size, const AxisSize y_size, const AxisSize z_size);

    /**
     * @brief A constructor for a 2D tissue
     *
     * @param name is the tissue name
     * @param mutants is the vector of mutants
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const std::string& name, const std::vector<MutantProperties>& mutants, const AxisSize x_size, const AxisSize y_size);

    /**
     * @brief Get the initial iterator for the tissue species
     *
     * @return the initial iterator for the tissue species
     */
    inline std::vector<Species>::const_iterator begin() const
    {
        return std::begin(species);
    }

    /**
     * @brief Get the final iterator for the tissue species
     *
     * @return the final iterator for the tissue species
     */
    inline std::vector<Species>::const_iterator end() const
    {
        return std::end(species);
    }

    /**
     * @brief Get the tissue species properties
     *
     * @return the vector of the tissue species properties
     */
    std::vector<SpeciesProperties> get_species_properties() const;

    /**
     * @brief Get a tissue species by identifier
     *
     * @param species_id is the species identifier
     * @return a constant reference to the tissue species
     */
    const Species& get_species(const SpeciesId& species_id) const;

    /**
     * @brief Get a tissue species by identifier
     *
     * @param species_id is the species identifier
     * @return a non-constant reference to the tissue species
     */
    Species& get_species(const SpeciesId& species_id);

    /**
     * @brief Get a tissue species by name
     *
     * @param species_name is the species name
     * @return a constant reference to the tissue species
     */
    const Species& get_species(const std::string& species_name) const;

    /**
     * @brief Get a tissue species by name
     *
     * @param species_name is the species name
     * @return a non-constant reference to the tissue species
     */
    Species& get_species(const std::string& species_name);

    /**
     * @brief Add a mutant to the tissue
     *
     * @param mutant_properties is the mutant properties of the mutant
     * @return a reference to the updated object
     */
    Tissue& add_mutant(const MutantProperties& mutant_properties);

    /**
     * @brief Test whether a position is valid in a tissue
     *
     * @param position is the position to be tested
     * @return `true` if and only if `position` is a valid
     *      position for the tissue
     */
    inline bool is_valid(const PositionInTissue& position) const
    {
        return (static_cast<size_t>(position.x)<space.size() &&
                static_cast<size_t>(position.y)<space[0].size() &&
                static_cast<size_t>(position.z)<space[0][0].size() &&
                position.x>=0 && position.y>=0 && position.z>=0);
    }

    /**
     * @brief Get the number of tissue dimensions
     *
     * @return the number of tissue dimensions
     */
    inline size_t num_of_species() const
    {
        return species.size();
    }

    /**
     * @brief Get the number of cells in the simulated tissue
     *
     * This method returns the total number of cells
     * in the simulated space. This number accounts for
     * both mutated cells and normal ones.
     *
     * @return the number of cells in the simulated tissue
     */
    inline size_t num_of_cells() const
    {
        return space.size() * space[0].size() * space[0][0].size();
    }

    /**
     * @brief Get the number of mutated cells in the tissue
     *
     * @return the number of mutated cells in the tissue
     */
    size_t num_of_mutated_cells() const;

    /**
     * @brief Get an iterator over the species having the same mutant
     *
     * @param mutant_id is the identifier of the mutant whose species are aimed
     * @return an interator over the tissue's species having `mutant_id` as
     *       mutant identifier
     */
    inline SpeciesView get_mutant_species(const MutantId& mutant_id) const
    {
        return SpeciesView(species, mutant_pos.at(mutant_id));
    }

    /**
     * @brief Get the cell in a position
     *
     * @param position is the position of the aimed cell
     * @return a cell in tissue proxy
     */
    CellInTissueProxy operator()(const PositionInTissue& position);

    /**
     * @brief Get the cell in a position
     *
     * @param position is the position of the aimed cell
     * @return a constant cell in tissue proxy
     */
    const CellInTissueConstantProxy operator()(const PositionInTissue& position) const;

    /**
     * @brief Get tissue name
     *
     * @return a constant reference to the tissue name
     */
    inline const std::string& get_name() const
    {
        return name;
    }

    /**
     * @brief Get number of dimensions
     *
     * @return a constant reference to the number of dimensions
     */
    inline const uint8_t& num_of_dimensions() const
    {
        return dimensions;
    }

    /**
     * @brief Count the contiguous mutated cells in a direction
     *
     * @param position is the position from which the cells are counted
     * @param direction is the counting direction
     * @return the number of contiguous mutated cells from
     *      `from_position` towards `directions`
     */
    size_t count_mutated_cells_from(const PositionInTissue position,
                                    const Direction& direction) const;

    /**
     * @brief Push contiguous mutated cells in a direction
     *
     * This method pushes the contiguous mutated cells from a position
     * towards a direction and returns the list of the cells that are
     * pushed outside the tissue.
     *
     * @param from_position is the position from which the cells are pushed
     * @param direction is the push direction
     * @return the list of the cells that have been pushed outside the
     *      tissue border
     */
    std::list<Cell> push_cells(const PositionInTissue from_position, const Direction& direction);

    /**
     * @brief Get the tissue size
     *
     * @return the tissue size for the 3 dimensions
     */
    std::vector<AxisSize> size() const;

    /**
     * @brief Get the logic variable of a species cardinality
     *
     * @param species_name is the name of a species
     * @return the logic variable associated to the cardinality of the species
     *      whose name is `species_name`
     */
    inline Logics::Variable get_cardinality_variable(const std::string& species_name) const
    {
        return Logics::Variable(get_species(species_name).get_id(), species_name);
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
        return Logics::Variable(event_type, get_species(species_name).get_id(), species_name);
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
        return Logics::Variable(species_id, get_species(species_id).get_name());
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
        return Logics::Variable(event_type, species_id, get_species(species_id).get_name());
    }

    /**
     * @brief Save a tissue in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & size()
                & name
                & species
                & mutant_pos;
    }

    /**
     * @brief Load a tissue from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded tissue
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static Tissue load(ARCHIVE& archive)
    {
        std::vector<AxisSize> sizes;

        archive & sizes;

        Tissue tissue(sizes);

        archive & tissue.name
                & tissue.species
                & tissue.mutant_pos;

        size_t i=0;
        for (auto& species : tissue.species) {
            tissue.id_pos[species.get_id()] = i;
            tissue.name_pos[species.get_name()] = i;
            ++i;
        }

        tissue.register_species_cells();

        return tissue;
    }

    template<typename TISSUE_TYPE>
    friend class BaseCellInTissueProxy;
    friend class Simulation;
};

/**
 * @brief Test whether two iterators differs
 *
 * @param a is the first iterator to compare
 * @param b is the second iterator to compare
 * @return `true` if and only if the two iterators
 *      do not refer to the same object
 */
inline bool operator!=(const Tissue::SpeciesView::const_iterator& a, const Tissue::SpeciesView::const_iterator& b)
{
    return !(a==b);
}

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_TISSUE__
