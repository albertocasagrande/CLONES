/**
 * @file tissue.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines tissue class
 * @version 1.7
 * @date 2026-06-20
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
#include <memory> // mutant_view::basic_iterator

#include "archive.hpp"
#include "time.hpp"
#include "species.hpp"
#include "cell.hpp"
#include "logics.hpp"

#include "error.hpp"

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

class TissueSimulation;

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

    std::set<std::string>   epistate_names;  //!< The names of the epigenetic states

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
     * @brief Count the number of cells in a neighborhood in a species
     *
     * @param position is the central position of the neighborhood
     * @param species_id is the identifier of the searched species
     * @return the number of cells in a neighborhood of `position` having
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
     * @brief Add a species
     *
     * @param species_properties are the properties of the species to
     *    be added
     */
    Tissue& add_species(const SpeciesProperties& species_properties);

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

    /**
     * @brief A class to go through the species of a mutant
     *
     * This class allows to have a view of the tissue's species of a
     * specific mutant.
     *
     * @tparam IS_CONSTANT is a Boolean flag to establish whether
     *      the view is constant
     */
    template <bool IS_CONSTANT>
    class basic_mutant_view
    {
        using vector_type = std::conditional_t<IS_CONSTANT, const std::vector<Species>,
                                               std::vector<Species>>;
        using species_type = std::conditional_t<IS_CONSTANT, const Species, Species>;

        vector_type* species;       //!< The species vector
        const std::vector<size_t>* species_pos;    //!< The species position vector

        /**
         * @brief A constructor
         *
         * @param species is a reference to the tissue's species vector
         * @param species_pos is a vector of valid positions for `vector`
         */
        basic_mutant_view(vector_type* species, const std::vector<size_t>& species_pos):
            species{species}, species_pos{&species_pos}
        {}

    public:
        /**
         * @brief An iterator over species views
         *
         * @tparam IS_CONSTANT_IT is a Boolean flag to establish whether
         *      the iterator is constant
         */
        template <bool IS_CONSTANT_IT>
        class basic_iterator {
            using vector_type = std::conditional_t<IS_CONSTANT_IT, const std::vector<Species>,
                                                   std::vector<Species>>;
            using species_type = std::conditional_t<IS_CONSTANT_IT, const Species, Species>;

            vector_type* species;    //!< The tissue's species vector
            std::shared_ptr<std::vector<size_t>::const_iterator> it; //!< A constant iterator over a vector of position in `species`

            /**
             * @brief A protected constructor
             *
             * @param species is the vector of the tissue's species
             * @param it is a constant iterator over a vector of the positions in `species`
             */
            basic_iterator(vector_type* species, const std::vector<size_t>::const_iterator it):
                species{species}, it{std::make_shared<std::vector<size_t>::const_iterator>(it)}
            {}
        public:
            using difference_type   =   std::ptrdiff_t;
            using value_type        =   Species;
            using pointer           =   species_type*;
            using reference         =   species_type&;
            using iterator_category =   std::random_access_iterator_tag;

            /**
             * @brief A copy constructor
             *
             * This copy constructor enables the implicit conversion from `basic_iterator<false>`
             * into `basic_iterator<true>`.
             *
             * @tparam PARAM_CONST is the template parameter of the copied object.
             * @param orig is the original basic iterator
             */
            template <bool PARAM_CONST>
                requires (IS_CONSTANT_IT)
            basic_iterator(const basic_iterator<PARAM_CONST>& orig) :
                species{orig.species}, it{orig.it}
            {}

            /**
             * @brief An empty constructor
             */
            basic_iterator():
                species{nullptr}, it{nullptr}
            {}

            /**
             * @brief Reference operator
             *
             * @return a reference to the species by the iterator
             */
            inline reference operator*() const
            {
                return (*species)[**it];
            }

            /**
             * @brief Pointer operator
             *
             * @return a pointer to the species by the iterator
             */
            inline pointer operator->() const
            {
                return &((*species)[**it]);
            }

            /**
             * @brief The prefix increment
             *
             * @return a reference to the updated object
             */
            inline basic_iterator<IS_CONSTANT_IT>& operator++()
            {
                ++(*it);

                return *this;
            }

            /**
             * @brief The postfix increment
             *
             * @return a copy of the original object
             */
            basic_iterator<IS_CONSTANT_IT> operator++(int)
            {
                basic_iterator<IS_CONSTANT_IT> copy(*this);

                this->operator++();

                return copy;
            }

            /**
             * @brief The prefix decrement
             *
             * @return a reference to the updated object
             */
            inline basic_iterator<IS_CONSTANT_IT>& operator--()
            {
                --(*it);

                return *this;
            }

            /**
             * @brief The postfix decrement
             *
             * @return a copy of the original object
             */
            basic_iterator<IS_CONSTANT_IT> operator--(int)
            {
                basic_iterator<IS_CONSTANT_IT> copy(*this);

                this->operator--();

                return copy;
            }

            /**
             * @brief Add operator
             *
             * @param delta is the value to add
             * @return a new iterator that points `delta` position ahead
             *      with respect to the original object
             */
            inline basic_iterator<IS_CONSTANT_IT> operator+(const int delta)
            {
                return basic_iterator<IS_CONSTANT_IT>(*species, *it + delta);
            }

            /**
             * @brief Subtract operator
             *
             * @param delta is the value to subtract
             * @return a new iterator that points `delta` position backwards
             *      with respect to the original object
             */
            inline basic_iterator<IS_CONSTANT_IT> operator-(const int delta)
            {
                return basic_iterator<IS_CONSTANT_IT>(*species, *it - delta);
            }

            /**
             * @brief Inplace add operator
             *
             * @param delta is the value to add
             * @return a reference to the update object
             */
            inline basic_iterator<IS_CONSTANT_IT>& operator+=(const int& delta)
            {
                (*it) += delta;

                return *this;
            }

            /**
             * @brief Inplace subtract operator
             *
             * @param delta is the value to subtract
             * @return a reference to the update object
             */
            inline basic_iterator<IS_CONSTANT_IT>& operator-=(const int& delta)
            {
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
            friend inline bool operator==(const basic_iterator<IS_CONSTANT_IT>& a,
                                          const basic_iterator<IS_CONSTANT_IT>& b)
            {
                return (*a.it == *b.it) && (a.species == b.species);
            }

            /**
             * @brief Test whether two iterators differs
             *
             * @param a is the first iterator to compare
             * @param b is the second iterator to compare
             * @return `true` if and only if the two iterators
             *      do not refer to the same object
             */
            friend inline bool operator!=(const basic_iterator<IS_CONSTANT_IT>& a,
                                          const basic_iterator<IS_CONSTANT_IT>& b)
            {
                return !(a==b);
            }

            template<bool IS_CONSTANT_IT2>  friend class basic_iterator;
            template<bool IS_CONSTANT2> friend class basic_mutant_view;
        };
    public:

        using const_iterator = basic_iterator<true>;
        using iterator = std::conditional_t<IS_CONSTANT, void, basic_iterator<false>>;

        /**
         * @brief The empty constructor
         */
        basic_mutant_view():
            species{nullptr}, species_pos{nullptr}
        {}

        /**
         * @brief A copy constructor
         *
         * This copy constructor enables the implicit conversion from `basic_mutant_view<false>`
         * into `basic_mutant_view<true>`.
         *
         * @tparam PARAM_CONST is the template parameter of the copied object.
         * @param orig is the original basic species view
         */
        template <bool PARAM_CONST>
            requires (IS_CONSTANT)
        basic_mutant_view(const basic_mutant_view<PARAM_CONST>& orig) :
            species{orig.species}, species_pos{orig.species_pos}
        {}

        /**
         * @brief Get the species by epigenetic state name
         *
         * @param epistate_name is the epigenetic state name of the aimed species
         * @return a reference to the species having `epistate_name` as
         *      epigenetic state name. If none of the species in the view have this
         *      epigenetic state name, an `std::out_of_range` is thrown.
         */
        const Species& get_species_by_epistate(const std::string& epistate_name) const
        {
            for (const auto& species: *this) {
                if (species.get_epistate_name()==epistate_name) {
                    return species;
                }
            }

            if (size()==0) {
                throw Error<std::out_of_range>("The mutant has no epigenetic states.");
            }

            throw Error<std::out_of_range>("\"" + (*this)[0].get_mutant_name()
                                        + "\" does not have the epigenetic state \""
                                        + epistate_name + "\".");
        }


        /**
         * @brief Get the species by epigenetic state name
         *
         * @param epistate_name is the epigenetic state name of the aimed species
         * @return a reference to the species having `epistate_name` as
         *      epigenetic state name. If none of the species in the view have this
         *      epigenetic state name, an `std::out_of_range` is thrown.
         */
        Species& get_species_by_epistate(const std::string& epistate_name)
            requires(!IS_CONSTANT)
        {
            for (Species& species: *this) {
                if (species.get_epistate_name()==epistate_name) {
                    return species;
                }
            }

            if (size()==0) {
                throw Error<std::out_of_range>("The mutant has no epigenetic states.");
            }

            throw Error<std::out_of_range>("\"" + (*this)[0].get_mutant_name()
                                        + "\" does not have the epigenetic state \""
                                        + epistate_name + "\".");
        }

        /**
         * @brief Index operator
         *
         * @param index is the index to access
         * @return a non-constant reference to the `index`-th species
         *      in the view
         */
        inline Species& operator[](const size_t& index)
            requires(!IS_CONSTANT)
        {
            return (*species)[(*species_pos)[index]];
        }

        /**
         * @brief Index operator
         *
         * @param index is the index to access
         * @return a constant reference to the `index`-th species
         *      in the view
         */
        inline const Species& operator[](const size_t& index) const
        {
            return (*species)[(*species_pos)[index]];
        }

        /**
         * @brief Index operator
         *
         * @param epistate_name is the epigenetic state name of the aimed species
         * @return a constant reference to the species having `epistate_name` as
         *      epigenetic state name. If none of the species in the view have this
         *      epigenetic state name, an `std::out_of_range` is thrown.
         */
        inline const Species& operator[](const std::string& epistate_name) const
        {
            return get_species_by_epistate(epistate_name);
        }

        /**
         * @brief Index operator
         *
         * @param epistate_name is the epigenetic state name of the aimed species
         * @return a non-constant reference to the species having `epistate_name` as
         *      epigenetic state name. If none of the species in the view have this
         *      epigenetic state name, an `std::out_of_range` is thrown.
         */
        inline Species& operator[](const std::string& epistate_name)
            requires(!IS_CONSTANT)
        {
            return get_species_by_epistate(epistate_name);
        }

        /**
         * @brief Get the mutant properties
         *
         * @return The mutant properties
         */
        inline const MutantProperties& get_mutant() const
        {
            return operator[](0).get_mutant_properties();
        }

        /**
         * @brief Get the view size
         *
         * @return the view size
         */
        inline size_t size() const
        {
            return species_pos->size();
        }

        /**
         * @brief Get the view begin
         *
         * @return a constant iterator to the first element
         *      in the view
         */
        inline const_iterator begin() const
        {
            return const_iterator(species, species_pos->begin());
        }

        /**
         * @brief Get the view end
         *
         * @return a constant iterator to the view end
         */
        inline const_iterator end() const
        {
            return const_iterator(species, species_pos->end());
        }

        /**
         * @brief Get the view begin
         *
         * @return a constant iterator to the first element
         *      in the view
         */
        inline iterator begin()
            requires (!IS_CONSTANT)
        {
            return iterator(species, species_pos->begin());
        }

        /**
         * @brief Get the view end
         *
         * @return a constant iterator to the view end
         */
        inline iterator end()
            requires (!IS_CONSTANT)
        {
            return iterator(species, species_pos->end());
        }

        /**
         * @brief Get the numer of cells in the species view
         *
         * @return the number of cells in the species view
         */
        size_t num_of_cells() const
        {
            size_t num_of_cells{0};
            for (auto species_it=begin(); species_it != end(); ++species_it) {
                num_of_cells += species_it->num_of_cells();
            }

            return num_of_cells;
        }

        friend class Tissue;
    };
public:

    using const_mutant_view = basic_mutant_view<true>;
    using mutant_view = basic_mutant_view<false>;

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

                oss << ") does not belong to the tissue.";

                throw Error<std::out_of_range>(oss.str());
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
         * @param event_type is the event type for which the availability is tested
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
         * @brief Get the valid neighborhood positions
         *
         * @return The list of valid neighborhood positions
         */
        inline std::list<PositionInTissue> get_neighborhood_positions() const
        {
            return tissue->get_neighborhood_positions(position);
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

            throw Error<std::runtime_error>("Wild-type cell.");
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
    Tissue(const std::string& name, const AxisSize x_size, const AxisSize y_size,
           const AxisSize z_size);

    /**
     * @brief A constructor for a 2D tissue
     *
     * @param name is the tissue name
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const std::string& name, const AxisSize x_size, const AxisSize y_size);

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
     * @brief Get a tissue species by name
     *
     * @param species_name is the species name
     * @return a constant reference to the tissue species
     */
    const Species& operator[](const std::string& species_name) const;

    /**
     * @brief Get a tissue species by name
     *
     * @param species_name is the species name
     * @return a non-constant reference to the tissue species
     */
    Species& operator[](const std::string& species_name);

    /**
     * @brief Add a mutant to the tissue
     *
     * @param mutant_properties is the mutant properties of the mutant
     * @return a reference to the updated object
     */
    Tissue& add_mutant(const MutantProperties& mutant_properties);

    /**
     * @brief Test whether the tissue knowns a species
     *
     * @param species_name is the species name
     * @return `true` if and only if the species has been added to the
     *      tissue
     */
    bool knowns_species(const std::string& species_name) const;

    /**
     * @brief Test whether the tissue knowns a mutant
     *
     * @param mutant_name is the mutant name
     * @return `true` if and only if the mutant has been added to the
     *      tissue
     */
    bool knowns_mutant(const std::string& mutant_name) const;

    /**
     * @brief Test whether the tissue knowns an epigenetic state
     *
     * @param epistate_name is the epigenetic state name
     * @return `true` if and only if the epigenetic state has been
     *      added to the tissue
     */
    bool knowns_epigenetic_state(const std::string& epistate_name) const;

    /**
     * @brief Add an epigenetic state to the tissue mutants
     *
     * @param epistate_name is the name of the epigenetic state
     * @return a reference to the updated tissue
     */
    Tissue& add_epigenetic_state(const std::string& epistate_name);

    /**
     * @brief Get the epigenetic state names
     *
     * @return The epigenetic state names
     */
    inline const std::set<std::string>& get_epigenetic_state_names() const
    {
        return epistate_names;
    }

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
     * @brief Get the valid neighborhood positions
     *
     * @param position is the the position whose neighbor positions
     *      are aimed
     * @return The list of valid neighborhood positions
     */
    std::list<PositionInTissue>
    get_neighborhood_positions(const PositionInTissue& position) const;

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
     * @brief Get a view of a mutant species
     *
     * @param mutant_name is the name of the mutant whose species are aimed
     * @return a const interator over the tissue's species having `mutant_name` as
     *       mutant name
     */
    const_mutant_view get_mutant_view(const std::string& mutant_name) const;

    /**
     * @brief Get a view of a mutant species
     *
     * @param mutant_id is the identifier of the mutant whose species are aimed
     * @return a const interator over the tissue's species having `mutant_id` as
     *       mutant identifier
     */
    const_mutant_view get_mutant_view(const MutantId mutant_id) const;


    /**
     * @brief Get a view of a mutant species
     *
     * @param mutant_name is the name of the mutant whose species are aimed
     * @return an interator over the tissue's species having `mutant_name` as
     *       mutant name
     */
    mutant_view get_mutant_view(const std::string& mutant_name);

    /**
     * @brief Get a view of a mutant species
     *
     * @param mutant_id is the identifier of the mutant whose species are aimed
     * @return an interator over the tissue's species having `mutant_id` as
     *       mutant identifier
     */
    mutant_view get_mutant_view(const MutantId mutant_id);

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
     * @param duplicate_internal_cells is a Boolean flag that, when set to
     *      be `false`, removes internal cells from the cells available for
     *      the duplication
     * @return the list of the cells that have been pushed outside the
     *      tissue border
     */
    std::list<Cell> push_cells(const PositionInTissue from_position,
                               const Direction& direction,
                               const bool duplicate_internal_cells);

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
    friend class TissueSimulation;
};

}   // Evolutions

}   // Mutants

}   // CLONES

#endif // __CLONES_TISSUE__
