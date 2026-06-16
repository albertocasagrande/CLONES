/**
 * @file mutant_properties.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines mutant properties
 * @version 1.4
 * @date 2026-06-16
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

#ifndef __CLONES_MUTANT_PROPERTIES__
#define __CLONES_MUTANT_PROPERTIES__

#include <map>
#include <string>
#include <sstream>
#include <ranges>

#include "species_name.hpp"

#include "archive.hpp"
#include "cell_event.hpp"

namespace CLONES
{

namespace Mutants
{

class MutantProperties;

namespace Evolutions
{

    class Species;

}   // Evolutions

class EventRate
{
    CellEventType type;   //!< Event type
    SpeciesId dst;        //!< Destination species
    double rate;          //!< Event rate

public:
    /**
     * @brief The empty constructor
     */
    EventRate():
        EventRate{CellEventType::ANY, 0, 0.0}
    {}

    /**
     * @brief A constructor
     *
     * This is a generic constructor for
     *
     * @param type is the event type
     * @param dst is the destination of the species
     * @param rate is the rate of the event
     */
    EventRate(const CellEventType type,
              const SpeciesId dst, const double rate):
        type{type}, dst{dst}, rate{rate}
    {}

    /**
     * @brief Get the event type
     *
     * @return The event type
     */
    inline const CellEventType& get_type() const
    {
        return type;
    }

    /**
     * @brief Get the event destination
     *
     * @return The event destination
     */
    inline const SpeciesId& get_dst() const
    {
        return dst;
    }

    /**
     * @brief Get the event rate
     *
     * @return The event rate
     */
    inline const double& get_rate() const
    {
        return rate;
    }
};

/**
 * @brief A class to represent species properties
 *
 * A mutant is characterized by a (potentially unknown) genotype.
 * The same mutant may have different rates depending on its
 * methylation signature. A species is a mutant with a methylation
 * signature.
 */
class SpeciesProperties
{
public:

    using RateType = double;                                    //!< The rate type
    using SpeciesIdRate = std::map<SpeciesId, RateType>;        //!< The event type rate map
    using EventRate = std::map<CellEventType, SpeciesIdRate>;   //!< The event rate map

private:

    static unsigned int counter;    //!< The total number of species along the computation

    SpeciesId id;               //!< The species identifier
    MutantId mutant_id;         //!< The mutant identifier

    SpeciesName species_name;   //!< The species name

    EventRate event_rates;  //!< Event rates

    /**
     * @brief The empty constructor
     */
    SpeciesProperties();

    /**
     * @brief A constructor
     *
     * @param mutant is the mutant of the new object
     * @param mutant is the name of the epigenetic state of the new object
     */
    SpeciesProperties(const MutantProperties& mutant,
                      const std::string& epistate_name);
public:
    /**
     * @brief Get the species identifier
     *
     * @return a constant reference to the identifier
     */
    inline const SpeciesId& get_id() const
    {
        return id;
    }

    /**
     * @brief Get the name of the species epigenetic state
     *
     * @return a constant reference to name of the epigenetic state
     */
    inline const std::string& get_epistate_name() const
    {
        return species_name.get_epistate_name();
    }

    /**
     * @brief Get the mutant name
     *
     * @return a constant reference to the mutant name
     */
    inline const std::string& get_mutant_name() const
    {
        return species_name.get_mutant_name();
    }

    /**
     * @brief Get the species name
     *
     * @return the species name
     */
    inline std::string get_name() const
    {
        return static_cast<std::string>(species_name);
    }

    /**
     * @brief Get the mutant identifier
     *
     * @return a constant reference to the mutant identifier
     */
    inline const MutantId& get_mutant_id() const
    {
        return mutant_id;
    }

    /**
     * @brief Get the rate of an event
     *
     * @param event is the event whose rate is required
     * @param dst_species is the destination species identifier
     *      of the event
     * @return if the rate of `event` has been set, then
     *       the rate of `event`. The value 0 otherwise
     */
    double get_rate(const CellEventType& event,
                    const SpeciesId& dst_species) const;

    /**
     * @brief Get event rates
     *
     * @return a constant reference to event rates
     */
    inline const EventRate& get_rates() const
    {
        return event_rates;
    }

    /**
     * @brief Set the rate of an event
     *
     * @param event is the event whose rate is set the
     *      event. If the event is neither `CellEventType::DEATH`
     *      nor `CellEventType::DUPLICATION`, a `std::domain_error`
     *      exception is thrown
     * @param rate is the new rate for the event
     * @return a non-constant reference to species properties
     */
    SpeciesProperties& set_rate(const CellEventType& event,
                                const double rate);

    /**
     * @brief Set the rate of an event
     *
     * @param event is the event whose rate is set the
     *      event. If the event is not among `CellEventType::DEATH`,
     *      `CellEventType::DUPLICATION`,
     *      and `CellEventType::DUP_AND_EPI_SWITCH` or `dst_species` and this
     *      object do not belong to the same mutant, a `std::domain_error`
     *      exception is thrown
     * @param dst_species is the destination species of the
     *      event
     * @param rate is the new rate for the event
     * @return a non-constant reference to species properties
     */
    SpeciesProperties& set_rate(const CellEventType& event,
                                const SpeciesProperties& dst_species,
                                const double rate);

    /**
     * @brief Save the species in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & id
                & mutant_id
                & species_name
                & event_rates;
    }

    /**
     * @brief Load a species from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded species
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SpeciesProperties load(ARCHIVE& archive)
    {
        SpeciesProperties species_properties;

        archive & species_properties.id
                & species_properties.mutant_id
                & species_properties.species_name
                & species_properties.event_rates;

        if (SpeciesProperties::counter < static_cast<unsigned int>(species_properties.id+1)) {
            SpeciesProperties::counter = species_properties.id+1;
        }

        return species_properties;
    }

    friend class MutantProperties;
    friend class Evolutions::Species;
};

/**
 * @brief A class representing mutant properties
 */
class MutantProperties
{
    static unsigned int counter; //!< The total number of mutant along the computation

    MutantId id;        //!< The mutant identifier
    std::string name;   //!< The mutant name

    std::map<std::string, SpeciesProperties> species;   //!< The species properties of this mutant
public:

    /**
     * @brief A constructor
     *
     * This constructor builds a mutant and its species.
     * The number of species is determined by the size of the
     * epigenetic state names.
     *
     * @param name is the name of the mutant
     * @param epi_states is a list of the epigenetic state names of this mutant
     */
    MutantProperties(const std::string& name, std::list<std::string>&& epi_states);

    /**
     * @brief A constructor
     *
     * This constructor builds a mutant and its species.
     * The number of species is determined by the size of the
     * epigenetic state names.
     *
     * @param name is the name of the mutant
     * @param epi_states is a list of the epigenetic state names of this mutant
     */
    MutantProperties(const std::string& name, const std::list<std::string>& epi_states);

    /**
     * @brief A constructor
     *
     * This constructor builds a mutant consisting in one species.
     *
     * @param name is the name of the mutant
     */
    explicit MutantProperties(const std::string& name);

    /**
     * @brief Get the species identifier
     *
     * @return a constant reference to the identifier
     */
    inline const MutantId& get_id() const
    {
        return id;
    }

    /**
     * @brief Get the mutant name
     *
     * @return a constant reference to the mutant name
     */
    inline const std::string& get_name() const
    {
        return name;
    }

    /**
     * @brief Get the species associated to this mutant
     *
     * @return a constant reference to the species
     */
    inline const std::map<std::string, SpeciesProperties>& get_species() const
    {
        return species;
    }

    /**
     * @brief Get the species associated to this mutant
     *
     * @return a non-constant reference to the species
     */
    inline std::map<std::string, SpeciesProperties>& get_species()
    {
        return species;
    }

    /**
     * @brief Get the species by epigenetic state name
     *
     * @param epistate_name is the epigenetic state name of the aimed species
     * @return A constant reference to the species having `epistate_name` as the
     *      epigenetic state name
     */
    const SpeciesProperties& operator[](const std::string& epistate_name) const;

    /**
     * @brief Get the species by epigenetic state name
     *
     * @param epistate_name is the epigenetic state name of the aimed species
     * @return A non-constant reference to the species having `epistate_name` as
     *      the epigenetic state name
     */
    SpeciesProperties& operator[](const std::string& epistate_name);

    /**
     * @brief Test whether two mutants have the same set of epigenetic states
     *
     * @param mutant is a `MutantProperties` object
     * @return `true` iff this object and `mutant` have the same set of
     *      epigenetic states
     */
    inline bool has_the_same_epistates(const MutantProperties& mutant) const
    {
        return have_the_same_epistates(*this, mutant);
    }

    /**
     * @brief Test whether two mutants have the same set of epigenetic states
     *
     * @param mutant_a is a `MutantProperties` object
     * @param mutant_b is a `MutantProperties` object
     * @return `true` iff `mutant_a` and `mutant_b` have the same set of
     *      epigenetic states
     */
    static bool have_the_same_epistates(const MutantProperties& mutant_a,
                                        const MutantProperties& mutant_b);
};

/**
 * @brief Test whether two species properties are the same
 *
 * @param a is a `SpeciesProperties` object
 * @param b is a `SpeciesProperties` object
 * @return `true` iff the two species properties are the same
 */
bool operator==(const SpeciesProperties& a, const SpeciesProperties& b);

/**
 * @brief Test whether two species properties differ
 *
 * @param a is a `SpeciesProperties` object
 * @param b is a `SpeciesProperties` object
 * @return `true` iff the two species properties differ
 */
inline bool operator!=(const SpeciesProperties& a, const SpeciesProperties& b)
{
    return !(a==b);
}

}   // Mutants

}   // CLONES

namespace std
{

/**
 * @brief Write information about a species in an output stream
 *
 * @param out is the output stream
 * @param species is the species to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const CLONES::Mutants::SpeciesProperties& species);

/**
 * @brief Write information about a mutant in an output stream
 *
 * @param out is the output stream
 * @param mutant is the mutant to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const CLONES::Mutants::MutantProperties& mutant);

}   // std


#endif // __CLONES_MUTANT_PROPERTIES__