/**
 * @file mutational_properties.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to represent the mutational properties
 * @version 1.6
 * @date 2025-11-24
 *
 * @copyright Copyright (c) 2023-2025
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

#ifndef __RACES_MUTATIONAL_PROPERTIES__
#define __RACES_MUTATIONAL_PROPERTIES__

#include <map>
#include <list>
#include <string>

#include "sid.hpp"
#include "cna.hpp"
#include "mutation_spec.hpp"

#include "mutation_list.hpp"
namespace RACES
{

namespace Mutations
{

/**
 * @brief Rates of the instant passenger mutations for a species
 */
struct PassengerRates
{
    double indel;   //!< The species indel rate
    double snv;     //!< The species SNV rate
    double cna;     //!< The species CNA rate

    /**
     * @brief The empty constructor
     */
    PassengerRates();

    /**
     * @brief A constructor
     *
     * @param indel_rate is the the SNV rate
     * @param SNV_rate is the the SNV rate
     * @param CNA_rate is the the CNA rate
     */
    PassengerRates(const double& indel_rate, const double& SNV_rate,
                          const double& CNA_rate);

    /**
     * @brief Save passenger rates
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & indel
                & snv
                & cna;
    }

    /**
     * @brief Load instant passenger rates
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load instant passenger rates object
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static PassengerRates load(ARCHIVE& archive)
    {
        PassengerRates p_rates;

        archive & p_rates.indel
                & p_rates.snv
                & p_rates.cna;

        return p_rates;
    }
};

/**
 * @brief A template to represent timed object
 *
 * This class represents objects whose values changes during time. The values of these
 * objects stand still during time intervals of the type [begin, end) and can change
 * as soon as a new time interval begin. The time intervals and the values are set
 * by users. The first time interval always starts at time 0 and initially holds from
 * time 0 up to the end of the time horizon.
 *
 * @tparam T is the type of timed objects
 */
template<typename T>
class Timed
{
    std::map<Time, T> time_map;  //!< a map from the switch time to the values

    /**
     * @brief Find an iterator to a value at a given time
     *
     * This method finds the (constant) iterator referring to the value at the
     * specified time. It is used to simplify the implementations of non-static
     * methods by handling both constant and non-constant time interval parameters.
     *
     * @tparam MAP is the type of `time_map`. It can be constant or non-constant
     * @param time_map is a map from `Time` to the values
     * @param time is the time for which the value is searched
     * @return the object in `time_map` whose key is the greatest key smaller
     *      or equal to `time`
     */
    template<typename MAP>
    static auto find_iterator_at(MAP& time_map, const Time& time)
    {
        if (time<0) {
            throw std::domain_error("The time must be non-negative.");
        }

        auto found = time_map.upper_bound(time);

        // actually this check is not required because by the constructor
        // add the default time zone from 0
        if (found != time_map.begin()) {
            --found;
        }

        return found;
    }
public:
    /**
     * @brief The empty constructor
     *
     * This constructor initialized the timed objects by adding a default
     * value from the time 0.
     */
    Timed():
        time_map{{0, T{}}}
    {}

    /**
     * @brief Get the value at a given time
     *
     * @param time is the time at which the object is evaluated
     * @return a constant reference to the value of the object at
     *      time `time`
     */
    const T& at(const Time& time) const
    {
        auto found = find_iterator_at(time_map, time);

        return found->second;
    }

    /**
     * @brief Add a new value from a specified time up the time horizon
     *
     * This method assigns a new value to the object from a specified time
     * up the time horizon. The values that the object had from the
     * specified time before the method call are deleted.
     *
     * @param time is the time from which the object assume the new value
     * @return a reference to the new value
     */
    inline T& from(const Time& time)
    {
        delete_down_to(time);

        return change_value_at(time);
    }

    /**
     * @brief Delete the values assumed after a specified time
     *
     * @param time is the greatest time for which the values are
     *      not deleted
     */
    void delete_down_to(const Time& time)
    {
        auto it = time_map.end();

        if (it != time_map.begin()) {
            --it;
        }
        while (it != time_map.begin()) {
            if (it->first <= time) {
                return;
            }

            it = --time_map.erase(it);
        }
    }

    /**
     * @brief Change value at a given time
     *
     * @param time is the time at which the value changes
     * @param value is the new value of the object
     * @return a reference to the value of the object from
     *      `time` up to the next value change
     */
    T& change_value_at(const Time& time, const T& value)
    {
        auto found = find_iterator_at(time_map, time);
        if (found->first == time) {
            found->second = value;
        } else {
            std::tie(found, std::ignore) = time_map.emplace(time, value);
        }

        return found->second;
    }

    /**
     * @brief Change value at a given time
     *
     * @param time is the time at which the value changes
     * @return a reference to the value of the object from
     *      `time` up to the next value change
     */
    T& change_value_at(const Time& time)
    {
        auto found = find_iterator_at(time_map, time);
        if (found->first != time) {
            std::tie(found, std::ignore) = time_map.emplace(time, found->second);
        }

        return found->second;
    }

    /**
     * @brief Delete value change
     *
     * @param time is the time at which the value change must be deleted
     */
    void delete_change_at(const Time& time)
    {
        if (time == 0) {
            throw std::domain_error("The time interval beginning at time 0 "
                                    "cannot be removed");
        }

        auto found = time_map.find(time);

        if (found == time_map.end()) {
            throw std::runtime_error("No is time interval beginning at time "
                                     + std::to_string(time));
        }

        time_map.remove(found);
    }

    /**
     * @brief Test whether the object value changes at a given time
     *
     * @param time is the time at which the value change is queried
     * @return `true` if and only if the value of the current object
     *      changes at time `time`
     */
    inline bool switches_at(const Time& time) const
    {
        return time_map.find(time) != time_map.end();
    }

    /**
     * @brief Get the constant iterator to the time-value map begin
     *
     * @return the constant iterator to the time-value map begin
     */
    inline std::map<Time, T>::const_iterator begin() const
    {
        return time_map.begin();
    }

    /**
     * @brief Get the constant iterator to the time-value map end
     *
     * @return the constant iterator to the time-value map end
     */
    inline std::map<Time, T>::const_iterator end() const
    {
        return time_map.end();
    }

    /**
     * @brief Get the number of different values of this object
     *
     * @return the number of different values of this object
     */
    inline size_t num_of_values() const
    {
        return time_map.size();
    }

    /**
     * @brief Save timed data
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & time_map;
    }

    /**
     * @brief Load timed data
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load timed data
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static Timed<T> load(ARCHIVE& archive)
    {
        Timed<T> timed_data;

        archive & timed_data.time_map;

        return timed_data;
    }
};

/**
 * @brief A class representing the mutational properties of all the species
 *
 */
class MutationalProperties
{
    std::map<std::string, Timed<PassengerRates>> passenger_rates;      //!< The species passenger rates
    std::map<std::string, DriverMutations> driver_mutations;    //!< The mutation per mutant

    /**
     * @brief Get the passenger rates for a species
     *
     * @tparam MAP is the type of `passenger_rates`. It can be constant or non-constant
     * @param passenger_rates is a map associating species names to their rates
     * @param species_name is the name of the species whose rates are aimed
     * @return a (constant) reference to the rates associated to `species_name` by
     *      `passenger_rates`
     */
    template<typename MAP>
    static auto& get_passenger_rates(MAP& passenger_rates, const std::string& species_name)
    {
        auto found = passenger_rates.find(species_name);
        if (found == passenger_rates.end()) {
            throw std::domain_error("\"" + species_name + "\" has no mutational rate");
        }

        return found->second;
    }

    /**
     * @brief Get the passenger rates for a species
     *
     * @tparam MAP is the type of `passenger_rates`. It can be constant or non-constant
     * @param passenger_rates is a map associating species names to their rates
     * @param mutant_name is the mutant of the species whose rates are aimed
     * @param epistate is the epigenetic state representation of the species whose
     *      rates are aimed
     * @return a (constant) reference to the rates associated by `passenger_rates` to
     *      species whose mutant and epigenetic state are `mutant_name` and `epistate`,
     *      respectively
     */
    template<typename T>
    inline static auto& get_passenger_rates(T& passenger_rates, const std::string& mutant_name,
                                            const std::string& epistate)
    {
        const auto species_name = mutant_name + epistate;

        return get_passenger_rates(passenger_rates, species_name);
    }

    /**
     * @brief Get the passenger rates for a species
     *
     * @tparam MAP is the type of `passenger_rates`. It can be constant or non-constant
     * @param passenger_rates is a map associating species names to their rates
     * @param mutant_name is the mutant of the species whose rates are aimed
     * @param signature is the epigenetic signature of the species whose rates are
     *      aimed
     * @return a (constant) reference to the rates associated by `passenger_rates` to
     *      species whose mutant and epigenetic signature are `mutant_name` and
     *      `signature`, respectively
     */
    template<typename T>
    static auto& get_passenger_rates(T& passenger_rates, const std::string& mutant_name,
                                     const Mutants::MethylationSignature& signature)
    {
        const auto epistate = Mutants::MutantProperties::signature_to_string(signature);

        return get_passenger_rates(passenger_rates, mutant_name, epistate);
    }
public:

    /**
     * @brief The empty constructor
     */
    MutationalProperties();

    /**
     * @brief Add the properties of a mutant
     *
     * @param mutant_name is the name of the mutant
     * @param epistate_passenger_rates is a map from epigenetic state to
     *          passenger rates
     * @param driver_SIDs is a list of driver SIDs
     * @param driver_CNAs is a list of driver CNAs
     * @param wg_doubling is a Boolean flag to enable whole genome doubling
     * @return a reference to the updated object
     */
    MutationalProperties& add_mutant(const std::string& mutant_name,
                                     const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                     const std::list<MutationSpec<SID>>& driver_SIDs={},
                                     const std::list<CNA>& driver_CNAs={},
                                     const bool& wg_doubling=false);

    /**
     * @brief Add the properties of a mutant
     *
     * @param mutant_name is the name of the mutant
     * @param epistate_passenger_rates is a map from epigenetic state to
     *          passenger rates
     * @param driver_SIDs is a list of driver SIDs
     * @param driver_CNAs is a list of driver CNAs
     * @param application_order is the list of mutation application order
     * @return a reference to the updated object
     */
    MutationalProperties& add_mutant(const std::string& mutant_name,
                                     const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                     const std::list<MutationSpec<SID>>& driver_SIDs,
                                     const std::list<CNA>& driver_CNAs,
                                     const std::list<DriverMutations::MutationType>& application_order);

    /**
     * @brief Change the passenger rates from a given time stamp
     *
     * This method changes the passenger rates for some species from a given timestamp
     * of the simulated time up to the time horizon. The previously recorder rates after
     * the provided timestamp are deleted.
     *
     * @param time is a time stamp
     * @param mutant_name is a mutant name
     * @param epistate_passenger_rates is a map associating epigenetic state representations
     *      to the corresponding rates
     * @return a reference to the updated object
     */
    MutationalProperties& change_rates_from(const Time& time, const std::string& mutant_name,
                                            const std::map<std::string, PassengerRates>& epistate_passenger_rates);

    /**
     * @brief Get the species-passenger rates map
     *
     * @return a constant reference to the species-passenger rates map
     */
    inline const std::map<std::string, Timed<PassengerRates>>& get_passenger_rates() const
    {
        return passenger_rates;
    }

    /**
     * @brief Get the timed passenger rates of a species
     *
     * @param mutant_name is the mutant of the species whose rates are aimed
     * @param epistate is the epigenetic state representation of the species whose
     *      rates are aimed
     * @return a constant reference to the timed rates associated by `passenger_rates`
     *      to the species whose mutant and epigenetic state are `mutant_name` and
     *      `epistate`, respectively
     */
    inline const Timed<PassengerRates>&
    get_passenger_rates(const std::string& mutant_name, const std::string& epistate) const
    {
        return get_passenger_rates(passenger_rates, mutant_name, epistate);
    }

    /**
     * @brief Get the timed passenger rates of a species
     *
     * @param mutant_name is the mutant of the species whose rates are aimed
     * @param epistate is the epigenetic state representation of the species whose
     *      rates are aimed
     * @return a reference to the timed rates associated by `passenger_rates`
     *      to the species whose mutant and epigenetic state are `mutant_name` and
     *      `epistate`, respectively
     */
    inline Timed<PassengerRates>&
    get_passenger_rates(const std::string& mutant_name, const std::string& epistate)
    {
        return get_passenger_rates(passenger_rates, mutant_name, epistate);
    }

    /**
     * @brief Get the timed passenger rates of a species
     *
     * @param mutant_name is the mutant of the species whose rates are aimed
     * @param signature is the epigenetic signature of the species whose rates are
     *      aimed
     * @return a constant reference to the timed rates associated by
     *      `passenger_rates` to the species whose mutant and epigenetic signature
     *      are `mutant_name` and `signature`, respectively
     */
    inline const Timed<PassengerRates>&
    get_passenger_rates(const std::string& mutant_name,
                        const Mutants::MethylationSignature& signature) const
    {
        return get_passenger_rates(passenger_rates, mutant_name, signature);
    }

    /**
     * @brief Get the timed passenger rates of a species
     *
     * @param mutant_name is the mutant of the species whose rates are aimed
     * @param signature is the epigenetic signature of the species whose rates are
     *      aimed
     * @return a reference to the timed rates associated by `passenger_rates` to
     *      the species whose mutant and epigenetic signature are `mutant_name`
     *      and `signature`, respectively
     */
    inline Timed<PassengerRates>&
    get_passenger_rates(const std::string& mutant_name,
                        const Mutants::MethylationSignature& signature)
    {
        return get_passenger_rates(passenger_rates, mutant_name, signature);
    }

    /**
     * @brief Get the timed passenger rates of a species
     *
     * @param species_name is the name of the species whose rates are aimed
     * @return a constant reference to the timed rates associated by
     *      `passenger_rates` to the species whose name is `species_name`
     */
    inline const Timed<PassengerRates>&
    get_passenger_rates(const std::string& species_name) const
    {
        return get_passenger_rates(passenger_rates, species_name);
    }

    /**
     * @brief Get the timed passenger rates of a species
     *
     * @param species_name is the name of the species whose rates are aimed
     * @return a reference to the timed rates associated by `passenger_rates` to
     *      the species whose name is `species_name`
     */
    inline Timed<PassengerRates>&
    get_passenger_rates(const std::string& species_name)
    {
        return get_passenger_rates(passenger_rates, species_name);
    }

    /**
     * @brief Get the driver mutation map
     *
     * @return a constant reference to map associating a mutant name to its
     *      driver mutation
     */
    inline const std::map<std::string, DriverMutations>& get_driver_mutations() const
    {
        return driver_mutations;
    }

    /**
     * @brief Set the driver mutation map
     *
     * @param mutant_drivers is the new list of driver mutations
     */
    void set_mutant_drivers(const DriverMutations& mutant_drivers);

    /**
     * @brief Save mutational properties
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & passenger_rates
                & driver_mutations;
    }

    /**
     * @brief Load mutational properties
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load mutational properties object
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static MutationalProperties load(ARCHIVE& archive)
    {
        MutationalProperties m_properties;

        archive & m_properties.passenger_rates
                & m_properties.driver_mutations;

        return m_properties;
    }
};

}   // Mutations

}   // RACES

#endif // __RACES_MUTATIONAL_PROPERTIES__
