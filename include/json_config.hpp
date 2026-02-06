/**
 * @file json_config.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for reading JSON configurations
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

#ifndef __CLONES_JSON_CONFIG__
#define __CLONES_JSON_CONFIG__

#include <map>
#include <vector>

#include <nlohmann/json.hpp>

#include "simulation.hpp"
#include "position_set.hpp"
#include "mutation_engine.hpp"
#include "timed_event.hpp"

namespace CLONES
{

class ConfigReader
{
    /**
     * @brief Get a sample region corner
     *
     * @param sampler_region_json is the JSON of the region to sample
     * @param corner_field_name is the corner field name
     * @return std::vector<CLONES::Mutants::Evolutions::AxisPosition>
     */
    static std::vector<CLONES::Mutants::Evolutions::AxisPosition>
    get_corner(const nlohmann::json& sampler_region_json, const std::string& corner_field_name);

    /**
     * @brief Get the passenger mutation rates of one epigenetic state
     *
     * @param epistate_rates_json is the JSON of the epistate passenger rates
     * @return a pair epigenetic state/passenger mutation rate
     */
    static std::pair<std::string, CLONES::Mutations::PassengerRates>
    get_epistate_passenger_rates(const nlohmann::json& epistate_rates_json);

    /**
     * @brief Get the CNA type
     *
     * @param CNA_json is the JSON of the CNA
     * @return the type of the CNA
     */
    static CLONES::Mutations::CNA::Type get_CNA_type(const nlohmann::json& CNA_json);

    /**
     * @brief Add a CNA to a list
     *
     * @param CNAs is the list of CNA in which the read CNA should be inserted
     * @param CNA_json is the JSON of the CNA to be inserted into the list
     */
    static void add_CNA(std::list<CLONES::Mutations::CNA>& CNAs,
                        const nlohmann::json& CNA_json);

    /**
     * @brief Add a SID mutation to a list
     *
     * @param mutant_name is the mutant name
     * @param SIDs is the list of SID in which the read SID should be inserted
     * @param SID_json is the JSON of the SID to be inserted into the list
     */
    static void
    add_SID(const std::string& mutant_name,
            std::list<CLONES::Mutations::MutationSpec<CLONES::Mutations::SID>>& SIDs,
            const nlohmann::json& SID_json);

    /**
     * @brief Get the passenger rates
     *
     * @param passenger_rates_json is the JSON of the passenger rates
     * @return a map associating an epigenetic status to its mutation rate
     */
    static std::map<std::string, CLONES::Mutations::PassengerRates>
    get_passenger_rates(const nlohmann::json& passenger_rates_json);

    /**
     * @brief Add the mutations to SID and CNA lists
     *
     * @param mutant_name is the mutant name
     * @param SIDs is the list of SIDs in which the new SIDs must be inserted
     * @param CNAs is the list of CNAs in which the new CNAs must be inserted
     * @param mutation_json is the JSON of the mutation
     */
    static void
    schedule_mutation(const std::string& mutant_name,
                      std::list<CLONES::Mutations::MutationSpec<CLONES::Mutations::SID>>& SIDs,
                      std::list<CLONES::Mutations::CNA>& CNAs,
                      const nlohmann::json& mutation_json);

    /**
     * @brief Get the fraction field
     *
     * @param fraction_json is the JSON of the faction field
     * @return the fraction value
     */
    static double get_fraction(const nlohmann::json& fraction_json);

    /**
     * @brief Add mutational properties
     *
     * @param mutational_properties are the mutational properties of all the species
     * @param mutational_properties_json is the JSON of the properties
     */
    static void
    add_mutational_properties(CLONES::Mutations::MutationalProperties& mutational_properties,
                              const nlohmann::json& mutational_properties_json);

public:
    /**
     * @brief Get the default exposure
     *
     * @param mutation_name is the mutation type name (i.e., either "indel" or "SBS")
     * @param exposures_json is the JSON of the exposures
     * @return a map associating to a set of SBS signatures their percentage
     *          in the default mutational configuration
     */
    static CLONES::Mutations::MutationalExposure
    get_default_exposure(const std::string& mutation_name,
                         const nlohmann::json& exposures_json);

    /**
     * @brief Add the timed exposures to a mutation engine
     *
     * @param mutation_name is the mutation type name (i.e., either "indel" or "SBS")
     * @param exposures_json is the JSON of the exposures
     */
    static std::map<double, CLONES::Mutations::MutationalExposure>
    get_timed_exposures(const std::string& mutation_name,
                        const nlohmann::json& exposures_json);

    /**
     * @brief Get a sample region
     *
     * @param sample_region_json is the sample region JSON
     * @param field_name is the name of the field
     * @return the rectangle position set
     */
    static Mutants::RectangleSet
    get_sample_region(const nlohmann::json& sample_region_json,
                      const std::string& field_name="sample region");

    /**
     * @brief Get a sample specification
     *
     * @param sample_specification_json is the sample specification JSON
     * @param field_name is the name of the field
     * @return a sample specification
     */
    static Mutants::Evolutions::SampleSpecification
    get_sample_specification(const nlohmann::json& sample_specification_json,
                             const std::string& field_name="sample");

    /**
     * @brief Collect the mutations
     *
     * @param[in] mutant_name is the mutant name
     * @param[in,out] SIDs is the list of mutant SIDs
     * @param[in,out] CNAs is the list of mutant CNAs
     * @param[in] mutations_json is the JSON of mutations
     */
    static void
    collect_mutations(const std::string& mutant_name,
                      std::list<CLONES::Mutations::MutationSpec<CLONES::Mutations::SID>>& SIDs,
                      std::list<CLONES::Mutations::CNA>& CNAs,
                      const nlohmann::json& mutations_json);

    /**
     * @brief Get the species mutational properties
     *
     * @param configuration_json is the JSON of the simulation configuration
     * @return the species mutational properties
     */
    static CLONES::Mutations::MutationalProperties
    get_mutational_properties(const nlohmann::json& configuration_json);

    /**
     * @brief Extract a timed event from a JSON object
     *
     * @param simulation is a simulation
     * @param name2mutant is the map from name to mutant
     * @param timed_event_json is the JSON of the timed event
     * @return the time event described in `timed_event_json`
     */
    static CLONES::Mutants::Evolutions::TimedEvent
    get_timed_event(const CLONES::Mutants::Evolutions::Simulation& simulation,
                    const std::map<std::string, CLONES::Mutants::MutantProperties> name2mutant,
                    const nlohmann::json& timed_event_json);

    /**
     * @brief Raise an exception if a field is not available in a JSON
     *
     * @param field_name is the aimed field name
     * @param json is the JSON that must contain the field
     * @param context is a description of the context
     */
    static
    void expecting(const std::string& field_name, const nlohmann::json& json,
                   const std::string& context="");

    /**
     * @brief Get a value of a JSON field
     *
     * @tparam TYPE is the type of the JSON field value
     * @param field_name is the aimed field name
     * @param json is the JSON that must contain the field
     * @param context is a description of the context
     * @return the value of the JSON field `field_name`
     */
    template<typename TYPE>
    static TYPE get_from(const std::string& field_name, const nlohmann::json& json,
                         const std::string& context="")
    {
        expecting(field_name, json, context);

        return json[field_name].template get<TYPE>();
    }

    /**
     * @brief Get the number of wild-type genome alleles from the configuration
     *
     * @param simulation_json is the JSON simulation configuration
     * @param chromosome_ids is the vector of the genome chromosome ids
     * @return a pair whose first element is the default number of allele per
     *          chromosome and the second element is a map reporting the exceptions
     *          The map associates the chromosome names to the corresponding
     *          number of alleles
     */
    static std::map<Mutations::ChromosomeId, size_t>
    get_number_of_alleles(const nlohmann::json& simulation_json,
                          const std::vector<Mutations::ChromosomeId>& chromosome_ids);

    /**
     * @brief Get the exposure
     *
     * @param exposure_json the JSON of the exposure
     * @return a map associating a set of SBS signature to their percentage
     *          in the exposure
     */
    static CLONES::Mutations::MutationalExposure
    get_exposure(const nlohmann::json& exposure_json);

    /**
     * @brief Get the number of germline mutations per kilobase
     *
     * @param simulation_json is the JSON simulation configuration
     * @return the number of germline mutations per kilobase
     */
    static double
    get_number_of_germline_mutations_per_kbase(const nlohmann::json& simulation_json)
    {
        return get_from<double>("germline mutations per kbase", simulation_json,
                                "The passengers simulation configuration");
    }

    /**
     * @brief Get the number of neoplastic mutations
     *
     * @param simulation_json is the JSON simulation configuration
     * @param type_name is the type name of the mutations (i.e., "SNV" or "indel")
     * @return the number of neoplastic mutations
     */
    static size_t
    get_number_of_pre_neoplastic_mutations(const nlohmann::json& simulation_json,
                                           const std::string& type_name)
    {
        std::string desc = "number of pre-neoplastic " + type_name + "s";

        return get_from<size_t>(desc, simulation_json,
                                "The passengers simulation configuration");
    }
};

}   // CLONES

#endif // __CLONES_JSON_CONFIG__
