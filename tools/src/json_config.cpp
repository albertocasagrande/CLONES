/**
 * @file json_config.cpp
 * @author Alberto Casagrande (alberto.casagrande@units.it)
 * @brief Implements classes and function for reading JSON configurations
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

#include "json_config.hpp"

#include "timed_event.hpp"

namespace CLONES
{

std::vector<CLONES::Mutants::Evolutions::AxisPosition>
ConfigReader::get_corner(const nlohmann::json& sampler_region_json, const std::string& corner_field_name)
{
    if (!sampler_region_json.contains(corner_field_name)) {
        throw std::runtime_error("The \"sampler region\" must contain a "
                                    "\"" + corner_field_name + "\" field");
    }

    auto corner_json = sampler_region_json[corner_field_name];

    if (!corner_json.is_array()) {
        throw std::runtime_error("The \"" + corner_field_name + "\" must be an array of natural values");
    }

    std::vector<CLONES::Mutants::Evolutions::AxisPosition> corner;

    for (const auto& axis_value : corner_json) {
        corner.push_back(axis_value.template get<CLONES::Mutants::Evolutions::AxisPosition>());
    }

    if (corner.size()<2 || corner.size()>3) {
        throw std::runtime_error("\""+corner_field_name+"\" must be either a 2D or 3D array");
    }

    return corner;
}

std::pair<std::string, CLONES::Mutations::PassengerRates>
ConfigReader::get_epistate_passenger_rates(const nlohmann::json& epistate_rates_json)
{
    if (!epistate_rates_json.is_object()) {
        throw std::runtime_error("All the elements in \"passenger rates\" "
                                    "must be objects");
    }

    double indel_rate{0}, SNV_rate{0}, CNA_rate{0};
    if (epistate_rates_json.contains("indel")) {
        indel_rate = epistate_rates_json["indel"].get<double>();
    }
    if (epistate_rates_json.contains("SNV")) {
        SNV_rate = epistate_rates_json["SNV"].get<double>();
    }
    if (epistate_rates_json.contains("CNA")) {
        CNA_rate = epistate_rates_json["CNA"].get<double>();
    }

    return {
        get_from<std::string>("epistate", epistate_rates_json,
                              "All the elements in \"passenger rates\""),
        CLONES::Mutations::PassengerRates(indel_rate, SNV_rate, CNA_rate)
    };
}

CLONES::Mutations::CNA::Type
ConfigReader::get_CNA_type(const nlohmann::json& CNA_json)
{
    auto subtype_str = get_from<std::string>("subtype", CNA_json, "All the CNAs");
    for (auto& subtype_chr: subtype_str) {
        subtype_chr = tolower(subtype_chr);
    }

    if (subtype_str == "amplification") {
        return CLONES::Mutations::CNA::Type::AMPLIFICATION;
    }

    if (subtype_str == "deletion") {
        return CLONES::Mutations::CNA::Type::DELETION;
    }

    throw std::runtime_error("Unknown CNA \""+
                                CNA_json["subtype"].template get<std::string>()+"\"");
}

void
ConfigReader::add_CNA(std::list<CLONES::Mutations::CNA>& CNAs,
                      const nlohmann::json& CNA_json)
{
    using namespace CLONES::Mutations;

    auto chr_str = get_from<std::string>("chromosome", CNA_json, "All the CNAs");
    auto position = get_from<ChrPosition>("position", CNA_json, "All the CNAs");
    auto length = get_from<ChrPosition>("length", CNA_json, "All the CNAs");
    auto allele = get_from<uint8_t>("allele", CNA_json, "All the CNAs");

    GenomicPosition genomic_position(GenomicPosition::stochr(chr_str), position);

    switch (get_CNA_type(CNA_json)) {
        case CNA::Type::AMPLIFICATION:
            CNAs.push_back(CNA::new_amplification(genomic_position, length,
                                                  allele, 0, Mutation::DRIVER));
            break;
        case CNA::Type::DELETION:
            CNAs.push_back(CNA::new_deletion(genomic_position, length, allele,
                                             Mutation::DRIVER));
            break;
        default:
            throw std::runtime_error("Unsupported CNA type");
    }
}

void
ConfigReader::add_SID(const std::string& mutant_name,
                      std::list<CLONES::Mutations::MutationSpec<CLONES::Mutations::SID>>& SIDs,
                      const nlohmann::json& SID_json)
{
    using namespace CLONES::Mutations;

    std::string ref = "?";

    if (SID_json.count("ref")>0) {
        ref = get_from<std::string>("ref", SID_json, "All the mutations");
    }
    auto chr_str = get_from<std::string>("chromosome", SID_json, "All the mutations");
    auto alt = get_from<std::string>("alt", SID_json, "All the mutations");
    auto position = get_from<ChrPosition>("position", SID_json, "All the mutations");
    AlleleId allele_id = RANDOM_ALLELE;
    if (SID_json.count("allele")>0) {
        allele_id = get_from<AlleleId>("allele", SID_json, "All the mutations");
    }

    GenomicPosition genomic_position(GenomicPosition::stochr(chr_str), position);

    SIDs.emplace_back(allele_id, genomic_position, ref, alt, mutant_name);
}

std::map<std::string, CLONES::Mutations::PassengerRates>
ConfigReader::get_passenger_rates(const nlohmann::json& passenger_rates_json)
{
    std::map<std::string, CLONES::Mutations::PassengerRates> passenger_rates;

    if (!passenger_rates_json.is_array()) {
        throw std::runtime_error("The \"passenger rates\" field must be an array "
                                "of passenger mutation rates");
    }

    for (const auto& passenger_rate_json : passenger_rates_json) {
        passenger_rates.insert(get_epistate_passenger_rates(passenger_rate_json));
    }

    return passenger_rates;
}

void
ConfigReader::schedule_mutation(const std::string& mutant_name,
                                std::list<CLONES::Mutations::MutationSpec<CLONES::Mutations::SID>>& SIDs,
                                std::list<CLONES::Mutations::CNA>& CNAs,
                                const nlohmann::json& mutation_json)
{
    if (!mutation_json.is_object()) {
        throw std::runtime_error("All the elements in \"mutations\" must be objects");
    }

    auto type = get_from<std::string>("type", mutation_json,
                                      "All the elements in \"mutations\"");
    if (type=="SNV" || type=="indel") {
        add_SID(mutant_name, SIDs, mutation_json);

        return;
    }

    if (type=="CNA") {
        add_CNA(CNAs, mutation_json);

        return;
    }

    throw std::runtime_error("Unsupported mutations type \""+type+"\"");
}

double
ConfigReader::get_fraction(const nlohmann::json& fraction_json)
{
    const double fraction = fraction_json.template get<double>();

    if (fraction<0) {
        throw std::runtime_error("The \"fraction\" field must be a non-negative value");
    }

    return fraction;
}

void
ConfigReader::add_mutational_properties(CLONES::Mutations::MutationalProperties& mutational_properties,
                                        const nlohmann::json& mutational_properties_json)
{
    if (!mutational_properties_json.is_object()) {
        throw std::runtime_error("All the elements in \"driver properties\" "
                                    "must be objects");
    }

    auto mutant_name = get_from<std::string>("name", mutational_properties_json,
                                             "All the elements in \"driver properties\"");

    if (!mutational_properties_json.contains("passenger rates")) {
        throw std::runtime_error("All the elements in \"driver properties\" "
                                "must contain a \"passenger rates\" field");
    }
    auto passenger_rates = get_passenger_rates(mutational_properties_json["passenger rates"]);

    using namespace CLONES::Mutations;
    std::list<MutationSpec<SID>> SIDs;
    std::list<CNA> CNAs;
    if (mutational_properties_json.contains("driver mutations")) {
        collect_mutations(mutant_name, SIDs, CNAs, mutational_properties_json["driver mutations"]);
    }

    mutational_properties.add_mutant(mutant_name, passenger_rates, SIDs, CNAs);
}

CLONES::Mutations::MutationalExposure
ConfigReader::get_default_exposure(const std::string& mutation_name, const nlohmann::json& exposures_json)
{
    if (!exposures_json.is_object()) {
        throw std::runtime_error("The \"" + mutation_name + "\" field must be an object");
    }

    expecting("default", exposures_json, "The \"" + mutation_name + "\" field");

    return get_exposure(exposures_json["default"]);
}

std::map<double, CLONES::Mutations::MutationalExposure>
ConfigReader::get_timed_exposures(const std::string& mutation_name,
                                  const nlohmann::json& exposures_json)
{
    std::map<double, CLONES::Mutations::MutationalExposure> timed_exposures;

    if (!exposures_json.is_object()) {
        throw std::runtime_error("The \"" + mutation_name + "\" field must be an object");
    }

    if (exposures_json.contains("timed")) {

        auto& timed_coeff_json = exposures_json["timed"];

        if (!timed_coeff_json.is_array()) {
            throw std::runtime_error("The optional \"timed\" field must be an array");
        }

        for (const auto& timed_json : timed_coeff_json) {
            if (!timed_json.is_object()) {
                throw std::runtime_error("The elements of the \"timed\" field must be objects");
            }

            double time = get_from<double>("time", timed_json,
                                            "The elements of the \"timed\" field");
            if (timed_exposures.count(time)>0) {
                std::ostringstream oss;

                oss << "Two elements of the \"timed\" field have time " << time << ".";
                throw std::runtime_error(oss.str());
            }

            expecting("exposure", timed_json, "The elements of the \"timed\" field");

            timed_exposures.insert({time, get_exposure(timed_json["exposure"])});
        }
    }

    return timed_exposures;
}

Mutants::Evolutions::SampleSpecification
ConfigReader::get_sample_specification(const nlohmann::json& sample_specification_json,
                                       const std::string& field_name)
{
    if (!sample_specification_json.is_object()) {
        throw std::runtime_error("The \""+field_name+"\" field must be objects");
    }

    auto region = ConfigReader::get_sample_region(sample_specification_json, field_name);

    if (sample_specification_json.contains("name")) {
        auto name = sample_specification_json["name"].template get<std::string>();

        return Mutants::Evolutions::SampleSpecification(name, region);
    }

    return Mutants::Evolutions::SampleSpecification(region);
}

Mutants::RectangleSet
ConfigReader::get_sample_region(const nlohmann::json& sample_region_json,
                                const std::string& field_name)
{
    if (!sample_region_json.is_object()) {
        throw std::runtime_error("The \""+field_name+"\" field must be objects");
    }

    ConfigReader::expecting("lower corner", sample_region_json, field_name.c_str());
    ConfigReader::expecting("upper corner", sample_region_json, field_name.c_str());

    auto lower_vector = get_corner(sample_region_json, "lower corner");
    auto upper_vector = get_corner(sample_region_json, "upper corner");

    if (lower_vector.size() != upper_vector.size()) {
        throw std::runtime_error("The \"lower corner\" and \"upper corner\" arrays must "
                                    "have the same size");
    }

    using namespace CLONES::Mutants::Evolutions;

    if (lower_vector.size()==2) {
        return {{lower_vector[0], lower_vector[1]},
                {upper_vector[0], upper_vector[1]}};
    }

    return {{lower_vector[0], lower_vector[1], lower_vector[2]},
            {upper_vector[0], upper_vector[1], upper_vector[2]}};
}

void
ConfigReader::collect_mutations(const std::string& mutant_name,
                                std::list<CLONES::Mutations::MutationSpec<CLONES::Mutations::SID>>& SIDs,
                                std::list<CLONES::Mutations::CNA>& CNAs,
                                const nlohmann::json& mutations_json)
{
    if (!mutations_json.is_array()) {
        throw std::runtime_error("The \"mutations\" field must be an array "
                                    "of mutations");
    }

    for (const auto& mutation_json : mutations_json) {
        schedule_mutation(mutant_name, SIDs, CNAs, mutation_json);
    }
}

CLONES::Mutations::MutationalProperties
ConfigReader::get_mutational_properties(const nlohmann::json& configuration_json)
{
    using namespace CLONES::Mutations;

    MutationalProperties mutational_properties;

    expecting("mutant properties", configuration_json, "The passengers simulation configuration");

    auto& mutational_properties_json = configuration_json["mutant properties"];

    if (!mutational_properties_json.is_array()) {
        throw std::runtime_error("The \"mutant properties\" field must be an array "
                                 "of species mutational properties");
    }

    for (const auto& species_properties_json : mutational_properties_json) {
        add_mutational_properties(mutational_properties, species_properties_json);
    }

    return mutational_properties;
}

CLONES::Mutants::Evolutions::TimedEvent
get_timed_mutation(const std::map<std::string, CLONES::Mutants::MutantProperties> mutants,
                   const nlohmann::json& timed_mutation_json)
{
    using namespace CLONES::Mutants::Evolutions;

    ConfigReader::expecting("time", timed_mutation_json, "Every timed mutation description");

    const auto time = timed_mutation_json["time"].template get<Time>();

    ConfigReader::expecting("original mutant", timed_mutation_json, "Every timed mutation description");

    const auto& orig_mutant = mutants.at(timed_mutation_json["original mutant"].template get<std::string>());

    ConfigReader::expecting("mutated mutant", timed_mutation_json, "Every timed mutation description");

    const auto& mutated_mutant = mutants.at(timed_mutation_json["mutated mutant"].template get<std::string>());

    SimulationEventWrapper mutation({orig_mutant, mutated_mutant});

    return {time, mutation};
}

double get_rate(const nlohmann::json& rate_json)
{
    double rate = rate_json.template get<double>();

    if (rate < 0) {
        throw std::domain_error("the rate must be non-negative");
    }

    return rate;
}

CLONES::Mutants::CellEventType
get_cell_event_type_by_name(const std::string& event_name)
{

    for (const auto& [type, name]: Mutants::cell_event_names) {
        if (name == event_name) {
            return type;
        }
    }

    throw std::runtime_error("Unknown cell event type \""+event_name+"\"");
}

const CLONES::Mutants::Evolutions::Species&
find_species_by_name(const CLONES::Mutants::Evolutions::Simulation& simulation,
                                const std::string& name)
{
    for (const auto& species : simulation.tissue()) {
        if (name == species.get_name()) {
            return species;
        }
    }

    throw std::runtime_error("Unknown species \""+name+"\"");
}

CLONES::Mutants::Evolutions::TimedEvent
get_timed_rate_update(const CLONES::Mutants::Evolutions::Simulation& simulation,
                      const nlohmann::json& timed_rate_update_json)
{
    using namespace CLONES::Mutants::Evolutions;

    const std::string descr("Every timed rate update description");

    std::vector<std::string> fields{"time", "mutant", "status", "rate name", "rate"};

    for (const auto& field: fields) {
        ConfigReader::expecting(field, timed_rate_update_json, descr);
    }
    const auto time = timed_rate_update_json["time"].template get<Time>();

    const auto name = (timed_rate_update_json["mutant"].template get<std::string>()
                        + timed_rate_update_json["status"].template get<std::string>());

    const Species& species = find_species_by_name(simulation, name);

    const auto rate_name = timed_rate_update_json["rate name"].template get<std::string>();
    CLONES::Mutants::CellEventType cell_event_type = get_cell_event_type_by_name(rate_name);

    SimulationEventWrapper rate_update({species.get_id(), cell_event_type,
                                        get_rate(timed_rate_update_json["rate"])});

    return {time, rate_update};
}

CLONES::Mutants::Evolutions::TimedEvent
get_timed_sampling(const nlohmann::json& timed_sampling_json)
{
    using namespace CLONES::Mutants::Evolutions;

    const std::string descr("Every timed sampling description");

    std::vector<std::string> fields{"time", "sample"};

    for (const auto& field: fields) {
        ConfigReader::expecting(field, timed_sampling_json, descr);
    }
    const auto time = timed_sampling_json["time"].template get<Time>();

    auto sample_specification = ConfigReader::get_sample_specification(timed_sampling_json["sample"]);

    Sampling sampling(sample_specification);
    return {time, SimulationEventWrapper(sampling)};
}

CLONES::Mutants::Evolutions::TimedEvent
ConfigReader::get_timed_event(const CLONES::Mutants::Evolutions::Simulation& simulation,
                              const std::map<std::string, CLONES::Mutants::MutantProperties> mutants,
                              const nlohmann::json& timed_event_json)
{
    expecting("type", timed_event_json, "The timed event description");

    std::string type_name =  timed_event_json["type"].template get<std::string>();

    if (type_name == "driver mutation") {
        return get_timed_mutation(mutants, timed_event_json);
    }

    if (type_name == "liveness rate update") {
        return get_timed_rate_update(simulation, timed_event_json);
    }

    if (type_name == "sampling") {
        return get_timed_sampling(timed_event_json);
    }

    throw std::domain_error("Unsupported timed event \""+type_name+"\"");
}

void ConfigReader::expecting(const std::string& field_name, const nlohmann::json& json,
                             const std::string& context)
{
    if (!json.contains(field_name)) {
        if (context.length()>0) {
            throw std::runtime_error(context+ " must contains a \""+field_name+"\" field");
        } else {
            throw std::runtime_error("Expected a missing \""+field_name+"\" field");
        }
    }
}

void read_num_of_alleles_exception(std::map<Mutations::ChromosomeId, size_t>& num_of_alleles,
                                   const nlohmann::json& exception_json)
{
    if (!exception_json.is_object()) {
        throw std::runtime_error("Each element of the \"exceptions\" field in \"number of "
                                 "alleles\" must be an object");
    }

    auto alleles_in_chr = ConfigReader::get_from<int>("number of alleles", exception_json,
                                                      "Each element of the \"exceptions\" field in \"number of "
                                                      "alleles\"");

    if (alleles_in_chr <= 0) {
        throw std::runtime_error("The default number of alleles must be positive");
    }

    auto chromosome_name = ConfigReader::get_from<std::string>("chromosome", exception_json,
                                                               "Each element of the \"exceptions\" field in \"number of "
                                                               "alleles\"");
    using namespace Mutations;
    num_of_alleles[GenomicPosition::stochr(chromosome_name)] = alleles_in_chr;
}

void read_num_of_alleles_exceptions(std::map<Mutations::ChromosomeId, size_t>& num_of_alleles,
                                    const nlohmann::json& num_of_alleles_json)
{
    if (num_of_alleles_json.count("exceptions")==0) {
        return;
    }

    auto exceptions_json = num_of_alleles_json["exceptions"];

    if (!exceptions_json.is_array()) {
        throw std::runtime_error("The \"exceptions\" field in \"number of "
                                    "alleles\" must be an array");
    }

    for (const auto& exception_json: exceptions_json) {
        read_num_of_alleles_exception(num_of_alleles, exception_json);
    }
}

std::map<Mutations::ChromosomeId, size_t>
ConfigReader::get_number_of_alleles(const nlohmann::json& simulation_json,
                                    const std::vector<Mutations::ChromosomeId>& chromosome_ids)
{
    expecting("number of alleles", simulation_json, "The passengers simulation configuration");

    auto num_of_alleles_json = simulation_json["number of alleles"];

    if (!num_of_alleles_json.is_object()) {
        throw std::runtime_error("The \"number of alleles\" field must be an object");
    }

    int default_number = get_from<int>("default", num_of_alleles_json,
                                       "The \"number of alleles\" field");

    if (default_number <= 0) {
        throw std::runtime_error("The default number of alleles must be positive");
    }

    std::map<Mutations::ChromosomeId, size_t> num_of_alleles;
    for (const auto& chr_id : chromosome_ids) {
        num_of_alleles[chr_id] = default_number;
    }

    read_num_of_alleles_exceptions(num_of_alleles, num_of_alleles_json);

    return num_of_alleles;
}

CLONES::Mutations::MutationalExposure
ConfigReader::get_exposure(const nlohmann::json& exposure_json)
{
    if (!exposure_json.is_array()) {
        throw std::runtime_error("Exposures must be an array of object");
    }

    double total=0;

    CLONES::Mutations::MutationalExposure exposure;
    for (const auto& exposures_json : exposure_json) {

        if (!exposures_json.is_object()) {
            throw std::runtime_error("An exposure must be an object");
        }

        const std::string name = get_from<std::string>("name", exposures_json,
                                                        "Every exposure");
        if (exposure.count(name)>0) {
            throw std::runtime_error("\""+name+"\" already among exposures");
        }

        double fraction = get_from<double>("fraction", exposures_json,
                                            "Every exposure");

        total += fraction;
        if (total>1) {
            std::ostringstream oss;

            oss << "The sum of the \"fraction\" fields must be 1. "
                << "The difference is " << (1-total);

            throw std::runtime_error(oss.str());
        }

        exposure[name] = fraction;
    }

    if (std::abs(1-total)>10*std::numeric_limits<double>::epsilon()) {
        std::ostringstream oss;

        oss << "The sum of the \"fraction\" fields must be 1. "
            << "The difference is " << (1-total);

        throw std::runtime_error(oss.str());
    }

    return exposure;
}

}   // CLONES
