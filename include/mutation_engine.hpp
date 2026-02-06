/**
 * @file mutation_engine.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to place mutations on a descendants forest
 * @version 1.26
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

#ifndef __CLONES_MUTATION_ENGINE__
#define __CLONES_MUTATION_ENGINE__

#include <map>
#include <stack>
#include <random>
#include <ostream>
#include <iterator>

#include "phylogenetic_forest.hpp"
#include "mutant_id.hpp"

#include "context_index.hpp"
#include "rs_index.hpp"
#include "genome_mutations.hpp"
#include "sbs_signature.hpp"
#include "id_signature.hpp"
#include "mutational_properties.hpp"
#include "driver_storage.hpp"

#include "filter.hpp"
#include "utils.hpp"
#include "warning.hpp"


namespace CLONES
{

namespace Mutations
{

/**
 * @brief Mutation Statistics
 */
class MutationStatistics
{
    /**
     * @brief SID statistics
     */
    struct SIDStatistics
    {
        size_t mutated_alleles;            //!< Number of mutated alleles
        size_t num_of_cells;               //!< Total number of cells containing the SID

        /**
         * @brief The empty constructor
         */
        SIDStatistics();
    };

    struct SampleMutationStatistics
    {
        size_t num_of_cells;          //!< Number of recorded cells

        std::map<std::shared_ptr<SID>, SIDStatistics> SIDs;      //!< SIDs
        std::list<std::shared_ptr<CNA>> CNAs;   //!< CNAs

        /**
         * @brief The empty constructor
         */
        SampleMutationStatistics();
    };

    SampleMutationStatistics overall_statistics;    //!< Overall statistics
    std::map<std::string, SampleMutationStatistics> sample_statistics;  //!< Statistics per sample

public:
    /**
     * @brief The empty constructor
     */
    MutationStatistics();

    /**
     * @brief Record the SIDs of a cell genome
     *
     * @param sample_name is the name of the sample from which the cell has
     *          been extracted
     * @param cell_mutations are the mutations of the cell whose statistics
     *          is about to be recorded
     * @return a reference to the updated object
     */
    MutationStatistics& record(const std::string& sample_name,
                               const CellGenomeMutations& cell_mutations);

    /**
     * @brief Record the SIDs in a list of sample genomic mutations
     *
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] progress_bar is a progress bar pointer
     * @return a reference to the updated object
     */
    MutationStatistics&
    record(const std::list<CLONES::Mutations::SampleGenomeMutations>& mutations_list,
           UI::ProgressBar* progress_bar=nullptr);

    /**
     * @brief Record the SIDs in a list of sample genomic mutations
     *
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] progress_bar is a progress bar pointer
     * @return a reference to the updated object
     */
    inline MutationStatistics&
    record(const std::list<CLONES::Mutations::SampleGenomeMutations>& mutations_list,
           UI::ProgressBar& progress_bar)
    {
        return record(mutations_list, &progress_bar);
    }

    /**
     * @brief Write in a stream a table summarizing SID statistics
     *
     * @param os is the output stream
     * @param separator is the column separator
     * @return a reference to the updated stream
     */
    std::ostream& write_SIDs_table(std::ostream& os, const char separator='\t');

    /**
     * @brief Write in a stream a table summarizing CNA statistics
     *
     * @param os is the output stream
     * @param separator is the column separator
     * @return a reference to the updated stream
     */
    std::ostream& write_CNAs_table(std::ostream& os, const char separator='\t');
};

/**
 * @brief The mutation exposure type
 *
 * A mutational signature is the probability distribution of mutation types.
 * The signature depends on the environmental context and, because of that,
 * more than one signature may be active at the same time with different
 * probabilities. An exposure is a discrete probability distribution among
 * mutational signatures.
 */
using MutationalExposure = std::map<std::string, double>;

/**
 * @brief The mutation engine
 *
 * The objects of this class place mutations on the cell genome according
 * to a descendants forest
 *
 * @tparam GENOME_WIDE_POSITION is the type used to represent genome-wise position
 * @tparam RANDOM_GENERATOR is the type of random generator
 */
template<typename GENOME_WIDE_POSITION = uint32_t, typename RANDOM_GENERATOR = std::mt19937_64>
class MutationEngine
{
    /**
     * @brief The inverse cumulative signature
     *
     * The mutational signatures are the probability distributions of mutation
     * types. If we arbitrary assign an order to mutation type, given a mutation
     * and a mutation type MUTATION_TYPE, we can compute the probability of
     * randomly choosing a mutation type smaller or equal to MUTATION_TYPE.
     * The inverse cumulative signature associates such a probability to
     * MUTATION_TYPE itself for all the mutation types MUTATION_TYPE.
     *
     * @tparam MUTATION_TYPE is a mutation type
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    using InverseCumulativeSignature = std::map<double, MUTATION_TYPE>;

    /**
     * @brief A collection of inverse cumulative signatures
     *
     * @tparam MUTATION_TYPE is a mutation type
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    using InverseCumulativeSignatures = std::map<std::string, InverseCumulativeSignature<MUTATION_TYPE>>;

    /**
     * @brief The SBS inverse cumulative signature type
     */
    using InverseCumulativeSBS = InverseCumulativeSignature<SBSType>;

    /**
     * @brief The ID inverse cumulative signature type
     */
    using InverseCumulativeID = InverseCumulativeSignature<IDType>;

    /**
     * @brief A collection of SBS inverse cumulative signatures
     */
    using InverseCumulativeSBSs = InverseCumulativeSignatures<SBSType>;

    /**
     * @brief A collection of ID inverse cumulative signatures
     */
    using InverseCumulativeIDs = InverseCumulativeSignatures<IDType>;

    /**
     * @brief The stack of contexts removed from the context index
     */
    using ContextStack = std::stack<std::pair<SBSContext, GENOME_WIDE_POSITION>>;

    /**
     * @brief The stack of the repetition removed from the repetition index
     */
    using IDTypeStack = std::stack<IDType>;

    RANDOM_GENERATOR generator; //!< the random generator

    ContextIndex<GENOME_WIDE_POSITION> context_index;  //!< the genome context index
    ContextStack context_stack;                        //!< the stack of the contexts removed from `context_index`

    RSIndex rs_index;       //!< the genome repetition index
    IDTypeStack rs_stack;   //!< the stack of the repetition removed from `rs_index`

    std::map<Time, MutationalExposure> timed_exposures[2];  //!< the timed exposures
    InverseCumulativeSBSs inv_cumulative_SBSs;              //!< the inverse cumulative SBS
    InverseCumulativeIDs inv_cumulative_IDs;                //!< the inverse cumulative INDEL

    MutationalProperties mutational_properties;  //!< the species mutational properties

    std::shared_ptr<GenomeMutations> germline_mutations; //!< the germline mutations

    GenomeMutations dm_genome;          //!< the driver mutations genome

    DriverStorage driver_storage;       //!< the driver storage

    std::vector<CNA> passenger_CNAs;    //!< the admissible passenger CNAs

    unsigned int driver_CNA_min_distance;         //!< the minimum distance between a driver mutation and a passenger CNA

    WarningFunction warning;            //!< the warning function

    /**
     * @brief Select a random value in a set
     *
     * @tparam T is the type of objects in the set
     * @param values is a list of values from which a random value must be selected
     * @return a random value among those in `values`
     */
    template<typename T>
    T select_random_value(const std::list<T>& values)
    {
        std::uniform_int_distribution<> u_dist(0,values.size()-1);
        size_t pos = u_dist(generator);

        size_t index=0;
        auto it = values.begin();
        while (index!=pos) {
            ++it;
            ++index;
        }

        return *it;
    }

    /**
     * @brief Compute the inverse cumulative signature

     * The mutational signatures are the probability distributions of mutation
     * types. If we arbitrary assign an order to mutation type, given a mutation
     * and a mutation type T, we can compute the probability of randomly choosing
     * a mutation type smaller or equal to T. The inverse cumulative signature
     * associates such a probability to T itself for all the mutation types T.
     *
     * @param signature is a mutation signature
     * @return the inverse cumulative signature of `signature`
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    static InverseCumulativeSignature<MUTATION_TYPE>
    get_inv_cumulative_signature(const Signature<MUTATION_TYPE>& signature)
    {
        InverseCumulativeSignature<MUTATION_TYPE> inv_cumulative;
        double cumulative_prob = 0;
        for (const auto& [m_type, prob]: signature) {
            if (prob != 0) {
                cumulative_prob += prob;

                inv_cumulative.insert({cumulative_prob, m_type});
            }
        }

        return inv_cumulative;
    }

    /**
     * @brief Get the exposure according to a cell birth time
     *
     * The active exposure depend on the simulation time. It must be selected in
     * agreement with cell birth times.
     *
     * @tparam MUTATION_TYPE is the type of the mutation whose exposures are aimed
     * @param mutation_type is the type of mutation affected by the exposure
     * @param cell is the cell whose birth time select the exposure
     * @return a constant reference to exposure that are associated
     *         to `cell` birth time
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    const MutationalExposure&
    get_active_exposure(const Mutants::Cell& cell) const
    {
        auto exposure_it = get_timed_exposures(MUTATION_TYPE::type()).upper_bound(cell.get_birth_time());

        return (--exposure_it)->second;
    }

    /**
     * @brief Select a random SBS
     *
     * This method selects a random SBS among whose having a specified
     * SBS type and available in a context index. The selected SBS is
     * extracted from the context index and inserted into a stack to
     * revert the selection.
     *
     * @param[in] m_type is the SBS type of the SBS to be selected
     * @return a SBS whose type is `m_type` and which was available in
     *          the context index
     */
    MutationSpec<SID> select_SID(const SBSType& m_type)
    {
        using namespace CLONES::Mutations;

        SBSContext context = m_type.get_context();
        SBSContext compl_context = context.get_reverse_complement();

        size_t total_pos = context_index[context].size();
        total_pos += context_index[compl_context].size();

        size_t index;
        {
            std::uniform_int_distribution<> u_dist(0,total_pos-1);
            index = u_dist(generator);
        }

        char alt_base = m_type.get_replace_base();
        if (index >= context_index[context].size()) {
            index -= context_index[context].size();
            context = compl_context;

            alt_base = GenomicSequence::get_complement(alt_base);
        }

        GENOME_WIDE_POSITION pos;

        if (infinite_sites_model) {
            pos = context_index.extract(context, index);

            context_stack.push({context, pos});
        } else {
            pos = context_index[context][index];
        }
        auto genomic_pos = context_index.get_genomic_position(pos);

        return {RANDOM_ALLELE, genomic_pos.chr_id, genomic_pos.position,
                context.get_central_nucleotide(), alt_base, Mutation::UNDEFINED};
    }

    /**
     * @brief Select a random indel
     *
     * This method selects a random indel among whose having a specified ID type
     * and available in the repetition index. The selected indel is extracted
     * from the repetition index and inserted into a stack to revert the selection.
     *
     * @param[in] id_type is the ID type of the indel to be selected
     * @return an indel whose type is `m_type` and which was available in the
     *      repeated sequence index
     */
    MutationSpec<SID> select_SID(const IDType& id_type)
    {
        using namespace CLONES::Mutations;

        if (infinite_sites_model) {
            rs_stack.push(id_type);
        }
        const auto& repetition = rs_index.select(id_type, infinite_sites_model);

        std::string alt(1, repetition.prev_base);
        std::string ref = alt;
        ref.append(repetition.unit);

        if (id_type.insertion) {
            std::swap(alt, ref);
        }

        GenomicPosition pos(repetition.g_position);
        --pos.position;

        return {RANDOM_ALLELE, pos, ref, alt, Mutation::UNDEFINED};
    }

    /**
     * @brief Randomly select the number of mutations according to a Poisson distribution
     *
     * @param genome_size is the genome size
     * @param passenger_mutation_rate is the rate of the passanger mutations
     * @return is the randomly selected number of mutations
     */
    size_t number_of_mutations(GenomeMutations::Length genome_size,
                               const double& passenger_mutation_rate)
    {
        auto mean = genome_size*passenger_mutation_rate;

        std::poisson_distribution<> p_dist(mean);

        return p_dist(generator);
    }

    /**
     * @brief Place the CNAs in a cell genome in the phylogenetic forest
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param num_of_mutations is the number of CNA to apply
     * @param nature is the nature of the mutations
     */
    void place_CNAs(PhylogeneticForest::node* node, GenomeMutations& cell_mutations,
                    const size_t& num_of_mutations, const Mutation::Nature& nature)
    {
        // while some of the requested CNAs have not been set
        size_t new_CNAs{0};
        size_t available{passenger_CNAs.size()};
        while (new_CNAs < num_of_mutations) {
            const size_t last_pos = available-1;
            std::uniform_int_distribution<size_t> u_dist(0, last_pos);

            // randomly select a CNA
            size_t index = u_dist(generator);

            bool to_be_placed = true;
            while (to_be_placed) {
                CNA cna = passenger_CNAs[index];

                cna.nature = nature;
                if (place_CNA(cna, cell_mutations)) {
                    // if the CNA has been successfully applied,
                    // then increase the number of new CNAs
                    ++new_CNAs;

                    if (node != nullptr) {
                        node->add_new_mutation(cna);
                    }

                    to_be_placed = false;
                } else {
                    if (--available == 0) {
                        throw std::runtime_error("None of the "
                                                + std::to_string(passenger_CNAs.size())
                                                + " admitted CNAs can be applied.");
                    }

                    if (index<available) {
                        // here passenger_CNAs is shuffled and it will not be
                        // resorted as at the beginning of the computation, but
                        // this is still acceptable because the CNAs are
                        // randomly selected
                        std::swap(passenger_CNAs[index], passenger_CNAs[available]);
                    } else {
                        index = available-1;
                    }
                }
            }
        }
    }

    /**
     * @brief Place the CNAs in a cell genome in the phylogenetic forest
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param CNA_rate is the rate of passegers CNAs
     */
    void place_CNAs(PhylogeneticForest::node& node, GenomeMutations& cell_mutations,
                    const double& CNA_rate)
    {
        const auto num_of_CNAs = number_of_mutations(cell_mutations.allelic_size(), CNA_rate);

        place_CNAs(&node, cell_mutations, num_of_CNAs, Mutation::PASSENGER);
    }

    /**
     * @brief Get the available alleles for placing a SID
     *
     * @param mutation is the mutation to be placed
     * @param cell_mutations are the non-germinal mutations of a cell
     * @return a list of the identifiers of the available alleles
     */
    std::list<AlleleId> get_available_alleles(const SID& mutation,
                                              const GenomeMutations& cell_mutations) const
    {
        std::list<AlleleId> available_alleles;

        if (infinite_sites_model) {
            const auto& chr_mutations = cell_mutations.get_chromosome(mutation.chr_id);
            if ((dm_genome.has_context_free(mutation)
                    || (mutation.nature == Mutation::DRIVER))
                  && chr_mutations.has_context_free(mutation)
                  && germline_mutations->has_context_free(mutation)) {

                for (const auto& [allele_id, allele]: chr_mutations.get_alleles()) {
                    if (allele.strictly_contains(mutation)) {
                        available_alleles.push_back(allele_id);
                    }
                }
            }

            return available_alleles;
        }

        auto germline_alleles = germline_mutations->get_alleles_with_context_free_for(mutation);
        auto driver_alleles = dm_genome.get_alleles_with_context_free_for(mutation);

        if (germline_alleles.size()==0 || driver_alleles.size()==0) {
            return available_alleles;
        }

        std::set<AlleleId> germline_set(germline_alleles.begin(), germline_alleles.end());
        std::set<AlleleId> driver_set(driver_alleles.begin(), driver_alleles.end());

        const auto& chr_mutations = cell_mutations.get_chromosome(mutation.chr_id);

        for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
            if (allele.has_context_free(mutation)
                    && germline_set.count(allele.get_history().front())>0
                    && (driver_set.count(allele.get_history().front())>0
                        || (mutation.nature == Mutation::DRIVER))) {
                available_alleles.push_back(allele_id);
            }
        }

        return available_alleles;
    }

    /**
     * @brief Place a passenger SID
     *
     * @param[in,out] mutation is the SID to place
     * @param[in,out] cell_mutations are the cell mutations
     * @return `true` if and only if the SID placement has succeed. This
     *          method returns `false` when the SID cannot be placed
     *          because its context is not free or no allele are available
     *          for it.
     */
    bool place_SID(MutationSpec<SID>& mutation, GenomeMutations& cell_mutations)
    {
        if (infinite_sites_model && mutation.nature != Mutation::DRIVER) {
            if (!dm_genome.has_context_free(mutation)) {
                return false;
            }
        }

        if (mutation.allele_id == RANDOM_ALLELE) {
            auto candidate_alleles = get_available_alleles(mutation, cell_mutations);

            if (candidate_alleles.size()==0) {
                return false;
            }

            mutation.allele_id = select_random_value(candidate_alleles);
        }

        // try to apply the selected SID
        return cell_mutations.apply(mutation);
    }

    /**
     * @brief Get the number of indiced SBSs of a given type
     *
     * @param mutation_type is the type of the aimed SBS
     * @return the number of indiced SBSs of type `mutation_type`
     */
    inline size_t count_available(const SBSType& mutation_type) const
    {
        return context_index[mutation_type.get_context()].size();
    }

    /**
     * @brief Get the number of indiced indels of a given type
     *
     * @param mutation_type is the type of the aimed SBS
     * @return the number of indiced SBSs of type `mutation_type`
     */
    inline size_t count_available(const IDType& mutation_type) const
    {
        return rs_index.count_available_for(mutation_type);
    }

    /**
     * @brief Get the inverse cumulative signatures of a type
     *
     * @tparam MUTATION_TYPE is the mutation type
     * @return the inverse cumulative signatures of the type
     *      `MUTATION_TYPE`
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    const InverseCumulativeSignatures<MUTATION_TYPE>&
    get_inverse_cumulative_signatures() const
    {
        if constexpr(std::is_base_of_v<MUTATION_TYPE, SBSType>) {
            return inv_cumulative_SBSs;
        }

        if constexpr(std::is_base_of_v<MUTATION_TYPE, IDType>) {
            return inv_cumulative_IDs;
        }

        throw std::runtime_error("Unsupported mutation type.");
    }

    /**
     * @brief Place the SIDs in a cell genome in the phylogenetic forest
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param signature_name is the name of the signature used to find the mutations
     * @param num_of_mutations is the number of SIDs to apply
     * @param nature is the nature of the mutations
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    void place_SIDs(PhylogeneticForest::node* node, GenomeMutations& cell_mutations,
                    const std::string& signature_name, const size_t& num_of_mutations,
                    const Mutation::Nature& nature)
    {
        // get the inverse cumulative signature
        const auto& inv_cum_signatures = get_inverse_cumulative_signatures<MUTATION_TYPE>();
        const auto& inv_cum_signature = inv_cum_signatures.at(signature_name);

        std::uniform_real_distribution<> u_dist(0,1);

        // while some of the requested mutations have not been set
        size_t new_mutations{0};
        while (new_mutations < num_of_mutations) {
            // randomly pick a mutation type according to the current signature
            auto it = inv_cum_signature.lower_bound(u_dist(generator));

            if ((it != inv_cum_signature.end())
                    && (count_available(it->second) > 0)) {

                // select a mutation among those having the selected mutation
                // type and that are available in the corresponding index
                auto sid = select_SID(it->second);

                sid.cause = signature_name;
                sid.nature = nature;

                if (place_SID(sid, cell_mutations)) {
                    // if the mutation has been successfully applied,
                    // then increase the number of new mutations
                    ++new_mutations;

                    if (node != nullptr) {
                        node->add_new_mutation(sid);
                    }
                }
            } else {
                std::ostringstream oss;
                if constexpr(std::is_base_of_v<SBSType, MUTATION_TYPE>) {
                    oss << "No mutations for context \"" << it->second
                        << "\".";
                    warning(WarningType::NO_MUT_FOR_CONTEXT, oss.str());
                }
                
                if constexpr(std::is_base_of_v<IDType, MUTATION_TYPE>) {
                    oss << "No mutations for repetition pattern \""
                        << it->second << "\".";
                    warning(WarningType::NO_MUT_FOR_RPATTERN, oss.str());
                }
            }
        }
    }

    /**
     * @brief Get the passenger rate of a specific mutation type
     *
     * @tparam MUTATION_TYPE is the mutation type of the aimed rate
     * @param rates are the passenger rates
     * @return the passenger rate of the mutation type `MUTATION_TYPE`
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    static const double& get_rate(const PassengerRates& rates)
    {
        if constexpr(std::is_base_of_v<MUTATION_TYPE, SBSType>) {
            return rates.snv;
        }

        if constexpr(std::is_base_of_v<MUTATION_TYPE, IDType>) {
            return rates.indel;
        }

        throw std::runtime_error("Unsupported mutation type.");
    }

    /**
     * @brief Place mutations on the genome of a phylogenetic forest cell
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param rates are the passenger mutation rates
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    void place_passengers(PhylogeneticForest::node& node, GenomeMutations& cell_mutations,
                          const PassengerRates& rates)
    {
        auto remaining_mutations = number_of_mutations(cell_mutations.allelic_size(),
                                                       get_rate<MUTATION_TYPE>(rates));

        const auto& exposure = get_active_exposure<MUTATION_TYPE>(node);

        auto exposure_it = exposure.begin();
        double missing_signature_prob{1.0};
        while (remaining_mutations > 0 && exposure_it != exposure.end()) {

            // get name and probability of the current signature
            const auto& signature_name = exposure_it->first;
            const auto& signature_prob = exposure_it->second;
            
            // update exposure pointer
            ++exposure_it;

            // evaluate the choice probability for the current signature
            double choice_prob;
            if (exposure_it != exposure.end()) {
                choice_prob = signature_prob/missing_signature_prob;
            } else {

                // this branch is to avoid numeric approximation errors
                choice_prob = 1.0;
            }

            // evaluate the number of mutation associated to the current signature
            std::binomial_distribution<size_t> b_dist(remaining_mutations,
                                                      choice_prob);
            const size_t num_of_signature_muts = b_dist(generator);

            place_SIDs<MUTATION_TYPE>(&node, cell_mutations, signature_name,
                                      num_of_signature_muts, Mutation::PASSENGER);

            remaining_mutations -= num_of_signature_muts;
            missing_signature_prob -= signature_prob;
        }
    }

    /**
     * @brief Seach for a random allele in a chromosome containing a region
     *
     * @param[out] allele_id is the identifier of the found allele
     * @param[in] chr_mutations is the chromosome in which the allele is searched
     * @param[in] region is the region that must be contained by the found allele
     * @param[in] no_driver_mutations is a Boolean flag to establish whether
     *          alleles containing driver mutations in the region are discharged
     * @param[in] multiple_alleles_required is a Boolean flag to require more
     *          than one allele covering the region
     * @return `true` if and only if one of the alleles of `chr_mutations`
     *          contains the region `region`. In this case, the parameter
     *          `allele_id` is set to the identifier of such an allele.
     *          When `no_driver_mutations` is set to `true`, the search
     *          must avoid alleles containing driver mutations in `region`.
     */
    bool select_allele_containing(AlleleId& allele_id,
                                  const ChromosomeMutations& chr_mutations,
                                  const GenomicRegion& region,
                                  const bool& no_driver_mutations,
                                  const bool& multiple_alleles_required)
    {
        const auto& chr_alleles = chr_mutations.get_alleles();

        if (chr_alleles.size()==0) {
            return false;
        }

        // choose a random element in the map
        auto it = chr_alleles.begin();
        std::advance(it, generator() % chr_alleles.size());

        size_t available_alleles{0};
        bool found{false};

        // cyclically iterate over the map from the choosen element
        for (size_t i=0; i<chr_alleles.size(); ++i, ++it) {

            // if the map end has been reached, restart from the begin
            if (it == chr_alleles.end()) {
                it = chr_alleles.begin();
            }

            const auto& allele = it->second;

            // if the allele contains the region
            if (allele.contains(region)) {
                ++available_alleles;

                if (found) {
                    // if an allele satisfying a the required conditions
                    // was found, this is not the first allele containing
                    // the region
                    return true;
                } else {
                    // if either the driver mutation removal is allowed
                    // or the region does not contain any driver mutations
                    if (!(no_driver_mutations
                            && allele.has_driver_mutations_in(region))) {

                        // if there are more than one allele containing
                        // the region, this is the selected allele
                        allele_id = it->first;
                        found = true;
                    }
                }
            }
        }

        return (found && (!multiple_alleles_required
                          || available_alleles>1));
    }

    /**
     * @brief Seach for an allele in which the CNA can be applied
     *
     * @param[in,out] cna is the CNA which should be applied
     * @param[in] chr_mutations is the chromosome in which the allele is searched
     * @param[in] no_driver_mutations is a Boolean flag to establish whether
     *          alleles containing driver mutations in the region are discharged
     * @return `true` if and only if one of the alleles of `chr_mutations`
     *          admits the application of `cna`.
     *          When `no_driver_mutations` is set to `true`, the search
     *          must avoid alleles containing driver mutations in `region`.
     *          If `cna` is an amplification and the source is `RANDOM_ALLELE` or
     *          the `cna` is a deletion and the destination is `RANDOM_ALLELE`,
     *          then, when `true` is returned the source or the destination of
     *          `cna`, respectively, are updated to the identifier of the found
     *          allele.
     */
    bool find_allele_for(CNA& cna, ChromosomeMutations& chr_mutations,
                         const bool& no_driver_mutations)
    {
        const auto cna_region = cna.get_region();
        switch(cna.type) {
            case CNA::Type::AMPLIFICATION:
                if (cna.source == RANDOM_ALLELE) {
                    return select_allele_containing(cna.source, chr_mutations,
                                                    cna_region, no_driver_mutations,
                                                    false);
                }
                return true;
            case CNA::Type::DELETION:
                if (cna.dest == RANDOM_ALLELE) {
                    return select_allele_containing(cna.dest, chr_mutations,
                                                    cna_region, no_driver_mutations,
                                                    avoid_homozygous_losses);
                }
                return true;
            default:
                throw std::runtime_error("Unsupported CNA type");
        }
    }

    /**
     * @brief Try to place a CNA
     *
     * @param[in,out] cna is the CNA to place
     * @param[in] cell_mutations are the cell mutations
     * @return `true` if and only if the CNA placement has succeed. This method
     *          returns `false` when the CNA cannot be placed because no allele
     *          are available for it.
     *          When `no_driver_mutations` is set to `true`, the search
     *          must avoid alleles containing driver mutations in `region`.
     *          If `cna` is an amplification and the source is `RANDOM_ALLELE` or
     *          the `cna` is a deletion and the destination is `RANDOM_ALLELE`,
     *          then, when `true` is returned the source or the destination of
     *          `cna`, respectively, are updated to the identifier of the found
     *          allele.
     */
    bool place_CNA(CNA& cna, GenomeMutations& cell_mutations)
    {
        auto& chr_mutations = cell_mutations.get_chromosome(cna.chr_id);
        if (!find_allele_for(cna, chr_mutations, (cna.nature != Mutation::DRIVER))) {
            return false;
        }

        const CNA& const_cna = cna;

        return chr_mutations.apply(const_cna);
    }

    /**
     * @brief Place a mutation
     *
     * @tparam MUTATION_TYPE is the type of mutation to be placed
     * @param mutation is the mutation to be placed
     * @param cell_mutations is the cell genome where the mutation must be placed
     * @return `true` if and only if the mutation was correctly placed
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<Mutation, MUTATION_TYPE>, bool> = true>
    bool place_mutation(MUTATION_TYPE& mutation, GenomeMutations& cell_mutations)
    {
        if constexpr (std::is_base_of_v<CNA, MUTATION_TYPE>) {
            return place_CNA(mutation, cell_mutations);
        }

        if constexpr (std::is_base_of_v<SID, MUTATION_TYPE>) {
            return place_SID(mutation, cell_mutations);
        }

        throw std::runtime_error("Unsupported mutation type.");
    }

    /**
     * @brief Place a driver mutation
     *
     * @tparam MUTATION_TYPE is the type of mutation
     * @param mutation is the placed mutation
     * @param node is the node in which the mutation first occurred
     * @param mutant_name is the name of the mutant which the mutation is
     *   a driver of
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<Mutation, MUTATION_TYPE>, bool> = true>
    void place_driver_mutation(MUTATION_TYPE& mutation, PhylogeneticForest::node& node,
                               GenomeMutations& cell_mutations,
                               const std::string& mutant_name)
    {
        mutation.nature = Mutation::DRIVER;
        mutation.cause = mutant_name;

        if (place_mutation(mutation, cell_mutations)) {
            node.add_new_mutation(mutation);
        } else {
            std::ostringstream oss;

            oss << mutant_name << "'s driver mutation "
                << mutation << " cannot be placed.";

            throw std::runtime_error(oss.str());
        }
    }

    /**
     * @brief Place the driver mutations
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param driver_mutations is the map associating mutant ids and mutations
     */
    void place_driver_mutations(PhylogeneticForest::node& node, GenomeMutations& cell_mutations,
                                std::map<Mutants::MutantId, DriverMutations>& driver_mutations)
    {
        if (node.is_root() || node.get_mutant_id()!=node.parent().get_mutant_id()) {
            auto& mutant_mp = driver_mutations.at(node.get_mutant_id());

            auto SID_it = mutant_mp.SIDs.begin();
            auto CNA_it = mutant_mp.CNAs.begin();

            for (auto order_it = mutant_mp.application_order.begin();
                 order_it != mutant_mp.application_order.end(); ++order_it) {

                switch(*order_it) {
                    case DriverMutations::SID_TURN:
                    {
                        place_driver_mutation(*SID_it, node, cell_mutations, mutant_mp.name);

                        ++SID_it;
                    }
                    break;
                    case DriverMutations::CNA_TURN:
                    {
                        place_driver_mutation(*CNA_it, node, cell_mutations, mutant_mp.name);

                        ++CNA_it;
                    }
                    break;
                    case DriverMutations::WGD_TURN:
                        cell_mutations.duplicate_alleles();
                        break;
                    default:
                        throw std::runtime_error("Unsupported driver mutation type.");
                }
            }
        }
    }

    /**
     * @brief Place the mutations on the genomes of a descendants forest node
     *
     * This method recursively places the mutations on the nodes of a descendants
     * forest. All the genome mutations associated to the forest leaves are
     * collected and saved in a container that partition them according
     * to the cell sample.
     *
     * @param[in] node is a phylogenetic forest node representing a cell
     * @param[in] ancestor_mutations is the genomic mutation of the ancestor
     * @param[in] passenger_rates is the map associating species ids to their passenger rates
     * @param[in] driver_mutations is the map associating mutant ids and mutations
     * @param[in] ancestor_mutations is the genomic mutation of the ancestor
     * @param[in,out] visited_nodes is the number of visited nodes
     * @param[in,out] progress_bar is a progress bar pointer
     */
    void place_mutations(PhylogeneticForest::node& node, const GenomeMutations& ancestor_mutations,
                         const std::map<Mutants::SpeciesId, PassengerRates>& passenger_rates,
                         std::map<Mutants::MutantId, DriverMutations>& driver_mutations,
                         size_t& visited_nodes, UI::ProgressBar *progress_bar)
    {
        using namespace CLONES::Mutations;

        const size_t context_stack_size = context_stack.size();
        const size_t rs_stack_size = rs_stack.size();

        CellGenomeMutations cell_mutations(static_cast<const Mutants::Cell&>(node),
                                           ancestor_mutations);

        place_driver_mutations(node, cell_mutations, driver_mutations);

        auto rates_it = passenger_rates.find(node.get_species_id());
        if (rates_it == passenger_rates.end()) {
            const auto& species_name = node.get_forest().get_species_name(node.get_species_id());

            throw std::runtime_error("Unknown species \""+ species_name +"\"");
        }
        place_passengers<SBSType>(node, cell_mutations, rates_it->second);
        place_passengers<IDType>(node, cell_mutations, rates_it->second);
        place_CNAs(node, cell_mutations, rates_it->second.cna);

        ++visited_nodes;
        if (progress_bar != nullptr) {

            size_t percentage = (100*visited_nodes)/(node.get_forest()).num_of_nodes();
            if (percentage>progress_bar->get_progress()) {
                progress_bar->set_progress(percentage);
            }
        }

        if (node.is_leaf()) {
            auto mut_ptr = std::make_shared<CellGenomeMutations>();

            std::swap(cell_mutations, *mut_ptr);

            auto& phylo_forest = node.get_forest();
            auto cell_id = static_cast<const Mutants::Cell&>(node).get_id();

            phylo_forest.leaves_mutations[cell_id] = mut_ptr;
        } else {
            for (auto child: node.children()) {
                place_mutations(child, cell_mutations, passenger_rates, driver_mutations,
                                visited_nodes, progress_bar);
            }
        }

        if (!infinite_sites_model) {
            // reverse context index extractions
            while (get_stack_size<SBSType>() > context_stack_size) {
                restore_last_extracted_from_index<SBSType>();
            }

            // reverse repetition index extractions
            while (get_stack_size<IDType>() > rs_stack_size) {
                restore_last_extracted_from_index<IDType>();
            }
        }
    }

    /**
     * @brief Get the size of the stack associated to a mutation type
     *
     * @tparam MUTATION_TYPE is the mutation type of the
     *      of the stack
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    size_t get_stack_size() const
    {
        if constexpr(std::is_base_of_v<MUTATION_TYPE, SBSType>) {
            return context_stack.size();
        }

        if constexpr(std::is_base_of_v<MUTATION_TYPE, IDType>) {
            return rs_stack.size();
        }

        throw std::runtime_error("Unsupported mutation type.");
    }

    /**
     * @brief Re-insert the last extracted from an index
     *
     * @tparam MUTATION_TYPE is the mutation type of the
     *      of the index
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    void restore_last_extracted_from_index()
    {
        if constexpr(std::is_base_of_v<MUTATION_TYPE, SBSType>) {
            const auto& top_stack = context_stack.top();
            context_index.insert(top_stack.first, top_stack.second);
            context_stack.pop();

            return;
        }

        if constexpr(std::is_base_of_v<MUTATION_TYPE, IDType>) {
            rs_index.restore(rs_stack.top());
            rs_stack.pop();

            return;
        }

        throw std::runtime_error("Unsupported mutation type.");
    }

    /**
     * @brief Get the species rate map
     *
     * This method returns a map from species identifier to the corresponding
     * passenger rates.
     *
     * @param descendants_forest is a descendent forest
     * @return a map associating a species in `descendants_forest` to its
     *          passenger rates
     */
    std::map<Mutants::SpeciesId, PassengerRates>
    get_species_rate_map(const Mutants::DescendantsForest& descendants_forest) const
    {
        using namespace CLONES::Mutants;

        std::map<SpeciesId, PassengerRates> species_rates;

        for (const auto& [species_id, species_data] : descendants_forest.get_species_data()) {
            const auto mutant_name = descendants_forest.get_mutant_name(species_data.mutant_id);
            const auto species_name = mutant_name
                                      + MutantProperties::signature_to_string(species_data.signature);

            auto passenger_it = mutational_properties.get_passenger_rates().find(species_name);

            if (passenger_it == mutational_properties.get_passenger_rates().end()) {
                throw std::runtime_error("\""+species_name+"\" has no mutational rate");
            }

            species_rates[species_id] = passenger_it->second;
        }

        return species_rates;
    }

    /**
     * @brief Get the driver mutation map
     *
     * This method returns a map from the mutant identifier to the correspoding
     * driver mutations.
     *
     * @param descendants_forest is a descendent forest
     * @return a map associating the mutants in `descendants_forest` to the
     *              correspoding driver mutations.
     */
    std::map<Mutants::MutantId, DriverMutations>
    get_driver_mutation_map(const Mutants::DescendantsForest& descendants_forest) const
    {
        using namespace CLONES::Mutants;

        std::map<MutantId, DriverMutations> driver_mutations;

        for (const auto& [species_id, species_data] : descendants_forest.get_species_data()) {
            const auto mutant_name = descendants_forest.get_mutant_name(species_data.mutant_id);

            auto driver_it = mutational_properties.get_driver_mutations().find(mutant_name);

            if (driver_it == mutational_properties.get_driver_mutations().end()) {
                throw std::runtime_error("\"" + mutant_name + "\" is unknown");
            }
            driver_mutations[species_data.mutant_id] = driver_it->second;
        }

        return driver_mutations;
    }

    /**
     * @brief Filter a list of CNAs according the chromosomes in the context index
     *
     * This method returns the vector of CNAs that are contained in the parameter
     * vector and lay in one of the chromosomes mentioned in the context index.
     *
     * @param CNAs is a list of CNAs
     * @return the vector of CNAs that are contained in `CNAs` and lay in one of
     *          the chromosomes mentioned in the context index.
     */
    std::vector<CNA>
    filter_CNA_by_chromosome_ids(const std::vector<CNA>& CNAs) const
    {
        auto chr_ids = context_index.get_chromosome_ids();

        std::set<ChromosomeId> chr_id_set(chr_ids.begin(),chr_ids.end());

        std::vector<CNA> filtered;
        for (const auto& cna : CNAs) {
            if (chr_id_set.count(cna.chr_id)>0) {
                filtered.push_back(cna);
            }
        }

        return filtered;
    }

    /**
     * @brief Check germline-reference sequence consistency
     *
     * This method checks whether germline contains the reference sequence chromosomes
     * and throws a `std::out_of_range` exception this is not the case.
     *
     * @param context_index is a context index
     * @param germline_mutations are the germline mutations
     */
    static void check_genomes_consistency(const ContextIndex<GENOME_WIDE_POSITION>& context_index,
                                          const GenomeMutations& germline_mutations)
    {
        const auto context_chr_ids = context_index.get_chromosome_ids();
        const auto& germline_chromosomes = germline_mutations.get_chromosomes();

        for (const auto chr_id : context_chr_ids) {
            if (germline_chromosomes.find(chr_id) == germline_chromosomes.end()) {
                throw std::out_of_range("The germline does not contain Chr. "
                                        + GenomicPosition::chrtos(chr_id)
                                        + " which is present in the reference.");
            }
        }
    }

    /**
     * @brief Get the exposures
     *
     * @return a constant reference to the exposures
     */
    inline std::map<Time, MutationalExposure>&
    get_timed_exposures(const MutationType::Type& mutation_type)
    {
        return timed_exposures[static_cast<size_t>(mutation_type)];
    }

    /**
     * @brief Split mutational exposure by mutation type
     *
     * @param exposure is an exposure
     * @return a map from mutation type to the corresponding mutational exposure
     */
    std::map<MutationType::Type, MutationalExposure>
    split_exposures(const MutationalExposure& exposure) const
    {
        std::map<MutationType::Type, MutationalExposure> exposures;

        for (const auto& [name, coeff]: exposure) {
            if (inv_cumulative_SBSs.count(name)>0) {
                exposures[MutationType::Type::SBS][name]=coeff;
            } else if (inv_cumulative_IDs.count(name)>0) {
                exposures[MutationType::Type::INDEL][name]=coeff;
            } else {
                throw std::runtime_error("Unknown signature " + name + ".");
            }
        }

        for (const auto& [m_type, s_exposure]: exposures) {
            std::string desc = std::string(m_type==MutationType::Type::SBS? "SNV": "indel")
                                + " exposure";
            validate_probability(s_exposure, desc);
        }

        return exposures;
    }

    /*
     * @brief Validate a probability map
     *
     * This static method validates a probability by testing whether its
     * values sum up to 1. If this is not the case, this method
     * throw a domain error exception.
     *
     * @param prob_map is a map whose values are the key probabilities
     * @param prob_name is the name of the probability map
     * @param approx_error is the approximation error
     */
    template<typename KEYS, typename VALUES>
    static void validate_probability(const std::map<KEYS, VALUES>& prob_map,
                                     const std::string& prob_name,
                                     const VALUES& approx_error=1e-13)
    {
        VALUES sum{0};
        for (const auto& [key, value]: prob_map) {
            sum += value;
        }

        if (std::abs(1-sum)>approx_error) {
            std::ostringstream oss;

            oss << "The " << prob_name << " must sum up to 1: " << prob_map
                << " sums up to " << sum << ".";

            throw std::domain_error(oss.str());
        }
    }

    /**
     * @brief Check whether any species has a mutation rate greater than 0
     *
     * @param mutation_type is a mutation type (i.e., INDEL or SBS)
     * @return `true` if and only if any specified species has a rate
     *  for `mutation_type` greater than 0
     */
    bool may_place_passengers(const MutationType::Type& mutation_type) const
    {
        for (const auto& [species, passenger_rates] :
                    mutational_properties.get_passenger_rates()) {
            switch (mutation_type) {
                case MutationType::Type::INDEL:
                    if (passenger_rates.indel>0) {
                        return true;
                    }
                    break;
                case MutationType::Type::SBS:
                    if (passenger_rates.snv>0) {
                        return true;
                    }
                    break;
                default:
                    throw std::domain_error("Unsupported mutation type "
                                            + std::to_string(static_cast<size_t>(mutation_type)));
            }
        }

        return false;
    }

    /**
     * @brief Place SNV and indels in normal cells
     *
     * @tparam MUTATION_TYPE is the type of the SIDs
     * @param cell_mutations is the cell genome
     * @param signature_name is the name of the signature of the SIDs
     * @param num_of_mutations_dist is a distribution for the number of
     *      SIDs to be placed
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    inline void place_SIDs_in_normal(GenomeMutations& cell_mutations,
                                     const std::string& signature_name,
                                     std::normal_distribution<double>& num_of_mutations_dist)
    {
        auto num_mutations = static_cast<size_t>(std::max(0.0, num_of_mutations_dist(generator)));
        place_SIDs<MUTATION_TYPE>(nullptr, cell_mutations, signature_name,
                                  num_mutations, Mutation::PRENEOPLASTIC);
    }

    /**
     * @brief Get the region list of a list of mutations
     *
     * This method returns the region list of a list of mutations being
     * either SIDs or CNAs. The returned regions are expanded by a shadow.
     *
     * @tparam MUTATION_TYPE is the type of mutations
     * @param mutations is the list of mutations whose regions are aimed
     * @param expansion_size is the region expansion size
     * @return a list of regions of `mutations` expanded by
     *      `expansion_size`
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<SID, MUTATION_TYPE>
                                || std::is_base_of_v<CNA, MUTATION_TYPE> , bool> = true>
    static std::list<GenomicRegion>
    get_region_list(const std::list<MUTATION_TYPE>& mutations,
                    const unsigned int& expansion_size=0)
    {
        std::list<GenomicRegion> region_list;

        for (const auto& mutation: mutations) {
            GenomicRegion::Length region_length(mutation.size_in_ref());
            ChrPosition region_begin;

            if (mutation.begin()<=expansion_size) {
                region_begin = 1;
                region_length += mutation.begin() + expansion_size;
            } else {
                region_begin = mutation.begin() - expansion_size;
                region_length += 2*expansion_size;
            }

            region_list.emplace_back(GenomicPosition(mutation.chr_id, region_begin),
                                     region_length);
        }

        return region_list;
    }

    /**
     * @brief Filter out passenger CNAs that overlaps one region in a list
     *
     * @param region_filter is the list of regions that does not overlap the
     *      valid passenger CNAs
     */
    inline void filter_passenger_CNAs(std::list<GenomicRegion>&& region_filter)
    {
        filter_passenger_CNAs(region_filter);
    }

    /**
     * @brief Filter out passenger CNAs that overlaps one region in a list
     *
     * @param region_filter is the list of regions that does not overlap the
     *      valid passenger CNAs
     */
    void filter_passenger_CNAs(const std::list<GenomicRegion>& region_filter)
    {
        auto passenger_it = passenger_CNAs.begin();

        while (passenger_it != passenger_CNAs.end()) {
            const auto cna_region{passenger_it->get_region()};
            if (cna_region.overlaps(region_filter)) {
                std::swap(*passenger_it, *passenger_CNAs.rbegin());
                passenger_CNAs.pop_back();
            } else {
                ++passenger_it;
            }
        }
    }
public:
    bool infinite_sites_model;   //!< a flag to enable/disable infinite sites model

    bool avoid_homozygous_losses;       //!< a flag to avoiding homozygous losses

    /**
     * @brief The empty constructor
     */
    MutationEngine():
        driver_CNA_min_distance(0), warning(CLONES::warning),
        infinite_sites_model(true), avoid_homozygous_losses(true)
    {}

    /**
     * @brief A constructor
     *
     * @param context_index is the genome context index
     * @param repetition_index is the genome repetition index
     * @param SBS_signatures is the map of the SBS signatures
     * @param ID_signatures is the map of the indel signatures
     * @param germline_mutations are the germline mutations
     * @param driver_storage is the storage of driver mutations
     * @param passenger_CNAs is the vector of the admissible passenger CNAs
     * @param driver_CNA_min_distance is the minimum distance between a driver
     *      mutation and a passenger CNA
     * @param warning is the warning function
     */
    MutationEngine(ContextIndex<GENOME_WIDE_POSITION>& context_index,
                   RSIndex& repetition_index,
                   const std::map<std::string, SBSSignature>& SBS_signatures,
                   const std::map<std::string, IDSignature>& ID_signatures,
                   const GenomeMutations& germline_mutations,
                   const DriverStorage& driver_storage,
                   const std::vector<CNA>& passenger_CNAs={},
                   const unsigned int& driver_CNA_min_distance=10000,
                   WarningFunction warning=CLONES::warning):
        MutationEngine(context_index, repetition_index, SBS_signatures,
                       ID_signatures, MutationalProperties(),
                       germline_mutations, driver_storage, passenger_CNAs,
                       driver_CNA_min_distance, warning)
    {}

    /**
     * @brief A constructor
     *
     * @param context_index is the genome context index
     * @param repetition_index is the genome repetition index
     * @param SBS_signatures is the map of the SBS signatures
     * @param ID_signatures is the map of the indel signatures
     * @param mutational_properties are the mutational properties of all the species
     * @param germline_mutations are the germline mutations
     * @param driver_storage is the storage of driver mutations
     * @param passenger_CNAs is the vector of the admissible passenger CNAs
     * @param driver_CNA_min_distance is the minimum distance between a driver
     *      mutation and a passenger CNA
     * @param warning is the warning function
     */
    MutationEngine(ContextIndex<GENOME_WIDE_POSITION>& context_index,
                   RSIndex& repetition_index,
                   const std::map<std::string, SBSSignature>& SBS_signatures,
                   const std::map<std::string, IDSignature>& ID_signatures,
                   const MutationalProperties& mutational_properties,
                   const GenomeMutations& germline_mutations,
                   const DriverStorage& driver_storage,
                   const std::vector<CNA>& passenger_CNAs={},
                   const unsigned int& driver_CNA_min_distance=10000,
                   WarningFunction warning=CLONES::warning):
        generator(), context_index(context_index), rs_index(repetition_index),
        mutational_properties(mutational_properties),
        germline_mutations(std::make_shared<GenomeMutations>(germline_mutations)),
        dm_genome(germline_mutations.copy_structure()),
        driver_storage(driver_storage),
        passenger_CNAs{filter_CNA_by_chromosome_ids(passenger_CNAs)},
        driver_CNA_min_distance{driver_CNA_min_distance},
        warning{warning}, infinite_sites_model{true}, avoid_homozygous_losses{true}
    {
        MutationEngine::check_genomes_consistency(context_index, germline_mutations);

        for (const auto& [name, SBS_signature]: SBS_signatures) {
            inv_cumulative_SBSs[name] = get_inv_cumulative_signature(SBS_signature);
        }

        for (const auto& [name, ID_signatures]: ID_signatures) {
            inv_cumulative_IDs[name] = get_inv_cumulative_signature(ID_signatures);
        }
    }

    /**
     * @brief Add the properties of a mutant
     *
     * This method add the properties of a mutant (i.e., its passenger rates,
     * its driver mutations) and all its species to the mutations engine.
     *
     * @param name is the name of the mutant
     * @param epistate_passenger_rates is a map from epigenomic status to
     *          passenger rate
     * @param driver_SIDs is a list of driver SIDs
     * @param driver_CNAs is a list of driver CNAs
     * @param wg_doubling is a Boolean flag to enable whole genome doubling
     * @return a reference to the updated object
     */
    MutationEngine& add_mutant(const std::string& name,
                               const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                               const std::list<MutationSpec<SID>>& driver_SIDs={},
                               const std::list<CNA>& driver_CNAs={},
                               const bool& wg_doubling=false)
    {
        auto application_order = DriverMutations::get_default_order(driver_SIDs, driver_CNAs,
                                                                    wg_doubling);

        return add_mutant(name, epistate_passenger_rates, driver_SIDs,
                          driver_CNAs, application_order);
    }

    /**
     * @brief Add the properties of a mutant
     *
     * This method add the properties of a mutant (i.e., its passenger rates,
     * its driver mutations) and all its species to the mutations engine.
     *
     * @param name is the name of the mutant
     * @param epistate_passenger_rates is a map from epigenomic status to
     *          passenger rate
     * @param driver_SIDs is a list of driver SIDs
     * @param driver_CNAs is a list of driver CNAs
     * @param application_order is the list of application order
     * @return a reference to the updated object
     */
    MutationEngine& add_mutant(const std::string& name,
                               const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                               const std::list<MutationSpec<SID>>& driver_SIDs,
                               const std::list<CNA>& driver_CNAs,
                               const std::list<DriverMutations::MutationType>& application_order)
    {
        mutational_properties.add_mutant(name, epistate_passenger_rates, driver_SIDs,
                                         driver_CNAs, application_order);

        filter_passenger_CNAs(get_region_list(driver_SIDs, driver_CNA_min_distance));
        filter_passenger_CNAs(get_region_list(driver_CNAs, driver_CNA_min_distance));

        for (auto sid: driver_SIDs) {
            sid.allele_id = RANDOM_ALLELE;
            sid.nature = Mutation::DRIVER;
            place_SID(sid, dm_genome);
        }

        return *this;
    }

    /**
     * @brief Validate the names of an exposure
     *
     * This method validates the names of an exposure by searching them among the keys
     * of a map.
     *
     * @tparam VALUE_TYPE is the map value type
     * @param signature_map is a map whose keys are the signature names
     * @param exposure is an exposure
     */
    template<typename VALUE_TYPE>
    static void validate_signature_names(const std::map<std::string, VALUE_TYPE>& signature_map,
                                         const MutationalExposure& exposure)
    {
        for (const auto& [name, coeff]: exposure) {
            if (signature_map.count(name)==0) {
                throw std::runtime_error("Unknown signature " + name + ".");
            }
        }
    }

    /**
     * @brief Add a timed exposure
     *
     * This method add a exposure that will be applied from the specified simulation
     * time on.
     *
     * @param time is the simulation time from which the exposure will be applied
     * @param exposure is a SNV-indel exposure
     * @return a reference to the updated object
     */
    MutationEngine& add(const Time& time, const MutationalExposure& exposure)
    {
        if (time<0) {
            throw std::domain_error("Simulation time is a non-negative value");
        }

        for (const auto& [m_type, mt_exposure] : split_exposures(exposure)) {
            auto& ttimed_exposures = get_timed_exposures(m_type);
            if (ttimed_exposures.count(time)>0) {
                throw std::runtime_error("Another exposure has been set for time "
                                         + std::to_string(time) + ".");
            }
            ttimed_exposures.emplace(Time(time), mt_exposure);
        }

        return *this;
    }

    /**
     * @brief Add a timed exposure
     *
     * This method add a exposure that will be applied from the specified simulation
     * time on.
     *
     * @param time is the simulation time from which the exposure will be applied
     * @param exposure is an SNV-indel exposure
     * @return a reference to the updated object
     */
    MutationEngine& add(const Time& time, MutationalExposure&& exposure)
    {
        return add(time, exposure);
    }

    /**
     * @brief Add the default exposure
     *
     * This method add a default exposure.
     *
     * @param default_exposure is the exposure at time 0
     * @return a reference to the updated object
     */
    MutationEngine& add(const MutationalExposure& default_exposure)
    {
        return add(0, default_exposure);
    }

    /**
     * @brief Place genomic mutations on a descendants forest
     *
     * @param descendants_forest is a descendants forest
     * @param num_of_preneoplatic_SNVs is the number of preneoplastic SNVs
     * @param num_of_preneoplatic_indels is the number of preneoplastic indels
     * @param seed is the random generator seed
     * @param preneoplatic_SNV_signature_name is the pre-neoplastic SNV signature name
     * @param preneoplatic_indel_signature_name is the pre-neoplastic indel signature name
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    inline
    PhylogeneticForest
    place_mutations(const Mutants::DescendantsForest& descendants_forest,
                    const size_t& num_of_preneoplatic_SNVs,
                    const size_t& num_of_preneoplatic_indels,
                    const int& seed=0,
                    const std::string& preneoplatic_SNV_signature_name="SBS1",
                    const std::string& preneoplatic_indel_signature_name="ID1")
    {
        return place_mutations(descendants_forest, num_of_preneoplatic_SNVs,
                               num_of_preneoplatic_indels, nullptr, seed,
                               preneoplatic_SNV_signature_name,
                               preneoplatic_indel_signature_name);
    }

    /**
     * @brief Place genomic mutations on a descendants forest
     *
     * @param descendants_forest is a descendants forest
     * @param num_of_preneoplatic_SNVs is the number of preneoplastic SNVs
     * @param num_of_preneoplatic_indels is the number of preneoplastic indels
     * @param progress_bar is a progress bar pointer
     * @param seed is the random generator seed
     * @param preneoplatic_SNV_signature_name is the pre-neoplastic SNV signature name
     * @param preneoplatic_indel_signature_name is the pre-neoplastic indel signature name
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    PhylogeneticForest
    place_mutations(const Mutants::DescendantsForest& descendants_forest,
                    const size_t& num_of_preneoplatic_SNVs,
                    const size_t& num_of_preneoplatic_indels,
                    UI::ProgressBar *progress_bar, const int& seed=0,
                    const std::string& preneoplatic_SNV_signature_name="SBS1",
                    const std::string& preneoplatic_indel_signature_name="ID1")
    {
        using namespace CLONES::Mutants;
        using namespace CLONES::Mutants::Evolutions;
        using namespace CLONES::Mutations;

        for (const auto& [name, type]:
                std::map<std::string, MutationType::Type>({{"SBS", MutationType::Type::SBS},
                                                           {"indel", MutationType::Type::INDEL}})) {
            if (may_place_passengers(type) && !has_default_exposure(type)) {
                throw std::runtime_error("The default exposure for " + name
                                         + "s has not been set yet.");
            }
        }

        generator.seed(seed);

        PhylogeneticForest forest;

        forest.germline_mutations = germline_mutations;
        forest.mutational_properties = mutational_properties;

        static_cast<DescendantsForest&>(forest) = descendants_forest;

        auto species_rates = get_species_rate_map(descendants_forest);
        auto driver_mutations = get_driver_mutation_map(descendants_forest);

        auto chr_regions = context_index.get_chromosome_regions();

        auto wild_type_structure = germline_mutations->copy_structure();

        size_t visited_node = 0;
        for (auto& root: forest.get_roots()) {
            GenomeMutations mutations = wild_type_structure;

            // place preneoplastic mutations
            place_SIDs<SBSType>(&root, mutations, preneoplatic_SNV_signature_name,
                                num_of_preneoplatic_SNVs, Mutation::PRENEOPLASTIC);
            place_SIDs<IDType>(&root, mutations, preneoplatic_indel_signature_name,
                               num_of_preneoplatic_indels, Mutation::PRENEOPLASTIC);

            place_mutations(root, mutations, species_rates, driver_mutations,
                            visited_node, progress_bar);
        }

        // save mutant driver back in mutational properties
        for (const auto& [mutant_id, mutant_drivers] : driver_mutations) {
            forest.mutational_properties.set_mutant_drivers(mutant_drivers);
        }

        // reverse index extractions
        reset_index<SBSType>();
        reset_index<IDType>();

        return forest;
    }

    /**
     * @brief Place genomic mutations on a descendants forest
     *
     * @param descendants_forest is a descendants forest
     * @param num_of_preneoplatic_SNVs is the number of preneoplastic SNVs
     * @param num_of_preneoplatic_indels is the number of preneoplastic indels
     * @param progress_bar is a progress bar pointer
     * @param seed is the random generator seed
     * @param preneoplatic_SNV_signature_name is the pre-neoplastic SNV signature name
     * @param preneoplatic_indel_signature_name is the pre-neoplastic indel signature name
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    PhylogeneticForest
    place_mutations(const Mutants::DescendantsForest& descendants_forest,
                    const size_t& num_of_preneoplatic_SNVs,
                    const size_t& num_of_preneoplatic_indels,
                    UI::ProgressBar& progress_bar, const int& seed=0,
                    const std::string& preneoplatic_SNV_signature_name="SBS1",
                    const std::string& preneoplatic_indel_signature_name="ID1")
    {
        return place_mutations(descendants_forest, num_of_preneoplatic_SNVs,
                               num_of_preneoplatic_indels, &progress_bar, seed,
                               preneoplatic_SNV_signature_name, preneoplatic_indel_signature_name);
    }

    /**
     * @brief Get the a sample of wild-type cells
     *
     * @param pn_mutations is a list of pre-neoplastic mutations
     * @return a sample of a single wild-type cell
     */
    inline SampleGenomeMutations
    get_wild_type_sample(const std::list<MutationSpec<SID>>& pn_mutations)
    {
         std::binomial_distribution<size_t> bin_dist(0, 0);

        return get_wild_type_sample(1, bin_dist, bin_dist, "SBS1", "ID1", pn_mutations);
    }

    /**
     * @brief Get the a sample of wild-type cells
     *
     * @param forest is a phylogenetic forest
     * @param num_of_cells is the number of cells in the wild-type samples
     * @param SNV_dist is the distribution of the number of cell specific SNVs
     * @param indel_dist is the distribution of the number of cell specific indels
     * @param SNV_signature_name is the SNV signature name
     * @param indel_signature_name is the indel signature name
     * @param forest_preneoplastic_prob is the probability for a pre-neoplastic
     *   mutation in the forest to be selected for a wild-type sample
     * @return a list of wild-type samples; each of them has size `num_of_cells`.
     */
    std::list<SampleGenomeMutations>
    get_wild_type_samples(const PhylogeneticForest& forest, const size_t& num_of_cells,
                          std::normal_distribution<double>& SNV_dist,
                          std::normal_distribution<double>& indel_dist,
                          std::normal_distribution<double>& common_SNV_dist,
                          std::normal_distribution<double>& common_indel_dist,
                          const std::string& SNV_signature_name="SBS1",
                          const std::string& indel_signature_name="ID1",
                          const double& forest_preneoplastic_prob=1)
    {
        using namespace Mutants::Evolutions;

        std::list<SampleGenomeMutations> samples;

        std::uniform_real_distribution<double> pnp_dist;
        for (const auto& [cell_id, pnp_mutations] : forest.get_preneoplastic_mutations()) {
            auto sample_common_mutations = germline_mutations->copy_structure();
            if (pnp_dist(generator) < forest_preneoplastic_prob) {
                for (auto mut_it = pnp_mutations.begin(); mut_it != pnp_mutations.end();
                    ++mut_it) {

                    switch(mut_it.get_type()) {
                        case MutationList::SID_TURN:
                        {
                            const auto& mutation = mut_it.get_last_SID();

                            sample_common_mutations.apply(mutation);
                            break;
                        }
                        case MutationList::CNA_TURN:
                        {
                            const auto& mutation = mut_it.get_last_CNA();

                            sample_common_mutations.apply(mutation);
                            break;
                        }
                        case MutationList::WGD_TURN:
                        {
                            sample_common_mutations.duplicate_alleles();
                        }
                        default:
                            break;
                    }
                }
            }

            place_SIDs_in_normal<SBSType>(sample_common_mutations, SNV_signature_name,
                                          common_SNV_dist);

            place_SIDs_in_normal<IDType>(sample_common_mutations, indel_signature_name,
                                         common_indel_dist);

            SampleGenomeMutations sample_mutations("wild-type_" + std::to_string(cell_id),
                                                   germline_mutations);
            for (size_t i=0; i<num_of_cells; ++i) {
                auto mutations = std::make_shared<CellGenomeMutations>(sample_common_mutations);

                place_SIDs_in_normal<SBSType>(sample_common_mutations, SNV_signature_name,
                                              SNV_dist);
                place_SIDs_in_normal<IDType>(sample_common_mutations, indel_signature_name,
                                             indel_dist);

                sample_mutations.mutations.push_back(mutations);
            }

            samples.push_back(std::move(sample_mutations));
        }

        return samples;
    }

    /**
     * @brief Check whether a default mutation type exposure has been set
     *
     * @param mutation_type is the type of mutation affected by the exposure
     * @return `true` if and only if a default mutation type exposure has been set
     */
    inline bool has_default_exposure(const MutationType::Type& mutation_type) const
    {
        const auto& timed_type_exposures = get_timed_exposures(mutation_type);

        return timed_type_exposures.find(0.0) != timed_type_exposures.end();
    }

    /**
     * @brief Check whether the default exposures have been set
     *
     * @return `true` if and only if the default exposures have been set
     */
    inline bool has_default_exposures() const
    {
        return has_default_exposure(MutationType::Type::SBS)
            && has_default_exposure(MutationType::Type::INDEL);
    }

    /**
     * @brief Get the mutational properties
     *
     * @return a constant reference to the mutational properties
     */
    inline const MutationalProperties& get_mutational_properties() const
    {
        return mutational_properties;
    }

    /**
     * @brief Get the exposures
     *
     * @return a constant reference to the exposures
     */
    inline const std::map<Time, MutationalExposure>&
    get_timed_exposures(const MutationType::Type& mutation_type) const
    {
        return timed_exposures[static_cast<size_t>(mutation_type)];
    }

    /**
     * @brief Get the driver storage
     *
     * @return a constant reference to the driver storage
     */
    inline const DriverStorage& get_driver_storage() const
    {
        return driver_storage;
    }

    /**
     * @brief Get the germline mutations
     *
     * @return the germline mutations
     */
    inline const GenomeMutations& get_germline_mutations() const
    {
        return *germline_mutations;
    }

    /**
     * @brief Reset the index to its original state
     *
     * @tparam MUTATION_TYPE is the mutation type of the
     *      of the index to be reset
     */
    template<typename MUTATION_TYPE,
             std::enable_if_t<std::is_base_of_v<MutationType, MUTATION_TYPE>, bool> = true>
    void reset_index()
    {
        while (get_stack_size<MUTATION_TYPE>()>0) {
            restore_last_extracted_from_index<MUTATION_TYPE>();
        }
    }
};

}   // Mutations

}   // CLONES

#endif // __CLONES_MUTATION_ENGINE__
