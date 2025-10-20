/**
 * @file read_simulator.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes to simulate sequencing
 * @version 1.27
 * @date 2025-10-20
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

#ifndef __RACES_READ_SIMULATOR__
#define __RACES_READ_SIMULATOR__

#include <list>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <cctype>
#include <random>
#include <ranges>
#include <memory>
#include <iomanip> // std::setw and std::setfill

#include "sid.hpp"
#include "genome_mutations.hpp"
#include "phylogenetic_forest.hpp"
#include "tissue_sample.hpp"
#include "read.hpp"

#include "fasta_chr_reader.hpp"

#include "progress_bar.hpp"
#include "variables.hpp"

#include "utils.hpp"

#include "sequencer.hpp"

#if WITH_MATPLOT
#include <matplot/matplot.h>
#endif // WITH_MATPLOT

namespace RACES
{

namespace Mutations
{

/**
 * @brief The sequencing simulation namespace
 */

namespace SequencingSimulations
{

using BaseCoverage = uint16_t;

/**
 * @brief Sequencing targets
 *
 * The objects of this class store for each sample the number
 * and the kind of the wild-cell contaminants. Moreover, their
 * maintains a copy of the wild-cell genome mutations.
 */
class SequencingTargets
{
public:

    /**
     * @brief The sample target
     *
     * The object of this class store the number and the kind
     * of the wild-cell contaminants of a sample.
     */
    class SampleTarget
    {
    public:
        std::string sample_name;    //!< the name of the sample
        Mutants::Evolutions::TissueSampleId sample_id;  //!< the identifier of the sample

        size_t num_of_tumour_cells; //!< the number of tumour cells in the sample
        std::map<Mutants::CellId, size_t> num_of_wild_types;    //!< the number of wild-type cells per forest root

        /**
         * @brief The empty constructor
         */
        SampleTarget();

        /**
         * @brief A constructor
         *
         * This constructor takes the number of sample cells and computes the number
         * of wild-type cells needed to achieve a specified purity. The wild-type cells
         * may or may be not added with the pre-neoplastic mutations. In the later case,
         * the number of wild-type cells having the same pre-neoplastic mutations of
         * each forest root is randomly chosen by using a multinomial distribution
         * whose probabilities respect the ratio of the sample cells in each tree.
         *
         * @tparam RANDOM_GENERATOR is a random number generator type
         * @param[in,out] generator is a random number generator
         * @param[in] sample is the sample for which the target is built
         * @param[in] sample_cells_per_tree is a map associating to each forest root
         *      the number of sample in the corresponding subtree
         * @param[in] purity is the aimed sequencing purity
         * @param[in] with_pre_neoplastic_mutations is a Boolean flag enabling/disabling
         *      pre-neoplastic mutations from target
         */
        template<typename RANDOM_GENERATOR>
        SampleTarget(RANDOM_GENERATOR& generator,
                     const Mutants::Evolutions::TissueSample& sample,
                     const std::map<Mutants::CellId, size_t>& sample_cells_per_tree,
                     const double& purity, const bool with_pre_neoplastic_mutations=true):
            sample_name{sample.get_name()},
            sample_id{sample.get_id()}
        {
            size_t num_of_contaminants;
            if (purity>0) {
                num_of_tumour_cells = sample.get_cell_ids().size();
                num_of_contaminants = static_cast<size_t>(num_of_tumour_cells*(1-purity)/purity);
            } else {
                num_of_tumour_cells = 0;

                if (with_pre_neoplastic_mutations) {
                    num_of_contaminants = sample.get_cell_ids().size();
                } else {
                    num_of_contaminants = 1;
                }
            }

            if (!with_pre_neoplastic_mutations) {
                num_of_wild_types.emplace(EMBRYO_CELL_ID, num_of_contaminants);

                return;
            }

            auto missing_cells = sample.get_cell_ids().size();
            for (const auto& [root_id, num_of_cells]: sample_cells_per_tree) {
                const double prob = static_cast<double>(num_of_cells)/missing_cells;
                missing_cells -= num_of_cells;

                std::binomial_distribution<size_t> dist(num_of_contaminants, prob);
                const size_t contaminant_from_root = dist(generator);

                if (contaminant_from_root>0) {
                    num_of_contaminants -= contaminant_from_root;
                    num_of_wild_types.emplace(root_id, contaminant_from_root);
                }
            }
        }

        /**
         * @brief Build the target for a normal sample
         *
         * @return the target for a normal sample
         */
        static SampleTarget build_normal_sample_target();

        /**
         * @brief Get the number of cells in the sample target
         *
         * @return the number of cells in the sample target
         */
        size_t num_of_cells() const;
    };
private:
    PhylogeneticForest const* phylo_forest;  //!< A pointer to the phylogenetic forest

    double seq_purity;    //!< The aimed sequencing purity
    bool pre_neoplastic;  //!< A Boolean flag for pre-neoplastic mutation sequencing
    bool germinal;        //!< A Boolean flag for germinal mutation sequencing

    std::map<Mutants::CellId, CellGenomeMutations> wt_genomes;  //!< The tree wild-type cell somatic mutations
    std::map<Mutants::Evolutions::TissueSampleId, SampleTarget> target_map;  //!< The sample id-target map

    /**
     * @brief Build the map sample identifier-target
     *
     * @tparam RANDOM_GENERATOR is a random number generator type
     * @param[in,out] generator is a random number generator
     * @param[in] forest is a phylogenetic forest
     * @param[in] produce_normal_sample is a Boolean flag to add/avoid a normal sample
     * @param[in] purity is the aimed sequencing purity
     * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
     *      mutations
     * @return the map sample identifier-target
     */
    template<typename RANDOM_GENERATOR>
    static std::map<Mutants::Evolutions::TissueSampleId, SampleTarget>
    build_target_map(RANDOM_GENERATOR& generator, const PhylogeneticForest& forest,
                     const bool& produce_normal_sample, const double& purity,
                     const bool& with_pre_neoplastic)
    {
        std::map<Mutants::Evolutions::TissueSampleId, SampleTarget> sample_targets;
        const auto cells_per_root = forest.count_cells_per_root();
        for (const auto& sample : forest.get_samples()) {
            sample_targets.emplace(sample.get_id(),
                                    SampleTarget{generator, sample,
                                                 cells_per_root.at(sample.get_id()),
                                                 purity, with_pre_neoplastic});
        }

        if (produce_normal_sample) {
            sample_targets.emplace(NORMAL_SAMPLE_ID,
                                   SampleTarget::build_normal_sample_target());
        }

        return sample_targets;
    }

    /**
     * @brief A constructor
     *
     * This method is not meant to be public because it assumes that wild-type genomes and
     * sample targets are consistent to the Boolean flag and the purity.
     *
     * @param[in] forest is the phylogenetic forest of the samples to be sequenced
     * @param[in] purity is the aimed purity
     * @param[in] with_pre_neoplastic is a Boolean flag to enable/avoid pre-neoplastic
     *      mutations
     * @param[in] with_germinal is a Boolean flag to enable/avoid pre-neoplastic
     *      mutations
     * @param[in,out] wild_type_genomes is a map associating each forest tree to the corresponding
     *      wild-type genome mutations
     * @param[in,out] sample_targets is a map associating each sample id to the corresponding
     *      sample target
     */
    SequencingTargets(const PhylogeneticForest& forest, const double& purity,
                      const bool& with_pre_neoplastic,
                      const bool& with_germinal,
                      std::map<Mutants::CellId, CellGenomeMutations>&& wild_type_genomes,
                      std::map<Mutants::Evolutions::TissueSampleId, SampleTarget>&& sample_targets);
public:

    /**
     * @brief Get the sequencing targets
     *
     * @tparam RANDOM_GENERATOR is a random number generator type
     * @param[in,out] generator is a random number generator
     * @param[in] forest is a phylogenetic forest
     * @param[in] produce_normal_sample is a Boolean flag to add/avoid a normal sample
     * @param[in] purity is the aimed sequencing purity
     * @param[in] with_pre_neoplastic is a Boolean flag to add/avoid pre-neoplastic
     *      mutations
     * @param[in] with_germinal is a Boolean flag to add/avoid germinal mutations
     * @return the sequencing target of the samples in the phylogenetic forest
     */
    template<typename RANDOM_GENERATOR>
    inline static SequencingTargets
    get_targets(RANDOM_GENERATOR& generator, const PhylogeneticForest& forest,
                const bool& produce_normal_sample, const double& purity,
                const bool with_pre_neoplastic_mutations=true,
                const bool with_germinal=true)
    {
        auto target_map = build_target_map(generator, forest, produce_normal_sample,
                                               purity, with_pre_neoplastic_mutations);

        return {forest, purity, with_pre_neoplastic_mutations, with_germinal,
                forest.get_wild_type_genomes(with_pre_neoplastic_mutations, false),
                std::move(target_map)};
    }

    /**
     * @brief The phylogenetic forest of the samples
     *
     * @return A constant reference to the phylogenetic forest of the samples
     */
    inline const PhylogeneticForest& forest() const
    {
        return *phylo_forest;
    }

    /**
     * @brief The aimed sequencing purity
     *
     * @return a constant reference to the aimed sequencing purity
     */
    inline const double& purity() const
    {
        return seq_purity;
    }

    /**
     * @brief Check whether pre-neoplastic mutations must be sequenced
     *
     * @return `true` if and only if pre-neoplastic mutations must be
     *      sequenced
     */
    inline const bool& with_pre_neoplastic() const
    {
        return pre_neoplastic;
    }
    /**
     * @brief Check whether germinal mutations must be sequenced
     *
     * @return `true` if and only if germinal mutations must be sequenced
     */
    inline const bool& with_germinal() const
    {
        return germinal;
    }

    /**
     * @brief Get the wild type genomes
     *
     * @return a constant reference to a map associating each tree to the
     *      corresponding wild-type mutations
     */
    inline const std::map<Mutants::CellId, CellGenomeMutations>& wild_type_genomes() const
    {
        return wt_genomes;
    }

    /**
     * @brief Get sample target map
     *
     * @return a constant reference to a map associating each sample id to the
     *      corresponding sample target
     */
    inline const std::map<Mutants::Evolutions::TissueSampleId, SampleTarget>&
    sample_targets() const
    {
        return target_map;
    }
};

/**
 * @brief Simulated sequencing chromosome statistics
 */
class ChrCoverage
{
    ChromosomeId chr_id;                //!< The chromosome id
    std::vector<BaseCoverage> coverage; //!< The chromosome base coverage
public:

    /**
     * @brief The empty constructor
     */
    ChrCoverage();

    /**
     * @brief A constructor
     *
     * @param[in] chromosome_id is the identifier of the chromosome whose coverage refer to
     * @param[in] size is the size of the chromosome whose coverage refer to
     */
    ChrCoverage(const ChromosomeId& chromosome_id, const GenomicRegion::Length& size);

    /**
     * @brief Increase the coverage of a region
     *
     * @param[in] begin_pos is the initial position of the region whose coverage must be increase
     * @param[in] read_size is the size of the region whose coverage must be increase
     */
    void increase_coverage(const ChrPosition& begin_pos, const size_t& read_size);

    /**
     * @brief Increase the coverage of a region
     *
     * @param[in] region is the region whose coverage must be increase
     */
    inline void increase_coverage(const GenomicRegion& region)
    {
        increase_coverage(region.begin(),
                          static_cast<const GenomicRegion::Length&>(region.size()));
    }

    /**
     * @brief Increase the coverage of a region
     *
     * @param[in,out] region is the region whose coverage must be increase
     */
    inline void increase_coverage(GenomicRegion&& region)
    {
        increase_coverage(region.begin(),
                          static_cast<const GenomicRegion::Length&>(region.size()));
    }

    /**
     * @brief Get the identifier of the chromosome whose coverage refer to
     *
     * @return a constant reference to the identifier of the chromosome whose
     *          coverage refer to
     */
    inline const ChromosomeId& get_chr_id() const
    {
        return chr_id;
    }

    /**
     * @brief Get the length of the chromosome whose coverage refer to
     *
     * @return the size of the chromosome whose coverage refer to
     */
    inline size_t size() const
    {
        return coverage.size();
    }

    /**
     * @brief Check whether the chromosome coverage is initialized
     *
     * @return `true` if and only if this object was created by using
     *      a non-empty constructor or was copied by an initialized
     *      object
     */
    inline bool is_initialized() const
    {
        return coverage.size()!=0;
    }

    /**
     * @brief Get the coverage of a genomic position
     *
     * @param[in] position is the genomic position whose coverage is aimed
     * @return a constant reference to the coverage of the aimed position
     */
    const BaseCoverage& get_coverage(const GenomicPosition& position) const;

    /**
     * @brief Get the chromosome coverage
     *
     * @return a constant reference to the chromosome coverage vector
     */
    inline const std::vector<BaseCoverage>& get_coverage_vector() const
    {
        return coverage;
    }

    /**
     * @brief Join a chromosome coverage to the current one
     *
     * This method updates the coverage of the current object by
     * adding the coverage of the parameter. Moreover, it adds to
     * the current object mutation occurrences the occurrences of the
     * parameter.
     *
     * @param[in,out] chr_coverage is the sequencing chromosome coverage to join
     * @return a reference to the updated coverage
     */
    ChrCoverage& operator+=(ChrCoverage&& chr_coverage);

    /**
     * @brief Join a chromosome coverage to the current one
     *
     * This method updates the coverage of the current object by
     * adding the coverage of the parameter. Moreover, it adds to
     * the current object mutation occurrences the occurrences of the
     * parameter.
     *
     * @param[in] chr_coverage is the sequencing chromosome coverage to join
     * @return a reference to the updated coverage
     */
    ChrCoverage& operator+=(const ChrCoverage& chr_coverage);

    /**
     * @brief Save sequencing chromosome coverage in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param[in,out] archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & chr_id
                & coverage;
    }

    /**
     * @brief Load sequencing chromosome coverage from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param[in,out] archive is the input archive
     * @return the sequencing chromosome coverage
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static ChrCoverage load(ARCHIVE& archive)
    {
        ChrCoverage chr_coverage;

        archive & chr_coverage.chr_id
                & chr_coverage.coverage;

        return chr_coverage;
    }
};

/**
 * @brief Sequencing data of any SID mutation
 */
struct SIDData
{
    BaseCoverage num_of_occurrences;        //!< The number of occurrences
    std::set<std::string> causes;           //!< The causes of the SID
    std::set<Mutation::Nature> nature_set;  //!< The nature of the SID

    /**
     * @brief The empty constructor
     */
    SIDData();

    /**
     * @brief Construct a new SIDData object
     *
     * @param[in] mutation is the mutation whose data refer to
     * @param[in] num_of_occurrences is the number of occurrences of `mutation`
     */
    SIDData(const SID& mutation, const BaseCoverage num_of_occurrences);

    /**
     * @brief Update the SID data
     *
     * This method increases the number of occurrences of the
     * SID and adds the cause and type of a SID to the
     * causes and types stored in the SID data.
     *
     * @param[in] mutation is the SID whose cause and type should be
     *          added to the causes and types stored in the
     *          data
     */
    void update(const SID& mutation);

    /**
     * @brief Update the current object
     *
     * This method adds all the data contained in an SIDData
     * object to the current object. It sums the number of
     * occurrences of the parameter to those of the current
     * object. Moreover, it also adds the causes and the
     * types of the parameter to the current object causes
     * and types.
     *
     * @param[in] data is the SIDData object whose data is added
     *          to the current object
     */
    void update(const SIDData& data);
};

/**
 * @brief Simulated sequencing chromosome statistics
 */
class ChrSampleStatistics : public ChrCoverage
{
    std::map<SID, SIDData> SID_data;    //!< The SID data about the simulated reads

    bool with_germline;     //!< A Boolean flag to account/avoid germinal mutations

    /**
     * @brief Add the mutation of a chromosome to the statistics
     *
     * @param[in] chromosome is the chromosome whose mutations must be
     *     added to the statistics
     */
    void account_for(const ChromosomeMutations& chromosome);
public:

    /**
     * @brief The empty constructor
     */
    ChrSampleStatistics();

    /**
     * @brief A copy constructor
     *
     * @param[in] original is the object to be copied
     */
    ChrSampleStatistics(const ChrSampleStatistics& original);

    /**
     * @brief A constructor
     *
     * @param[in] chromosome_id is the identifier of the chromosome whose sequencing
     *          statistics are collected
     * @param[in] size is the size of the chromosome whose sequencing
     *          statistics are collected
     * @param[in] forest is the phylogenetic forest of the samples to be sequenced
     * @param[in] with_missed_SID is a Boolean flag to account/avoid SIDs that are not
     *      sequenced
     * @param[in] with_germline is a Boolean flag to account/avoid germinal mutations
     */
    ChrSampleStatistics(const ChromosomeId& chromosome_id, const GenomicRegion::Length& size,
                        const PhylogeneticForest& forest, const bool& with_missed_SID=false,
                        const bool& with_germline=false);

    /**
     * @brief Update the data associated to an SID
     *
     * @param[in] mutation is the SID whose data must be updated
     */
    void update(const SID& mutation);

    /**
     * @brief Collect the data from a read
     *
     * @param[in] read is a read
     */
    void add(const Read& read);

    /**
     * @brief Get the collected SID data
     *
     * @return a constant reference to the collected SID data
     */
    inline const std::map<SID, SIDData>& get_data() const
    {
        return SID_data;
    }

    /**
     * @brief Get the collected SID data
     *
     * @return a constant reference to the collected SID data
     */
    inline std::map<SID, SIDData>& get_data()
    {
        return SID_data;
    }

    /**
     * @brief Get the collected data of an SID
     *
     * @param[in] mutation is an SID
     * @return the data of `mutation`
     */
    SIDData get_data(const SID& mutation) const;

    /**
     * @brief Join other chromosome statistics to the current one
     *
     * This method updates the coverage of the current object by
     * adding the coverage of the parameter. Moreover, it adds to
     * the current object SID occurrences the occurrences of the
     * parameter.
     *
     * @param[in,out] chr_stats is the sequencing chromosome statistics to join
     * @return a reference to the updated objects
     */
    ChrSampleStatistics& operator+=(ChrSampleStatistics&& chr_stats);

    /**
     * @brief Join other chromosome statistics to the current one
     *
     * This method updates the coverage of the current object by
     * adding the coverage of the parameter. Moreover, it adds to
     * the current object SID occurrences the occurrences of the
     * parameter.
     *
     * @param[in] chr_stats is the sequencing chromosome statistics to join
     * @return a reference to the updated objects
     */
    ChrSampleStatistics& operator+=(const ChrSampleStatistics& chr_stats);
};

/**
 * @brief Statistics about simulated sequencing on a set of tissue samples
 */
class SampleSetStatistics;

/**
 * @brief Statistics about simulated sequencing of a tissue sample
 */
class SampleStatistics
{
    std::filesystem::path data_dir;     //!< The directory in which the data is saved
    std::string sample_name;            //!< The name of the sample

    std::set<ChromosomeId> chr_ids;     //!< The identifiers of the chromosomes included in the statistics

    std::map<SID, SIDData> SID_data;    //!< The SID data in the sample
    std::map<GenomicPosition, BaseCoverage> locus_coverage;   //!< The coverage of any position hosting an SID

    bool save_coverage; //!< A flag to enable/disable storage of coverage data

    /**
     * @brief Get the filename of the chromosome coverage file
     *
     * @param[in] chr_id is the identifier of the chromosome whose coverage is
     *          going to be saved/load
     * @return the filename of the chromosome coverage file
     */
    std::filesystem::path get_coverage_filename(const ChromosomeId& chr_id) const;


    /**
     * @brief Add chromosome statistics coverage to the genome statistics
     *
     * @param[in] chr_stats is the chromosome statistics
     */
    void add_coverage(const ChrSampleStatistics& chr_stats);
public:

    /**
     * @brief A constructor
     *
     * @param[in] data_directory is the path of the directory containing the data files
     * @param[in] sample_name is the name of the sample whose statistics refer to
     * @param[in] save_coverage is a flag to enable/disable storage of coverage data
     */
    SampleStatistics(const std::filesystem::path& data_directory, const std::string& sample_name,
                     const bool save_coverage=false);

    /**
     * @brief Add chromosome statistics to the genome statistics
     *
     * @param[in,out] chr_stats is the chromosome statistics
     */
    inline void add_chr_statistics(ChrSampleStatistics&& chr_stats);

    /**
     * @brief Add chromosome statistics to the genome statistics
     *
     * @param[in] chr_stats is the chromosome statistics
     */
    void add_chr_statistics(const ChrSampleStatistics& chr_stats);

    /**
     * @brief Get the chromosome coverage by its identifier
     *
     * @param[in] chr_id is the id of the chromosome whose sequencing coverage is aimed
     * @return the sequencing coverage of the chromosome whose identifier is `chr_id`
     */
    ChrCoverage get_chr_coverage(const ChromosomeId& chr_id) const;

    /**
     * @brief Get the chromosome statistics
     *
     * @return a constant reference to the chromosome statistics
     */
    inline const std::set<ChromosomeId>& get_chromosome_ids() const
    {
        return chr_ids;
    }

    /**
     * @brief Get the sample name
     *
     * @return the name of the sample whose data refer to
     */
    inline const std::string& get_sample_name() const
    {
        return sample_name;
    }

    /**
     * @brief Get the SID data
     *
     * @return a constant reference to the SID data
     */
    inline const std::map<SID, SIDData>& get_data() const
    {
        return SID_data;
    }

    /**
     * @brief Get the collected data of an SID
     *
     * @param[in] mutation is an SID
     * @return the data of `mutation`
     */
    SIDData get_data(const SID& mutation) const;

    /**
     * @brief Get the coverage of the positions in which SIDs occur
     *
     * @return a constant reference to the coverage of the positions in which SIDs
     *      occur
     */
    inline const std::map<GenomicPosition, BaseCoverage>& get_coverage() const
    {
        return locus_coverage;
    }

    /**
     * @brief Get the coverage of a position in which an SID occurs
     *
     * @return a constant reference to the coverage of `pos` assuming that
     *      `pos` is a position in which an SID occurs
     */
    const BaseCoverage& get_coverage(const GenomicPosition& pos) const;

    /**
     * @brief Get the data directory path
     *
     * @return as constant reference to the path of the directory
     *          which contained the object saved data
     */
    inline const std::filesystem::path& get_data_directory() const
    {
        return data_dir;
    }

    /**
     * @brief Save sequencing sample statistics in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param[in,out] archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & data_dir
                & sample_name
                & save_coverage
                & chr_ids
                & SID_data
                & locus_coverage;
    }

    /**
     * @brief Load sequencing sample statistics from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param[in,out] archive is the input archive
     * @return the sequencing sample statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static SampleStatistics load(ARCHIVE& archive)
    {
        std::string data_dir, sample_name;
        bool save_coverage;

        archive & data_dir
                & sample_name
                & save_coverage;

        SampleStatistics sample_stats(data_dir, sample_name, save_coverage);

        archive & sample_stats.chr_ids
                & sample_stats.SID_data
                & sample_stats.locus_coverage;

        return sample_stats;
    }

    friend class SampleSetStatistics;
};

/**
 * @brief A read simulator
 *
 * @tparam RANDOM_GENERATOR is the random generator type used to sample genome
 */
template<typename RANDOM_GENERATOR=std::mt19937_64>
class ReadSimulator;

class SampleSetStatistics
{
    std::unique_ptr<std::map<std::string, SampleStatistics>> stats_map;  //!< The map of the sample statistics

    std::set<ChromosomeId> chr_ids; //!< The identifiers of the chromosomes whose statistics have been collected

    std::filesystem::path data_dir; //!< The directory in which the data is saved

    /**
     * @brief Stream the header of a VAF csv
     *
     * @param[in,out] os is the output stream
     * @param[in] separator is the column separator symbol
     * @return the updated stream
     */
    std::ofstream& stream_VAF_csv_header(std::ofstream& os, const char& separator) const;
public:
    /**
     * @brief A constructor
     *
     * @param[in] data_directory is the path of the directory containing the data files
     */
    SampleSetStatistics(const std::filesystem::path& data_directory);

    /**
     * @brief Check whether this object contains statistics about a sample
     *
     * @param[in] sample_name is the sample name
     * @return `true` if and only if this object contains statistics about
     *          a sample named `sample_name`
     */
    inline std::map<std::string, SampleStatistics>::const_iterator find(const std::string& sample_name) const
    {
        return stats_map->find(sample_name);
    }

    /**
     * @brief Get the first constant iterator of a map from the sample names to the corresponding statistics
     *
     * @return the first constant iterator of a map from the sample names to the corresponding statistics
     */
    inline std::map<std::string, SampleStatistics>::const_iterator begin() const
    {
        return stats_map->begin();
    }

    /**
     * @brief Get the last constant iterator of a map from the sample names to the corresponding statistics
     *
     * @return the last constant iterator of a map from the sample names to the corresponding statistics
     */
    inline std::map<std::string, SampleStatistics>::const_iterator end() const
    {
        return stats_map->end();
    }

    /**
     * @brief Get the identifiers of the represented chromosome
     *
     * @return a constant reference to the identifiers of the chromosomes whose statistics have been
     *          collected
     */
    std::set<ChromosomeId> get_represented_chromosomes() const;

    /**
     * @brief Add a new chromosome to the statistics of a sample
     *
     * @param[in] sample_name is the sample name
     * @param[in] chr_statistics are the statistics about a chromosome
     * @param[in] save_coverage is a flag to enable/disable storage of coverage data
     * @return a constant reference to the current statistics of the sample whose name is
     *          `sample_name`
     */
    const SampleStatistics& add_chr_statistics(const std::string& sample_name,
                                               ChrSampleStatistics&& chr_statistics,
                                               const bool save_coverage=false);

    /**
     * @brief Add a new chromosome to the statistics of a sample
     *
     * @param[in] sample_name is the sample name
     * @param[in] chr_statistics are the statistics about a chromosome
     * @param[in] save_coverage is a flag to enable/disable storage of coverage data
     * @return a constant reference to the current statistics of the sample whose name is
     *          `sample_name`
     */
    const SampleStatistics& add_chr_statistics(const std::string& sample_name,
                                               const ChrSampleStatistics& chr_statistics,
                                               const bool save_coverage=false);

    /**
     * @brief Add a new sample statistics
     *
     * @param[in] sample_statistics are the sample statistics to add
     * @return a constant reference to the added sample statistics
     */
    const SampleStatistics& add_statistics(const SampleStatistics& sample_statistics);

    /**
     * @brief Get a sample statistics
     *
     * @param[in] sample_name is the sample name
     * @return a constant reference to the sample statistics
     */
    const SampleStatistics& operator[](const std::string& sample_name) const;

    /**
     * @brief Get the data directory path
     *
     * @return as constant reference to the path of the directory
     *          which contained the object saved data
     */
    inline const std::filesystem::path& get_data_directory() const
    {
        return data_dir;
    }

    /**
     * @brief Get the list of the sample names
     *
     * @return the list of the sample names
     */
    std::list<std::string> get_sample_names() const;

    /**
     * @brief Get the number of samples in the set
     *
     * @return the number of samples in the set
     */
    inline size_t num_of_samples() const
    {
        return stats_map->size();
    }

    /**
     * @brief Check whether the statistics is canonical
     *
     * A statistics is canonical when all its samples stores the occurrences of
     * the very same SIDs.
     *
     * @return `true` if and only if current statistics is canonical
     */
    bool is_canonical() const;

    /**
     * @brief Canonize the statistics
     *
     * A statistics is canonical when all its samples stores the occurrences of
     * the very same SIDs.
     * This method alters the sample statistics of this object and canonizes it.
     */
    void canonize();

    /**
     * @brief Save CSV files reporting the SID VAFs
     *
     * This method saves one CSV for each chromosome in the reference
     * sequence. Each CSV reports the SID VAFs in the corresponding
     * chromosome.
     *
     * @param[in] base_name is the prefix of the name of the CSV files.
     *          The filename format is "<base_name><chromosome_name>.csv"
     * @param[in,out] progress_bar_stream is the output stream for the progress bar
     * @param[in] quiet is a Boolean flag to enable/disable a progress bar
     */
    void save_VAF_CSVs(const std::string& base_name,
                       std::ostream& progress_bar_stream=std::cout,
                       const bool& quiet = false) const;


    /**
     * @brief Save CSV files reporting the SID VAFs
     *
     * This method saves one CSV for each chromosome in the reference
     * sequence. Each CSV reports the SID VAFs in the corresponding
     * chromosome. The file names will have the format
     * "chr_<chromosome_name>.csv".
     *
     * @param[in,out] progress_bar_stream is the output stream for the progress bar
     * @param[in] quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_VAF_CSVs(std::ostream& progress_bar_stream=std::cout,
                              const bool& quiet = false) const
    {
        save_VAF_CSVs("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a CSV file reporting the SID VAFs in a chromosome
     *
     * @param[in] filename is the image filename
     * @param[in] chromosome_id is the identifier of the chromosome whose
     *          coverage must be represented
     */
    void save_VAF_CSV(const std::filesystem::path& filename,
                      const ChromosomeId& chromosome_id) const;

#if WITH_MATPLOT
    /**
     * @brief Save images representing the chromosome coverages
     *
     * This method saves one image for each chromosome in the reference
     * sequence. Each image depicts the coverage of the corresponding
     * chromosome.
     *
     * @param[in] base_name is the prefix of the name of the CSV files.
     *          The filename format is "<base_name><chromosome_name>_coverage.jpg"
     * @param[in,out] progress_bar_stream is the output stream for the progress bar
     * @param[in] quiet is a Boolean flag to enable/disable a progress bar
     */
    void save_coverage_images(const std::string& base_name,
                              std::ostream& progress_bar_stream=std::cout,
                              const bool& quiet = false) const;

    /**
     * @brief Save images representing the chromosome coverages
     *
     * This method saves one image for each chromosome in the reference
     * sequence. Each image depicts the coverage of the corresponding
     * chromosome. The file names will have the format
     * "chr_<chromosome_name>_coverage.jpg".
     *
     * @param[in,out] progress_bar_stream is the output stream for the progress bar
     * @param[in] quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_coverage_images(std::ostream& progress_bar_stream=std::cout,
                                     const bool& quiet = false) const
    {
        save_coverage_images("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a image representing the coverage of a chromosome
     *
     * @param[in] filename is the image filename
     * @param[in] chromosome_id is the identifier of the chromosome whose
     *          coverage must be represented
     */
    void save_coverage_image(const std::filesystem::path& filename,
                             const ChromosomeId& chromosome_id) const;

    /**
     * @brief Save images representing the histogram of SIDs
     *
     * This method saves one image for each chromosome in the reference
     * sequence. Each image depicts the SID histogram of the corresponding
     * chromosome.
     *
     * @param[in] base_name is the prefix of the name of the CSV files.
     *          The filename syntax is "<base_name><chromosome_name>_hist.jpg"
     * @param[in,out] progress_bar_stream is the output stream for the progress bar
     * @param[in] quiet is a Boolean flag to enable/disable a progress bar
     */
    void save_SID_histograms(const std::string& base_name,
                             std::ostream& progress_bar_stream=std::cout,
                             const bool& quiet = false) const;


    /**
     * @brief Save images representing the histogram of SIDs
     *
     * This method saves one image for each chromosome in the reference
     * sequence. Each image depicts the SID histogram of the corresponding
     * chromosome. The file names will have the format
     * "chr_<chromosome_name>_hist.jpg".
     *
     * @param[in,out] progress_bar_stream is the output stream for the progress bar
     * @param[in] quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_SID_histograms(std::ostream& progress_bar_stream=std::cout,
                                    const bool& quiet = false) const
    {
        save_SID_histograms("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a image representing the histogram of a chromosome SIDs
     *
     * @param[in] filename is the image filename
     * @param[in] chromosome_id is the identifier of the chromosome whose
     *          SIDs must be represented
     */
    void save_SID_histogram(const std::filesystem::path& filename,
                            const ChromosomeId& chromosome_id) const;
#endif // WITH_MATPLOT

    /**
     * @brief Save sequencing sample set statistics in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param[in,out] archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & data_dir
                & *stats_map;
    }

    /**
     * @brief Load sequencing sample set statistics from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param[in,out] archive is the input archive
     * @return the sequencing sample set statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static SampleSetStatistics load(ARCHIVE& archive)
    {
        std::string data_dir;

        archive & data_dir;

        SampleSetStatistics statistics(data_dir);

        archive & *(statistics.stats_map);

        for (const auto& [sample_name, sample_stats] : *(statistics.stats_map) ) {
            for (const auto& chr_id : sample_stats.get_chromosome_ids()) {
                statistics.chr_ids.insert(chr_id);
            }
        }

        return statistics;
    }

    template<typename RANDOM_GENERATOR>
    friend class ReadSimulator;
};

template<typename RANDOM_GENERATOR>
class ReadSimulator
{
public:
    /**
     * @brief The type of simulated reads
     */
    enum class ReadType {
        SINGLE_READ,
        PAIRED_READ
    };

    /**
     * @brief SAM generator mode
     */
    enum class Mode {
        CREATE,     // Create
        OVERWRITE,  // Overwrite files
        UPDATE      // Add new files
    };

    ReadType read_type;     //!< The type of the produced reads

    size_t read_size;       //!< The produced-read size
    std::binomial_distribution<uint32_t> insert_size;     //!< The paired-read insert size distribution

    bool write_SAM;         //!< A Boolean flag to write SAM files

private:
    bool save_coverage;     //!< A flag to enable/disable storage of coverage data

    RANDOM_GENERATOR random_generator;  //!< The random generator

    std::filesystem::path output_directory;         //!< The output directory
    std::filesystem::path ref_genome_filename;      //!< The reference genome FASTA filename

    size_t num_of_reads;    //!< Number of already placed reads

    size_t hamming_distance_threshold;  //!< Hamming distance threshold for reads

    std::string template_name_prefix;   //!< The template name prefix

    /**
     * @brief A structure to maintain read simulation data
     */
    struct ReadSimulationData
    {
        size_t non_covered_allelic_size;   //!< Allelic size to cover yet
        size_t missing_templates;          //!< Number of templates to be placed

        /**
         * @brief The empty constructor
         */
        ReadSimulationData():
            ReadSimulationData{0,0}
        {}

        /**
         * @brief A constructor
         *
         * @param[in] non_covered_allelic_size is the allelic size to cover yet
         * @param[in] missing_templates is the number of templates to be placed
         */
        ReadSimulationData(const size_t& non_covered_allelic_size,
                           const size_t& missing_templates):
            non_covered_allelic_size{non_covered_allelic_size},
            missing_templates{missing_templates}
        {}
    };

    /**
     * @brief Get the template read data
     *
     * @param[in] template_first_position is the template first position
     * @param[in] template_size is the template size
     * @return a vector containing a pair flag-first-base position for
     *       each read in the template
     */
    std::vector<std::pair<int, ChrPosition>>
    get_template_read_data(const ChrPosition& template_first_position,
                           const size_t& template_size)
    {
        std::uniform_int_distribution<int> c_dist(0,1);

        bool rev_comp = c_dist(random_generator)==0;

        std::vector<std::pair<int, ChrPosition>> template_read_data;

        if (read_type == ReadType::PAIRED_READ) {
            int first_flag = 0x1 | 0x2 | 0x20;
            int second_flag = 0x1 | 0x2 | 0x10;

            if (rev_comp) {
                first_flag |= 0x80;
                second_flag |= 0x40;
            } else {
                first_flag |= 0x40;
                second_flag |= 0x80;
            }
            template_read_data.push_back({first_flag, template_first_position});

            template_read_data.push_back({second_flag,
                                          template_first_position+template_size-read_size});

            return template_read_data;
        }

        int flag = 0;
        if (rev_comp) {
            flag = 0x10;
        }

        return {{flag, template_first_position}};
    }

    /**
     * @brief Get the name of a template
     *
     * @param[in] chr_id is the chromosome identifier from which the
     *      template comes from
     * @param[in] template_id is the identifier of the template
     * @return the name of the template whose identifier is
     *      `template_id`.
     */
    std::string get_template_name(const ChromosomeId& chr_id,
                                  const size_t& template_id) const
    {
        std::ostringstream oss_name;

        oss_name << template_name_prefix
                 << static_cast<size_t>(chr_id)
                 << std::setfill('0') << std::setw(11)
                 << template_id;

        return oss_name.str();
    }

    /**
     * @brief Encode cell identifier in a string
     *
     * @param[in] cell_id is the cell identifier to be encoded
     * @param[in] code_length is the cell code length
     * @param[in] alphabet is the code alphabet
     * @return a string that encodes the cell id by using
     *      `alphabet` as alphabet
     */
    std::string
    encode_cell_id(const RACES::Mutants::CellId& cell_id,
                   const uint8_t code_length=8,
                   const std::string alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    {
        std::string code(code_length, '0');

        uint64_t cell_id_value = static_cast<uint64_t>(cell_id);
        for (auto code_it=code.rbegin(); code_it!=code.rend(); ++code_it) {
            const uint64_t remainder = cell_id_value % alphabet.size();
            *code_it = alphabet[remainder];
            cell_id_value /= alphabet.size();
        }

        return code;
    }

    /**
     * @brief Write the SAM alignment line corresponding to a simulated read
     *
     * This method extracts the template from the chromosome sequence, obtains the
     * template reads, applies germinal and somatic mutations to it, simulates
     * sequencing errors, collects sequencing statistics, and, whenever it is
     * required, writes the read alignments in the SAM stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in,out] chr_statistics are the chromosome statistics about a tissue sample
     * @param[in] chr_data is the data about the chromosome from which the simulated
     *   template come from
     * @param[in] cell_id is the identifier of the fragment cell
     * @param[in] germline_mutations is a the map from genomic positions to the
     *   corresponding germline SIDs to be inserted into the read
     * @param[in] somatic_mutations is a map from genomic positions to the corresponding
     *   somatic SIDs
     * @param[in] template_begin_pos is the first position of the simulated template
     * @param[in] template_size is the size of the template
     * @param[out] SAM_stream is a pointer to a SAM stream. If `SAM_stream` is
     *   different from `nullptr` and the hamming distance between the produced read
     *   and the reference fragment is below the threshold, the read alignments is
     *   written in the SAM stream
     * @param[in] sample_name is the tissue sample name
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void process_template(SEQUENCER& sequencer, ChrSampleStatistics& chr_statistics,
                          const RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>& chr_data,
                          const RACES::Mutants::CellId& cell_id,
                          const std::map<GenomicPosition, std::shared_ptr<SID>>& germline_mutations,
                          const std::map<GenomicPosition, std::shared_ptr<SID>>& somatic_mutations,
                          const ChrPosition& template_begin_pos, const size_t& template_size,
                          std::ostream* SAM_stream, const std::string& sample_name="")
    {
        const auto template_read_data = get_template_read_data(template_begin_pos, template_size);

        const size_t template_id = num_of_reads/template_read_data.size();
        const std::string template_name = get_template_name(chr_data.chr_id, template_id);
        Read read[2];
        std::string qual[2];
        size_t hamming_dist[2];
        size_t over_threshold{0};
        for (size_t i=0; i<template_read_data.size(); ++i) {
            const auto& read_first_position = template_read_data[i].second;

            const GenomicPosition genomic_position{chr_data.chr_id, read_first_position};

            read[i] = Read{chr_data.nucleotides, germline_mutations, somatic_mutations,
                           genomic_position, read_size};

            qual[i] = sequencer.simulate_seq(read[i], genomic_position, i==1);
            hamming_dist[i] = read[i].Hamming_distance();

            if (hamming_dist[i] < hamming_distance_threshold) {
                ++over_threshold;

                chr_statistics.add(read[i]);
            }
        }


        if (SAM_stream != nullptr) {
            for (size_t i=0; i<template_read_data.size(); ++i) {
                if (hamming_dist[i] < hamming_distance_threshold) {
                    const auto& genomic_position = read[i].get_genomic_position();

                    int flag = template_read_data[i].first;

                    const auto mapq = 33 + 60;

                    if (over_threshold==1) {
                        flag = 0x10 & flag;
                    }

                    *SAM_stream << template_name                        // QNAME
                                << '\t' << flag                         // FLAG
                                << '\t' << chr_data.name                // RNAME
                                << '\t' << genomic_position.position    // POS
                                << '\t' << mapq                         // MAPQ
                                << '\t' << read[i].get_CIGAR();         // CIGAR
                    if (read_type == ReadType::PAIRED_READ && over_threshold > 1) {
                        const auto paired_pos = read[1-i].get_genomic_position().position;

                        *SAM_stream << '\t' << '='               // RNEXT
                                    << '\t' << paired_pos        // PNEXT
                                    << '\t' << (i==0?"":"-")
                                    << template_size;            // TLEN
                    } else {
                        *SAM_stream << '\t' << '*'   // RNEXT
                                    << '\t' << '0'   // PNEXT
                                    << '\t' << '0';  // TLEN
                    }
                    *SAM_stream << '\t' << read[i].get_sequence()  // SEQ
                                << '\t' << qual[i];                // QUAL

                    // TAGS
                    if (sample_name.size()!=0) {
                        *SAM_stream << "\tRG:Z:" << sample_name;  // The read group
                    }
                    *SAM_stream << "\tNM:i:" << hamming_dist[i];  // The Hamming distance
                    *SAM_stream << "\tCB:Z:"
                                << encode_cell_id(cell_id);       // The cell identifier

                    *SAM_stream << std::endl;
                }
            }
        }
    }

    /**
     * @brief Generate the reads in an allelic fragment
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in,out] sample_simulation_data are the read simulation data relative to the considered `sample`
     * @param[in,out] chr_statistics are the chromosome statistics of the tissue sample
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] sample_name is the name of the considered tissue sample
     * @param[in] cell_id is the identifier of the fragment cell
     * @param[in] germline_fragment is a pointer to the germinal the allele fragment.
     *          No germinal mutation is considered when it is set to `nullptr`
     * @param[in] fragment is the allele fragment for which reads must be generated
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void generate_fragment_reads(SEQUENCER& sequencer,
                                 ReadSimulationData& sample_simulation_data,
                                 ChrSampleStatistics& chr_statistics,
                                 const RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>& chr_data,
                                 const std::string& sample_name,
                                 const RACES::Mutants::CellId& cell_id,
                                 const AlleleFragment* germline_fragment,
                                 const AlleleFragment& fragment,
                                 const size_t& total_steps, size_t& steps,
                                 RACES::UI::ProgressBar& progress_bar, std::ostream* SAM_stream)
    {
        auto insert_size_mean = insert_size.p()*insert_size.t();
        auto insert_size_stddev = (1-insert_size.p())*insert_size_mean;
        auto threshold_template_size = ((read_type==ReadType::PAIRED_READ)
                                          ?2*read_size+insert_size_mean-insert_size_stddev:
                                          read_size);

        if (fragment.size() >= threshold_template_size && sample_simulation_data.missing_templates > 0) {
            const double hit_probability = static_cast<double>(fragment.size())/
                                        sample_simulation_data.non_covered_allelic_size;
            std::binomial_distribution<size_t> b_dist(sample_simulation_data.missing_templates,
                                                        hit_probability);

            const size_t num_of_frag_templates = std::min(b_dist(random_generator),
                                                          sample_simulation_data.missing_templates);

            const double fragment_progress_ratio = (100*static_cast<double>(fragment.size())/num_of_frag_templates)/total_steps;

            const double current_progress = 100*static_cast<double>(steps)/total_steps;

            const auto& somatic_mutations = fragment.get_mutations();

            auto first_possible_begin = fragment.begin();
            auto last_possible_begin = static_cast<ChrPosition>(fragment.end())
                                            -read_size+1;

            std::uniform_int_distribution<ChrPosition> dist(first_possible_begin, last_possible_begin);
            size_t simulated_templates = 0;
            while (simulated_templates < num_of_frag_templates) {
                auto begin_pos = dist(random_generator);

                auto template_size = ((read_type==ReadType::PAIRED_READ)
                                        ?2*read_size+insert_size(random_generator):read_size);

                if (begin_pos+template_size<=fragment.end()+1) {
                    if (germline_fragment != nullptr) {
                        process_template(sequencer, chr_statistics, chr_data, cell_id,
                                         germline_fragment->get_mutations(), somatic_mutations,
                                         begin_pos, template_size, SAM_stream, sample_name);
                    } else {
                        std::map<GenomicPosition, std::shared_ptr<SID>> germline_mutations;

                        process_template(sequencer, chr_statistics, chr_data, cell_id,
                                         germline_mutations, somatic_mutations,
                                         begin_pos, template_size, SAM_stream, sample_name);
                    }

                    num_of_reads += ((read_type == ReadType::PAIRED_READ)?2:1);
                    ++simulated_templates;
                    --sample_simulation_data.missing_templates;
                }
                progress_bar.set_progress(current_progress+simulated_templates*fragment_progress_ratio);
            }
        }

        sample_simulation_data.non_covered_allelic_size -= fragment.size();
        steps += fragment.size();
    }

    /**
     * @brief Generate simulated reads on a cell chromosome and write their SAM alignments
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in,out] sample_simulation_data is the read simulation data of the sample
     * @param[in,out] chr_statistics are the chromosome statistics of the tissue sample
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] sample_name is the name of the sample
     * @param[in] cell_id is the identifier of the considered cell
     * @param[in] germline_chr_mut is a pointer to the germinal chromosome mutations of the
     *          considered cell. No germinal mutation is considered when it is set to `nullptr`
     * @param[in] chr_mutations are the somatic genome mutations of the considered cell
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     * @return the sample statistics of the chromosome
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void
    generate_chromosome_reads(SEQUENCER& sequencer,
                              ReadSimulationData& sample_simulation_data,
                              ChrSampleStatistics& chr_statistics,
                              const RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>& chr_data,
                              const std::string& sample_name,
                              const Mutants::CellId& cell_id,
                              const ChromosomeMutations* germline_chr_mut,
                              const ChromosomeMutations& chr_mutations,
                              const size_t& total_steps, size_t& steps,
                              RACES::UI::ProgressBar& progress_bar,
                              std::ostream* SAM_stream)
    {
        Allele const* germline_allele{nullptr};
        for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
            if (germline_chr_mut != nullptr) {
                const auto& germline_allele_id = allele.get_history().front();
                germline_allele = &(germline_chr_mut->get_allele(germline_allele_id));
            }

            AlleleFragment const* germline_fragment{nullptr};
            for (const auto& [position, fragment] : allele.get_fragments()) {
                if (germline_allele != nullptr) {
                    // searching for the germline fragment AFTER `fragment`
                    auto germline_it = germline_allele->get_fragments().upper_bound(position);

                    // update `germline_mutations` so to the last germline fragment starting
                    // BEFORE or in the SAME POSITION of `fragment`
                    --germline_it;

                    germline_fragment = &(germline_it->second);
                }

                // since `cell_mutations` derives from `germline_mutations`, `germline_allele`
                // fully contains `allele`
                generate_fragment_reads(sequencer, sample_simulation_data, chr_statistics,
                                        chr_data, sample_name, cell_id, germline_fragment,
                                        fragment, total_steps, steps, progress_bar,
                                        SAM_stream);
            }

            progress_bar.set_progress(100*steps/total_steps);
        }
    }

    /**
     * @brief Generate simulated tumour reads on a chromosome and write their SAM alignments
     *
     * This method takes a sample target and the corresponding sample forest, it generates
     * simulated reads over the sample tumour cells up to a specified coverage, and write the
     * corresponding SAM alignments in a stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in,out] sample_simulation_data is the read simulation data of the sample
     * @param[in] chr_statistics are the chromosome statistics of the tissue sample
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] sample_target is the sample target
     * @param[in] sample_forest is the sample forest
     * @param[in] wild_type_genomes are the wild type genomes
     * @param[in] germline_chr_mut is a pointer to the germinal chromosome mutations of the
     *          considered cell. No germinal mutation is considered when it is set to `nullptr`
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     * @return the sample statistics of the chromosome
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    ChrSampleStatistics&
    generate_chr_tumour_reads(SEQUENCER& sequencer, ReadSimulationData& sample_simulation_data,
                              ChrSampleStatistics& chr_statistics,
                              const RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>& chr_data,
                              const SequencingTargets::SampleTarget& sample_target,
                              const PhylogeneticForest& sample_forest,
                              const ChromosomeMutations* germline_chr_mut, const size_t& total_steps,
                              size_t& steps, RACES::UI::ProgressBar& progress_bar,
                              std::ostream* SAM_stream=nullptr)
    {
        for (const auto& [cell_id, chr_mutations]:
                sample_forest.get_leaf_mutation_tour(chr_data.chr_id)) {

            generate_chromosome_reads(sequencer, sample_simulation_data, chr_statistics,
                                      chr_data, sample_target.sample_name, cell_id,
                                      germline_chr_mut, chr_mutations, total_steps,
                                      steps, progress_bar, SAM_stream);
        }

        return chr_statistics;
    }


    /**
     * @brief Generate simulated wild type reads on a chromosome and write their SAM alignments
     *
     * This method takes a sample target and the corresponding sample forest, it generates
     * simulated reads over the sample wild type cells up to a specified coverage, and write the
     * corresponding SAM alignments in a stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in,out] sample_simulation_data is the read simulation data of the sample
     * @param[in] chr_statistics are the chromosome statistics of the tissue sample
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] sample_target is the sample target
     * @param[in] sample_forest is the sample forest
     * @param[in] wild_type_genomes are the wild type genomes
     * @param[in] germline_chr_mut is a pointer to the germinal chromosome mutations of the
     *          considered cell. No germinal mutation is considered when it is set to `nullptr`
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     * @return the sample statistics of the chromosome
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    ChrSampleStatistics&
    generate_chr_wild_type_reads(SEQUENCER& sequencer, ReadSimulationData& sample_simulation_data,
                                 ChrSampleStatistics& chr_statistics,
                                 const RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>& chr_data,
                                 const SequencingTargets::SampleTarget& sample_target,
                                 const std::map<Mutants::CellId, CellGenomeMutations>& wild_type_genomes,
                                 const ChromosomeMutations* germline_chr_mut, const size_t& total_steps,
                                 size_t& steps, RACES::UI::ProgressBar& progress_bar,
                                 std::ostream* SAM_stream=nullptr)
    {
        for (const auto& [cell_id, multiplicity] : sample_target.num_of_wild_types) {
            const auto& wild_cell_mutations = wild_type_genomes.at(cell_id);

            for (size_t i=0; i<multiplicity; ++i) {
                const auto& chr_mutations = wild_cell_mutations.get_chromosome(chr_data.chr_id);

                generate_chromosome_reads(sequencer, sample_simulation_data, chr_statistics,
                                          chr_data, sample_target.sample_name, cell_id,
                                          germline_chr_mut, chr_mutations, total_steps,
                                          steps, progress_bar, SAM_stream);
            }
        }

        return chr_statistics;
    }

    /**
     * @brief Generate simulated reads on a chromosome and write their SAM alignments
     *
     * This method takes a list of genome mutations and a reference genome chromosome, it
     * generates simulated reads over it up to a specified coverage of the mutated genomes,
     * and write the corresponding SAM alignments in a stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in,out] sample_simulation_data is the read simulation data of the sample
     * @param[in] chr_statistics are the chromosome statistics of the tissue sample
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] sample_target is the sample target
     * @param[in] forest is the phylogenetic forest
     * @param[in] wild_type_genomes are the wild type genomes
     * @param[in] with_germinal_mutations is a Boolean flag to enable/disable germinal mutations
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     * @return the sample statistics of the chromosome
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    ChrSampleStatistics
    generate_chromosome_reads(SEQUENCER& sequencer, ReadSimulationData& sample_simulation_data,
                              ChrSampleStatistics chr_statistics,
                              const RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>& chr_data,
                              const SequencingTargets::SampleTarget& sample_target,
                              const PhylogeneticForest& forest,
                              const std::map<Mutants::CellId, CellGenomeMutations>& wild_type_genomes,
                              const bool& with_germinal_mutations, const size_t& total_steps,
                              size_t& steps, RACES::UI::ProgressBar& progress_bar,
                              std::ostream* SAM_stream=nullptr)
    {
        ChromosomeMutations const* germline_chr_mut{nullptr};
        if (with_germinal_mutations) {
            const auto& germline_mutations = forest.get_germline_mutations();

            germline_chr_mut = &(germline_mutations.get_chromosome(chr_data.chr_id));
        }

        if (sample_target.sample_id != NORMAL_SAMPLE_ID) {
            const auto sample_forest = forest.get_subforest_for({sample_target.sample_name});
            generate_chr_tumour_reads(sequencer, sample_simulation_data, chr_statistics,
                                      chr_data, sample_target, sample_forest, germline_chr_mut,
                                      total_steps, steps, progress_bar, SAM_stream);
        }

        generate_chr_wild_type_reads(sequencer, sample_simulation_data, chr_statistics,
                                     chr_data, sample_target, wild_type_genomes,
                                     germline_chr_mut, total_steps, steps, progress_bar,
                                     SAM_stream);

        return chr_statistics;
    }

    /**
     * @brief Generate simulated reads on a chromosome and write their SAM alignments
     *
     * This method takes a list of genome mutations and a reference genome chromosome, it
     * generates simulated reads over it up to a specified coverage of the mutated genomes,
     * and write the corresponding SAM alignments in a stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in, out] statistics are the sample statistics that will be updated with the chromosome read data
     * @param[in,out] simulation_data is a list of simulated data corresponding to `mutations_list` elements
     * @param[in] targets are the sequencing targets
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in] missed_SID_statistics is a Boolean flag that enable/disable statistics about SIDs that
     *      never occurred in the chromosome reads
     * @param[in] germinal_statistics is a Boolean flag that enable/disable statistics about germinal mutations
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void generate_chromosome_reads(SEQUENCER& sequencer, SampleSetStatistics& statistics,
                                   std::map<Mutants::Evolutions::TissueSampleId, ReadSimulationData>& simulation_data,
                                   const RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>& chr_data,
                                   const SequencingTargets& targets, const size_t& total_steps, size_t& steps,
                                   const bool& missed_SID_statistics, const bool& germinal_statistics,
                                   RACES::UI::ProgressBar& progress_bar,
                                   std::ostream* SAM_stream=nullptr)
    {
        auto chr_name = GenomicPosition::chrtos(chr_data.chr_id);
        progress_bar.set_message("Processing chr. " + chr_name);

        ChrSampleStatistics basic_chr_stats{chr_data.chr_id,
                                            static_cast<GenomicRegion::Length>(chr_data.length),
                                            targets.forest(), missed_SID_statistics, germinal_statistics};

        for (const auto& [sample_id, sample_target] : targets.sample_targets()) {
            auto chr_stats = generate_chromosome_reads(sequencer, simulation_data.at(sample_id), basic_chr_stats,
                                                       chr_data, sample_target, targets.forest(),
                                                       targets.wild_type_genomes(), targets.with_germinal(),
                                                       total_steps, steps, progress_bar, SAM_stream);

            statistics.add_chr_statistics(sample_target.sample_name, std::move(chr_stats),
                                          save_coverage);
        }
    }

    /**
     * @brief Get the SAM stream
     *
     * @param[in] chr_data in the information available about the chromosome to be considered
     * @param[in] targets are the sequencing targets
     * @param[in] mutations_list is a list of sample mutations
     * @param[in] base_name is the prefix of the filename
     * @return the SAM stream
     */
    std::ofstream
    get_SAM_stream(const RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>& chr_data,
                   const SequencingTargets& targets,
                   const std::string& base_name="chr_",
                   const std::string& platform="ILLUMINA") const
    {
        auto chr_str = GenomicPosition::chrtos(chr_data.chr_id);

        std::string SAM_filename = base_name + chr_str + ".sam";
        auto SAM_file_path = output_directory/SAM_filename;

        size_t i{0};
        while (std::filesystem::exists(SAM_file_path)) {
            if (!std::filesystem::is_regular_file(SAM_file_path)) {
                throw std::runtime_error("\"" + to_string(SAM_file_path)
                                         + "\" is not a regular file");
            }
            SAM_filename = base_name + chr_str + "_"
                            + std::to_string(i++) + ".sam";
            SAM_file_path = output_directory/SAM_filename;
        }

        auto SAM_stream = std::ofstream(SAM_file_path);

        // write the SAM header
        SAM_stream << "@HD\tVN:1.6\tSO:unknown" << std::endl;
        SAM_stream << "@SQ\tSN:" << chr_data.name << "\tLN:" << chr_data.length
                << "\tAN:chr"<< chr_str << ", chromosome" << chr_str
                << ",chromosome_" << chr_str << ",chr_" << chr_str << std::endl;
        for (const auto& [sample_id, sample_target] : targets.sample_targets()) {
            SAM_stream << "@RG\tID:" << sample_target.sample_name
                       << "\tSM:" << sample_target.sample_name
                       << "\tPL:" << platform << std::endl;
        }

        return SAM_stream;
    }

    /**
     * @brief Compute the initial read simulation data
     *
     * @param[in] targets are the sequencing targets
     * @param[in] chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param[in] coverage is the aimed coverage
     * @return the list of the initial sample read simulation data. The order
     *       of the returned list matches that of `mutations_list`
     */
    std::map<Mutants::Evolutions::TissueSampleId, ReadSimulationData>
    get_initial_data(const SequencingTargets& targets,
                     const std::set<ChromosomeId>& chromosome_ids,
                     const double& coverage)
    {
        // compute complete and non relevant allelic sizes
        std::map<Mutants::Evolutions::TissueSampleId, ReadSimulationData> simulation_data;
        std::map<Mutants::Evolutions::TissueSampleId, size_t> non_relevant;

        const auto& forest = targets.forest();
        for (const auto& [leaf_id, leaf_mutations]: forest.get_leaf_mutation_tour()) {
            const auto leaf = forest.get_node(leaf_id);
            const auto sample_id = leaf.get_sample().get_id();

            size_t& non_relevant_in_sample = non_relevant[sample_id];
            ReadSimulationData& sample_data = simulation_data[sample_id];
            for (const auto& [chr_id, cell_chr]: leaf_mutations.get_chromosomes()) {
                const size_t allelic_size = cell_chr.allelic_size();

                if (chromosome_ids.count(chr_id)) {
                    sample_data.non_covered_allelic_size += allelic_size;
                } else {
                    non_relevant_in_sample += allelic_size;
                }
            }
        }

        const auto& wild_type_genomes = targets.wild_type_genomes();
        for (const auto& [sample_id, sample_target]: targets.sample_targets()) {
            size_t& non_relevant_in_sample = non_relevant[sample_id];
            ReadSimulationData& sample_data = simulation_data[sample_id];
            for (const auto& [cell_id, num_of_cells]: sample_target.num_of_wild_types) {
                const auto& wild_type_genome = wild_type_genomes.at(cell_id);
                for (const auto& [chr_id, cell_chr]: wild_type_genome.get_chromosomes()) {
                    const size_t allelic_size = num_of_cells * cell_chr.allelic_size();

                    if (chromosome_ids.count(chr_id)) {
                        sample_data.non_covered_allelic_size += allelic_size;
                    } else {
                        non_relevant_in_sample += allelic_size;
                    }
                }
            }
        }

        // compute missing templates
        const size_t total_read_size = (read_type==ReadType::PAIRED_READ?2:1)*read_size;
        const size_t num_of_templates =  static_cast<size_t>((forest.get_germline_mutations().size()
                                                               *coverage)/total_read_size);
        for (auto& [sample_id, sample_data] : simulation_data) {

            const size_t& non_relevant_in_sample = non_relevant[sample_id];

            // compute the probability for a template to cover a nucleotide in an allele of
            // a relevant chromosome
            const size_t total_allelic_size = sample_data.non_covered_allelic_size + non_relevant_in_sample;
            const double hit_probability = 1-static_cast<double>(non_relevant_in_sample)/total_allelic_size;

            // get the number of templates that cover a relevant chromosome
            std::binomial_distribution<size_t> b_dist(num_of_templates, hit_probability);
            sample_data.missing_templates = b_dist(random_generator);
        }

        return simulation_data;
    }

    /**
     * @brief Get the identifiers of the chromosomes in the genome
     *
     * This method deduces the identifiers of the chromosomes in the
     * considered genome.
     *
     * @param[in] genome_mutations is a genome mutations
     * @return the set of identifiers of the genome chromosomes
     */
    static std::set<ChromosomeId>
    get_genome_chromosome_ids(const GenomeMutations& genome)
    {
        std::set<ChromosomeId> chr_ids;
        for (const auto& [chr_id, chr_mutations] : genome.get_chromosomes()) {
            chr_ids.insert(chr_id);
        }

        return chr_ids;

    }

    /**
     * @brief Generate simulated reads
     *
     * This method takes a list of sample mutations and it generates simulated reads up to a
     * specified coverage of the mutated genomes. The base-coverage and SID occurrences are
     * collected and, upon request, the SAM alignments corresponding to the produced reads
     * are saved.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in,out] targets are the sequencing targets
     * @param[in] chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param[in] coverage is the aimed coverage
     * @param[in] base_name is the prefix of the filename
     * @param[in] missed_SID_statistics is a Boolean flag that enable/disable statistics about SIDs that
     *      never occurred in the chromosome reads
     * @param[in] germinal_statistics is a Boolean flag that enable/disable statistics about germinal mutations
     * @param[in] progress_bar is the progress bar
     * @return the statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleSetStatistics generate_reads(SEQUENCER& sequencer, SequencingTargets& targets,
                                       const std::set<ChromosomeId>& chromosome_ids,
                                       const double& coverage, const std::string& base_name,
                                       const bool& missed_SID_statistics, const bool& germinal_statistics,
                                       UI::ProgressBar& progress_bar)
    {
        auto read_simulation_data = get_initial_data(targets, chromosome_ids, coverage);

        using namespace RACES::IO::FASTA;

        IndexedReader<ChromosomeData<Sequence>> chr_reader(ref_genome_filename);

        const auto relevant_chr_names = get_relevant_chr_names(chr_reader, targets.forest(),
                                                               chromosome_ids);

        size_t total_steps = 2*relevant_chr_names.size()*targets.sample_targets().size();
        for (const auto& [sample_id, sample_data]: read_simulation_data) {
            total_steps += sample_data.non_covered_allelic_size;
        }

        SampleSetStatistics statistics(output_directory);

        size_t steps{0};
        ChromosomeData<Sequence> chr_data;
        for (const auto& chr_name: relevant_chr_names) {
            chr_reader.read(chr_data, chr_name, progress_bar);
            progress_bar.set_progress((100*(++steps))/total_steps);

            if (write_SAM) {
                std::ofstream SAM_stream = get_SAM_stream(chr_data, targets, base_name,
                                                          sequencer.get_platform_name());

                generate_chromosome_reads(sequencer, statistics, read_simulation_data,
                                          chr_data, targets, total_steps, steps,
                                          missed_SID_statistics, germinal_statistics,
                                          progress_bar, &SAM_stream);
            } else {
                generate_chromosome_reads(sequencer, statistics, read_simulation_data,
                                          chr_data, targets, total_steps, steps,
                                          missed_SID_statistics, germinal_statistics,
                                          progress_bar);
            }

            statistics.chr_ids.insert(chr_data.chr_id);
        }

        progress_bar.set_progress(100, "Reads simulated");

        return statistics;
    }

    /**
     * @brief A constructor
     *
     * @param[in] output_directory is the output directory
     * @param[in] ref_genome_filename is reference genome filename
     * @param[in] read_type is the type of the produced-read, i.e., single or paired-read
     * @param[in] read_size is the size of the output reads
     * @param[in] insert_size_distribution is the insert size distribution
     * @param[in] mode is the SAM generator output mode
     * @param[in] save_coverage is a flag to enable/disable storage of coverage data
     * @param[in] template_name_prefix is the template name prefix
     * @param[in] seed is the random generator seed
     */
    ReadSimulator(const std::filesystem::path& output_directory,
                  const std::filesystem::path& ref_genome_filename,
                  const ReadType read_type, const size_t& read_size,
                  const std::binomial_distribution<uint32_t>& insert_size_distribution,
                  const Mode mode, const bool save_coverage,
                  const std::string& template_name_prefix, const int& seed):
        read_type(read_type), read_size(read_size), insert_size(insert_size_distribution),
        write_SAM(false), save_coverage(save_coverage), random_generator(seed),
        output_directory(output_directory), ref_genome_filename(ref_genome_filename),
        num_of_reads(0), hamming_distance_threshold(30),
        template_name_prefix(template_name_prefix)
    {
        namespace fs = std::filesystem;

        if (fs::exists(this->output_directory)) {
            if (mode==Mode::CREATE) {
                throw std::domain_error("\"" + to_string(output_directory)
                                        + "\" already exists");
            }

            if (!fs::is_directory(this->output_directory)) {
                throw std::domain_error("\"" + to_string(output_directory)
                                        + "\" is not a directory");
            }

            if (mode==Mode::OVERWRITE) {
                fs::remove_all(this->output_directory);
            }
        }

        fs::create_directory(this->output_directory);

        if (!fs::exists(this->ref_genome_filename)) {
            throw std::domain_error("\"" + to_string(ref_genome_filename)
                                    + "\" does not exist");
        }

        if (!fs::is_regular_file(this->ref_genome_filename)) {
            throw std::domain_error("\"" + to_string(ref_genome_filename)
                                    + "\" is not a regular file");
        }

        std::ifstream ref_stream(ref_genome_filename);
        int c = ref_stream.get();

        if (c != EOF && c != '>') {
            throw std::domain_error("\"" + to_string(ref_genome_filename)
                                    + "\" is not a FASTA file");
        }
    }

    /**
     * @brief Get the sorted list of relevant names
     *
     * This method identifies the chromosomes whose identifiers appear
     * in both `mutations_list` and `chromosome_ids` and build the
     * list of their names sorted according their position in the
     * file indexed by `chr_reader`.
     *
     * @param[in] chr_reader is the indexed reader of the FASTA file
     *      containing the chromosome sequences
     * @param[in] forest is a phylogenetic forest
     * @param[in] chromosome_ids is a list of chromosome identifiers
     * @return the list of their names sorted according their
     *      position in the file indexed by `chr_reader`
     */
    std::list<std::string>
    get_relevant_chr_names(const RACES::IO::FASTA::IndexedReader<RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>>& chr_reader,
                           const PhylogeneticForest& forest,
                           const std::set<ChromosomeId>& chromosome_ids)
    {
        const auto chr_ids = get_genome_chromosome_ids(forest.get_germline_mutations());
        std::set<ChromosomeId> relevant_ids;
        std::set_intersection(chr_ids.begin(), chr_ids.end(),
                              chromosome_ids.begin(), chromosome_ids.end(),
                              std::inserter(relevant_ids, relevant_ids.begin()));

        std::list<std::string> relevant_names;
        for (const auto& chr_id : relevant_ids) {
            relevant_names.push_back(GenomicPosition::chrtos(chr_id));
        }

        return chr_reader.get_index().sort_according_offset(relevant_names);
    }

public:

    /**
     * @brief The empty constructor
     */
    ReadSimulator():
        read_type(ReadType::SINGLE_READ), read_size(0),
        insert_size(std::binomial_distribution<uint32_t>(0,0)), write_SAM(false),
        save_coverage(false), random_generator(), output_directory(""), ref_genome_filename(""),
        num_of_reads(0), hamming_distance_threshold(30), template_name_prefix("r")
    {}

    /**
     * @brief Create a simulator that produces single reads
     *
     * @param[in] output_directory is the output directory
     * @param[in] ref_genome_filename is reference genome filename
     * @param[in] read_size is the size of the output reads
     * @param[in] mode is the SAM generator output mode
     * @param[in] save_coverage is a flag to enable/disable storage of coverage data
     * @param[in] template_name_prefix is the template name prefix
     * @param[in] seed is the random generator seed
     */
    ReadSimulator(const std::filesystem::path& output_directory,
                  const std::filesystem::path& ref_genome_filename,
                  const size_t& read_size, const Mode mode=Mode::CREATE,
                  const bool& save_coverage=false,
                  const std::string& template_name_prefix="r", const int& seed=0):
        ReadSimulator(output_directory, ref_genome_filename, ReadType::SINGLE_READ, read_size,
                      std::binomial_distribution<uint32_t>(0, 0), mode, save_coverage,
                      template_name_prefix, seed)
    {}

    /**
     * @brief Create a simulator that produces paired reads
     *
     * @param[in] output_directory is the SAM output directory
     * @param[in] ref_genome_filename is reference genome filename
     * @param[in] read_size is the size of the output reads
     * @param[in] insert_size_distribution is the insert size distribution
     * @param[in] mode is the SAM generator output mode
     * @param[in] save_coverage is a flag to enable/disable storage of coverage data
     * @param[in] template_name_prefix is the template name prefix
     * @param[in] seed is the random generator seed
     */
    ReadSimulator(const std::filesystem::path& output_directory,
                  const std::filesystem::path& ref_genome_filename,
                  const size_t& read_size,
                  const std::binomial_distribution<uint32_t>& insert_size, const Mode mode=Mode::CREATE,
                  const bool& save_coverage=false,
                  const std::string& template_name_prefix="r", const int& seed=0):
        ReadSimulator(output_directory, ref_genome_filename, ReadType::PAIRED_READ, read_size,
                      insert_size, mode, save_coverage, template_name_prefix, seed)
    {}

    /**
     * @brief Generate simulated reads for a list of sample genome mutations
     *
     * This method generates simulated reads for a list of sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] forest is a phylogenetic forest
     * @param[in] chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param[in] coverage is the aimed coverage
     * @param[in] produce_normal_sample is a Boolean flag to produce/avoid a normal sample
     * @param[in] purity is ratio between the number of sampled tumour cells, which are
     *              represented in the mutation list, and that of the overall sampled
     *              cells which contains normal cells too
     * @param[in] with_pre_neoplastic is a Boolean flag to consider/avoid pre-neoplastic
     *              mutations (default: true)
     * @param[in] with_germinal is a Boolean flag to consider/avoid germinal mutations
     *              (default: true)
     * @param[in] base_name is the prefix of the filename (default: "chr_")
     * @param[in] missed_SID_statistics is a Boolean flag that enable/disable statistics
     *      about SIDs that never occurred in the chromosome reads
     * @param[in] germinal_statistics is a Boolean flag that enable/disable statistics
     *      about germinal mutations
     * @param[in] progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param[in] quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleSetStatistics operator()(SEQUENCER& sequencer, const PhylogeneticForest& forest,
                                   const std::set<ChromosomeId>& chromosome_ids,
                                   const double& coverage, const bool& produce_normal_sample,
                                   const double purity, const bool& with_pre_neoplastic=true,
                                   const bool& with_germinal=true, const std::string& base_name="chr_",
                                   const bool& missed_SID_statistics=false,
                                   const bool& germinal_statistics=false,
                                   std::ostream& progress_bar_stream=std::cout,
                                   const bool quiet=false)
    {
        using namespace RACES;

        if (coverage < 0) {
            throw std::domain_error("The coverage value must be non-negative: got "
                                    + std::to_string(coverage) + ".");
        }

        if (purity < 0 || purity > 1) {
            throw std::domain_error("The purity value must be in [0,1]: got "
                                    + std::to_string(purity) + ".");

        }

        if (read_type == ReadType::PAIRED_READ
                && !sequencer.supports_paired_reads()) {
            throw std::domain_error(sequencer.get_model_name()
                                    + " does not support paired reads.");
        }

        auto targets = SequencingTargets::get_targets(random_generator, forest,
                                                      produce_normal_sample,
                                                      purity,
                                                      with_pre_neoplastic,
                                                      with_germinal);

        UI::ProgressBar progress_bar(progress_bar_stream, quiet);

        return generate_reads<SEQUENCER>(sequencer, targets, chromosome_ids, coverage,
                                         base_name, missed_SID_statistics,
                                         germinal_statistics, progress_bar);
    }

    /**
     * @brief Generate simulated reads for a list of sample genome mutations
     *
     * This method generates simulated reads for a list of sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] forest is a phylogenetic forest
     * @param[in] chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param[in] coverage is the aimed coverage
     * @param[in] with_pre_neoplastic is a Boolean flag to consider/avoid pre-neoplastic
     *              mutations (default: true)
     * @param[in] with_germinal is a Boolean flag to consider/avoid germinal mutations
     *              (default: true)
     * @param[in] base_name is the prefix of the filename (default: "chr_")
     * @param[in] missed_SID_statistics is a Boolean flag that enable/disable statistics
     *      about SIDs that never occurred in the chromosome reads
     * @param[in] germinal_statistics is a Boolean flag that enable/disable statistics
     *      about germinal mutations
     * @param[in] progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param[in] quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    inline SampleSetStatistics operator()(SEQUENCER& sequencer, const PhylogeneticForest& forest,
                                          const std::set<ChromosomeId>& chromosome_ids,
                                          const double& coverage, const bool& with_pre_neoplastic=true,
                                          const bool& with_germinal=true, const std::string& base_name="chr_",
                                          const bool& missed_SID_statistics=false,
                                          const bool& germinal_statistics=false,
                                          std::ostream& progress_bar_stream=std::cout,
                                          const bool quiet=false)
    {
        return operator()(sequencer, forest, chromosome_ids, coverage, false, 1.0, with_pre_neoplastic,
                          with_germinal, base_name, missed_SID_statistics, germinal_statistics,
                          progress_bar_stream, quiet);
    }

    /**
     * @brief Generate simulated reads for a list of sample genome mutations
     *
     * This method generates simulated reads for a list of sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] forest is a phylogenetic forest
     * @param[in] coverage is the aimed coverage
     * @param[in] produce_normal_sample is a Boolean flag to produce/avoid a normal sample
     * @param[in] purity is ratio between the number of sampled tumour cells, which are
     *              represented in the mutation list, and that of the overall sampled
     *              cells which contains normal cells too
     * @param[in] with_pre_neoplastic is a Boolean flag to consider/avoid pre-neoplastic
     *              mutations (default: true)
     * @param[in] with_germinal is a Boolean flag to consider/avoid germinal mutations
     *              (default: true)
     * @param[in] base_name is the prefix of the filename (default: "chr_")
     * @param[in] missed_SID_statistics is a Boolean flag that enable/disable statistics
     *      about SIDs that never occurred in the chromosome reads
     * @param[in] germinal_statistics is a Boolean flag that enable/disable statistics
     *      about germinal mutations
     * @param[in] progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param[in] quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    inline SampleSetStatistics operator()(SEQUENCER& sequencer, const PhylogeneticForest& forest,
                                          const double& coverage, const bool& produce_normal_sample,
                                          const double purity, const bool& with_pre_neoplastic=true,
                                          const bool& with_germinal=true, const std::string& base_name="chr_",
                                          const bool& missed_SID_statistics=false,
                                          const bool& germinal_statistics=false,
                                          std::ostream& progress_bar_stream=std::cout,
                                          const bool quiet=false)
    {
        const auto chromosome_ids = get_genome_chromosome_ids(forest.get_germline_mutations());

        return operator()(sequencer, forest, chromosome_ids, coverage, produce_normal_sample,
                          purity, with_pre_neoplastic, with_germinal, base_name, missed_SID_statistics,
                          germinal_statistics, progress_bar_stream, quiet);
    }

    /**
     * @brief Generate simulated reads for a list of sample genome mutations
     *
     * This method generates simulated reads for a list of sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] forest is a phylogenetic forest
     * @param[in] coverage is the aimed coverage
     * @param[in] with_pre_neoplastic is a Boolean flag to consider/avoid pre-neoplastic
     *              mutations (default: true)
     * @param[in] with_germinal is a Boolean flag to consider/avoid germinal mutations
     *              (default: true)
     * @param[in] base_name is the prefix of the filename (default: "chr_")
     * @param[in] missed_SID_statistics is a Boolean flag that enable/disable statistics
     *      about SIDs that never occurred in the chromosome reads
     * @param[in] germinal_statistics is a Boolean flag that enable/disable statistics
     *      about germinal mutations
     * @param[in] progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param[in] quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    inline SampleSetStatistics operator()(SEQUENCER& sequencer, const PhylogeneticForest& forest,
                                          const double& coverage, const bool& with_pre_neoplastic=true,
                                          const bool& with_germinal=true, const std::string& base_name="chr_",
                                          const bool& missed_SID_statistics=false,
                                          const bool& germinal_statistics=false,
                                          std::ostream& progress_bar_stream=std::cout,
                                          const bool quiet=false)
    {
        return operator()(sequencer, forest, coverage, false, 1.0, with_pre_neoplastic, with_germinal,
                          base_name, missed_SID_statistics, germinal_statistics, progress_bar_stream,
                          quiet);
    }

    /**
     * @brief Generate simulated reads for a list of sample genome mutations
     *
     * This method generates simulated reads for a list of sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] forest is a phylogenetic forest
     * @param[in] coverage is the aimed coverage
     * @param[in] with_pre_neoplastic is a Boolean flag to consider/avoid pre-neoplastic
     *              mutations (default: true)
     * @param[in] with_germinal is a Boolean flag to consider/avoid germinal mutations
     *              (default: true)
     * @param[in] base_name is the prefix of the filename (default: "chr_")
     * @param[in] missed_SID_statistics is a Boolean flag that enable/disable statistics
     *      about SIDs that never occurred in the chromosome reads
     * @param[in] germinal_statistics is a Boolean flag that enable/disable statistics
     *      about germinal mutations
     * @param[in] progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param[in] quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<RACES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    inline SampleSetStatistics operator()(SEQUENCER& sequencer, const PhylogeneticForest& forest,
                                          const double& coverage, const double& purity,
                                          const bool& with_pre_neoplastic=true,
                                          const bool& with_germinal=true, const std::string& base_name="chr_",
                                          const bool& missed_SID_statistics=false,
                                          const bool& germinal_statistics=false,
                                          std::ostream& progress_bar_stream=std::cout,
                                          const bool quiet=false)
    {
        return operator()(sequencer, forest, coverage, false, purity, with_pre_neoplastic, with_germinal,
                          base_name, missed_SID_statistics, germinal_statistics, progress_bar_stream, quiet);
    }

    /**
     * @brief Get the simulated read size
     *
     * @return the simulated read size
     */
    inline size_t get_read_size() const
    {
        return read_size;
    }

    /**
     * @brief Enable SAM file writing
     *
     * @param[in] enable is a flag to enable(true)/disable(false) SAM file writing
     */
    inline void enable_SAM_writing(const bool& enable)
    {
        write_SAM = enable;
    }
};

}   // SequencingSimulations

}   // Mutations

}   // RACES

#endif // __RACES_READ_SIMULATOR__
