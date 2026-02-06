/**
 * @file read_simulator.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes to simulate sequencing
 * @version 1.24
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

#ifndef __CLONES_READ_SIMULATOR__
#define __CLONES_READ_SIMULATOR__

#include <list>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <cctype>
#include <random>
#include <iomanip> // std::setw and std::setfill

#include "sid.hpp"
#include "genome_mutations.hpp"
#include "read.hpp"

#include "fasta_chr_reader.hpp"

#include "progress_bar.hpp"
#include "variables.hpp"

#include "utils.hpp"

#include "sequencer.hpp"

#if WITH_MATPLOT
#include <matplot/matplot.h>
#endif // WITH_MATPLOT

namespace CLONES
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
     * @param chromosome_id is the identifier of the chromosome whose coverage refer to
     * @param size is the size of the chromosome whose coverage refer to
     */
    ChrCoverage(const ChromosomeId& chromosome_id, const GenomicRegion::Length& size);

    /**
     * @brief Increase the coverage of a region
     *
     * @param begin_pos is the initial position of the region whose coverage must be increase
     * @param read_size is the size of the region whose coverage must be increase
     */
    void increase_coverage(const ChrPosition& begin_pos, const size_t& read_size);

    /**
     * @brief Increase the coverage of a region
     *
     * @param region is the region whose coverage must be increase
     */
    inline void increase_coverage(const GenomicRegion& region)
    {
        increase_coverage(region.begin(),
                          static_cast<const GenomicRegion::Length&>(region.size()));
    }

    /**
     * @brief Increase the coverage of a region
     *
     * @param region is the region whose coverage must be increase
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
     * @param position is the genomic position whose coverage is aimed
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
     * @param chr_coverage is the sequencing chromosome coverage to join
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
     * @param chr_coverage is the sequencing chromosome coverage to join
     * @return a reference to the updated coverage
     */
    ChrCoverage& operator+=(const ChrCoverage& chr_coverage);

    /**
     * @brief Save sequencing chromosome coverage in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
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
     * @param archive is the input archive
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
     * @param mutation is the mutation whose data refer to
     * @param num_of_occurrences is the number of occurrences of `mutation`
     */
    SIDData(const SID& mutation, const BaseCoverage num_of_occurrences);

    /**
     * @brief Update the SID data
     *
     * This method increases the number of occurrences of the
     * SID and adds the cause and type of a SID to the
     * causes and types stored in the SID data.
     *
     * @param mutation is the SID whose cause and type should be
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
     * @param data is the SIDData object whose data is added
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

    /**
     * @brief Add the mutation of a chromosome to the statistics
     *
     * @param chromosome is the chromosome whose mutations must be
     *     added to the statistics
     */
    void account_for(const ChromosomeMutations& chromosome);
public:

    /**
     * @brief The empty constructor
     */
    ChrSampleStatistics();

    /**
     * @brief A constructor
     *
     * @param[in] chromosome_id is the identifier of the chromosome whose sequencing
     *          statistics are collected
     * @param[in] size is the size of the chromosome whose sequencing
     *          statistics are collected
     * @param[in] mutations_list is a list of sample mutations
     */
    ChrSampleStatistics(const ChromosomeId& chromosome_id, const GenomicRegion::Length& size,
                        const std::list<SampleGenomeMutations>& mutations_list={});

    /**
     * @brief Update the data associated to an SID
     *
     * @param mutation is the SID whose data must be updated
     */
    void update(const SID& mutation);

    /**
     * @brief Collect the data from a read
     *
     * @param read is a read
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
     * @brief Get the collected data of an SID
     *
     * @param mutation is an SID
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
     * @param chr_stats is the sequencing chromosome statistics to join
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
     * @param chr_stats is the sequencing chromosome statistics to join
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
     * @param chr_id is the identifier of the chromosome whose coverage is
     *          going to be saved/load
     * @return the filename of the chromosome coverage file
     */
    std::filesystem::path get_coverage_filename(const ChromosomeId& chr_id) const;
public:

    /**
     * @brief A constructor
     *
     * @param data_directory is the path of the directory containing the data files
     * @param sample_name is the name of the sample whose statistics refer to
     * @param save_coverage is a flag to enable/disable storage of coverage data
     */
    SampleStatistics(const std::filesystem::path& data_directory, const std::string& sample_name,
                     const bool save_coverage=false);

    /**
     * @brief Add chromosome statistics to the genome statistics
     *
     * @param chr_stats is the chromosome statistics
     */
    inline void add_chr_statistics(ChrSampleStatistics&& chr_stats)
    {
        add_chr_statistics(chr_stats);
    }

    /**
     * @brief Add chromosome statistics to the genome statistics
     *
     * @param chr_stats is the chromosome statistics
     */
    void add_chr_statistics(const ChrSampleStatistics& chr_stats);

    /**
     * @brief Get the chromosome coverage by its identifier
     *
     * @param chr_id is the id of the chromosome whose sequencing coverage is aimed
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
     * @param mutation is an SID
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
     * @param archive is the output archive
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
     * @param archive is the input archive
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
    std::map<std::string, SampleStatistics> stats_map;  //!< The map of the sample statistics

    std::set<ChromosomeId> repr_chr_ids; //!< The identifiers of the chromosomes whose statistics have been collected

    std::filesystem::path data_dir; //!< The directory in which the data is saved

    /**
     * @brief Stream the header of a VAF csv
     *
     * @param os is the output stream
     * @param separator is the column separator symbol
     * @return the updated stream
     */
    std::ofstream& stream_VAF_csv_header(std::ofstream& os, const char& separator) const;
public:
    /**
     * @brief A constructor
     *
     * @param data_directory is the path of the directory containing the data files
     */
    SampleSetStatistics(const std::filesystem::path& data_directory);

    /**
     * @brief Check whether this object contains statistics about a sample
     *
     * @param sample_name is the sample name
     * @return `true` if and only if this object contains statistics about
     *          a sample named `sample_name`
     */
    inline std::map<std::string, SampleStatistics>::const_iterator find(const std::string& sample_name) const
    {
        return stats_map.find(sample_name);
    }

    /**
     * @brief Get the first constant iterator of a map from the sample names to the corresponding statistics
     *
     * @return the first constant iterator of a map from the sample names to the corresponding statistics
     */
    inline std::map<std::string, SampleStatistics>::const_iterator begin() const
    {
        return stats_map.begin();
    }

    /**
     * @brief Get the last constant iterator of a map from the sample names to the corresponding statistics
     *
     * @return the last constant iterator of a map from the sample names to the corresponding statistics
     */
    inline std::map<std::string, SampleStatistics>::const_iterator end() const
    {
        return stats_map.end();
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
     * @param sample_name is the sample name
     * @param chr_statistics are the statistics about a chromosome
     * @param save_coverage is a flag to enable/disable storage of coverage data
     * @return a constant reference to the current statistics of the sample whose name is
     *          `sample_name`
     */
    const SampleStatistics& add_chr_statistics(const std::string& sample_name,
                                               const ChrSampleStatistics& chr_statistics,
                                               const bool save_coverage=false);

    /**
     * @brief Add a new sample statistics
     *
     * @param sample_statistics are the sample statistics to add
     * @return a constant reference to the added sample statistics
     */
    const SampleStatistics& add_statistics(const SampleStatistics& sample_statistics);

    /**
     * @brief Get a sample statistics
     *
     * @param sample_name is the sample name
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
        return stats_map.size();
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
     * @param base_name is the prefix of the name of the CSV files.
     *          The filename format is "<base_name><chromosome_name>.csv"
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
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
     * @param base_name is the prefix of the name of the CSV files.
     *          The file name format is "<base_name><chromosome_name>.csv"
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_VAF_CSVs(std::ostream& progress_bar_stream=std::cout,
                              const bool& quiet = false) const
    {
        save_VAF_CSVs("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a CSV file reporting the SID VAFs in a chromosome
     *
     * @param filename is the image filename
     * @param chromosome_id is the identifier of the chromosome whose
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
     * @param base_name is the prefix of the name of the CSV files.
     *          The filename format is "<base_name><chromosome_name>_coverage.jpg"
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
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
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_coverage_images(std::ostream& progress_bar_stream=std::cout,
                                     const bool& quiet = false) const
    {
        save_coverage_images("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a image representing the coverage of a chromosome
     *
     * @param filename is the image filename
     * @param chromosome_id is the identifier of the chromosome whose
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
     * @param base_name is the prefix of the name of the CSV files.
     *          The filename syntax is "<base_name><chromosome_name>_hist.jpg"
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
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
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_SID_histograms(std::ostream& progress_bar_stream=std::cout,
                                    const bool& quiet = false) const
    {
        save_SID_histograms("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a image representing the histogram of a chromosome SIDs
     *
     * @param filename is the image filename
     * @param chromosome_id is the identifier of the chromosome whose
     *          SIDs must be represented
     */
    void save_SID_histogram(const std::filesystem::path& filename,
                            const ChromosomeId& chromosome_id) const;
#endif // WITH_MATPLOT

    /**
     * @brief Save sequencing sample set statistics in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & data_dir
                & stats_map;
    }

    /**
     * @brief Load sequencing sample set statistics from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the sequencing sample set statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static SampleSetStatistics load(ARCHIVE& archive)
    {
        std::string data_dir;

        archive & data_dir;

        SampleSetStatistics statistics(data_dir);

        archive & statistics.stats_map;

        for (const auto& [sample_name, sample_stats] : statistics.stats_map ) {
            for (const auto& chr_id : sample_stats.get_chromosome_ids()) {
                statistics.repr_chr_ids.insert(chr_id);
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
            ReadSimulationData(0,0)
        {}

        /**
         * @brief A constructor
         *
         * @param non_covered_allelic_size is the allelic size to cover yet
         * @param missing_templates is the number of templates to be placed
         */
        ReadSimulationData(const size_t& non_covered_allelic_size,
                           const size_t& missing_templates):
            non_covered_allelic_size(non_covered_allelic_size),
            missing_templates(missing_templates)
        {}
    };

    /**
     * @brief Get the template read data
     *
     * @param template_first_position is the template first position
     * @param template_size is the template size
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
     * @param chr_id is the chromosome identifier from which the
     *      template comes from
     * @param template_id is the identifier of the template
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
     * @param cell_id is the cell identifier to be encoded
     * @param code_length is the cell code length
     * @param alphabet is the code alphabet
     * @return a string that encodes the cell id by using
     *      `alphabet` as alphabet
     */
    std::string
    encode_cell_id(const CLONES::Mutants::CellId& cell_id,
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
     * @param[in] chr_data is the data about the chromosome from which the simulated
     *   template come from
     * @param[in] cell_id is the identifier of the fragment cell
     * @param[in] germlines is a map from genomic positions to the corresponding
     *   germline SIDs
     * @param[in] passengers is a map from genomic positions to the corresponding
     *   passenger SIDs
     * @param[in] template_begin_pos is the first position of the simulated template
     * @param[in] template_size is the size of the template
     * @param[in,out] chr_statistics are the chromosome statistics about a tissue sample
     * @param[out] SAM_stream is a pointer to a SAM stream. If `SAM_stream` is
     *   different from `nullptr` and the hamming distance between the produced read
     *   and the reference fragment is below the threshold, the read alignments is
     *   written in the SAM stream
     * @param[in] sample_name is the tissue sample name
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void process_template(SEQUENCER& sequencer,
                          const CLONES::IO::FASTA::ChromosomeData<CLONES::IO::FASTA::Sequence>& chr_data,
                          const CLONES::Mutants::CellId& cell_id,
                          const std::map<GenomicPosition, std::shared_ptr<SID>>& germlines,
                          const std::map<GenomicPosition, std::shared_ptr<SID>>& passengers,
                          const ChrPosition& template_begin_pos,
                          const size_t& template_size, ChrSampleStatistics& chr_statistics,
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

            read[i] = Read{chr_data.nucleotides, germlines, passengers,
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
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] cell_id is the identifier of the fragment cell
     * @param[in] fragment is the allele fragment for which reads must be generated
     * @param[in,out] sample_simulation_data are the read simulation data relative to the considered `sample`
     * @param[in] sample_name is the name of the considered tissue sample
     * @param[in,out] chr_statistics are the chromosome statistics of the tissue sample
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void generate_fragment_reads(SEQUENCER& sequencer,
                                 const CLONES::IO::FASTA::ChromosomeData<CLONES::IO::FASTA::Sequence>& chr_data,
                                 const CLONES::Mutants::CellId& cell_id,
                                 const AlleleFragment& germline_fragment,
                                 const AlleleFragment& fragment,
                                 ReadSimulationData& sample_simulation_data, const std::string& sample_name,
                                 ChrSampleStatistics& chr_statistics,
                                 const size_t& total_steps, size_t& steps,
                                 CLONES::UI::ProgressBar& progress_bar, std::ostream* SAM_stream)
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

            const auto& germlines = germline_fragment.get_mutations();
            const auto& passengers = fragment.get_mutations();

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

                    process_template(sequencer, chr_data, cell_id, germlines, passengers, begin_pos,
                                    template_size, chr_statistics, SAM_stream, sample_name);

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
     * @brief Generate simulated reads on a chromosome and write their SAM alignments
     *
     * This method takes a list of genome mutations and a reference genome chromosome, it
     * generates simulated reads over it up to a specified coverage of the mutated genomes,
     * and write the corresponding SAM alignments in a stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] sample_genome_mutations is a sample genome mutations
     * @param[in,out] sample_simulation_data is the read simulation data of the sample
     * @param[in] chr_statistics are the chromosome statistics of the tissue sample
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     * @return the sample statistics of the chromosome
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    ChrSampleStatistics
    generate_chromosome_reads(SEQUENCER& sequencer,
                              const CLONES::IO::FASTA::ChromosomeData<CLONES::IO::FASTA::Sequence>& chr_data,
                              const SampleGenomeMutations& sample_genome_mutations,
                              ReadSimulationData& sample_simulation_data,
                              ChrSampleStatistics chr_statistics,
                              const size_t& total_steps, size_t& steps,
                              CLONES::UI::ProgressBar& progress_bar,
                              std::ostream* SAM_stream=nullptr)
    {
        for (const auto& cell_mutations: sample_genome_mutations.mutations) {
            const auto& chr_mutations = cell_mutations->get_chromosome(chr_data.chr_id);
            const auto& germline_chr_mut = sample_genome_mutations.germline_mutations->get_chromosome(chr_data.chr_id);

            for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
                const auto& germline_allele_id = allele.get_history().front();
                const auto& germline_allele = germline_chr_mut.get_allele(germline_allele_id);

                for (const auto& [position, fragment] : allele.get_fragments()) {
                    // searching for the germline fragment AFTER `fragment`
                    auto germline_it = germline_allele.get_fragments().upper_bound(position);

                    // update `germline_mutations` so to the last germline fragment starting
                    // BEFORE or in the SAME POSITION of `fragment`
                    --germline_it;

                    // since `cell_mutations` derives from `germline_mutations`, `germline_allele`
                    // fully contains `allele`
                    generate_fragment_reads(sequencer, chr_data, cell_mutations->get_id(),
                                            germline_it->second, fragment, sample_simulation_data,
                                            sample_genome_mutations.name, chr_statistics,
                                            total_steps, steps, progress_bar, SAM_stream);
                }

                progress_bar.set_progress(100*steps/total_steps);
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
     * @param[out] statistics are the statistics about the generated reads
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] simulation_data is a list of simulated data corresponding to `mutations_list` elements
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void generate_chromosome_reads(SEQUENCER& sequencer,
                                   SampleSetStatistics& statistics,
                                   const CLONES::IO::FASTA::ChromosomeData<CLONES::IO::FASTA::Sequence>& chr_data,
                                   const std::list<SampleGenomeMutations>& mutations_list,
                                   std::list<ReadSimulationData>& simulation_data,
                                   const size_t& total_steps, size_t& steps,
                                   CLONES::UI::ProgressBar& progress_bar,
                                   std::ostream* SAM_stream=nullptr)
    {
        auto chr_name = GenomicPosition::chrtos(chr_data.chr_id);
        progress_bar.set_message("Processing chr. " + chr_name);

        auto simulation_data_it = simulation_data.begin();

        ChrSampleStatistics basic_chr_stats{chr_data.chr_id,
				                            static_cast<GenomicRegion::Length>(chr_data.length),
                                            mutations_list};

        std::list<ChrSampleStatistics> chr_stats_list;
        for (const auto& mutations : mutations_list) {
            chr_stats_list.push_back(generate_chromosome_reads(sequencer, chr_data, mutations,
                                                               *simulation_data_it, basic_chr_stats,
                                                               total_steps, steps, progress_bar,
                                                               SAM_stream));
            ++simulation_data_it;
        }

        auto chr_stats_it = chr_stats_list.cbegin();
        for (const auto& mutations : mutations_list) {
            statistics.add_chr_statistics(mutations.name, *chr_stats_it, save_coverage);

            ++chr_stats_it;
        }
    }

    /**
     * @brief Get the SAM stream
     *
     * @param chr_data in the information available about the chromosome to be considered
     * @param mutations_list is a list of sample mutations
     * @param base_name is the prefix of the filename
     * @return the SAM stream
     */
    std::ofstream
    get_SAM_stream(const CLONES::IO::FASTA::ChromosomeData<CLONES::IO::FASTA::Sequence>& chr_data,
                   const std::list<SampleGenomeMutations>& mutations_list,
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
        for (const auto& sample_mutations : mutations_list) {
            SAM_stream << "@RG\tID:" << sample_mutations.name
                       << "\tSM:" << sample_mutations.name
                       << "\tPL:" << platform << std::endl;
        }

        return SAM_stream;
    }

    /**
     * @brief Compute the initial read simulation data
     *
     * @param mutations_list is a list of sample mutations
     * @param chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param coverage is the aimed coverage
     * @return the list of the initial sample read simulation data. The order
     *       of the returned list matches that of `mutations_list`
     */
    std::list<ReadSimulationData>
    get_initial_read_simulation_data(const std::list<SampleGenomeMutations>& mutations_list,
                                     const std::set<ChromosomeId>& chromosome_ids,
                                     const double& coverage)
    {
        std::list<ReadSimulationData> read_simulation_data;

        size_t total_read_size = (read_type==ReadType::PAIRED_READ?2:1)*read_size;
        for (const auto& sample_mutations : mutations_list) {
            size_t missing_templates{0}, non_covered_allelic_size{0}, non_relevant_allelic_size{0};

            // collect the data from all the chromosomes to evaluate the overall hit probability per
            // allelic nucleotide
            for (const auto& [chr_id, germline_chr]: sample_mutations.germline_mutations->get_chromosomes()) {
                missing_templates += static_cast<size_t>((germline_chr.size()
                                                            *coverage)/total_read_size);

                const bool chr_is_non_relevant = !chromosome_ids.count(chr_id);

                for (const auto& cell_mutations: sample_mutations.mutations) {
                    const auto& chr_mutations = cell_mutations->get_chromosome(chr_id);
                    const size_t allelic_size = chr_mutations.allelic_size();

                    non_covered_allelic_size += allelic_size;
                    if (chr_is_non_relevant) {
                        non_relevant_allelic_size += allelic_size;
                    }
                }
            }

            // compute the probability for a template to cover a nucleotide in an allele of
            // a relevant chromosome
            const double hit_probability = 1-static_cast<double>(non_relevant_allelic_size)/
                                                                 non_covered_allelic_size;

            // get the number of templates that cover a relevant chromosome
            std::binomial_distribution<size_t> b_dist(missing_templates, hit_probability);
            missing_templates = b_dist(random_generator);

            // remove the allelic size of non-relevant chromsomes
            non_covered_allelic_size -= non_relevant_allelic_size;

            // build a read simulation data accounting for the relevant non-covered alleles
            // and the missing templates
            read_simulation_data.emplace_back(non_covered_allelic_size, missing_templates);
        }

        return read_simulation_data;
    }

    /**
     * @brief Get the identifiers of the chromosomes in the genome
     *
     * This method deduces the identifiers of the chromosomes in the
     * considered genome. The method assumes all the samples in the
     * mutations list to have the same chromosomes.
     *
     * @param mutations_list is a list of sample mutations
     * @return the set of identifiers of the genome chromosomes
     */
    static std::set<ChromosomeId>
    get_genome_chromosome_ids(const std::list<SampleGenomeMutations>& mutations_list)
    {
        for (const auto& sample_mutations : mutations_list) {
            for (const auto& mutations : sample_mutations.mutations) {
                std::set<ChromosomeId> ids;
                for (const auto& [chr_id, chr_mutations] : mutations->get_chromosomes()) {
                    ids.insert(chr_id);
                }

                return ids;
            }
        }

        return {};
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
     * @param mutations_list is a list of sample mutations
     * @param chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param coverage is the aimed coverage
     * @param base_name is the prefix of the filename
     * @param progress_bar is the progress bar
     * @return the statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleSetStatistics generate_reads(SEQUENCER& sequencer,
                                       const std::list<SampleGenomeMutations>& mutations_list,
                                       const std::set<ChromosomeId>& chromosome_ids,
                                       const double& coverage, const std::string& base_name,
                                       UI::ProgressBar& progress_bar)
    {
        auto read_simulation_data = get_initial_read_simulation_data(mutations_list, chromosome_ids,
                                                                     coverage);


        using namespace CLONES::IO::FASTA;

        IndexedReader<ChromosomeData<Sequence>> chr_reader(ref_genome_filename);

        const auto relevant_chr_names = get_relevant_chr_names(chr_reader, mutations_list,
                                                               chromosome_ids);

        size_t total_steps = 2*relevant_chr_names.size()*mutations_list.size();
        for (const auto& sample_data: read_simulation_data) {
            total_steps += sample_data.non_covered_allelic_size;
        }

        SampleSetStatistics statistics(output_directory);

        size_t steps{0};
        ChromosomeData<Sequence> chr_data;
        for (const auto& chr_name: relevant_chr_names) {
            chr_reader.read(chr_data, chr_name, progress_bar);
            progress_bar.set_progress((100*(++steps))/total_steps);

            if (write_SAM) {
                std::ofstream SAM_stream = get_SAM_stream(chr_data, mutations_list,
                                                          base_name,
                                                          sequencer.get_platform_name());

                generate_chromosome_reads(sequencer, statistics, chr_data, mutations_list,
                                          read_simulation_data, total_steps, steps,
                                          progress_bar, &SAM_stream);
            } else {
                generate_chromosome_reads(sequencer, statistics, chr_data, mutations_list,
                                          read_simulation_data, total_steps,
                                          steps, progress_bar);
            }

            statistics.repr_chr_ids.insert(chr_data.chr_id);
        }

        progress_bar.set_progress(100, "Reads simulated");

        return statistics;
    }

    /**
     * @brief A constructor
     *
     * @param output_directory is the output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_type is the type of the produced-read, i.e., single or paired-read
     * @param read_size is the size of the output reads
     * @param insert_size_distribution is the insert size distribution
     * @param mode is the SAM generator output mode
     * @param save_coverage is a flag to enable/disable storage of coverage data
     * @param template_name_prefix is the template name prefix
     * @param seed is the random generator seed
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
     * @param chr_reader is the indexed reader of the FASTA file
     *      containing the chromosome sequences
     * @param mutations_list is a list of sample mutations
     * @param chromosome_ids is a list of chromosome identifiers
     * @return the list of their names sorted according their
     *      position in the file indexed by `chr_reader`
     */
    std::list<std::string>
    get_relevant_chr_names(const CLONES::IO::FASTA::IndexedReader<CLONES::IO::FASTA::ChromosomeData<CLONES::IO::FASTA::Sequence>>& chr_reader,
                           const std::list<SampleGenomeMutations>& mutations_list,
                           const std::set<ChromosomeId>& chromosome_ids)
    {
        const auto chr_ids = get_genome_chromosome_ids(mutations_list);
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
     * @param output_directory is the output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_size is the size of the output reads
     * @param mode is the SAM generator output mode
     * @param save_coverage is a flag to enable/disable storage of coverage data
     * @param template_name_prefix is the template name prefix
     * @param seed is the random generator seed
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
     * @param output_directory is the SAM output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_size is the size of the output reads
     * @param insert_size_distribution is the insert size distribution
     * @param mode is the SAM generator output mode
     * @param save_coverage is a flag to enable/disable storage of coverage data
     * @param template_name_prefix is the template name prefix
     * @param seed is the random generator seed
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
     * @param sequencer is the sequencer
     * @param mutations_list is a list of sample mutations
     * @param chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param coverage is the aimed coverage
     * @param normal_sample is a normal sample
     * @param purity is ratio between the number of sampled tumour cells, which are
     *              represented in the mutation list, and that of the overall sampled
     *              cells which contains normal cells too
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleSetStatistics operator()(SEQUENCER& sequencer,
                                   std::list<CLONES::Mutations::SampleGenomeMutations> mutations_list,
                                   const std::set<ChromosomeId>& chromosome_ids,
                                   const double& coverage,
                                   const CLONES::Mutations::SampleGenomeMutations& normal_sample,
                                   const double purity, const std::string& base_name="chr_",
                                   std::ostream& progress_bar_stream=std::cout,
                                   const bool quiet=false)
    {
        using namespace CLONES;

        if (coverage < 0) {
            throw std::domain_error("The coverage value must be non-negative: got "
                                    + std::to_string(coverage) + ".");
        }

        if (purity < 0 || purity > 1) {
            throw std::domain_error("The purity value must be in [0,1]: got "
                                    + std::to_string(purity) + ".");

        }

        if (purity < 1 && normal_sample.mutations.size()==0) {
            throw std::domain_error("The purity value is lower than 1 and "
                                    "the normal sample does not contain cells.");

        }

        if (read_type == ReadType::PAIRED_READ && !sequencer.supports_paired_reads()) {
            throw std::domain_error(sequencer.get_model_name()
                                    + " does not support paired reads.");
        }

        if (mutations_list.size()==0) {
            return SampleSetStatistics(output_directory);
        }

        std::vector<std::shared_ptr<CellGenomeMutations>> normal_cells(normal_sample.mutations.begin(),
                                                                       normal_sample.mutations.end());

        std::uniform_int_distribution<size_t> selector(0, normal_cells.size()-1);

        for (auto& mutation_list : mutations_list) {
            if (purity>0) {
                size_t tumour_number = mutation_list.mutations.size();
                size_t normal_number = static_cast<size_t>(tumour_number*(1-purity)/purity);

                for (size_t i=0; i<normal_number; ++i) {
                    mutation_list.mutations.push_back(normal_cells[selector(random_generator)]);
                }
            } else {
                mutation_list.mutations = normal_sample.mutations;
            }
        }

        UI::ProgressBar progress_bar(progress_bar_stream, quiet);

        return generate_reads<SEQUENCER>(sequencer, mutations_list, chromosome_ids,
                                         coverage, base_name, progress_bar);
    }

    /**
     * @brief Generate simulated reads for a list of sample genome mutations
     *
     * This method generates simulated reads for a list of sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param sequencer is the sequencer
     * @param mutations_list is a list of sample mutations
     * @param chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param coverage is the aimed coverage
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    inline SampleSetStatistics operator()(SEQUENCER& sequencer,
                                          std::list<CLONES::Mutations::SampleGenomeMutations> mutations_list,
                                          const std::set<ChromosomeId>& chromosome_ids,
                                          const double& coverage, const std::string& base_name="chr_",
                                          std::ostream& progress_bar_stream=std::cout,
                                          const bool quiet=false)
    {
        return operator()(sequencer, mutations_list, chromosome_ids, coverage,
                          {"normal_sample", {}}, 1.0, base_name, progress_bar_stream,
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
     * @param sequencer is the sequencer
     * @param mutations_list is a list of sample mutations
     * @param coverage is the aimed coverage
     * @param normal_sample is a normal sample
     * @param purity is ratio between the number of sampled tumour cells, which are
     *              represented in the mutation list, and that of the overall sampled
     *              cells which contains normal cells too
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleSetStatistics operator()(SEQUENCER& sequencer,
                                   std::list<CLONES::Mutations::SampleGenomeMutations> mutations_list,
                                   const double& coverage,
                                   const CLONES::Mutations::SampleGenomeMutations& normal_sample,
                                   const double purity,
                                   const std::string& base_name="chr_",
                                   std::ostream& progress_bar_stream=std::cout,
                                   const bool quiet=false)
    {
        const auto chromosome_ids = get_genome_chromosome_ids(mutations_list);

        return operator()(sequencer, mutations_list, chromosome_ids, coverage,
                          normal_sample, purity, base_name, progress_bar_stream,
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
     * @param sequencer is the sequencer
     * @param mutations_list is a list of sample mutations
     * @param coverage is the aimed coverage
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleSetStatistics operator()(SEQUENCER& sequencer,
                                   std::list<CLONES::Mutations::SampleGenomeMutations> mutations_list,
                                   const double& coverage,
                                   const std::string& base_name="chr_",
                                   std::ostream& progress_bar_stream=std::cout,
                                   const bool quiet=false)
    {
        const auto chromosome_ids = get_genome_chromosome_ids(mutations_list);

        return operator()(sequencer, mutations_list, chromosome_ids, coverage,
                          {"normal_sample", {}}, 1.0, base_name,
                          progress_bar_stream, quiet);
    }

    /**
     * @brief Generate simulated reads for a sample genome mutations
     *
     * This method generates simulated reads for a sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param sequencer is the sequencer
     * @param mutations is a sample genome mutations
     * @param coverage is the aimed coverage
     * @param normal_sample is a normal sample
     * @param purity is ratio between the number of sampled tumour cells, which are
     *              represented in the mutation list, and that of the overall sampled
     *              cells which contains normal cells too
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleStatistics operator()(SEQUENCER& sequencer,
                                const CLONES::Mutations::SampleGenomeMutations& mutations,
                                const double& coverage,
                                const CLONES::Mutations::SampleGenomeMutations& normal_sample,
                                const double purity,
                                const std::string& base_name="chr_",
                                std::ostream& progress_bar_stream=std::cout,
                                const bool quiet=false)
    {
        using namespace CLONES::Mutations;
        std::list<SampleGenomeMutations> mutations_list{mutations};

        auto statistics = operator()(sequencer, mutations_list, coverage, normal_sample,
                                     purity, base_name, progress_bar_stream, quiet);

        return statistics[mutations.name];
    }


    /**
     * @brief Generate simulated reads for a sample genome mutations
     *
     * This method generates simulated reads for a sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param sequencer is the sequencer
     * @param mutations is a sample genome mutations
     * @param coverage is the aimed coverage
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    inline SampleStatistics operator()(SEQUENCER& sequencer,
                                       const CLONES::Mutations::SampleGenomeMutations& mutations,
                                       const double& coverage,
                                       const std::string& base_name="chr_",
                                       std::ostream& progress_bar_stream=std::cout,
                                       const bool quiet=false)
    {
        return operator()(sequencer, mutations, coverage, {"normal_sample", {}},
                          1.0, base_name, progress_bar_stream, quiet);
    }

    /**
     * @brief Generate simulated reads for a sample genome mutations
     *
     * This method generates simulated reads for a sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param sequencer is the sequencer
     * @param mutations is a sample genome mutations
     * @param chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param coverage is the aimed coverage
     * @param normal_sample is a normal sample
     * @param purity is ratio between the number of sampled tumour cells, which are
     *              represented in the mutation list, and that of the overall sampled
     *              cells which contains normal cells too
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleStatistics operator()(SEQUENCER& sequencer,
                                const CLONES::Mutations::SampleGenomeMutations& mutations,
                                const std::set<ChromosomeId>& chromosome_ids,
                                const double& coverage,
                                const CLONES::Mutations::SampleGenomeMutations& normal_sample,
                                const double purity,
                                const std::string& base_name="chr_",
                                std::ostream& progress_bar_stream=std::cout,
                                const bool quiet=false)
    {
        using namespace CLONES::Mutations;
        std::list<SampleGenomeMutations> mutations_list{mutations};

        auto statistics = operator()(sequencer, mutations_list, chromosome_ids, coverage,
                                     normal_sample, purity, base_name, progress_bar_stream,
                                     quiet);

        return statistics[mutations.name];
    }

    /**
     * @brief Generate simulated reads for a sample genome mutations
     *
     * This method generates simulated reads for a sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param sequencer is the sequencer
     * @param mutations is a sample genome mutations
     * @param chromosome_ids is the set of chromosome identifiers whose reads will be
     *      simulated
     * @param coverage is the aimed coverage
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<CLONES::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    inline SampleStatistics operator()(SEQUENCER& sequencer,
                                       const CLONES::Mutations::SampleGenomeMutations& mutations,
                                       const std::set<ChromosomeId>& chromosome_ids,
                                       const double& coverage,
                                       const std::string& base_name="chr_",
                                       std::ostream& progress_bar_stream=std::cout,
                                       const bool quiet=false)
    {
        return operator()(sequencer, mutations, chromosome_ids, coverage,
                          {"normal_sample", {}}, 1.0, base_name, progress_bar_stream,
                          quiet);
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
     * @param enable is a flag to enable(true)/disable(false) SAM file writing
     */
    inline void enable_SAM_writing(const bool& enable)
    {
        write_SAM = enable;
    }
};

}   // SequencingSimulations

}   // Mutations

}   // CLONES

#endif // __CLONES_READ_SIMULATOR__
