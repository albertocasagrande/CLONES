/**
 * @file germline.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines the functions to generate and load germline mutations
 * @version 1.1
 * @date 2025-10-11
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

#ifndef __RACES_GERMLINE__
#define __RACES_GERMLINE__

#include <map>
#include <string>
#include <filesystem>
#include <random>

#include "genome_mutations.hpp"
#include "fasta_chr_reader.hpp"

#include "driver_storage.hpp"

#include "progress_bar.hpp"

namespace RACES
{

namespace Mutations
{

/**
 * @brief A class that generate germline mutations
 *
 */
class GermlineMutations
{
    using Chromosome = RACES::IO::FASTA::ChromosomeData<RACES::IO::FASTA::Sequence>;

    const std::vector<char> bases;  //!< the base vector
    std::map<char, uint8_t> bases_pos;  //!< the position of the bases in the base vector
    std::uniform_int_distribution<uint8_t> base_dist; //!< the base distribution probability

    std::mt19937_64 rand_gen;   //!< the random generator

    const size_t genome_size;   //!< the genome size
    const size_t expected_mutations;    //!< the expected number of geminal mutations to place

    size_t not_processed_size;  //!< the length of the regions in which mutations has not been placed
    size_t mutations_not_placed;  //!< the number of mutations not placed yet

    GermlineMutations(const size_t& genome_size, const double& mutations_per_ref_kbase,
                      const int& seed=0);

    /**
     * @brief Get the percentage of progresses
     *
     * @return the percentage of progresses
     */
    inline unsigned int get_process_percentage() const
    {
        return 100-(100*(not_processed_size+mutations_not_placed))/(genome_size+expected_mutations);
    }

    /**
     * @brief Build a candidate germline SNV
     *
     * @param chromosome is the chromosome
     * @param position is the candidate germline SNV position
     * @return a candidate germline SNV on `sequence` centered on `position`
     */
    SID get_candidate_germline_SNV(const Chromosome& chromosome, const ChrPosition& position);

    /**
     * @brief Get the number of mutations in a chromosome
     *
     * @param chr_mutations are the chromosome mutations
     * @return the number of mutations to be placed in a chromosome having
     *          length `chr_size`
     */
    size_t get_mutations_in(const ChromosomeMutations& chr_mutations);

public:
    /**
     * @brief Generate germline mutations
     *
     * @param reference_fasta_filename is the name of the reference genome fasta file
     * @param chromosome_regions is the vector of chromosome regions
     * @param alleles_per_chromosome is the number of alleles in wild-type cells
     * @param driver_storage is the storage of driver mutations
     * @param mutations_per_ref_kbase is the number of germline mutations
     *              per kilobase of reference genome
     * @param seed is the seed of the random generator
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is Boolean flag to disable/enable a progress bar
     * @return the germline mutations
     */
    static GenomeMutations
    generate(const std::filesystem::path& reference_fasta_filename,
             const std::list<GenomicRegion>& chromosome_regions,
             const DriverStorage& driver_storage,
             const std::map<ChromosomeId, size_t>& alleles_per_chromosome,
             const double& mutations_per_ref_kbase,
             const int& seed=0, std::ostream& progress_bar_stream=std::cout,
             const bool quiet=true);

    /**
     * @brief Load germline mutations
     *
     * This method loads germline mutations from a _germline data file_.
     * These files are "\t"-separated CSV files having two columns: "chr" and
     * "file". Each row of these files corresponds to a chromosome and reports,
     * in order, the chromosome name and the relative path of a VCF file in
     * IGSR format (see https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40)
     * that describes the mutations in the chrosomome.
     *
     * @param germline_data_file is a germline data file
     * @param alleles_per_chromosome is the number of alleles in wild-type cells
     * @param subject is the subject whose germline mutations must be loaded
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable progress messages
     * @return The genome mutations of `subject` if they are contained in the
     *       VCF files
     */
    static GenomeMutations load(const std::filesystem::path& germline_data_file,
                                const std::map<ChromosomeId, size_t>& alleles_per_chromosome,
                                const std::string& subject,
                                std::ostream& progress_bar_stream=std::cout,
                                const bool quiet=false);
};

}   // Mutations

}   // RACES

#endif // __RACES_GERMLINE__