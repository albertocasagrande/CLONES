/**
 * @file id_context_index.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines ID context index
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

#ifndef __CLONES_ID_CONTEXT_INDEX__
#define __CLONES_ID_CONTEXT_INDEX__

#include <vector>
#include <random>
#include <utility>

#include <algorithm>    // std::fill

#include "genomic_position.hpp"
#include "genomic_region.hpp"
#include "utils.hpp"
#include "progress_bar.hpp"
#include "archive.hpp"

#include "index.hpp"

#include "id_context.hpp"
#include "fasta_chr_reader.hpp"

namespace CLONES
{

/**
 * @brief A `partition` specialization for `IDContext`
 *
 */
template<>
struct partition<Mutations::IDContext>
{
    /**
     * @brief Return a list of the values in a class
     *
     * @param[in] value is a representant of the class that is
     *      aimed
     * @return a list of the values in the class including
     *      `value`
     */
    static std::list<Mutations::IDContext>
    get_class_of(const Mutations::IDContext& context);
};

namespace Mutations
{

/**
 * @brief Repetition reference
 *
 * This class represents repetition references. The objects of this class maintain
 * the reference genomic position and the length of the repeated unit.
 */
struct RepetitionReference
{
    using RepetitionType = IDContext::FirstLevelType;

    GenomicPosition position;  //!< The genomic position of the repeated sequence
    RepetitionType unit_size;    //!< The unit size

    /**
     * @brief The empty constructor
     */
    RepetitionReference();

    /**
     * @brief A constructor
     *
     * @param chr_id is the identifier of the chromosome containing
     *      the repeated sequence
     * @param begin is the position of the repeated sequence first base
     *      in the chromosome
     * @param unit_size is the unit size
     */
    RepetitionReference(const ChromosomeId& chr_id, const ChrPosition begin,
                        const RepetitionType unit_size);

    /**
     * @brief Save a repetition in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & position
                & unit_size;
    }

    /**
     * @brief Load a repetition from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded repetition
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static RepetitionReference load(ARCHIVE& archive)
    {
        RepetitionReference repetition;

        archive & repetition.position
                & repetition.unit_size;

        return repetition;
    }

    CHECK_CONSTANT_SPACE_ON_DISK(position, unit_size)
};

}   // Mutations

}   //  CLONES


namespace std
{
    std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::RepetitionReference& repetition);
}

namespace CLONES
{

namespace Mutations
{

template<class RANDOM_GENERATOR = std::mt19937_64>
class IDContextIndex;

/**
 * @brief Builders for indel context indices
 *
 * The objects of this class build indel context indices.
 */
class IDContextIndexBuilder : public Archive::IndexBuilder<IDContext, RepetitionReference>
{
    using RepetitionType = IDContext::SecondLevelType;

    RepetitionType max_unit_size;   //!< The maximum considered size of the repetition unit

    /**
     * @brief Get the default cache size
     *
     * @return the default cache size
     */
    static constexpr size_t default_cache_size()
    {
        return 1000;
    }

    /**
     * @brief Initialize the suffix array vector
     *
     * @param s is the sequence whose suffix array is aimed
     * @param suffix_array is the vector that will contain the suffix array
     * @param classes is the vector that labels each sequence position with the
     *      class of its first nucleotide
     * @return the number of the different nucleotides in the sequence
     */
    static size_t init_suffix_array(const char* s,
                                    std::vector<ChrPosition>& suffix_array,
                                    std::vector<ChrPosition>& classes);
    /**
     * @brief Upgrade the (h-1)-suffix array to a (h)-suffix array
     *
     * @param h is the length of the prefixes sorted in the (h)-suffix array
     * @param h_suffix_array is a (h-1)-suffix array
     * @param h_classes is the vector of the classes of the (h-1)-suffix array
     * @param num_of_classes is the number of classes of the (h-1)-suffix array
     * @param tmp_a is a temporary vector having the same size of `h_suffix_array`
     * @param tmp_b is a temporary vector having the same size of `h_suffix_array`
     */
    static void update_suffix_array(const size_t h,
                                    std::vector<ChrPosition>& h_suffix_array,
                                    std::vector<ChrPosition>& h_classes,
                                    size_t& num_of_classes,
                                    std::vector<ChrPosition>& tmp_a,
                                    std::vector<ChrPosition>& tmp_b);

    /**
     * @brief Add a repeated sequence to the index
     *
     * @param[in] chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param[in] seq is the considered genomic sequence
     * @param[in] begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param[in] unit_size is the size of the repetition unit
     * @param[in] r_begin is the position of the repeated sequence first base
     *      on the considered genomic sequence
     * @param[in] r_end is the position of the repeated sequence last base
     *      on the considered genomic sequence
     * @param[in] covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void add_repetition(const ChromosomeId& chr_id, const char* seq,
                        const ChrPosition& begin, const size_t unit_size,
                        const ChrPosition& r_begin, const ChrPosition& r_end,
                        std::vector<bool>& covered);

    //! @private
    void add_null_heteropolymer(const ChromosomeId& chr_id, const size_t unit_size,
                                const ChrPosition& begin, const ChrPosition& r_begin);

    //! @private
    void add_null_homopolymer(const size_t nucleotide_index, const char* seq,
                            const ChromosomeId& chr_id, const ChrPosition& begin,
                            const ChrPosition& r_begin);

    /**
     * @brief Collect the candidate repeated sequences whose unit size is in [h, 2*h-1]
     *
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param h is the order of the suffix array
     * @param h_suffix_array is the (h)-suffix array
     * @param classes is the class vector of the (h)-suffix array
     * @return the candidate repeated sequences. The candidates are encoded in a map
     *      from the position of the candidate first base in the considered sequence
     *      to a map from the candidate unit size to the position of the final base.
     *      The returned repeated sequence are not added to the index yet because
     *      some of them may be fully included in other candidates
     */
    static std::map<ChrPosition, std::map<size_t, ChrPosition>>
    collect_candidates(const ChrPosition& begin, const size_t& h,
                    std::vector<ChrPosition>& h_suffix_array,
                    std::vector<ChrPosition>& classes);

    /**
     * @brief Collect the repeated sequences whose unit size is in [h, 2*h-1]
     *
     * @param[in] chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param[in] seq is the considered genomic sequence
     * @param[in] begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param[in] h is the order of the suffix array
     * @param[in,out] h_suffix_array is the (h)-suffix array
     * @param[in,out] classes is the class vector of the (h)-suffix array
     * @param[in,out] covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void add_repetitions(const ChromosomeId& chr_id, const char* seq,
                        const ChrPosition& begin, const size_t& h,
                        std::vector<ChrPosition>& h_suffix_array,
                        std::vector<ChrPosition>& classes,
                        std::vector<bool>& covered);

    /**
     * @brief Collect micro-homologies
     *
     * @param[in] chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param[in] seq is the considered genomic sequence
     * @param[in] begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param[in] covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void add_microhomologies(const ChromosomeId& chr_id, const char* seq,
                            const ChrPosition& begin, std::vector<bool>& covered);

    /**
     * @brief Collect non-repeated sequence
     *
     * @param[in] chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param[in] seq is the considered genomic sequence
     * @param[in] begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param[in] covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void add_non_repeated_seq(const ChromosomeId& chr_id, const char* seq,
                            const ChrPosition& begin, std::vector<bool>& covered);

    /**
     * @brief Add the repeated sequences of a genomic sequence
     *
     * @param[in] chr_id is the identifier of the chromosome containing the
     *      genomic sequence
     * @param[in] seq is the considered genomic sequence
     * @param[in] begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param[in] length is the length of the considered sequence
     * @param[in,out] progress_bar is the progress bar
     */
    std::vector<bool>
    add_repetitions(const ChromosomeId& chr_id, const char* seq,
                    const ChrPosition begin, const size_t& length,
                    UI::ProgressBar& progress_bar);

    /**
     * @brief Collect the data from a genomic sequence
     *
     * @param[in,out] random_generator is a random generator
     * @param[in] sequence is a genomic sequence
     * @param[in] chr_id is the identifier of the chromosome containing the
     *      genomic sequence
     * @param[in] begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param[in] length is the length of the considered sequence
     * @param[in,out] progress_bar is the progress bar
     */
    void add_contexts_from(const ChromosomeId& chr_id, const std::string& sequence,
                           const ChrPosition begin, size_t length,
                           UI::ProgressBar& progress_bar);

    //!< @private
    static IDContext::FirstLevelType get_unit_size_code(const size_t& unit_size);

    //!< @private
    static IDContext::SecondLevelType get_num_of_repetitions_code(const size_t& num_of_repetitions);

    //!< @private
    static IDContext::SecondLevelType get_homology_size_code(const size_t& homology_size);

    /**
     * @brief Add a polymeric repetition to the index
     *
     * @param genomic_position is the genomic position of the repetition's
     *      initial position
     * @param num_of_repetitions is the number of repetition of the unit
     * @param unit is the pointer to the unit first base in the sequence
     * @param unit_size is the unit size
     */
    void add_polymer(const GenomicPosition& genomic_position,
                    const size_t& num_of_repetitions,
                    const char* unit, const size_t& unit_size);

    /**
     * @brief Add to the index all the ID context in a chromosome
     *
     * @param[in] chr_id is the identifier of a chromosome
     * @param[in] chr_sequence is the nucleotide sequence of a chromosome
     * @param[in,out] progress_bar is the progress bar
     */
    inline void add_contexts_from(const ChromosomeId& chr_id, const std::string& chr_sequence,
                                  UI::ProgressBar& progress_bar)
    {
        std::set<GenomicRegion> no_regions_to_avoid;

        add_contexts_from(chr_id, chr_sequence, no_regions_to_avoid, progress_bar);
    }

    /**
     * @brief Add to the index all the ID context in a chromosome
     *
     * @param[in] chr_id is the identifier of a chromosome
     * @param[in] chr_sequence is the nucleotide sequence of a chromosome
     * @param[in] regions_to_avoid is the set of regions in the chromosome to be avoided
     * @param[in,out] progress_bar is the progress bar
     */
    void add_contexts_from(const ChromosomeId& chr_id, const std::string& chr_sequence,
                           const std::set<GenomicRegion>& regions_to_avoid,
                           UI::ProgressBar& progress_bar);

    /**
     * @brief A constructor
     *
     * This constructor loads an already built repeated sequence index from
     * the directory in which it is stored.
     *
     * @param index_path is the path to the repeated sequence index directory
     * @param max_unit_size is the maximum unit size
     * @param cache_size is the index read cache size
     */
    IDContextIndexBuilder(const std::filesystem::path index_path,
                          const RepetitionType max_unit_size = 50,
                          const size_t cache_size = default_cache_size());
public:

    /**
     * @brief Build an ID context index
     *
     * This method works as a wrapper. It is exclusively called when the
     * last parameter is not a progress bar and add a quiet progress bar
     * as the new last parameter.
     *
     * @tparam RANDOM_GENERATOR is the random generator type
     * @tparam ARGS is the pack of the parameter type
     * @param[in,out] args are the parameters
     * @return An repeated sequence index according `args`
     */
    template <typename RANDOM_GENERATOR, typename... ARGS>
    requires last_is_not<UI::ProgressBar, ARGS...>
    static inline IDContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator, ARGS... args)
    {
        UI::ProgressBar progress_bar;

        return build(random_generator, args..., progress_bar);
    }

    /**
     * @brief Build an ID context index
     *
     * @tparam RANDOM_GENERATOR is the random generator type
     * @param[in,out] random_generator is a random generator
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] sampling_delta is the number of repeated sequences to be found in
     *          order to record a repeated sequence in the index
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`
     */
    template<typename RANDOM_GENERATOR>
    static inline IDContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator, const std::filesystem::path& genome_fasta,
          const RepetitionType max_unit_size, const uint8_t sampling_delta,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, genome_fasta, {}, max_unit_size, sampling_delta,
                     progress_bar);
    }

    /**
     * @brief Build an ID context index
     *
     * @tparam RANDOM_GENERATOR is the random generator type
     * @param[in,out] random_generator is a random generator
     * @param[in] index_path is the path to the directory in which
     *         the index is stored
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] tmp_dir is the directory that will contains the
     *         temporary files
     * @param[in] cache_size is the size of the cache used to build the sample
     * @param[in] sampling_delta is the number of repeated sequences to be found in
     *          order to record a repeated sequence in the index
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`,
     *          but that are located outside the regions in `regions_to_avoid`
     */
    template<typename RANDOM_GENERATOR>
    static IDContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator, const std::filesystem::path& index_path,
          const std::filesystem::path& genome_fasta,
          const std::set<GenomicRegion>& regions_to_avoid, const RepetitionType max_unit_size,
          const std::filesystem::path& tmp_dir,
          const size_t cache_size, const uint8_t sampling_delta,
          UI::ProgressBar& progress_bar)
    {

        (void)sampling_delta;

        using namespace IO::FASTA;

        Reader<ChromosomeData<Sequence>> chr_reader(genome_fasta);

        auto regions_to_avoid_by_chr = split_by_chromosome_id(regions_to_avoid);
        const auto streamsize = chr_reader.get_stream_size();

        ChromosomeData<Sequence> chr;

        std::map<ChromosomeId, GenomicRegion::Length> chr_lengths;
        IDContextIndexBuilder builder{index_path, max_unit_size, cache_size};
        while (chr_reader.read(chr, progress_bar)) {
            progress_bar.set_progress(static_cast<uint8_t>(100*chr_reader.get_position()/streamsize),
                                      "Processing chr. " + GenomicPosition::chrtos(chr.chr_id));

            builder.add_contexts_from(chr.chr_id, chr.nucleotides,
                                      regions_to_avoid_by_chr[chr.chr_id], progress_bar);
            chr_lengths.emplace(chr.chr_id, chr.nucleotides.size());
        }

        progress_bar.set_progress(100, "Index initialised");
        progress_bar.init_new();

        builder.shuffle(random_generator, tmp_dir, progress_bar);

        builder.save_map_on_disk();

        {
            Archive::Binary::Out archive(index_path / get_ID_context_data_filename());

            archive & chr_lengths
                    & max_unit_size;
        }

        return IDContextIndex<RANDOM_GENERATOR>{index_path, cache_size};
    }

    /**
     * @brief Get the ID context index specific data filename
     *
     * @return the ID context index specific data filename
     */
    inline static constexpr std::string get_ID_context_data_filename()
    {
        return "ID_context_index_data.bin";
    }
};

/**
 * @brief Indices for indel contexts
 *
 * This class represents indices for indel contexts. Its objects maps
 * any indel type in a bucket containing the positions of the corresponding
 * repeated sequence in the genome. The index buckets are shuffled to
 * randomize the access order.
 *
 * @tparam RANDOM_GENERATOR is a random number generator type
 */
template<class RANDOM_GENERATOR>
class IDContextIndex : public Archive::IndexReader<IDContext, RepetitionReference,
                                                   RANDOM_GENERATOR>
{
public:
    using RepetitionType = IDContext::SecondLevelType;

private:
    RepetitionType max_unit_size;   //!< The maximum considered size of the repetition unit

    std::map<ChromosomeId, GenomicRegion::Length> chr_lengths;  //!< The lengths of the genome chromosomes

    /**
     * @brief Get the default cache size
     *
     * @return the default cache size
     */
    static constexpr size_t default_cache_size()
    {
        return 1000;
    }
public:
    /**
     * @brief The empty constructor
     */
    IDContextIndex():
        Archive::IndexReader<IDContext, RepetitionReference, RANDOM_GENERATOR>{}
    {}

    /**
     * @brief A constructor
     *
     * This constructor loads an already built repeated sequence index from
     * the directory in which it is stored.
     *
     * @param index_path is the path to the repeated sequence index directory
     * @param cache_size is the index read cache size
     */
    IDContextIndex(const std::filesystem::path index_path,
                   const size_t cache_size = default_cache_size()):
        Archive::IndexReader<IDContext, RepetitionReference, RANDOM_GENERATOR>{index_path, cache_size}
    {
        Archive::Binary::In archive(index_path / get_ID_context_data_filename());

        archive & chr_lengths
                & max_unit_size;
    }

    /**
     * @brief Get the chromosome lengths
     *
     * @return a map associating each indexed chromosome to its length
     */
    inline const std::map<ChromosomeId, GenomicRegion::Length>& get_chromosome_lengths() const
    {
        return chr_lengths;
    }

    /**
     * @brief Get the maximal unit size
     *
     * @return a constant reference to the maximal unit size
     */
    inline const RepetitionType& get_max_unit_size() const
    {
        return max_unit_size;
    }

    /**
     * @brief Build an ID context index
     *
     * This method works as a wrapper. It is exclusively called when the
     * last parameter is not a progress bar and add a quiet progress bar
     * as the new last parameter.
     *
     * @tparam ARGS is the pack of the parameter type
     * @param[in,out] args are the parameters
     * @return An repeated sequence index according `args`
     */
    template <typename... ARGS>
    requires last_is_not<UI::ProgressBar, ARGS...>
    static inline IDContextIndex<RANDOM_GENERATOR>
    build(ARGS... args)
    {
        UI::ProgressBar progress_bar;

        return build(args..., progress_bar);
    }

    /**
     * @brief Build an ID context index
     *
     * @param[in,out] random_generator is a random generator
     * @param[in] index_path is the path to the directory in which
     *         the index is stored
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] sampling_delta is the number of repeated sequences to be found in
     *          order to record a repeated sequence in the index
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`
     */
    static inline IDContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator, const std::filesystem::path& index_path,
          const std::filesystem::path& genome_fasta,
          const RepetitionType max_unit_size, const uint8_t sampling_delta,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta, {}, max_unit_size,
                     sampling_delta, progress_bar);
    }

    /**
     * @brief Build an ID context index
     *
     * @param[in,out] random_generator is a random generator
     * @param[in] index_path is the path to the directory in which
     *         the index is stored
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] tmp_dir is the directory that will contains the
     *         temporary files
     * @param[in] cache_size is the size of the cache used to build the sample
     * @param[in] sampling_delta is the number of repeated sequences to be found in
     *          order to record a repeated sequence in the index
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`,
     *          but that are located outside the regions in `regions_to_avoid`
     */
    static IDContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator, const std::filesystem::path& index_path,
          const std::filesystem::path& genome_fasta,
          const std::set<GenomicRegion>& regions_to_avoid,
          const RepetitionType max_unit_size, const std::filesystem::path& tmp_dir,
          const size_t cache_size, const uint8_t sampling_delta,
          UI::ProgressBar& progress_bar)
    {
        return IDContextIndexBuilder::build(random_generator, index_path, genome_fasta,
                                            regions_to_avoid, max_unit_size, tmp_dir,
                                            cache_size, sampling_delta, progress_bar);
    }


    /**
     * @brief Build an ID context index
     *
     * @param[in,out] random_generator is a random generator
     * @param[in] index_path is the path to the directory in which
     *         the index is stored
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`
     */
    static inline IDContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator, const std::filesystem::path& index_path,
          const std::filesystem::path& genome_fasta, const RepetitionType max_unit_size,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta, {}, max_unit_size,
                     progress_bar);
    }

    /**
     * @brief Build an ID context index
     *
     * @param[in,out] random_generator is a random generator
     * @param[in] index_path is the path to the directory in which
     *         the index is stored
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] cache_size is the size of the cache used to build the sample
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`,
     *          but that are located outside the regions in `regions_to_avoid`
     */
    static inline IDContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator, const std::filesystem::path& index_path,
          const std::filesystem::path& genome_fasta,
          const std::set<GenomicRegion>& regions_to_avoid,
          const RepetitionType max_unit_size,
          const size_t cache_size, UI::ProgressBar& progress_bar)
    {
        return IDContextIndexBuilder::build(random_generator, index_path, genome_fasta,
                                            regions_to_avoid, max_unit_size,
                                            std::filesystem::temp_directory_path(),
                                            cache_size, 1, progress_bar);
    }

    /**
     * @brief Build an ID context index
     *
     * @param[in,out] random_generator is a random generator
     * @param[in] index_path is the path to the directory in which
     *         the index is stored
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] tmp_dir is the directory that will contains the
     *         temporary files
     * @param[in] cache_size is the size of the cache used to build the sample
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`,
     *          but that are located outside the regions in `regions_to_avoid`
     */
    static inline IDContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator, const std::filesystem::path& index_path,
          const std::filesystem::path& genome_fasta,
          const std::set<GenomicRegion>& regions_to_avoid,
          const RepetitionType max_unit_size, const std::filesystem::path& tmp_dir,
          const size_t cache_size, UI::ProgressBar& progress_bar)
    {
        return IDContextIndexBuilder::build(random_generator, index_path, genome_fasta,
                                            regions_to_avoid, max_unit_size, tmp_dir,
                                            cache_size, 1, progress_bar);
    }

    /**
     * @brief Get the ID context index specific data filename
     *
     * @return the ID context index specific data filename
     */
    inline static constexpr std::string get_ID_context_data_filename()
    {
        return IDContextIndexBuilder::get_ID_context_data_filename();
    }
};

}   // Mutations

}   // CLONES

#endif // __CLONES_RS_INDEX__