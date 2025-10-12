/**
 * @file sbs_context_index.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements SBS context index and its builder
 * @version 1.0
 * @date 2025-10-12
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

#ifndef __RACES_SBS_CONTEXT_INDEX__
#define __RACES_SBS_CONTEXT_INDEX__

#include <filesystem>
#include <limits>
#include <ranges>
#include <concepts>

#include "index.hpp"
#include "sbs_context.hpp"
#include "genomic_position.hpp"
#include "genomic_sequence.hpp"
#include "fasta_utils.hpp"      // IO::FASTA::is_chromosome_header
#include "genomic_region.hpp"   // Mutations::GenomicRegion
#include "basic_IO.hpp"         // IO::get_stream_size

#include "utils.hpp"

#include "progress_bar.hpp"

namespace RACES
{

namespace Archive
{

/**
 * @brief A `partition` specialization for `SBSContext`
 *
 */
template<>
struct partition<Mutations::SBSContext>
{
    /**
     * @brief Return a list of the values in a class
     *
     * @param value is a representant of the class that is
     *      aimed
     * @return a list of the values in the class including
     *      `value`
     */
    static inline std::list<Mutations::SBSContext>
    get_class_of(const Mutations::SBSContext& context)
    {
        return std::list<Mutations::SBSContext>{context,
                                                context.get_reverse_complement()};
    }
};

}   // Archive

namespace Mutations
{

/**
 * @brief Indices for SBS context positions
 *
 * This class represents indices for SBS contexts. Its objects maps
 * any SBS context in a bucket containing the positions of the context
 * in the genome. The index buckets are shuffled to randomize the access
 * order.
 *
 * @tparam RANDOM_GENERATOR is a random number generator type
 */
template<class RANDOM_GENERATOR = std::mt19937_64>
class SBSContextIndex : public Archive::IndexReader<SBSContext, GenomicPosition, RANDOM_GENERATOR>
{
    /**
     * @brief Get the number of contexts
     *
     * @return the number of contexts
     */
    inline static constexpr size_t possible_context_code()
    {
        return std::numeric_limits<SBSContext::CodeType>::max() + 1;
    }

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
     * @brief A type to account the number of skipped contexts
     *
     * The indices can be sampled during their construction by providing
     * the number of contexts to skip before sampling one context. This
     * type is meant to record the number of contexts skipped during the
     * index construction.
     */
    using SkippedContexts = std::array<size_t, possible_context_code()>;

    std::map<ChromosomeId, GenomicRegion::Length> chr_lengths;  //!< The lengths of the genome chromosomes

    /**
     * @brief Initialize a skipped contexts array
     *
     * The skipped contexts array counts how many time a context
     * has not been inserted into the index since the last
     * insertion. This method creates a skipped contexts array
     * whose values are zeros.
     *
     * @return a skipped contexts array whose values are zeros
     */
    static SkippedContexts init_skipped_contexts()
    {
        SkippedContexts skipped_contexts;

        for (auto& value: skipped_contexts) {
            value = 0;
        }

        return skipped_contexts;
    }

    /**
     * @brief Update skipped context
     *
     * @param[in,out] skipped_contexts is an array counting how many time a context has not been
     *          inserted into the index since the last insertion
     * @param[in] context_code is the current mutational context code
     * @param[in] sampling_delta is the number of contexts to be found in order to record a context
     *          in the index
     * @return `true` if and only if either no context corresponding to the current `c_automaton`'s
     *          state has been found or the context has not been inserted  into the index
     *          `sampling_rate` times.
     */
    static bool update_skipped_contexts(SkippedContexts& skipped_contexts,
                                        const SBSContext::CodeType& context_code,
                                        const uint8_t sampling_delta)
    {
        if ((++(skipped_contexts[context_code]))==sampling_delta) {
            skipped_contexts[context_code] = 0;

            return true;
        }
        return false;
    }

    /**
     * @brief Skip all the remaining nucleotides of the current sequence in a FASTA stream
     *
     * This method discharge the remaining part of the current sequence, but it does consider
     * it as part of the genome.
     *
     * @param[in,out] fasta_stream is the input FASTA stream
     */
    static void skip_to_next_seq(std::ifstream& fasta_stream)
    {
        char in_char{'N'};
        while (in_char != '>' && !fasta_stream.eof())
        {
            fasta_stream.get(in_char);
        }

        if (!fasta_stream.eof()) {
            fasta_stream.unget();
        }
    }

    /**
     * @brief Find the mutational contexts in parts of a FASTA sequence and save their positions in the map
     *
     * This method finds the mutational contexts that lay the chromosome read from a FASTA stream and
     * outside a specified set of genomic regions.
     *
     * @param[in,out] index_builder is the index builder
     * @param[in,out] fasta_stream is the FASTA stream pointing at the first nucleotide of the considered sequence
     * @param[in] streamsize is the size of the FASTA stream in bytes
     * @param[in] chr_id is the identifier of the chromosome in the stream
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in,out] skipped_contexts is an array counting how many time a context has not been
     *          inserted into the index since the last insertion
     * @param[in] sampling_delta is the number of contexts to be found in order to record a context
     *          in the index
     * @param[in,out] progress_bar is the progress bar
     */
    static GenomicRegion::Length
    build_index_in_seq(Archive::IndexBuilder<SBSContext, GenomicPosition>& index_builder,
                       std::ifstream& fasta_stream, const size_t& streamsize,
                       const ChromosomeId& chr_id, const std::set<GenomicRegion>& regions_to_avoid,
                       SkippedContexts& skipped_contexts, const uint8_t sampling_delta,
                       UI::ProgressBar& progress_bar)
    {
        progress_bar.set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()
                                                        /streamsize));

        ExtendedContextAutomaton c_automata;

        GenomicPosition g_pos{chr_id, 0};
        auto region_it = regions_to_avoid.lower_bound({g_pos, 1});

        char last_char{'N'};
        while (last_char != '>' && !fasta_stream.eof()) {
            fasta_stream.get(last_char);
            if (GenomicSequence::is_a_DNA_base(last_char) || last_char == 'N') {
                ++g_pos.position;

                if (region_it != regions_to_avoid.end()) {
                    if (region_it->ends_before(g_pos)) {
                        ++region_it;
                    }
                }

                if (region_it != regions_to_avoid.end()) {
                    if (region_it->contains(g_pos)) {
                        last_char = 'N';
                    }
                }

                c_automata.update_state(last_char);

                if (c_automata.read_a_context()) {
                    if (update_skipped_contexts(skipped_contexts, c_automata.get_state(),
                                                sampling_delta)) {
                        index_builder.insert(c_automata.get_context(),
                                             {chr_id, g_pos.position-2});
                    }
                }
            }

            // update progress bar once every 2^22-1 nucleotides
            if ((g_pos.position & 0x3FFFFF)==0) {
                progress_bar.set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()
                                                                 /streamsize));
            }
        }

        if (!fasta_stream.eof()) {
            fasta_stream.unget();
        }

        return g_pos.position;
    }

    /**
     * @brief Split a set of genomic regions by chromosome id
     *
     * @param[in] genomic_regions is the set of genomic region to be split
     * @return a map that associates a chromosome id to the the set of genomic regions
     *     laying in the corresponding chromosome
     */
    static std::map<ChromosomeId, std::set<GenomicRegion> > split_by_chromosome_id(const std::set<GenomicRegion>& genomic_regions)
    {
        std::map<ChromosomeId, std::set<GenomicRegion> > split;

        for (const auto& genomic_region: genomic_regions) {
            split[genomic_region.get_chromosome_id()].insert(genomic_region);
        }

        return split;
    }

    /**
     * @brief Get the SBS context index specific data filename
     *
     * @return the SBS context index specific data filename
     */
    inline static constexpr std::string get_SBS_context_data_filename()
    {
        return "SBS_context_index_data.bin";
    }
public:
    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_stream is the stream of the chromosomes to be indexed
     * @param regions_to_avoid is the set of the regions that should
     *         not be considered by the index
     * @param tmp_dir is the directory that will contains the
     *         temporary files
     * @param cache_size is the size of the cache used to build the sample
     * @param sampling_delta is the number of loci of a specific context
     *         to be found before inserting one locus in the index
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     *          not laying in `regions_to_avoid` sampled according
     *          `sampling_delta`
     */
    static SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          std::ifstream& genome_stream,
          const std::set<GenomicRegion>& regions_to_avoid,
          const std::filesystem::path& tmp_dir,
          const size_t cache_size,
          const uint8_t sampling_delta,
          UI::ProgressBar& progress_bar)
    {
        if (!genome_stream.good()) {
            throw std::runtime_error("The genome FASTA stream is not readable");
        }

        auto regions_to_avoid_by_chr = split_by_chromosome_id(regions_to_avoid);

        auto streamsize = static_cast<size_t>(RACES::IO::get_stream_size(genome_stream));

        skip_to_next_seq(genome_stream);

        auto skipped_contexts = init_skipped_contexts();

        std::map<ChromosomeId, GenomicRegion::Length> chr_lengths;
        Archive::IndexBuilder<SBSContext, GenomicPosition> builder{index_path, cache_size};
        while (!genome_stream.eof()) {
            std::string sequence_title;
            getline(genome_stream, sequence_title);

            ChromosomeId chr_id;

            // if the sequence is a chromosome
            if (IO::FASTA::is_chromosome_header(sequence_title, chr_id)) {

                progress_bar.set_message("Processing chr. "
                                         + GenomicPosition::chrtos(chr_id));

                auto chr_length = build_index_in_seq(builder, genome_stream, streamsize,
                                                     chr_id, regions_to_avoid_by_chr[chr_id],
                                                     skipped_contexts, sampling_delta,
                                                     progress_bar);

                chr_lengths[chr_id] = chr_length;
            } else {
                progress_bar.set_progress(static_cast<uint8_t>(100*genome_stream.tellg()
                                            /streamsize), "Skipping a sequence");
                skip_to_next_seq(genome_stream);
            }
        }

        genome_stream.close();

        progress_bar.set_progress(100, "Index initialised");
        progress_bar.init_new();

        builder.shuffle(random_generator, tmp_dir, progress_bar);

        builder.save_map_on_disk();

        {
            Archive::Binary::Out archive(index_path / get_SBS_context_data_filename());

            archive & chr_lengths;
        }

        return {index_path, cache_size};
    }

    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_fasta is the path of the FASTA file containing the
     *         sequences of the chromosomes to be indexed
     * @param regions_to_avoid is the set of the regions that should not
     *         be considered by the index
     * @param tmp_dir is the directory that will contains the
     *         temporary files
     * @param cache_size is the size of the cache used to build the sample
     * @param sampling_delta is the number of loci of a specific context
     *         to be found before inserting one locus in the index
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     *          not laying in `regions_to_avoid`
     */
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          const std::filesystem::path& genome_fasta,
          const std::set<GenomicRegion>& regions_to_avoid,
          const std::filesystem::path& tmp_dir,
          const size_t cache_size,
          const uint8_t sampling_delta,
          UI::ProgressBar& progress_bar)
    {
        std::ifstream genome_fasta_stream(genome_fasta);

        if (!genome_fasta_stream.good()) {
            throw std::runtime_error("\"" + to_string(genome_fasta)
                                     + "\" does not exist");
        }

        return build(random_generator, index_path, genome_fasta_stream, regions_to_avoid,
                     tmp_dir, cache_size, sampling_delta, progress_bar);
    }

    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_fasta is the path of the FASTA file containing the
     *         sequences of the chromosomes to be indexed
     * @param regions_to_avoid is the set of the regions that should not
     *         be considered by the index
     * @param cache_size is the size of the cache used to build the sample
     * @param sampling_delta is the number of loci of a specific context
     *         to be found before inserting one locus in the index
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     *          not laying in `regions_to_avoid`
     */
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          const std::filesystem::path& genome_fasta,
          const std::set<GenomicRegion>& regions_to_avoid,
          const size_t cache_size,
          const uint8_t sampling_delta,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta, regions_to_avoid,
                     std::filesystem::temp_directory_path(), cache_size, sampling_delta,
                     progress_bar);
    }

    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_fasta is the path of the FASTA file containing the
     *         sequences of the chromosomes to be indexed
     * @param regions_to_avoid is the set of the regions that should not
     *         be considered by the index
     * @param cache_size is the size of the cache used to build the sample
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     *          not laying in `regions_to_avoid`
     */
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          const std::filesystem::path& genome_fasta,
          const std::set<GenomicRegion>& regions_to_avoid,
          const size_t cache_size,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta, regions_to_avoid,
                     cache_size, 1, progress_bar);
    }

    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_fasta is the path of the FASTA file containing the
     *         sequences of the chromosomes to be indexed
     * @param regions_to_avoid is the set of the regions that should not
     *         be considered by the index
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     *          not laying in `regions_to_avoid`
     */
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          const std::filesystem::path& genome_fasta,
          const std::set<GenomicRegion>& regions_to_avoid,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta, regions_to_avoid,
                     default_cache_size(), progress_bar);
    }

    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_fasta is the path of the FASTA file containing the
     *         sequences of the chromosomes to be indexed
     * @param tmp_dir is the directory that will contains the
     *         temporary files
     * @param cache_size is the size of the cache used to build the sample
     * @param sampling_delta is the number of loci of a specific context
     *         to be found before inserting one locus in the index
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     */
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          const std::filesystem::path& genome_fasta,
          const std::filesystem::path& tmp_dir,
          const size_t cache_size,
          const uint8_t sampling_delta,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta, {},
                     tmp_dir, cache_size, sampling_delta, progress_bar);
    }

    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_fasta is the path of the FASTA file containing the
     *         sequences of the chromosomes to be indexed
     * @param tmp_dir is the directory that will contains the
     *         temporary files
     * @param cache_size is the size of the cache used to build the sample
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     */
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          const std::filesystem::path& genome_fasta,
          const std::filesystem::path& tmp_dir,
          const size_t cache_size,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta,
                     tmp_dir, cache_size, 1, progress_bar);
    }

    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_fasta is the path of the FASTA file containing the
     *         sequences of the chromosomes to be indexed
     * @param cache_size is the size of the cache used to build the sample
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     */
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          const std::filesystem::path& genome_fasta,
          const size_t cache_size,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta,
                     std::filesystem::temp_directory_path(), cache_size,
                     progress_bar);
    }

    /**
     * @brief Buid a SBS context index
     *
     * This method builds a SBS context index.
     *
     * @param random_generator is a random generator
     * @param index_path is the path to the directory in which
     *         the index is stored
     * @param genome_fasta is the path of the FASTA file containing the
     *         sequences of the chromosomes to be indexed
     * @param progress_bar is a progress bar
     * @return An SBS context index containing the loci in `genome_stream`
     */
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(RANDOM_GENERATOR& random_generator,
          const std::filesystem::path index_path,
          const std::filesystem::path& genome_fasta,
          UI::ProgressBar& progress_bar)
    {
        return build(random_generator, index_path, genome_fasta,
                     default_cache_size(), progress_bar);
    }


    /**
     * @brief Buid a SBS context index
     *
     * This method works as a wrapper. It is exclusively called when the
     * last parameter is not a progress bar and add a quiet progress bar
     * as the new last parameter.
     *
     * @tparam ARGS is the pack of the parameter type
     * @param args are the parameters
     * @return An SBS context index according `args`
     */
    template <typename... ARGS>
    requires last_is_not<UI::ProgressBar, ARGS...>
    static inline SBSContextIndex<RANDOM_GENERATOR>
    build(ARGS... args)
    {
        UI::ProgressBar progress_bar;

        return build(args..., progress_bar);
    }

    /**
     * @brief The empty constructor
     */
    SBSContextIndex():
        Archive::IndexReader<SBSContext, GenomicPosition, RANDOM_GENERATOR>{}
    {}

    /**
     * @brief A constructor
     *
     * This constructor loads an already built SBS context index from
     * the directory in which it is stored.
     *
     * @param index_path is the path to the SBS context index directory
     * @param cache_size is the index read cache size
     */
    SBSContextIndex(const std::filesystem::path index_path,
                    const size_t cache_size = default_cache_size()):
        Archive::IndexReader<SBSContext, GenomicPosition, RANDOM_GENERATOR>{index_path, cache_size}
    {
        Archive::Binary::In archive(index_path / get_SBS_context_data_filename());

        archive & chr_lengths;
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
     * @brief Get the regions of the chromosomes
     *
     * @return the regions of the chromosomes
     */
    std::list<GenomicRegion> get_chromosome_regions() const
    {
        std::list<GenomicRegion> gr_list;

        for (const auto& [chr_id, chr_length] : chr_lengths) {
            gr_list.emplace_back(chr_id, chr_length);
        }

        return gr_list;
    }

    /**
     * @brief Get the chromosome identifiers of the reference sequence
     *
     * @return a vector containing the chromosome identifiers of the
     *      reference sequence
     */
    std::vector<ChromosomeId> get_chromosome_ids() const
    {
        std::vector<ChromosomeId> chromosome_ids(chr_lengths.size());

        auto it = chromosome_ids.begin();
        for (const auto& chr_id : std::views::keys(chr_lengths)) {
            *it = chr_id;

            ++it;
        }

        return chromosome_ids;
    }
};

}

}
#endif // __RACES_SBS_CONTEXT_INDEX__