/**
 * @file context_index.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to build a context index
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

#ifndef __CLONES_CONTEXT_INDEX__
#define __CLONES_CONTEXT_INDEX__

#include <map>
#include <set>
#include <memory>
#include <vector>
#include <fstream>
#include <type_traits>

#include "archive.hpp"
#include "utils.hpp"
#include "genomic_sequence.hpp"
#include "fasta_utils.hpp"      // IO::FASTA::is_chromosome_header
#include "genomic_region.hpp"   // Mutations::GenomicRegion
#include "basic_IO.hpp"         // IO::get_stream_size
#include "sbs_context.hpp"          // Mutations::ExtendedContextAutomaton

#include "progress_bar.hpp"


namespace CLONES
{

namespace Mutations
{

/**
 * @brief Genome-wise mutational context indices
 *
 * @tparam GENOME_WIDE_POSITION is the type used to represent genome-wise position
 */
template<typename GENOME_WIDE_POSITION = uint32_t,
         std::enable_if_t<std::is_integral_v<GENOME_WIDE_POSITION>
                          && std::is_unsigned_v<GENOME_WIDE_POSITION>, bool> = true>
struct ContextIndex
{
    /**
     * @brief The genome-wide position type
     */
    using GenomeWidePosition = GENOME_WIDE_POSITION;
protected:
    /**
     * @brief Maps associating mutational contexts to their position in the genome
     */
    using ContextPositionMap = std::map<SBSContext, std::vector<GenomeWidePosition> >;

    std::shared_ptr<ContextPositionMap> context2pos;             //!< the context-genomic positions map
    std::map<GenomeWidePosition, ChromosomeId> abs_pos2chr;      //!< the absolute genomic position-chromosome id map

    GenomeWidePosition genome_size;  //!< the genome size

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
    static std::array<size_t, 125> init_skipped_contexts()
    {
        std::array<size_t, 125> skipped_contexts;

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
     * @param[in] sampling_rate is the number of contexts to be found in order to record a context
     *          in the index
     * @return `true` if and only if either no context corresponding to the current `c_automaton`'s
     *          state has been found or the context has not been inserted  into the index
     *          `sampling_rate` times.
     */
    static bool update_skipped_contexts(std::array<size_t, 125>& skipped_contexts,
                                        const SBSContext::CodeType& context_code,
                                        const size_t& sampling_rate)
    {
        if ((++(skipped_contexts[context_code]))==sampling_rate) {
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
    void skip_to_next_seq(std::ifstream& fasta_stream)
    {
        char in_char{'N'};
        while (in_char != '>' && !fasta_stream.eof())
        {
            fasta_stream.get(in_char);

            if (ExtendedContextAutomaton::is_a_nucleotide(in_char)) {
                ++genome_size;
            }
        }

        if (!fasta_stream.eof()) {
            fasta_stream.unget();
        }
    }

    /**
     * @brief Skip all the remaining nucleotides of the current sequence in a FASTA stream
     *
     * This method discharge the current sequence and does not consider it as part of the
     * genome.
     *
     * @param[in,out] fasta_stream is the input FASTA stream
     */
    static void discharge_seq(std::ifstream& fasta_stream)
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
     * @param[in,out] fasta_stream is the FASTA stream pointing at the first nucleotide of the considered sequence
     * @param[in] streamsize is the size of the FASTA stream in bytes
     * @param[in] chr_id is the identifier of the chromosome in the stream
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in,out] skipped_contexts is an array counting how many time a context has not been
     *          inserted into the index since the last insertion
     * @param[in] sampling_rate is the number of contexts to be found in order to record a context
     *          in the index
     * @param[in,out] progress_bar is the progress bar
     */
    void build_index_in_seq(std::ifstream& fasta_stream, const std::streampos& streamsize,
                            const ChromosomeId& chr_id,
                            const std::set<GenomicRegion>& regions_to_avoid,
                            std::array<size_t, 125>& skipped_contexts,
                            const size_t& sampling_rate, UI::ProgressBar* progress_bar)
    {
        if (progress_bar != nullptr) {
            progress_bar->set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()/streamsize));
        }

        ExtendedContextAutomaton c_automata;

        GenomicPosition pos{chr_id, 0};
        auto region_it = regions_to_avoid.lower_bound({pos, 1});

        char last_char{'N'};
        while (last_char != '>' && !fasta_stream.eof()) {
            fasta_stream.get(last_char);
            if (GenomicSequence::is_a_DNA_base(last_char) || last_char == 'N') {
                ++pos.position;

                if (region_it != regions_to_avoid.end() && region_it->ends_before(pos)) {
                    ++region_it;
                }

                if (region_it != regions_to_avoid.end() && region_it->contains(pos)) {
                    last_char = 'N';
                }

                c_automata.update_state(last_char);

                if (c_automata.read_a_context()) {
                    if (update_skipped_contexts(skipped_contexts, c_automata.get_state(),
                                                sampling_rate)) {
                        (*context2pos)[c_automata.get_context()].emplace_back(genome_size);
                    }
                }
                ++genome_size;
            }

            // update progress bar once every 2^22-1 nucleotides
            if ((genome_size&0x3FFFFF)==0 && progress_bar != nullptr) {
                progress_bar->set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()/streamsize));
            }
        }

        if (!fasta_stream.eof()) {
            fasta_stream.unget();
        }
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
     * @brief Initialize `context2pos` map
     */
    void initialize_context2pos()
    {
        context2pos = std::make_shared<ContextPositionMap>();

        for (const auto& nucleotide1 : GenomicSequence::DNA_bases) {
            for (const auto& nucleotide2 : GenomicSequence::DNA_bases) {
                for (const auto& nucleotide3 : GenomicSequence::DNA_bases) {
                    SBSContext context(std::string{nucleotide1, nucleotide2, nucleotide3});
                    (*context2pos)[context] = std::vector<GENOME_WIDE_POSITION>();
                }
            }
        }
    }

    /**
     * @brief Reset the object and find the mutational contexts in a FASTA stream
     *
     * This method resets the current object, finds the mutational contexts in the chromosome sequences
     * contained in a FASTA stream, and saves their absolute positions in the object member
     * `ContextIndex::context2pos`. The chromosome sequences are identified by their names according
     * to the method template parameter.
     * At the same time, this method stores the absolute positions of each chromosome in the member
     * `ContextIndex::abs_pos2chr`.
     *
     * @param[in,out] fasta_stream is the genome FASTA stream
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in] sampling_rate is the number of contexts to be found in order to record a context
     *          in the index
     * @param[in,out] progress_bar is the progress bar
     */
    void reset_with(std::ifstream& fasta_stream, const std::set<GenomicRegion>& regions_to_avoid,
                    const size_t& sampling_rate, UI::ProgressBar* progress_bar)
    {
        if (!fasta_stream.good()) {
            throw std::runtime_error("the stream is not readable");
        }

        if (sampling_rate==0) {
            throw std::domain_error("The sampling rate must be positive");
        }

        auto streamsize = CLONES::IO::get_stream_size(fasta_stream);

        initialize_context2pos();
        abs_pos2chr.clear();
        genome_size = 0;
        auto skipped_contexts = init_skipped_contexts();

        auto regions_to_avoid_by_chr = split_by_chromosome_id(regions_to_avoid);

        skip_to_next_seq(fasta_stream);

        while (!fasta_stream.eof()) {
            std::string sequence_title;
            getline(fasta_stream, sequence_title);

            ChromosomeId chr_id;

            // if the sequence is a chromosome
            if (IO::FASTA::is_chromosome_header(sequence_title, chr_id)) {

                if (progress_bar != nullptr) {
                    progress_bar->set_message("Processing chr. " + GenomicPosition::chrtos(chr_id));
                }
                abs_pos2chr[genome_size+1] = chr_id;

                build_index_in_seq(fasta_stream, streamsize, chr_id, regions_to_avoid_by_chr[chr_id],
                                   skipped_contexts, sampling_rate, progress_bar);
            } else {
                if (progress_bar != nullptr) {
                    progress_bar->set_message("Discharging a sequence");
                    progress_bar->set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()/streamsize));
                }
                discharge_seq(fasta_stream);
            }
        }

        if (progress_bar != nullptr) {
            progress_bar->set_progress(100, "Context index built");
        }
    }

    /**
     * @brief Reset the object and find the mutational contexts in a FASTA file
     *
     * This method resets the current object, finds the mutational contexts in the chromosome sequences
     * contained in a FASTA file, and saves their absolute positions in the object member
     * `ContextIndex::context2pos`. The chromosome sequences are identified by their names according
     * to the method template parameter.
     * At the same time, this method stores the absolute positions of each chromosome in the member
     * `ContextIndex::abs_pos2chr`.
     *
     * @param[in,out] fasta_filename is the FASTA filename
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in] sampling_rate is the number of contexts to be found in order to record a context
     *          in the index
     * @param[in,out] progress_bar is the progress bar
     */
    inline void reset_with(const std::filesystem::path& fasta_filename, const std::set<GenomicRegion>& regions_to_avoid,
                           const size_t& sampling_rate, UI::ProgressBar* progress_bar)
    {

        std::ifstream genome_fasta_stream(fasta_filename);

        if (!genome_fasta_stream.good()) {
            std::ostringstream oss;

            oss << "\"" << to_string(fasta_filename) << "\" does not exist";
            throw std::runtime_error(oss.str());
        }

        reset_with(genome_fasta_stream, regions_to_avoid, sampling_rate, progress_bar);
    }

public:

    /**
     * @brief The empty constructor
     */
    ContextIndex():
        context2pos(std::make_shared<ContextPositionMap>()), genome_size(0)
    {}

    /**
     * @brief Find the context positions in a FASTA file
     *
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the contexts that lay in the sequences corresponding to a
     *          chromosome according to `CLONES::IO::FASTA::seq_name_decoders`
     */
    static inline ContextIndex build_index(const std::filesystem::path& genome_fasta,
                                           UI::ProgressBar* progress_bar=nullptr)
    {
        return build_index(genome_fasta, {}, 1, progress_bar);
    }

    /**
     * @brief Find the context positions in a FASTA file
     *
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] sampling_rate is the number of contexts to be found in order to record a context
     *          in the index
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the contexts that lay in the sequences corresponding to a
     *          chromosome according to `CLONES::IO::FASTA::seq_name_decoders`
     */
    static inline ContextIndex build_index(const std::filesystem::path& genome_fasta,
                                           const size_t& sampling_rate,
                                           UI::ProgressBar* progress_bar=nullptr)
    {
        return build_index(genome_fasta, {}, sampling_rate, progress_bar);
    }

    /**
     * @brief Find the context positions in some genomic fragments of a FASTA file
     *
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the contexts that lay in the sequences corresponding
     *      to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`,
     *      but that are located outside the regions in `regions_to_avoid`
     */
    static inline ContextIndex build_index(const std::filesystem::path& genome_fasta,
                                           const std::set<GenomicRegion>& regions_to_avoid,
                                           UI::ProgressBar* progress_bar=nullptr)
    {
        return build_index(genome_fasta, regions_to_avoid, 1, progress_bar);
    }

    /**
     * @brief Find the context positions in some genomic fragments of a FASTA file
     *
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in] sampling_rate is the number of contexts to be found in order to record a context
     *          in the index
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the contexts that lay in the sequences corresponding
     *      to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`,
     *      but that are located outside the regions in `regions_to_avoid`
     */
    static ContextIndex build_index(const std::filesystem::path& genome_fasta,
                                    const std::set<GenomicRegion>& regions_to_avoid,
                                    const size_t& sampling_rate,
                                    UI::ProgressBar* progress_bar=nullptr)
    {
        ContextIndex context_index;

        context_index.reset_with(genome_fasta, regions_to_avoid, sampling_rate, progress_bar);

        return context_index;
    }

    /**
     * @brief Get the context positions
     *
     * @return a constant reference to simulator context positions
     */
    inline const std::map<SBSContext, std::vector<GENOME_WIDE_POSITION> >& get_context_positions() const
    {
        return *context2pos;
    }

    /**
     * @brief Get the absolute genomic positions of a context
     *
     * @param context is the searched context
     * @return the vector of the absolute genomic positions
     */
    inline const std::vector<GENOME_WIDE_POSITION>& operator[](const SBSContext& context) const
    {
        return context2pos->at(context);
    }

    /**
     * @brief Extract a genome position from the index
     *
     * @param context is the context whose position will be extracted
     * @param index is the index of the position that will be extracted
     * @return the extracted absolute position for `context`
     */
    GENOME_WIDE_POSITION extract(const SBSContext& context, const size_t index)
    {
        auto& context_pos = context2pos->at(context);

        if (context_pos.size() < index+1) {
            throw std::out_of_range("The context does not have so many positions");
        }

        if (context_pos.size() > index+1) {
            std::swap(context_pos[index], context_pos.back());
        }

        auto extracted = context_pos.back();

        context_pos.pop_back();

        return extracted;
    }

    /**
     * @brief Insert a genomic position for a context
     *
     * This method inserts an absolute genomic position in the
     * last position of the position vector of a context
     *
     * @param context is the context whose position will be inserted
     * @param absolute_position is the absolute position to insert
     */
    void insert(const SBSContext& context, const GENOME_WIDE_POSITION& absolute_position)
    {
        context2pos->at(context).push_back(absolute_position);
    }

    /**
     * @brief Insert a genomic position for a context
     *
     * @param context is the context whose position will be inserted
     * @param absolute_position is the absolute position to insert
     * @param index is the index in the position vector of `context`
     *          in which `absolute_position` will inserted
     */
    void insert(const SBSContext& context, const GENOME_WIDE_POSITION& absolute_position,
                const size_t index)
    {
        auto& context_pos = context2pos->at(context);

        if (context_pos.size() < index) {
            throw std::out_of_range("The context does not have so many positions");
        }

        context_pos.push_back(absolute_position);

        if (context_pos.size() > index+1) {
            std::swap(context_pos[index], context_pos.back());
        }
    }

    /**
     * @brief Get the regions of the chromosomes
     *
     * @return the regions of the chromosomes
     */
    std::list<GenomicRegion> get_chromosome_regions() const
    {
        std::list<GenomicRegion> chr_regions;

        if (abs_pos2chr.size()>0) {
            auto it = abs_pos2chr.begin();
            auto next = it;

            while (++next != abs_pos2chr.end()) {
                chr_regions.emplace_back(it->second, static_cast<GenomicRegion::Length>(next->first-it->first));
                it = next;
            }

            chr_regions.emplace_back(it->second, static_cast<GenomicRegion::Length>(genome_size-it->first+1));
        }

        return chr_regions;
    }

    /**
     * @brief Get the chromosome identifiers of the reference sequence
     *
     * @return a vector containing the chromosome identifiers of the
     *      reference sequence
     */
    std::vector<ChromosomeId> get_chromosome_ids() const
    {
        std::vector<ChromosomeId> chromosome_ids;

        chromosome_ids.reserve(abs_pos2chr.size());
        for (const auto& [abs_pos, chr_id]: abs_pos2chr) {
            chromosome_ids.push_back(chr_id);
        }

        return chromosome_ids;
    }

    /**
     * @brief Get the number of chromosomes in the reference genome
     *
     * @return the number of chromosomes in the reference genome
     */
    inline size_t num_of_chromosomes() const
    {
        return abs_pos2chr.size();
    }

    /**
     * @brief Turn an absolute genomic position into the corresponding genomic position
     *
     * @param abs_position is an absolute genomic position
     * @return the genomic position corresponding to `abs_position`
     */
    GenomicPosition get_genomic_position(const GENOME_WIDE_POSITION& abs_position) const
    {
        auto it = abs_pos2chr.upper_bound(abs_position);
        --it;

        return GenomicPosition(it->second, static_cast<ChrPosition>(abs_position-it->first+1));
    }

    /**
     * @brief Get the genome size
     *
     * @return get the genome size
     */
    inline const GENOME_WIDE_POSITION& get_genome_size() const
    {
        return genome_size;
    }

    /**
     * @brief Save a simulator in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "CLONES Context Index", 0);

        uint8_t abs_pos_size = sizeof(GENOME_WIDE_POSITION);

        archive & abs_pos_size
                & *context2pos
                & abs_pos2chr
                & genome_size;
    }

    /**
     * @brief Load a simulator from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded simulator
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static ContextIndex load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "CLONES Context Index", 0);

        ContextIndex context_index;
        uint8_t abs_pos_size;

        archive & abs_pos_size;

        if (abs_pos_size != sizeof(GENOME_WIDE_POSITION)) {
            std::ostringstream oss;

            oss << "Absolute position size (" << sizeof(GENOME_WIDE_POSITION)
                << " bytes) does not match file object one (" << static_cast<size_t>(abs_pos_size)
                << " bytes)." ;

            throw std::runtime_error(oss.str());
        }

        context_index.context2pos = std::make_shared<ContextPositionMap>();

        archive & *(context_index.context2pos)
                & context_index.abs_pos2chr
                & context_index.genome_size;

        return context_index;
    }

    /**
     * @brief Read the number of bytes per absolute position in an input archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the archive from which the number of bytes per absolute
     *          position must be read
     * @return the number of bytes per absolute position
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static size_t read_bytes_per_absolute_position(ARCHIVE& archive)
    {
        uint8_t abs_pos_size;

        auto orig_pos = archive.tellg();

        archive & abs_pos_size;

        archive.seekg(orig_pos);

        return abs_pos_size;
    }
};

}   // CLONES

}   // Mutations


#endif // __CLONES_CONTEXT_INDEX__