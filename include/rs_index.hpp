/**
 * @file rs_index.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to compute the repeated substring index
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

#ifndef __CLONES_RS_INDEX__
#define __CLONES_RS_INDEX__

#include <string>
#include <vector>
#include <tuple>
#include <ostream>
#include <random>
#include <memory>

#include "genomic_position.hpp"
#include "genomic_region.hpp"
#include "utils.hpp"
#include "progress_bar.hpp"
#include "archive.hpp"

#include "id_signature.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief Repeated sequences
 *
 * This class represents repeated sequences.
 *
 * @tparam REPETITION_TYPE is the number of repetitions type
 */
template<typename REPETITION_TYPE>
struct Repetition
{
    GenomicPosition g_position;  //!< The genomic position of the repeated sequence
    REPETITION_TYPE num_of_repetitions; //!< The number of repetitions of the unit
    std::string unit;   //!< The repeated subsequence
    char prev_base;     //!< The base preceding the first repetition of the unit

    /**
     * @brief The empty constructor
     */
    Repetition():
        g_position(), num_of_repetitions(0), unit(""), prev_base('N')
    {}

    /**
     * @brief A constructor
     *
     * @param chr_id is the identifier of the chromosome containing
     *      the repeated sequence
     * @param begin is the position of the repeated sequence first base
     *      in the chromosome
     * @param num_of_repetitions is the number of repetition of the unit
     * @param unit is the pointer to the unit first base in the sequence
     * @param unit_size is the unit size
     * @param previous_base is base preceding the first repetition of
     *      the unit
     */
    Repetition(const ChromosomeId& chr_id, const ChrPosition begin,
               const REPETITION_TYPE num_of_repetitions,
               const char* unit, const uint8_t unit_size,
               const char previous_base):
        g_position(chr_id, begin), num_of_repetitions(num_of_repetitions),
        unit(unit, unit+unit_size), prev_base(previous_base)
    {
        if (unit_size==0) {
            throw std::domain_error("Unit size must be greater than 0.");
        }
    }

    /**
     * @brief Save a repetition in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & g_position
                & num_of_repetitions
                & unit
                & prev_base;
    }

    /**
     * @brief Load a repetition from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded repetition
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static Repetition<REPETITION_TYPE> load(ARCHIVE& archive)
    {
        Repetition<REPETITION_TYPE> repetition;

        archive & repetition.g_position
                & repetition.num_of_repetitions
                & repetition.unit
                & repetition.prev_base;

        return repetition;
    }
};

}   // Mutations

}   //  CLONES


namespace std
{
    template<typename RepetitionType>
    std::ostream& operator<<(std::ostream& os, const CLONES::Mutations::Repetition<RepetitionType>& repetition)
    {
        os << "[" << std::string(1, repetition.prev_base) << "]" << repetition.unit << " ("
            << repetition.begin << ":" << static_cast<size_t>(repetition.num_of_repetitions) << ")";

        return os;
    }
}

namespace CLONES
{

namespace Mutations
{


/**
 * @brief An index of the genome repeated sequences
 *
 * This index stores the some of the repeated sequences of genomic sequences.
 */
struct RSIndex
{
    using RepetitionType = IDType::RepetitionType;

    /**
     * @brief A repeated sequence storage
     *
     * The objects of this class store a sample of the identified
     * repeated sequences.
     */
    struct RepetitionStorage
    {
        size_t total_number;        //!< The total number of the identified repeated sequences
        size_t stored_repetitions;  //!< The number of repetitions stored in storage
        std::vector<Repetition<RepetitionType>> repetitions;  //!< The samples repeated sequences

        /**
         * @brief The empty constructor
         */
        RepetitionStorage();

        /**
         * @brief Extract a repetition
         *
         * @param pos is the position of the repetition to be removed
         * @return a constant reference to the extracted repetition
         */
        const Repetition<RepetitionType>& extract(const size_t& pos);

        /**
         * @brief Insert a repetition
         *
         * @param repetition is the repetition to insert
         */
        void push_back(Repetition<RepetitionType>&& repetition);

        /**
         * @brief Get the `n`-th repetition in the storage
         *
         * @param pos is the position of the aimed repetition in the storage
         * @return a constant reference to the `pos`-th repetition in the
         *      storage
         */
        inline const Repetition<RepetitionType>& operator[](const size_t& pos) const
        {
            return repetitions[pos];
        }

        /**
         * @brief Get the `n`-th repetition in the storage
         *
         * @param pos is the position of the aimed repetition in the storage
         * @return a reference to the `pos`-th repetition in the
         *      storage
         */
        inline Repetition<RepetitionType>& operator[](const size_t& pos)
        {
            return repetitions[pos];
        }

        /**
         * @brief Get the number of repetitions currently in the storage
         *
         * @return a constant reference to the number of repetitions currently
         *      in the storage
         */
        inline const size_t& size() const
        {
            return stored_repetitions;
        }

        /**
         * @brief Get the initial number of repetitions in the storage
         *
         * @return the initial number of repetitions in the storage
         */
        inline size_t full_size() const
        {
            return repetitions.size();
        }

        /**
         * @brief Restore the last extracted repetitions
         *
         * @param num_of_repetitions is the number of the last repetitions to
         *      be restored
         */
        inline void restore_extracted(const size_t num_of_repetitions=1)
        {
            stored_repetitions += num_of_repetitions;
        }

        /**
         * @brief Restore the extracted repetitions
         */
        inline void restore()
        {
            stored_repetitions = full_size();
        }

        /**
         * @brief Save a repetition storage in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & total_number
                    & stored_repetitions
                    & repetitions;
        }

        /**
         * @brief Load a repetition storage from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the loaded repetition storage
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
        inline static RepetitionStorage load(ARCHIVE& archive)
        {
            RepetitionStorage storage;

            archive & storage.total_number
                    & storage.stored_repetitions
                    & storage.repetitions;

            return storage;
        }
    };

    /**
     * @brief The second level index to repetition-storage map
     */
    using RepetitionMap = std::map<uint8_t, RepetitionStorage>;

    std::shared_ptr<std::map<uint8_t, RepetitionMap>> hetero_map;    //!< The heteropolymer map unit-size to repetition-map
    std::shared_ptr<std::map<char, RepetitionMap>> homo_map;         //!< Tha homopolymer map nucleotide to repetition-map
    std::shared_ptr<std::map<uint8_t, RepetitionMap>> micro_map;     //!< The micro-homology map deletion size to repetition-map

    std::mt19937_64 random_gen;     //!< The index random generator

    size_t max_unit_size;           //!< The maximum considered size of the repetition unit
    size_t max_stored_repetitions;  //!< The maximum number of stored repetitions

    /**
     * @brief Initialize the suffix array vector
     *
     * @param s is the sequence whose suffix array is aimed
     * @param suffix_array is the vector that will contain the suffix array
     * @param classes is the vector that labels each sequence position with the
     *      class of its first nucleotide
     * @return the number of the different nucleotides in the sequence
     */
    static size_t init_suffix_array(const char* s, std::vector<ChrPosition>& suffix_array,
                                    std::vector<ChrPosition>& classes);
    /**
     * @brief Upgrade the (h-1)-suffix array to a (h)-suffix array
     *
     * @param h is the length of the prefices sorted in the (h)-suffix array
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
     * @brief Check whether a nucleotide is "A"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is "A"
     */
    inline static bool is_A(const char& nucleotide)
    {
        return nucleotide=='A' || nucleotide=='a';
    }

    /**
     * @brief Check whether a nucleotide is "T"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is "T"
     */
    inline static bool is_T(const char& nucleotide)
    {
        return nucleotide=='T' || nucleotide=='t';
    }

    /**
     * @brief Check whether a nucleotide is "C"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is "C"
     */
    inline static bool is_C(const char& nucleotide)
    {
        return nucleotide=='C' || nucleotide=='c';
    }

    /**
     * @brief Check whether a nucleotide is "G"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is "G"
     */
    inline static bool is_G(const char& nucleotide)
    {
        return nucleotide=='G' || nucleotide=='g';
    }

    /**
     * @brief Check whether a nucleotide is either "A" or "T"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is either
     *      "A" or "T"
     */
    inline static bool is_AT(const char& nucleotide)
    {
        return is_A(nucleotide) || is_T(nucleotide);
    }

    /**
     * @brief Check whether a nucleotide is either "C" or "G"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is either
     *      "C" or "G"
     */
    inline static bool is_CG(const char& nucleotide)
    {
        return is_C(nucleotide) || is_G(nucleotide);
    }

    /**
     * @brief Add a repeated sequence to a repetition map
     *
     * @param r_map is the repetition map in which the repetition should
     *      be added
     * @param chr_id is the identifier of the chromosome containing
     *      the repeated sequence
     * @param begin is the position of the repeated sequence first base
     *      in the chromosome
     * @param index is the index of the repetition in the repetion map
     * @param num_of_repetitions is the number of repetition of the unit
     * @param unit is the pointer to the unit first base in the sequence
     * @param unit_size is the unit size
     * @param previous_base is base preceding the first repetition of
     *      the unit
     * @return `true` if and only if the repeated sequence has been
     *      sampled and saved in the storage
     */
    bool add(RepetitionMap& r_map, const ChromosomeId& chr_id,
             const ChrPosition& begin,
             const uint8_t index,
             const RepetitionType& num_of_repetitions,
             const char*& unit, const uint8_t& unit_size,
             const char& previous_base);

    /**
     * @brief Add a repeated sequence to the index
     *
     * @param s is the considered genomic sequence
     * @param chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param unit_size is the size of the repetition unit
     * @param r_begin is the position of the repeated sequence first base
     *      on the considered genomic sequence
     * @param r_end is the position of the repeated sequence last base
     *      on the considered genomic sequence
     * @param covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void add_repetition(const char* s, const ChromosomeId& chr_id,
                        const ChrPosition& begin, const size_t unit_size,
                        const ChrPosition& r_begin, const ChrPosition& r_end,
                        std::vector<bool>& covered);

    /**
     * @brief Build a random k-mer
     *
     * @param k is the size of the k-mer to build
     * @return a random k-mer
     */
    std::string build_k_mer(const uint8_t k);

    //! @private
    void add_null_heteropolymer(const char* s, const ChromosomeId& chr_id,
                                const uint8_t unit_size,
                                const ChrPosition& begin, const ChrPosition& r_begin);


    //! @private
    void add_null_homopolymer(const size_t nucleotide_index, const char* s,
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
     * @param s is the considered genomic sequence
     * @param chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param h is the order of the suffix array
     * @param h_suffix_array is the (h)-suffix array
     * @param classes is the class vector of the (h)-suffix array
     * @param covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void collect_repetitions(const char* s, const ChromosomeId& chr_id,
                             const ChrPosition& begin, const size_t& h,
                             std::vector<ChrPosition>& h_suffix_array,
                             std::vector<ChrPosition>& classes,
                             std::vector<bool>& covered);

    /**
     * @brief Collect micro-homologies
     *
     * @param s is the considered genomic sequence
     * @param chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void collect_microhomologies(const char* s, const ChromosomeId& chr_id,
                                 const ChrPosition& begin, std::vector<bool>& covered);

    /**
     * @brief Analyze non-repeated sequence
     *
     * @param s is the considered genomic sequence
     * @param chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void collected_non_repeated_seq(const char* s, const ChromosomeId& chr_id,
                                    const ChrPosition& begin, std::vector<bool>& covered);

    /**
     * @brief Add the repeated sequences of a genomic sequence
     *
     * @param s is the considered genomic sequence
     * @param[in] chr_id is the identifier of the chromosome containing the
     *      genomic sequence
     * @param[in] begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param[in] length is the length of the considered sequence
     * @param[in,out] progress_bar is the progress bar
     */
    std::vector<bool>
    collect_repetitions(const char* s, const ChromosomeId& chr_id,
                        const ChrPosition begin, const size_t& length,
                        UI::ProgressBar* progress_bar = nullptr);

    /**
     * @brief Collect the data from a genomic sequence
     *
     * @param[in] sequence is a genomic sequence
     * @param[in] chr_id is the identifier of the chromosome containing the
     *      genomic sequence
     * @param[in] begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param[in] length is the length of the considered sequence
     * @param[in,out] progress_bar is the progress bar
     */
    void collect_data_from(const std::string& sequence,
                           const ChromosomeId& chr_id,
                           const ChrPosition begin, size_t length,
                           UI::ProgressBar* progress_bar = nullptr);

    //!< @private
    static inline uint8_t get_first_index(const uint8_t& unit_size)
    {
        return static_cast<uint8_t>(unit_size>5?5:unit_size);
    }

    //!< @private
    static inline uint8_t get_second_index(const RSIndex::RepetitionType& num_of_repetitions,
                                           const bool& for_insertion)
    {
        return static_cast<uint8_t>(std::min(static_cast<RSIndex::RepetitionType>(for_insertion?5:6),
                                             num_of_repetitions));
    }

    /**
     * @brief Find the storages matching the parameters
     *
     * @tparam INDEX_TYPE is the polymer map index type
     * @param polymer_map is the polymer map (i.e., `hetero_map` or `homo_map`)
     * @param index is the polymer map index (i.e., unit size or homopolymer nucleotide)
     * @param num_of_repetitions is the number of repetitions of the searched polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion
     * @return the vector of the pointers of the storages matching the parameter
     *      specification
     */
    template<typename INDEX_TYPE>
    static std::vector<RepetitionStorage*>
    find_storages(std::map<INDEX_TYPE, RepetitionMap>& polymer_map, const INDEX_TYPE index,
                  const uint8_t& num_of_repetitions, const bool& for_insertion)
    {
        const uint8_t repeated_unit = get_second_index(num_of_repetitions, for_insertion);

        auto p_it = polymer_map.find(index);

        std::vector<RepetitionStorage*> r_vectors;

        if (p_it != polymer_map.end()) {
            auto &p_map = p_it->second;
            auto r_it = p_map.find(repeated_unit);
            if (r_it != p_map.end()) {
                r_vectors.push_back(&(r_it->second));

                if (for_insertion && repeated_unit==5) {
                    r_it = p_map.find(6);

                    if (r_it != p_map.end()) {
                        r_vectors.push_back(&(r_it->second));
                    }
                }
            }
        }

        return r_vectors;
    }

    /**
     * @brief Find the storages matching the parameters
     *
     * @tparam INDEX_TYPE is the polymer map index type
     * @param polymer_map is the polymer map (i.e., `hetero_map` or `homo_map`)
     * @param index is the polymer map index (i.e., unit size or homopolymer nucleotide)
     * @param num_of_repetitions is the number of repetitions of the searched polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion
     * @return the vector of the pointers of the storages matching the parameter
     *      specification
     */
    template<typename INDEX_TYPE>
    static std::vector<RepetitionStorage const*>
    find_storages(const std::map<INDEX_TYPE, RepetitionMap>& polymer_map, const INDEX_TYPE index,
                  const uint8_t& num_of_repetitions, const bool& for_insertion)
    {
        const uint8_t repeated_unit = get_second_index(num_of_repetitions, for_insertion);

        auto p_it = polymer_map.find(index);

        std::vector<RepetitionStorage const*> r_vectors;

        if (p_it != polymer_map.end()) {
            auto &p_map = p_it->second;
            auto r_it = p_map.find(repeated_unit);
            if (r_it != p_map.end()) {
                r_vectors.push_back(&(r_it->second));

                if (for_insertion && repeated_unit==5) {
                    r_it = p_map.find(6);

                    if (r_it != p_map.end()) {
                        r_vectors.push_back(&(r_it->second));
                    }
                }
            }
        }

        return r_vectors;
    }

    /**
     * @brief Get a polymer reference
     *
     * @tparam INDEX_TYPE is the polymer map index type
     * @param polymer_map is the polymer map (i.e., `hetero_map` or `homo_map`)
     * @param polymer_type is a description of the polymer map
     * @param index is the polymer map index (i.e., unit size or homopolymer nucleotide)
     * @param num_of_repetitions is the number of repetitions of the searched polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion
     * @return a pair repetition storage-position in the storage refering to a polymer
     *      satisfing the parameters
     */
    template<typename INDEX_TYPE>
    std::pair<RepetitionStorage*, size_t>
    find_a_polymer(std::map<INDEX_TYPE, RepetitionMap>& polymer_map, const std::string& polymer_type,
                   const INDEX_TYPE index, const uint8_t& num_of_repetitions, const bool& for_insertion)
    {
        auto s_pointers = find_storages(polymer_map, index, num_of_repetitions, for_insertion);

        size_t total_size{0};
        for (const auto& s_pointer: s_pointers) {
            total_size += s_pointer->size();
        }

        if (total_size>0) {

            std::uniform_int_distribution<size_t> dist(0, total_size-1);

            auto pos = dist(random_gen);

            for (const auto& r_pointer : s_pointers) {
                if (pos < r_pointer->size()) {
                    return {r_pointer, pos};
                }

                pos -= r_pointer->size();
            }
        }
        throw std::runtime_error("No such a " + polymer_type + " available.");
    }

    /**
     * @brief Find a heteropolymer
     *
     * @param unit_size is the unit size of the aimed heteropolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion
     * @return a pair repetition storage-position refering to a heteropolymer having
     *      unit size `unit_size` and whose unit is repeated `num_of_repetitions` times,
     *      when `for_insertion` is `false` or `num_of_repetitions` is smaller than 6,
     *      or either 5 or 6 times, when `for_insertion` is `true` and
     *      `num_of_repetitions` is 5
     */
    inline std::pair<RepetitionStorage*, size_t>
    find_a_heteropolymer(const uint8_t& unit_size, const uint8_t& num_of_repetitions,
                         const bool& for_insertion)
    {
        return find_a_polymer(*hetero_map, "heteropolymer", get_first_index(unit_size),
                              num_of_repetitions, for_insertion);
    }

    /**
     * @brief Find a homopolymer
     *
     * @param unit_nucleotide is the unit nucleotide of the aimed homopolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion or deletion
     * @return a pair repetition storage-position refering to a homopolymer whose
     *      unit is `unit_nucleotide` and it is repeated `num_of_repetitions`
     *      times, when `for_insertion` is `false` or `num_of_repetitions` is
     *      smaller than 6, or either 5 or 6 times, when `for_insertion` is `true`
     *      and `num_of_repetitions is 5
     */
    inline std::pair<RepetitionStorage*, size_t>
    find_a_homopolymer(const char& unit_nucleotide, const uint8_t& num_of_repetitions,
                       const bool& for_insertion)
    {
        return find_a_polymer(*homo_map, "homopolymer", canonize(unit_nucleotide),
                              num_of_repetitions, for_insertion);
    }

    /**
     * @brief Find a microhomolory
     *
     * @param homology_distance is the aimed micro-homology distance
     * @param homology_size is the aimed micro-homology size
     * @return a pair repetition storage-position refering to a micro-homology
     *      between sequences whose distance is `homology_distance` and whose
     *      size if `homology_size`.
     */
    inline std::pair<RepetitionStorage*, size_t>
    find_a_microhomology(const uint8_t& homology_distance, const uint8_t& homology_size)
    {
        return find_a_polymer(*micro_map, "microhomolory", get_first_index(homology_distance),
                              homology_size, false);
    }

    /**
     * @brief Canonize a nucleotide
     *
     * @param nucleotide is a nucleotide among 'A', 'C', 'G', and 'T'
     * @return if `nucleotide` is either 'A' or 'T', then 'T'; if
     *      `nucleotide` is either 'A' or 'T', then 'C'
     */
    static char canonize(const char& nucleotide);

    /**
     * @brief Add a polymeric repetition to the index
     *
     * @param chr_id is the identifier of the chromosome containing
     *      the repeated sequence
     * @param begin is the position of the repeated sequence first base
     *      in the chromosome
     * @param num_of_repetitions is the number of repetition of the unit
     * @param unit is the pointer to the unit first base in the sequence
     * @param unit_size is the unit size
     * @param previous_base is base preceding the first repetition of
     *      the unit
     */
    void add_polymer(const ChromosomeId& chr_id, const ChrPosition& begin,
                     const RepetitionType& num_of_repetitions,
                     const char* unit, const size_t& unit_size,
                     const char& previous_base);

    /**
     * @brief Add to the index all the repetitions in a chromosome
     *
     * @param[in] chr_id is the identifier of a chrosomome
     * @param[in] chr_sequence is the nucleotide sequence of a chrosomome
     * @param[in,out] progress_bar is the progress bar
     */
    void collect_data_from(const ChromosomeId& chr_id, const std::string& chr_sequence,
                           UI::ProgressBar* progress_bar = nullptr);
public:
    /**
     * @brief A constructor
     *
     * @param max_unit_size is the maximum considered size of the repetition unit
     * @param max_stored_repetitions is the maximum number of stored repetitions
     * @param seed is the random generator seed
     */
    RSIndex(const size_t max_unit_size = std::numeric_limits<size_t>::max(),
            const size_t max_stored_repetitions = std::numeric_limits<size_t>::max(),
            const int seed = 0);

    /**
     * @brief Randomly select a repetition according to an ID type
     *
     * @param id_type is an ID type
     * @param remove is a Boolean flag that establish whether the selected repetition
     *      must be removed from the index
     * @return a constant reference to a repetition that is consistent with `id_type`
     */
    const Repetition<RepetitionType>&
    select(const IDType& id_type, const bool remove=false);

    /**
     * @brief Restore the last extracted repetition of a given type
     *
     * @param id_type is the ID type of the repetition to be restored in the index
     */
    void restore(const IDType& id_type);

    /**
     * @brief Deep clone the repeated sequence index
     *
     * This method performs a complete copy of the repeated sequence index.
     *
     * @return a deep clone the repeated sequence index
     */
    RSIndex clone() const;

    /**
     * @brief Get the number of repetitions matching an indel type
     *
     * @param id_type is an ID type
     * @return the number of repetitions matching the indel type `id_type`
     */
    size_t count_available_for(const IDType& id_type) const;

    /**
     * @brief Find the context positions in a FASTA file
     *
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] max_stored_repetitions is the maximum number of stored repetitions
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`
     */
    static inline RSIndex build_index(const std::filesystem::path& genome_fasta,
                                      const size_t max_unit_size,
                                      const size_t max_stored_repetitions,
                                      UI::ProgressBar* progress_bar=nullptr)
    {
        return build_index(genome_fasta, max_unit_size, max_stored_repetitions,
                           0, progress_bar);
    }

    /**
     * @brief Find the context positions in a FASTA file
     *
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] max_stored_repetitions is the maximum number of stored repetitions
     * @param[in] seed is the random generator seed
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`
     */
    static inline RSIndex build_index(const std::filesystem::path& genome_fasta,
                                      const size_t max_unit_size,
                                      const size_t max_stored_repetitions,
                                      const int seed,
                                      UI::ProgressBar* progress_bar=nullptr)
    {
        return build_index(genome_fasta, {}, max_unit_size,
                           max_stored_repetitions, seed, progress_bar);
    }

    /**
     * @brief Find the context positions in some genomic fragments of a FASTA file
     *
     * @param[in] genome_fasta is the path of a FASTA file
     * @param[in] regions_to_avoid is a set of regions to avoid
     * @param[in] max_unit_size is the maximum considered size of the repetition unit
     * @param[in] max_stored_repetitions is the maximum number of stored repetitions
     * @param[in] seed is the random generator seed
     * @param[in,out] progress_bar is the progress bar
     * @return the index of the repetitions that lay in the sequences corresponding
     *          to a chromosome according to `CLONES::IO::FASTA::seq_name_decoders`,
     *          but that are located outside the regions in `regions_to_avoid`
     */
    static RSIndex build_index(const std::filesystem::path& genome_fasta,
                               const std::set<GenomicRegion>& regions_to_avoid,
                               const size_t max_unit_size,
                               const size_t max_stored_repetitions,
                               const int seed, UI::ProgressBar* progress_bar=nullptr);

    /**
     * @brief Restore the extracted polymers
     */
    void restore();

    /**
     * @brief Save a repeated sequence index in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "CLONES RS Index", 0);

        archive & sizeof(RepetitionType)
                & *hetero_map
                & *homo_map
                & *micro_map
                & max_unit_size
                & max_stored_repetitions;
    }

    /**
     * @brief Load a repeated sequence index from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded repeated sequence index
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static RSIndex load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "CLONES RS Index", 0);

        RSIndex rs_index;

        size_t size_of_repetition_type;

        archive & size_of_repetition_type;

        if (size_of_repetition_type != sizeof(RepetitionType)) {
            std::ostringstream oss;

            oss << "Unsupported archive type: the archive contains a RSIndex whose repetition type size is "
                << size_of_repetition_type << " whereas the repetition type size of this RSIndex is "
                << sizeof(RepetitionType) << ".";
            throw std::runtime_error(oss.str());
        }

        archive & *(rs_index.hetero_map)
                & *(rs_index.homo_map)
                & *(rs_index.micro_map)
                & rs_index.max_unit_size
                & rs_index.max_stored_repetitions;

        return rs_index;
    }
};

}   // Mutations

}   // CLONES

#endif // __CLONES_RS_INDEX__