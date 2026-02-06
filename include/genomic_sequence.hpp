/**
 * @file genomic_sequence.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a structure to handle genomic sequence
 * @version 1.3
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

#ifndef __CLONES_GENOMIC_SEQUENCE__
#define __CLONES_GENOMIC_SEQUENCE__

#include <string>
#include <vector>
#include <algorithm>

namespace CLONES
{

/**
 * @brief A class to handle genomic sequences
 */
struct GenomicSequence
{
    /**
     * @brief The vector of DNA bases
     */
    static std::vector<char> DNA_bases;

    /**
     * @brief Get the base index
     *
     * @param base is a nucleotide
     * @return is the index of `base` in
     *      `GenomicSequence::DNA_bases`
     */
    static size_t get_base_index(const char& base);

    /**
     * @brief Get the complemented base
     *
     * @param base is the base whose complement is request
     * @return the complemented base of `base`
     */
    static char get_complement(const char& base);

    /**
     * @brief Complement a sequence
     *
     * This method complements a sequence in-place.
     */
    static void complement(std::string& sequence);

    /**
     * @brief Get the complemented sequence
     *
     * @param sequence is the sequence whose complement is request
     * @return the complemented sequence of `sequence`
     */
    static std::string get_complement(const std::string& sequence);

    /**
     * @brief Reverse a sequence
     *
     * This method reveres a sequence in-place.
     */
    static inline void reverse(std::string& sequence)
    {
        std::reverse(sequence.begin(), sequence.end());
    }

    /**
     * @brief Get the reverse sequence
     *
     * @param sequence is the sequence whose reverse is request
     * @return the reverse sequence of `sequence`
     */
    static std::string get_reverse(const std::string& sequence);

    /**
     * @brief Reverse and complement a sequence
     *
     * This method reverses and complements a sequence in-place.
     */
    static inline void reverse_complement(std::string& sequence)
    {
        GenomicSequence::reverse(sequence);
        GenomicSequence::complement(sequence);
    }

    /**
     * @brief Get the reverse-complemented sequence
     *
     * @param sequence is the sequence whose reverse-complement is request
     * @return the reverse-complemented sequence of `sequence`
     */
    static std::string get_reverse_complement(const std::string& sequence);

    /**
     * @brief Check whether a symbol represents a DNA nucleic base
     *
     * @param symbol is the symbol to be tested
     * @return `true` if and only if symbol represents a base, i.e.,
     *       it is one among 'A', 'C', 'G', 'T', 'a', 'c', 'g', or 't'.
     */
    static bool is_a_DNA_base(const char& symbol);
};

}   // CLONES

#endif // __CLONES_GENOMIC_SEQUENCE__