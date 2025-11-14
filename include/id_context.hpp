/**
 * @file id_context.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines ID contexts
 * @version 1.0
 * @date 2025-11-14
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

#ifndef __RACES_ID_CONTEXT__
#define __RACES_ID_CONTEXT__

#include <cstdint>
#include <functional>   // std::less
#include <string>
#include <iostream>
#include <limits>

#include "archive.hpp"

namespace RACES
{

namespace Mutations
{

/**
 * @brief A class to represent ID context
 *
 * An ID context is a repeated sequence type. It consists in a fragment type
 * (either homopolymer, heteropolymer, or microhomology), a first level code
 * that is either a base (for homopolymer), or the size of the unit (for
 * heteropolymer and microhomology), and a second level code that is either
 * the number of repetitions (for homopolymer and heteropolymer) or the
 * homology size (for microhomology).
 */
struct IDContext
{
public:
    using FirstLevelType = uint8_t;
    using SecondLevelType = uint8_t;

    /**
     * @brief The type of the fragment involved in the mutation
     */
    enum class FragmentType
    {
        HOMOPOLYMER,    //!< A repeated sequence whose nucleotides are the same
        HETEROPOLYMER,  //!< A repeated sequence whose nucleotides may differ
        MICROHOMOLOGY   //!< A fragment followed by a sequence that matches its prefix
    };

protected:
    FragmentType ftype;         //!< The type of fragment
    FirstLevelType fl_code;     //!< The first level code
    SecondLevelType sl_code;    //!< The second level code
public:

    /**
     * @brief The empty constructor
     */
    IDContext();

    /**
     * @brief A constructor
     *
     * @param fragment_type is the fragment type
     * @param first_level_code is the first level code
     * @param first_level_code is the second level code
     */
    IDContext(const FragmentType& fragment_type,
              const FirstLevelType& first_level_code,
              const SecondLevelType& second_level_code);

    /**
     * @brief A constructor
     *
     * An indel context is conventionally represented by a string in the
     * form `{number}{'A','C','G','T','R','M'}{number}`. The first number
     * is the unit size for both homopolymer and heteropolymer, and the
     * fragment size for microhomology. When the string represents a
     * heteropolymer or a microhomology it is the first level code. The
     * character denotes the fragment type: 'R' stands for heteropolymer,
     * 'M' for microhomology, and 'A', 'C', 'G', and 'T' is the unit base
     * of a homopolymer. Dealing with homopolymer, the character is first
     * level code. The last number is the size of the microhomology (for
     * microhomologies) and the number of repetitions in the sequence
     * (for homo/heteropolymer). This number is used as the second level
     * code.
     *
     * @param context is the textual representation of an ID context
     */
    explicit IDContext(const std::string& context);

    /**
     * @brief Build the ID context of an homopolymer
     *
     * @param unit_base is the nucleotide of the homopolymer
     * @param num_of_repetitions is the number of repetitions
     * @return the ID context corresponding to the homopolymer having
     *      `unit_base` as the unit base and `num_of_repetitions`
     *      repetitions
     */
    static inline
    IDContext build_for_homopolymer(const char unit_base,
                                    const uint8_t num_of_repetitions)
    {
        return {FragmentType::HOMOPOLYMER,
                static_cast<const FirstLevelType &>(unit_base),
                num_of_repetitions};
    }

    /**
     * @brief Build the ID context of an heteropolymer
     *
     * @param unit_size is the length of the heteropolymer's repeated unit
     * @param num_of_repetitions is the number of repetitions
     * @return the ID context corresponding to the heteropolymer having
     *      a repeated unit with size `unit_size` and `num_of_repetitions`
     *      repetitions
     */
    static inline
    IDContext build_for_heteropolymer(const uint8_t unit_size,
                                      const uint8_t num_of_repetitions)
    {
        return {FragmentType::HETEROPOLYMER, unit_size, num_of_repetitions};
    }

    /**
     * @brief Build the ID context of an microhomology
     *
     * @param unit_size is the length of the microhomology's repeated unit
     * @param homology_size is the number of repetitions
     * @return the ID context corresponding to the heteropolymer having
     *      a repeated unit with size `unit_size` and `num_of_repetitions`
     *      repetitions
     */
    static inline
    IDContext build_for_microhomology(const uint8_t unit_size,
                                      const uint8_t homology_size)
    {
        return {FragmentType::MICROHOMOLOGY, unit_size, homology_size};
    }

    /**
     * @brief Check whether the ID context is defined
     *
     * @return `true` if and only if the current object has been created
     *      by using the empty constructor
     */
    inline bool is_defined() const
    {
        return (sl_code != std::numeric_limits<SecondLevelType>::max());
    }

    /**
     * @brief Get the fragment type
     *
     * @return the fragment type
     */
    inline const FragmentType& fragment_type() const
    {
        return ftype;
    }

    /**
     * @brief Get the first level code
     *
     * @return a constant reference to the first level code
     */
    inline const FirstLevelType& get_first_level_code() const
    {
        return fl_code;
    }

    /**
     * @brief Get the homopolymer unit base
     *
     * @return a constant reference to the homopolymer unit base
     */
    char unit_base() const;

    /**
     * @brief Get the heteropolymer or microhomology unit size
     *
     * @return the heteropolymer or microhomology unit size
     */
    const FirstLevelType& unit_size() const
    {
        if (ftype == FragmentType::HOMOPOLYMER) {
            throw std::runtime_error("IDContext::unit_size(): \"this\" is "
                                     "a homopolymer.");
        }

        return get_first_level_code();
    }

    /**
     * @brief Get the second level code
     *
     * @return a constant reference to the second level code
     */
    inline const SecondLevelType& get_second_level_code() const
    {
        return sl_code;
    }

    /**
     * @brief Get the number of repetitions of the (homo/hetero)-polymer
     *
     * @return the number of repetitions of either the homopolymer or
     *      the heteropolymer
     */
    const SecondLevelType& num_of_repetitions() const
    {
        if (ftype == FragmentType::MICROHOMOLOGY) {
            throw std::runtime_error("IDContext::num_of_repetitions(): \"this\" is "
                                     "a microhomology.");
        }

        return get_second_level_code();
    }

    /**
     * @brief Get the microhomology size
     *
     * @return the microhomology size
     */
    const SecondLevelType& microhomology_size() const
    {
        if (ftype != FragmentType::MICROHOMOLOGY) {
            throw std::runtime_error("IDContext::microhomology_size(): \"this\" is "
                                     "a (homo/hetero)-polymer.");
        }

        return get_second_level_code();
    }

    /**
     * @brief Test whether two ID contexts are equivalent
     *
     * @param context is the ID context to compare
     * @return `true` if and only if the two ID contexts represent
     *      the same repeated sequence type
     */
    inline bool operator==(const IDContext& context) const
    {
        return ((ftype == context.ftype)
                && (fl_code == context.fl_code)
                && (sl_code == context.sl_code));
    }

    /**
     * @brief Test whether two ID contexts differ
     *
     * @param context is the ID context to compare
     * @return `true` if and only if the two ID contexts represent
     *      different repeated sequence types
     */
    inline bool operator!=(const IDContext& context) const
    {
        return !(this->operator==(context));
    }

    /**
     * @brief Save a ID context in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<uint8_t>(ftype)
                & fl_code
                & sl_code;
    }

    /**
     * @brief Load a ID context from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded ID context
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static IDContext load(ARCHIVE& archive)
    {
        IDContext context;

        uint8_t ftype;
        archive & ftype
                & context.fl_code
                & context.sl_code;

        context.ftype = static_cast<FragmentType>(ftype);

        return context;
    }
};

}   // Mutations

}   // RACES

/**
 * @brief A specialization of `std::less` for `IDContext`
 */
template<>
struct std::less<RACES::Mutations::IDContext>
{
    bool operator()(const RACES::Mutations::IDContext &lhs,
                    const RACES::Mutations::IDContext &rhs) const;
};

namespace std
{

/**
 * @brief Stream the ID context in a stream
 *
 * Indel contexts are conventionally represented by a string in the
 * form `{number}{'A','C','G','T','R','M'}{number}`. The first number
 * is the unit size for both homopolymer and heteropolymer, and the
 * fragment size for microhomology. When the string represents a
 * heteropolymer or a microhomology it is the first level code. The
 * character denotes the fragment type: 'R' stands for heteropolymer,
 * 'M' for microhomology, and 'A', 'C', 'G', and 'T' is the unit base
 * of a homopolymer. Dealing with homopolymer, the character is first
 * level code. The last number is the size of the microhomology (for
 * microhomologies) and the number of repetitions in the sequence
 * (for polymer).
 *
 * @param out is the output stream
 * @param type is the ID context to stream
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const RACES::Mutations::IDContext& type);

/**
 * @brief Stream the ID context from a stream
 *
 * Indel contexts are conventionally represented by a string in the
 * form `{number}{'A','C','G','T','R','M'}{number}`. The first number
 * is the unit size for both homopolymer and heteropolymer, and the
 * fragment size for microhomology. When the string represents a
 * heteropolymer or a microhomology it is the first level code. The
 * character denotes the fragment type: 'R' stands for heteropolymer,
 * 'M' for microhomology, and 'A', 'C', 'G', and 'T' is the unit base
 * of a homopolymer. Dealing with homopolymer, the character is first
 * level code. The last number is the size of the microhomology (for
 * microhomologies) and the number of repetitions in the sequence
 * (for polymer).
 *
 * @param in is the input stream
 * @param type is the object where the streamed ID context will be placed
 * @return a reference to the input stream
 */
std::istream& operator>>(std::istream& in, RACES::Mutations::IDContext& type);

}   // std

#endif // __RACES_ID_CONTEXT__