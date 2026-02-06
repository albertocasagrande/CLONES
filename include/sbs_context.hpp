/**
 * @file sbs_context.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines SBS contexts and extended context automata
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

#ifndef __CLONES_SBS_CONTEXT__
#define __CLONES_SBS_CONTEXT__

#include <cstdint>
#include <functional>   // std::less
#include <string>
#include <iostream>
#include <array>
#include <limits>

#include "archive.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief A class to represent SBS context
 *
 * A context is a nucleic triplet where a SBS may
 * occur.
 * Every context univocally corresponds to a code that
 * represents the triplet. This class guarantees the
 * SBS context code to be in the interval
 * natural [0,63]: 5', 3', and the central nucleotide
 * have 4 possible values, i.e., A, C, G, and T.
 * Moreover, all the codes for contexts having either
 * a 'C' and a 'T' as central nucleotide are guaranteed
 * to be in the interval [0,31].
 */
struct SBSContext
{
    /**
     * @brief The type of SBS context type
     */
    using CodeType = uint8_t;

private:
    CodeType code;   //!< the code associated to the context
public:
    /**
     * @brief The empty constructor
     */
    SBSContext();

    /**
     * @brief A constructor
     *
     * @param[in] nucleic_triplet is the string of the nucleic triplet
     */
    SBSContext(const char* nucleic_triplet);

    /**
     * @brief A constructor
     *
     * @param[in] nucleic_triplet is the string of the nucleic triplet
     */
    SBSContext(const std::string& nucleic_triplet);

    /**
     * @brief Get the SBS context code
     *
     * @return a constant reference to the SBS context code
     */
    inline const CodeType& get_code() const
    {
        return code;
    }

    /**
     * @brief Check whether the SBS context is defined
     *
     * @return `true` if and only if the nucleic triplet was
     *          specified during the object creation
     */
    inline bool is_defined() const
    {
        return (code != std::numeric_limits<CodeType>::max());
    }

    /**
     * @brief Get the nucleic sequence of the context
     *
     * @return the nucleic sequence of the context
     */
    std::string get_sequence() const;

    /**
     * @brief Get the context central nucleotide
     *
     * @return the context central nucleotide
     */
    char get_central_nucleotide() const;

    /**
     * @brief Get the complement SBS context
     *
     * @return the complement SBS context
     */
    SBSContext get_complement() const;

    /**
     * @brief Get the reverse complement SBS context
     *
     * @return the reverse complement SBS context
     */
    SBSContext get_reverse_complement() const;

    /**
     * @brief Test whether two SBS contexts are equivalent
     *
     * @param[in] context is the SBS context to compare
     * @return `true` if and only if the two SBS contexts represent
     *      the same nucleic triplet
     */
    inline bool operator==(const SBSContext& context) const
    {
        return get_code() == context.get_code();
    }

    /**
     * @brief Test whether two SBS contexts differ
     *
     * @param[in] context is the SBS context to compare
     * @return `true` if and only if the two SBS contexts represent
     *      different nucleic triplets
     */
    inline bool operator!=(const SBSContext& context) const
    {
        return get_code() != context.get_code();
    }

    /**
     * @brief Save a SBS context in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param[in,out] archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & code;
    }

    /**
     * @brief Load a SBS context from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param[in,out] archive is the input archive
     * @return the loaded SBS context
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SBSContext load(ARCHIVE& archive)
    {
        SBSContext mcode;

        archive & mcode.code;

        return mcode;
    }

    /**
     * @brief Get the code of the complement SBS context
     *
     * @param[in] code is the SBS context code whose complement is request
     * @return the code of the complement context of the
     *      context whose code is `code`
     */
    static CodeType get_complement(const CodeType& code);

    /**
     * @brief Get the code of the reverse complement SBS context
     *
     * @param[in] code is the SBS context code whose reverse complement
     *      is request
     * @return the code of the reverse complement context of the
     *      context whose code is `code`
     */
    static CodeType get_reverse_complement(const CodeType& code);
};


/**
 * @brief An automaton for parsing extended context
 *
 * An extended context is a triplet in the alphabet
 * 'A', 'C', 'G', 'T', and 'N'.
 * An object of this class is an automaton whose state
 * correspond to an extended context. The automaton state
 * always matches the extended context of the last read
 * nucleotides assuming that "NNN" is the initial
 * configuration.
 */
struct ExtendedContextAutomaton
{
    /**
     * @brief The type of code for bases
     */
    using BaseCodeType = uint8_t;

    /**
     * @brief The type of the edges from a node
     *
     * This type represents edges leaving a node: the key in the
     * map is the code of the newly read base; the value is the
     * new context code.
     */
    using FromNodeEdgeType = std::map<BaseCodeType, SBSContext::CodeType>;

    std::array<bool, 125> is_a_context;         //!< a Boolean vector to flag proper context
    std::array<SBSContext, 125> contexts;       //!< the state corresponding contexts
    std::array<FromNodeEdgeType, 125> edges;    //!< the automata edges

    SBSContext::CodeType state;  //!< the current automaton state

    /**
     * @brief Get a code for a character
     *
     * @param[in] character is a possible base
     * @return a number in the interval [0,5]. The
     *      value 5 is returned for all the characters
     *      different from  'A', 'C', 'G', 'T', and 'N'.
     */
    static BaseCodeType base2code(const char& character);

    /**
     * @brief Get the state corresponding to an extended context
     *
     * @param[in] first is the first character of the extended context
     * @param[in] second is the second character of the extended context
     * @param[in] third is the third character of the extended context
     * @return the automaton state corresponding to the extended
     *      context formed by the three characters.
     */
    static SBSContext::CodeType get_state_for(const char& first, const char& second,
                                              const char& third);

public:
    /**
     * @brief Build an extended context automaton
     */
    ExtendedContextAutomaton();

    /**
     * @brief Get the current automaton state
     *
     * @return return the current automaton state
     */
    inline const SBSContext::CodeType& get_state() const
    {
        return state;
    }

    /**
     * @brief Check whether the automaton state corresponds to a proper context
     *
     * @return `true` if and only if the last three nucleotides read by the
     *      automaton are among 'A', 'C', 'G', 'T'
     */
    inline const bool& read_a_context() const
    {
        return is_a_context[state];
    }

    /**
     * @brief Get the context corresponding to the current state
     *
     * @return the context corresponding to the current state
     */
    inline const SBSContext& get_context() const
    {
        return contexts[state];
    }

    /**
     * @brief Read a character and update the automaton state
     *
     * @param[in] character is the read character
     * @return `true` if and only if the read character is a valid
     *      nucleotide character, i.e., it is among 'A', 'C', 'G',
     *      'T', and 'N'.
     */
    bool update_state(const char& character);

    /**
     * @brief Reset the automaton state
     *
     * @return a reference to the updated automaton
     */
    ExtendedContextAutomaton& reset();

    /**
     * @brief Check whether the character is a nucleotide
     *
     * @param[in] character is a character
     * @return `true` if and only if the read character is a valid
     *      nucleotide character, i.e., it is among 'A', 'C', 'G',
     *      'T', and 'N'.
     */
    static inline bool is_a_nucleotide(const char& character)
    {
        return base2code(character)<5;
    }
};

}   // Mutations

}   // CLONES

template<>
struct std::less<CLONES::Mutations::SBSContext>
{
    inline bool operator()(const CLONES::Mutations::SBSContext &lhs,
                           const CLONES::Mutations::SBSContext &rhs) const
    {
        return lhs.get_code() < rhs.get_code();
    }
};

namespace std
{

/**
 * @brief Stream the SBS context in a stream
 *
 * @param[in,out] out is the output stream
 * @param[in] context is the SBS context to stream
 * @return a reference to the output stream
 */
inline std::ostream& operator<<(std::ostream& out, const CLONES::Mutations::SBSContext& context)
{
    return (out << context.get_sequence());
}

}   // std

#endif // __CLONES_SBS_CONTEXT__