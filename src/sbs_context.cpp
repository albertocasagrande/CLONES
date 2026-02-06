/**
 * @file sbs_context.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements SBS contexts and extended context automata
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

#include <sstream>
#include <cstring>

#include "sbs_context.hpp"

#include "archive.hpp"

namespace CLONES
{

namespace Mutations
{


char decode_base(const uint8_t& code)
{
    switch(code) {
        case 0:
            return 'C';
        case 1:
            return 'T';
        case 2:
            return 'G';
        case 3:
            return 'A';
        default:
            throw std::domain_error("Unknown code");
    }
}

uint8_t encode_base(const char& base)
{
    switch(base) {
        case 'C':
        case 'c':
            return 0;
        case 'T':
        case 't':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'A':
        case 'a':
            return 3;
        default:
            throw std::domain_error("Unknown base");
    }
}

SBSContext::SBSContext():
    code(std::numeric_limits<CodeType>::max())
{}

SBSContext::SBSContext(const char* nucleic_triplet):
    code(0)
{
    if (nucleic_triplet[0]=='\0' || nucleic_triplet[1]=='\0'
            || nucleic_triplet[2]=='\0' || nucleic_triplet[3]!='\0') {

        throw std::domain_error("Expected a nucleic triplet");
    }

    // this is to have the central nucleotide associated
    // to the most significant bits
    const std::array<uint8_t, 3> shifts{0,4,2};
    for (unsigned int i=0; i<3; ++i) {
        code = code | (encode_base(nucleic_triplet[i]) << shifts[i]);
    }
}

SBSContext::SBSContext(const std::string& nucleic_triplet):
    code(0)
{
    if (nucleic_triplet.size()!=3) {
        throw std::domain_error("Expected a nucleic triplet. Got \""
                                + nucleic_triplet +"\"");
    }

    // this is to have the central nucleotide associated
    // to the most significant bits
    const std::array<uint8_t, 3> shifts{0,4,2};
    for (unsigned int i=0; i<3; ++i) {
        code = code | (encode_base(nucleic_triplet[i]) << shifts[i]);
    }
}

std::string SBSContext::get_sequence() const
{
    std::string sequence;

    if (!is_defined()) {
        return "???";
    }

    const uint8_t nucleotide_mask = 0x03;

    // this is to have the central nucleotide associated
    // to the most significant bits
    const std::array<uint8_t, 3> shifts{0,4,2};
    for (unsigned int i=0; i<3; ++i) {
        auto base = decode_base((code >> shifts[i])&nucleotide_mask);
        sequence.push_back(static_cast<char>(base));
    }

    return sequence;
}

char SBSContext::get_central_nucleotide() const
{
    if (!is_defined()) {
        return '?';
    }

    return static_cast<char>(decode_base((code >> 4)&0x03));
}

uint8_t SBSContext::get_complement(const uint8_t& code)
{
    if (code == std::numeric_limits<CodeType>::max()) {
        return code;
    }

    const uint8_t nucleotide_mask = 0x03;

    // this is to have the central nucleotide associated
    // to the most significant bits
    const std::array<uint8_t, 3> shifts{0,4,2};

    uint8_t complementary_code{0};
    for (unsigned int i=0; i<3; ++i) {
        auto base_code = (code >> shifts[i])&nucleotide_mask;

        // Complement base code by using this tricky feature
        // of base encoding
        base_code =  (base_code + 2) % 4;

        complementary_code = complementary_code | (base_code << shifts[i]);
    }

    return complementary_code;
}

uint8_t SBSContext::get_reverse_complement(const uint8_t& code)
{
    if (code == std::numeric_limits<CodeType>::max()) {
        return code;
    }

    const uint8_t nucleotide_mask = 0x03;

    // this is to have the central nucleotide associated
    // to the most significant bits
    const std::array<uint8_t, 3> shifts{0,4,2};

    uint8_t revcomp_code{0};
    for (unsigned int i=0; i<3; ++i) {
        auto base_code = (code >> shifts[2-i])&nucleotide_mask;

        // Complement base code by using this tricky feature
        // of base encoding
        base_code =  (base_code + 2) % 4;

        revcomp_code = revcomp_code | (base_code << shifts[i]);
    }

    return revcomp_code;
}

SBSContext SBSContext::get_complement() const
{
    SBSContext complement;

    complement.code = get_complement(code);

    return complement;
}


SBSContext SBSContext::get_reverse_complement() const
{
    SBSContext complement;

    complement.code = get_reverse_complement(code);

    return complement;
}

ExtendedContextAutomaton::BaseCodeType ExtendedContextAutomaton::base2code(const char& character)
{
    switch(character) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        case 'N':
        case 'n':
            return 4;
        default:
            return 5;
    }
}

SBSContext::CodeType ExtendedContextAutomaton::get_state_for(const char& first, const char& second, const char& third)
{
    return 5*(5*base2code(first)+base2code(second))+base2code(third);
}

ExtendedContextAutomaton::ExtendedContextAutomaton():
    state(get_state_for('N','N','N'))
{
    std::string context{"NNN"};
    const std::vector<char> alphabet{'A', 'C', 'G', 'T', 'N'};

    for (const auto& nucleotide2: alphabet) {
        context[2] = nucleotide2;
        for (const auto& nucleotide1: alphabet) {
            context[1] = nucleotide1;
            for (const auto& nucleotide0: alphabet) {
                context[0] = nucleotide0;
                size_t state = get_state_for(nucleotide0, nucleotide1, nucleotide2);

                if (context[0] != 'N' && context[1] != 'N' && context[2] != 'N') {
                    is_a_context[state] = true;

                    contexts[state] = SBSContext(context);
                } else {
                    is_a_context[state] = false;
                }
                auto& out_edges = edges[state];

                for (const auto& read_base: alphabet) {
                    out_edges[base2code(read_base)] = get_state_for(nucleotide1, nucleotide2, read_base);
                }
            }
        }
    }
}

bool ExtendedContextAutomaton::update_state(const char& base)
{
    uint8_t base_code = base2code(base);

    if (base_code<5) {
        state = edges[state][base_code];

        return true;
    }

    return false;
}

ExtendedContextAutomaton& ExtendedContextAutomaton::reset()
{
    state = get_state_for('N','N','N');

    return *this;
}

}   // Mutations

}   // CLONES

