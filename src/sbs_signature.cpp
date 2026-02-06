/**
 * @file sbs_signature.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements SBS signature
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

#include <array>
#include <cctype>  // toupper
#include <list>
#include <cmath>
#include <limits>

#include "sbs_signature.hpp"


namespace CLONES
{

namespace Mutations
{

char read_a_base(std::istream& in)
{
    char symbol;

    in >> symbol;
    switch(symbol) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return symbol;
        case 'a':
        case 'c':
        case 'g':
        case 't':
            return toupper(symbol);
        default:
        {
            std::ostringstream oss;

            oss << "'" << symbol << "' is not a base.";
            throw std::runtime_error(oss.str());
        }
    }
}

std::istream& read_symbol(std::istream& in, const char& symbol)
{
    char read_symbol;

    in >> read_symbol;

    if (read_symbol != symbol) {
        std::ostringstream oss;

        oss << " Expected '" << symbol << "'; got '" << read_symbol <<"'.";
        throw std::runtime_error(oss.str());
    }

    return in;
}

SBSType::SBSType():
    context(), replace_base('A')
{}

SBSType::SBSType(const SBSContext& context, const char& replace_base)
{
    auto central_nucleotide = context.get_central_nucleotide();

    if (central_nucleotide == replace_base) {
        std::ostringstream oss;

        oss << "Expected a replace base different from "
            << " the second nucleotide in the context. Got \""
            << context.get_sequence() +"\" and '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (!GenomicSequence::is_a_DNA_base(replace_base)) {
        std::ostringstream oss;

        oss << "Expected a replace base. Got '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (central_nucleotide != 'C' && central_nucleotide != 'T') {
        this->context = context.get_reverse_complement();
        this->replace_base = GenomicSequence::get_complement(replace_base);
    } else {
        this->context = context;
        this->replace_base = toupper(replace_base);
    }
}

SBSType::SBSType(const std::string& context, const char& replace_base):
    context(context), replace_base(toupper(replace_base))
{
    if (context[1] == replace_base) {
        std::ostringstream oss;

        oss << "Expected a replace base different from "
            << " the second nucleotide in the context. Got \""
            << context +"\" and '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (!GenomicSequence::is_a_DNA_base(replace_base)) {
        std::ostringstream oss;

        oss << "Expected a replace base. Got '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (context[1] != 'C' &&  context[1] != 'T' && context[1] != 'c' &&  context[1] != 't') {
        this->context = this->context.get_reverse_complement();
        this->replace_base = GenomicSequence::get_complement(replace_base);
    }
}

SBSType::SBSType(const std::string& type)
{
    std::istringstream iss(type);

    iss >> *this;
}

}  // Mutations

}  // CLONES

namespace std
{

bool less<CLONES::Mutations::SBSType>::operator()(const CLONES::Mutations::SBSType &lhs,
                                                         const CLONES::Mutations::SBSType &rhs) const
{
    const auto& lhs_code = lhs.get_context().get_code();
    const auto& rhs_code = rhs.get_context().get_code();

    return ((lhs_code < rhs_code) ||
            ((lhs_code == rhs_code) && (lhs.get_replace_base()<rhs.get_replace_base())));
}

std::ostream& operator<<(std::ostream& out, const CLONES::Mutations::SBSType& type)
{
    std::string type_sequence = type.get_context().get_sequence();

    out << type_sequence[0] << "["
        << type_sequence[1] << ">" << type.get_replace_base()
        << "]" << type_sequence[2];

    return out;
}

std::istream& operator>>(std::istream& in, CLONES::Mutations::SBSType& type)
{
    using namespace CLONES::Mutations;

    std::string seq;
    char replace_base;

    seq.push_back(read_a_base(in));
    read_symbol(in, '[');
    seq.push_back(read_a_base(in));
    read_symbol(in, '>');
    in >> replace_base;
    read_symbol(in, ']');
    seq.push_back(read_a_base(in));

    try {
        type = SBSType(seq, replace_base);
    } catch (std::domain_error& err) {
        throw std::runtime_error(err.what());
    }
    return in;
}

}   // std
