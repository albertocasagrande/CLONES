/**
 * @file id_context.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements ID context
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

#include "id_context.hpp"

#include "genomic_sequence.hpp"

namespace RACES
{

namespace Mutations
{

IDContext::IDContext():
    ftype{FragmentType::HOMOPOLYMER},
    fl_code{0}, sl_code{0}
{}

IDContext::IDContext(const FragmentType& fragment_type,
                     const FirstLevelType& first_level_code,
                     const SecondLevelType& second_level_code):
    ftype{fragment_type}, fl_code{first_level_code}, sl_code{second_level_code}
{
    if (fragment_type == FragmentType::HOMOPOLYMER) {
        if (!GenomicSequence::is_a_DNA_base(static_cast<const char&>(first_level_code))) {
            throw std::runtime_error("IDContext: Unknown base '" +
                    std::to_string(static_cast<const char&>(first_level_code)) + "'.");
        }
    }

}

IDContext::IDContext(const std::string& context)
{
    std::istringstream istype{context};
    char c;
    uint8_t num1, num2;

    istype >> num1 >> c >> num2;

    if (!istype) {
        throw std::domain_error("\"" + context + "\" does not represent an ID context: "
                                + "It does not have the form {number}{charater}{number} "
                                + "with {character} in {'A','C','G','T','M','R'}.");
    }

    switch (c) {
        case 'A':
        case 'a':
            fl_code = 'A';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'C':
        case 'c':
            fl_code = 'C';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'G':
        case 'g':
            fl_code = 'G';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'T':
        case 't':
            fl_code = 'T';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'R':
            fl_code = num1;
            ftype = FragmentType::HETEROPOLYMER;
            break;
        case 'M':
            fl_code = num1;
            ftype = FragmentType::MICROHOMOLOGY;
            break;
        default:
            throw std::domain_error("\"" + context + "\" does not represent an ID context.");
    }

    if (ftype == FragmentType::MICROHOMOLOGY) {
        sl_code = num2;
    } else {
        sl_code = num2;
    }
}

char IDContext::unit_base() const
{
    if (ftype != FragmentType::HOMOPOLYMER) {
        throw std::runtime_error("IDContext::unit_base(): \"this\" is not "
                                 "a homopolymer.");
    }

    return static_cast<char>(get_first_level_code());
}

}  // Mutations

}  // RACES

namespace std
{

bool less<RACES::Mutations::IDContext>::operator()(const RACES::Mutations::IDContext &lhs,
                                                   const RACES::Mutations::IDContext &rhs) const
{
    if (lhs.fragment_type() > rhs.fragment_type()) {
        return false;
    }

    if (lhs.fragment_type() < rhs.fragment_type()) {
        return true;
    }

    if (lhs.get_first_level_code() > rhs.get_first_level_code()) {
        return false;
    }

    if (lhs.get_first_level_code() < rhs.get_first_level_code()) {
        return true;
    }

    return lhs.get_second_level_code() < rhs.get_second_level_code();
}


std::ostream& operator<<(std::ostream& out, const RACES::Mutations::IDContext& context)
{
    using namespace RACES::Mutations;
    using FragmentType = IDContext::FragmentType;

    if (context.fragment_type() == FragmentType::HOMOPOLYMER) {
        out << "1" << context.unit_base();
    } else {
        out << static_cast<unsigned int>(context.unit_size())
            << (context.fragment_type() == FragmentType::HETEROPOLYMER?"R":"M");
    }

    out << static_cast<unsigned int>(context.get_second_level_code());

    return out;
}

std::istream& operator>>(std::istream& in, RACES::Mutations::IDContext& type)
{
    std::string str_type;

    in >> str_type;

    type = RACES::Mutations::IDContext(str_type);

    return in;
}

}   // std
