/**
 * @file id_signature.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements SBS signature
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

#include <limits>

#include "id_signature.hpp"


namespace CLONES
{

namespace Mutations
{

IDType::IDType():
    ftype{FragmentType::HOMOPOLYMER},
    fl_index{0}, sl_index{0}, insertion{true}
{}

IDType::IDType(const FragmentType& fragment_type, const uint8_t& first_level_index,
               const RepetitionType& second_level_index, const bool& insertion):
    ftype{fragment_type}, fl_index{first_level_index}, sl_index{second_level_index},
    insertion{insertion}
{}

IDType::IDType(const std::string& type):
    insertion{true}
{
    std::istringstream istype{type};
    std::vector<std::string> fields;
    std::string field;

    if (type.back()==':') {
        throw std::domain_error("\"" + type + "\" does not represent an ID type: "
                                + "it should contain 4 field separated by ':'.");
    }

    while (std::getline(istype, field, ':'))
    {
        fields.push_back(field);

        if (fields.size()>4) {
            throw std::domain_error("\"" + type + "\" does not represent an ID type: "
                                    + "it should contain 4 field separated by ':'.");
        }
    }

    if (fields.size()<4) {
        throw std::domain_error("\"" + type + "\" does not represent an ID type: "
                                + "it should contain 4 field separated by ':'.");
    }

    if (fields[2].size() != 1) {
        throw std::domain_error("\"" + type + "\" does not represent an ID type: "
                                + "\"" +fields[2] + "\" should be a character among "
                                + "'C', 'T', 'R', or 'M'.");
    }

    switch (fields[2][0]) {
        case 'C':
            fl_index = 'C';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'T':
            fl_index = 'T';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'R':
            fl_index = read_size<uint8_t>(fields[0], type);
            ftype = FragmentType::HETEROPOLYMER;
            break;
        case 'M':
            fl_index = read_size<uint8_t>(fields[0], type);
            ftype = FragmentType::MICROHOMOLOGY;
            break;
        default:
            throw std::domain_error("\""+type+"\" does not represent an ID type.");
    }

    if (ftype == FragmentType::MICROHOMOLOGY) {
        sl_index = read_size<uint8_t>(fields[3], type);
    } else {
        sl_index = read_size<RepetitionType>(fields[3], type);
    }

    if (fields[1] == "Del") {
        insertion = false;
        if (ftype != FragmentType::MICROHOMOLOGY) {
            ++sl_index;
        }
    } else if (fields[1] != "Ins") {
        throw std::domain_error("\"" + type + "\" does not represent an ID type: "
                                + "\"" + fields[1] + "\" should be either \"Ins\""
                                + "\"Del\".");
    }

}


}  // Mutations

}  // CLONES

namespace std
{

bool less<CLONES::Mutations::IDType>::operator()(const CLONES::Mutations::IDType &lhs,
                                                const CLONES::Mutations::IDType &rhs) const
{
    using namespace CLONES::Mutations;

    if (!lhs.insertion && rhs.insertion) {
        return true;
    }

    if (lhs.insertion && !rhs.insertion) {
        return false;
    }

    if (lhs.fl_index<rhs.fl_index) {
        return true;
    }

    if (lhs.fl_index>rhs.fl_index) {
        return false;
    }

    if (lhs.sl_index<rhs.sl_index) {
        return true;
    }

    if (lhs.sl_index>rhs.sl_index) {
        return false;
    }

    using Type = CLONES::Mutations::IDType::FragmentType;

    return (lhs.ftype == Type::HOMOPOLYMER && rhs.ftype != Type::HOMOPOLYMER)
            || (lhs.ftype == Type::HETEROPOLYMER && rhs.ftype == Type::MICROHOMOLOGY);
}


std::ostream& operator<<(std::ostream& out, const CLONES::Mutations::IDType& type)
{
    using namespace CLONES::Mutations;
    using FragmentType = IDType::FragmentType;

    if (type.ftype == FragmentType::HOMOPOLYMER) {
        out << "1:" << (type.insertion?"Ins":"Del") << ":"
            << (type.fl_index=='C'?"C":"T");
    } else {
        if (!type.insertion) {
            out << static_cast<unsigned int>(type.fl_index-1) << ":Del:";
        } else {
            out << static_cast<unsigned int>(type.fl_index) << ":Ins:";
        }

        out << (type.ftype == FragmentType::HETEROPOLYMER?"R":"M");
    }

    out << ":" << type.sl_index;

    return out;
}

std::istream& operator>>(std::istream& in, CLONES::Mutations::IDType& type)
{
    std::string str_type;

    in >> str_type;

    type = CLONES::Mutations::IDType(str_type);

    return in;
}

}   // std
