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
    IDContext{}
{}

IDType::IDType(const FragmentType& fragment_type, const FirstLevelType& first_level_code,
               const SecondLevelType& second_level_code, const bool& insertion):
    IDContext{fragment_type, first_level_code, second_level_code},
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
                                + "\"" + fields[2] + "\" should be a character among "
                                + "'A', 'C', 'G', 'T', 'R', or 'M'.");
    }

    FragmentType ftype;
    FirstLevelType fl_code;
    SecondLevelType sl_code;

    switch (fields[2][0]) {
        case 'A':
            fl_code = 'A';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'C':
            fl_code = 'C';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'G':
            fl_code = 'G';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'T':
            fl_code = 'T';
            ftype = FragmentType::HOMOPOLYMER;
            break;
        case 'R':
            fl_code = read_size<uint8_t>(fields[0], type);
            ftype = FragmentType::HETEROPOLYMER;
            break;
        case 'M':
            fl_code = read_size<uint8_t>(fields[0], type);
            ftype = FragmentType::MICROHOMOLOGY;
            break;
        default:
            throw std::domain_error("\""+type+"\" does not represent an ID type.");
    }

    sl_code = read_size<SecondLevelType>(fields[3], type);

    if (fields[1] == "Del") {
        insertion = false;
        if (ftype != FragmentType::MICROHOMOLOGY) {
            ++sl_code;
        }
    } else if (fields[1] != "Ins") {
        throw std::domain_error("\"" + type + "\" does not represent an ID type: "
                                + "\"" + fields[1] + "\" should be either \"Ins\""
                                + "\"Del\".");
    }

    static_cast<IDContext &>(*this) = IDContext{ftype, fl_code, sl_code};
}


}  // Mutations

}  // CLONES

namespace std
{

bool less<CLONES::Mutations::IDType>::operator()(const CLONES::Mutations::IDType &lhs,
                                                const CLONES::Mutations::IDType &rhs) const
{
    using namespace CLONES::Mutations;

    if (lhs.is_deletion() && rhs.is_insertion()) {
        return true;
    }

    if (lhs.is_insertion() && rhs.is_deletion()) {
        return false;
    }

    std::less<CLONES::Mutations::IDContext> less_context;

    return less_context(lhs, rhs);
}


std::ostream& operator<<(std::ostream& out, const CLONES::Mutations::IDType& type)
{
    using namespace CLONES::Mutations;
    using FragmentType = IDType::FragmentType;

    if (type.fragment_type() == FragmentType::HOMOPOLYMER) {
        out << "1" << (type.is_insertion()?":Ins:":":Del:")
            << type.unit_base();
    } else {
        out << static_cast<unsigned int>(type.unit_size())
            << (type.is_insertion()?":Ins:":":Del:")
            << (type.fragment_type() == FragmentType::HETEROPOLYMER?"R":"M");
    }

    out << ":" << type.get_second_level_code();

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
