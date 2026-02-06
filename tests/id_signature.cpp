/**
 * @file id_signature.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Some ID example
 * @version 1.1
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE id_signature

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>   // std::transform
#include <sstream>
#include <set>

#include <boost/test/unit_test.hpp>

#include "id_signature.hpp"


BOOST_AUTO_TEST_CASE(ID_type_create)
{
    using namespace CLONES::Mutations;

    BOOST_CHECK_NO_THROW(IDType());

    for (const IDType::FragmentType ftype : {IDType::FragmentType::HETEROPOLYMER,
                                             IDType::FragmentType::MICROHOMOLOGY}) {
        for (const IDType::FirstLevelType first_level_code: {1, 2, 3, 4, 5}) {
            for (const IDType::SecondLevelType second_level_code: {0, 1, 2, 3, 4, 5}) {
                for (const bool insertion: {true, false}) {
                    BOOST_CHECK_NO_THROW(IDType(ftype, first_level_code,
                                                second_level_code, insertion));
                }
            }
        }
    }

    for (const IDType::FirstLevelType first_level_code: {'a', 'A', 'c', 'C',
                                                         'g', 'G', 't', 'T'}) {
        for (const IDType::SecondLevelType second_level_code: {0, 1, 2, 3, 4, 5}) {
            for (const bool insertion: {true, false}) {
                BOOST_CHECK_NO_THROW(IDType(IDType::FragmentType::HOMOPOLYMER,
                                            first_level_code, second_level_code,
                                            insertion));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(ID_type_read)
{
    using namespace CLONES::Mutations;

    struct IDTypeData
    {
        IDType::FragmentType ftype;
        IDType::FirstLevelType first_level_code;
        IDType::SecondLevelType second_level_code;
        bool insertion;
    };

    std::list<std::pair<std::string, IDTypeData>> tests{
        {"2:Del:R:0",{IDType::FragmentType::HETEROPOLYMER, 2, 1, false}},
        {"3:Ins:R:0",{IDType::FragmentType::HETEROPOLYMER, 3, 0, true}},
        {"1:Del:C:3",{IDType::FragmentType::HOMOPOLYMER, 'C', 4, false}},
        {"1:Del:T:3",{IDType::FragmentType::HOMOPOLYMER, 'T', 4, false}},
        {"1:Ins:C:3",{IDType::FragmentType::HOMOPOLYMER, 'C', 3, true}},
        {"1:Ins:T:3",{IDType::FragmentType::HOMOPOLYMER, 'T', 3, true}},
        {"3:Ins:R:0",{IDType::FragmentType::HETEROPOLYMER, 3, 0, true}},
        {"3:Del:R:1",{IDType::FragmentType::HETEROPOLYMER, 3, 2, false}},
        {"3:Del:M:1",{IDType::FragmentType::MICROHOMOLOGY, 3, 1, false}}
    };

    for (const auto& [input, results]: tests) {
        std::istringstream is(input);
        IDType type;

        is >> type;

        BOOST_CHECK(type.fragment_type()==results.ftype);
        BOOST_CHECK_EQUAL(type.get_first_level_code(), results.first_level_code);
        BOOST_CHECK_EQUAL(type.get_second_level_code(), results.second_level_code);
        BOOST_CHECK_EQUAL(type.is_insertion(), results.insertion);
    }
}

BOOST_AUTO_TEST_CASE(ID_type_read_error)
{
    using namespace CLONES::Mutations;

    std::list<std::string> errors{
        "2:Del:R:0:", "2:Dela:R:0", "-2:Del:R:0", "2:Del:R:-10",
        "2:Del:S:0", "2:Del:R:", "2:Del:R", "2:Del:R:0:A"
    };

    for (const auto& error: errors) {
        std::istringstream is(error);
        IDType type;

        BOOST_CHECK_THROW(is >> type, std::domain_error);
    }
}

namespace std
{

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::set<T>& S)
{
    out << "{";
    if (S.size()>0) {
        auto it = S.begin();
        out << *it;

        while (++it != S.end()) {
            out << "," << *it;
        }
    }

    out << "}";

    return out;
}

}

BOOST_AUTO_TEST_CASE(ID_signature_load)
{
    using namespace CLONES::Mutations;

    std::set<std::string> signature_names;

    for (size_t i=1; i<24; ++i) {
        std::ostringstream oss;

        oss << "ID" << i;
        signature_names.insert(oss.str());
    }

    std::map<std::string, IDSignature> example;
    {
        std::ifstream in(ID_EXAMPLE, std::ios_base::in);
        BOOST_CHECK_NO_THROW(example = IDSignature::read_from_stream(in));
    }

    std::set<std::string> example_set;
    std::transform(example.begin(), example.end(),  std::inserter(example_set, example_set.end()),
                   [](auto pair){ return pair.first; });

    BOOST_CHECK_EQUAL(signature_names,example_set);
}


BOOST_AUTO_TEST_CASE(selective_ID_signature_load)
{
    using namespace CLONES::Mutations;

    std::set<std::string> signature_names{"ID3","ID20","ID1"};

    std::map<std::string, IDSignature> example;
    {
        std::ifstream in(ID_EXAMPLE, std::ios_base::in);
        BOOST_CHECK_NO_THROW(example = IDSignature::read_from_stream(in, signature_names));
    }

    std::set<std::string> example_set;
    std::transform(example.begin(), example.end(),  std::inserter(example_set, example_set.end()),
                   [](auto pair){ return pair.first; });

    BOOST_CHECK_EQUAL(signature_names,example_set);
}

BOOST_AUTO_TEST_CASE(ID_signature_expression)
{
    using namespace CLONES::Mutations;

    std::map<std::string, IDSignature> signatures;
    {
        std::ifstream in(ID_EXAMPLE, std::ios_base::in);
        signatures = IDSignature::read_from_stream(in);
    }

    double alpha = 1.0/signatures.size();

    SignatureExprResult<IDType> expr_result;
    for (const auto& [key, signature]: signatures) {
        expr_result = expr_result + alpha * signature;
    }

    IDSignature test = static_cast<IDSignature>(expr_result);
}

