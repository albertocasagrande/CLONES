/**
 * @file sbs_signature.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Some SBS example
 * @version 1.0
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
#define BOOST_TEST_MODULE sbs_signature

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>   // std::transform
#include <sstream>
#include <set>

#include <boost/test/unit_test.hpp>

#include "sbs_signature.hpp"


BOOST_AUTO_TEST_CASE(sbs_context_creation)
{
    using namespace CLONES::Mutations;

    BOOST_CHECK_NO_THROW(SBSContext());
    BOOST_CHECK_NO_THROW(SBSContext("AAA"));
    BOOST_CHECK_NO_THROW(SBSContext("ACA"));
    BOOST_CHECK_THROW(SBSContext("AA"), std::domain_error);
    BOOST_CHECK_THROW(SBSContext("AAAA"), std::domain_error);
    BOOST_CHECK_THROW(SBSContext("AZA"), std::domain_error);
}

struct ContextFixture
{
    std::list<CLONES::Mutations::SBSContext> contexts;

    ContextFixture()
    {
        std::vector<char> bases{'A', 'C', 'G', 'T'};

        for (const auto& seq0: bases) {
            for (const auto& seq1: bases) {
                for (const auto& seq2: bases) {
                    std::string seq{seq0,seq1,seq2};

                    contexts.emplace_back(seq);
                }
            }
        }
    }
};

BOOST_FIXTURE_TEST_SUITE( context_test, ContextFixture )

BOOST_AUTO_TEST_CASE(sbs_context_lower_case)
{
    using namespace CLONES::Mutations;
    std::vector<char> bases{'a', 'c', 'g', 't'};

    auto it=contexts.cbegin();

    for (const auto& seq0: bases) {
        for (const auto& seq1: bases) {
            for (const auto& seq2: bases) {
                std::string seq{seq0,seq1,seq2};
                SBSContext mc(seq);

                BOOST_CHECK_EQUAL(mc, *(it++));
            }
        }
    }

}

BOOST_AUTO_TEST_CASE(sbs_context_relation)
{
    using namespace CLONES::Mutations;

    for (const auto& context: contexts) {
        BOOST_CHECK_EQUAL(context, context);
    }

    for (auto it=contexts.cbegin(); it!=contexts.cend(); ++it) {
        for (auto first_it=contexts.cbegin(); first_it!=it; ++first_it) {
            BOOST_CHECK_NE(*first_it, *it);
            BOOST_CHECK_NE(*it, *first_it);
        }
        auto second_it=it;
        for (++second_it; second_it!=contexts.cend(); ++second_it) {
            BOOST_CHECK_NE(*second_it, *it);
            BOOST_CHECK_NE(*it, *second_it);
        }
    }
}

BOOST_AUTO_TEST_CASE(sbs_context_copy_by_sequence)
{
    using namespace CLONES::Mutations;

    for (const auto& context: contexts) {
        SBSContext copy(context.get_sequence());

        BOOST_CHECK_EQUAL(context, copy);
    }
}

BOOST_AUTO_TEST_CASE(sbs_context_sequence)
{
    using namespace CLONES::Mutations;

    for (const auto& context: contexts) {
        SBSContext copy(context.get_sequence());

        BOOST_CHECK_EQUAL(context.get_sequence(), copy.get_sequence());
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(sbs_context_complement)
{
    using namespace CLONES::Mutations;

    std::list<std::pair<std::string, std::string>> tests{
            {"AAA","TTT"},
            {"CAG","GTC"},
            {"TGG","ACC"},
            {"ACA","TGT"},
            {"GCC","CGG"},
        };

    for (const auto& [seq, c_seq]: tests) {
        {
            SBSContext context(seq);

            auto c_context = context.get_complement();

            BOOST_CHECK_EQUAL(c_context.get_sequence(), c_seq);
        }

        {
            SBSContext c_context(c_seq);

            auto context = c_context.get_complement();

            BOOST_CHECK_EQUAL(context.get_sequence(), seq);
        }
    }
}

BOOST_AUTO_TEST_CASE(sbs_context_reverse_complement)
{
    using namespace CLONES::Mutations;

    std::list<std::pair<std::string, std::string>> tests{
            {"AAA","TTT"},
            {"CAG","CTG"},
            {"TGG","CCA"},
            {"ACA","TGT"},
            {"GCC","GGC"},
        };

    for (const auto& [seq, c_seq]: tests) {
        {
            SBSContext context(seq);

            auto c_context = context.get_reverse_complement();

            BOOST_CHECK_EQUAL(c_context.get_sequence(), c_seq);
        }

        {
            SBSContext c_context(c_seq);

            auto context = c_context.get_reverse_complement();

            BOOST_CHECK_EQUAL(context.get_sequence(), seq);
        }
    }
}

BOOST_AUTO_TEST_CASE(SBS_type_create)
{
    using namespace CLONES::Mutations;

    BOOST_CHECK_NO_THROW(SBSType());

    for (const auto base: {'A', 'C', 'G', 'T'}){
        for (const auto seq0: {'A', 'C', 'G', 'T'}){
            for (const auto seq1: {'A', 'C', 'G', 'T'}){
                for (const auto seq2: {'A', 'C', 'G', 'T'}){
                    std::string seq{seq0,seq1,seq2};

                    SBSContext context(seq);

                    if (seq1!=base) {
                        BOOST_CHECK_NO_THROW(SBSType(seq, base));
                        BOOST_CHECK_NO_THROW(SBSType(context, base));
                    } else {
                        BOOST_CHECK_THROW(SBSType(seq, base), std::domain_error);
                        BOOST_CHECK_THROW(SBSType(context, base), std::domain_error);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(SBS_type_get_replace_base)
{
    using namespace CLONES::Mutations;

    for (const auto base: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
        for (const auto seq1: {'C', 'T', 'c', 't'}){
            if (seq1!=base) {
                for (const auto seq0: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
                    for (const auto seq2: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
                        std::string seq{seq0,seq1,seq2};
                        SBSType type(seq, base);

                        BOOST_CHECK_EQUAL(type.get_replace_base(), toupper(base));
                    }
                }
            }
        }
    }

    for (const auto base: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
        for (const auto seq1: {'A', 'G', 'a', 'g'}){
            if (seq1!=base) {
                for (const auto seq0: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
                    for (const auto seq2: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
                        std::string seq{seq0,seq1,seq2};
                        SBSType type(seq, base);

                        char c_base;
                        switch(base) {
                            case 'A':
                            case 'a':
                                c_base = 'T';
                                break;
                            case 'G':
                            case 'g':
                                c_base = 'C';
                                break;
                            case 'C':
                            case 'c':
                                c_base = 'G';
                                break;
                            default:
                            case 'T':
                            case 't':
                                c_base = 'A';
                                break;
                        }
                        BOOST_CHECK_EQUAL(type.get_replace_base(), c_base);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(SBS_type_get_complement_replace_base)
{
    using namespace CLONES::Mutations;

    for (const auto base: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
        for (const auto seq1: {'C', 'T', 'c', 't'}){
            if (seq1!=base) {
                for (const auto seq0: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
                    for (const auto seq2: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
                        std::string seq{seq0,seq1,seq2};
                        SBSType type(seq, base);

                        char c_base;
                        switch(base) {
                            case 'A':
                            case 'a':
                                c_base = 'T';
                                break;
                            case 'G':
                            case 'g':
                                c_base = 'C';
                                break;
                            case 'C':
                            case 'c':
                                c_base = 'G';
                                break;
                            default:
                            case 'T':
                            case 't':
                                c_base = 'A';
                                break;
                        }
                        BOOST_CHECK_EQUAL(type.get_complement_replace_base(), c_base);
                    }
                }
            }
        }
    }

    for (const auto base: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
        for (const auto seq1: {'A', 'G', 'a', 'g'}){
            if (seq1!=base) {
                for (const auto seq0: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
                    for (const auto seq2: {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}){
                        std::string seq{seq0,seq1,seq2};
                        SBSType type(seq, base);

                        BOOST_CHECK_EQUAL(type.get_complement_replace_base(), toupper(base));
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(SBS_type_read)
{
    using namespace CLONES::Mutations;

    std::list<std::pair<std::string, std::pair<std::string, char>>> tests{
        {"A[A>C]C",{"GTT", 'G'}},
        {"a[G>T]G",{"CCT", 'A'}}
    };

    for (const auto& [input, results]: tests) {
        std::istringstream is(input);
        SBSType type;

        is >> type;

        BOOST_CHECK_EQUAL(type.get_context().get_sequence(), results.first);
        BOOST_CHECK_EQUAL(type.get_replace_base(), results.second);
    }
}

BOOST_AUTO_TEST_CASE(SBS_type_read_error)
{
    using namespace CLONES::Mutations;

    std::list<std::string> errors{
        "A[A>A]T", "A", "A[A<C]T", "A<A>C]T"
    };

    for (const auto& error: errors) {
        std::istringstream is(error);
        SBSType type;

        BOOST_CHECK_THROW(is >> type, std::runtime_error);
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

BOOST_AUTO_TEST_CASE(SBS_signature_load)
{
    using namespace CLONES::Mutations;

    std::set<std::string> signature_names{"SBS3_GRCh37","SBS3_GRCh38","SBS3_mm9","SBS3_mm10","SBS3_rn6"};

    std::map<std::string, SBSSignature> example;
    {
        std::ifstream in(SBS_EXAMPLE, std::ios_base::in);
        BOOST_CHECK_NO_THROW(example = SBSSignature::read_from_stream(in));
    }

    std::set<std::string> example_set;

    std::transform(example.begin(), example.end(),  std::inserter(example_set, example_set.end()),
                   [](auto pair){ return pair.first; });

    BOOST_CHECK_EQUAL(signature_names,example_set);
}

BOOST_AUTO_TEST_CASE(selective_SBS_signature_load)
{
    using namespace CLONES::Mutations;

    std::set<std::string> signature_names{"SBS3_GRCh38","SBS3_mm10","SBS3_rn6"};

    std::map<std::string, SBSSignature> example;
    {
        std::ifstream in(SBS_EXAMPLE, std::ios_base::in);
        BOOST_CHECK_NO_THROW(example = SBSSignature::read_from_stream(in, signature_names));
    }

    std::set<std::string> example_set;
    std::transform(example.begin(), example.end(),  std::inserter(example_set, example_set.end()),
                   [](auto pair){ return pair.first; });

    BOOST_CHECK_EQUAL(signature_names,example_set);
}

BOOST_AUTO_TEST_CASE(SBS_signature_expression)
{
    using namespace CLONES::Mutations;

    std::map<std::string, SBSSignature> signatures;
    {
        std::ifstream in(SBS_EXAMPLE, std::ios_base::in);
        signatures = SBSSignature::read_from_stream(in);
    }

    double alpha = 1.0/signatures.size();

    SignatureExprResult<SBSType> expr_result;
    for (const auto& [key, signature]: signatures) {
        expr_result = expr_result + alpha * signature;
    }

    SBSSignature test = static_cast<SBSSignature>(expr_result);
}