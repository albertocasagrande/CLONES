/**
 * @file context_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing CLONES::Mutations::ContextIndex class
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
#define BOOST_TEST_MODULE context_index

#include <boost/test/unit_test.hpp>

#include <fstream>

#include "context_index.hpp"


BOOST_AUTO_TEST_CASE(context_index_creation)
{
    using namespace CLONES::Mutations;

    BOOST_CHECK_NO_THROW(ContextIndex());

    BOOST_CHECK_NO_THROW(ContextIndex<>::build_index(FASTA_FILE));

    std::set<GenomicRegion> regions{{{2,115}, 20},
                                    {{1,5}, 73},
                                    {{2,247}, 11}};

    BOOST_CHECK_NO_THROW(ContextIndex<>::build_index(FASTA_FILE, regions));

    BOOST_CHECK_THROW(ContextIndex<>::build_index("/TEST-ERROR"), std::runtime_error);
}

template<typename GENOME_WIDE_POSITION>
std::set<CLONES::Mutations::GenomicPosition> get_genomic_positions(const CLONES::Mutations::ContextIndex<GENOME_WIDE_POSITION>& context_index,
                                                                   const CLONES::Mutations::SBSContext& mutational_context)
{
    std::set<CLONES::Mutations::GenomicPosition> positions;

    for (const auto& abs_pos: context_index[mutational_context]) {
        positions.insert(context_index.get_genomic_position(abs_pos));
    }

    return positions;
}

struct ContextFixture
{
    using SBSContext = CLONES::Mutations::SBSContext;
    using GenomicPosition = CLONES::Mutations::GenomicPosition;

    std::map<SBSContext, std::set<GenomicPosition> > test_positions;

    ContextFixture():
        test_positions{
            {"ACT",{{1,77},{2,264},{3,6}}},
            {"GCG",{{1,31},{3,9}}},
            {"TCC",{{1,84},{2,296}}},
            {"TCT",{{1,62},{1,108},{2,164},{2,166}}},
            {"GCT",{{1,82},{2,128},{2,171},{2,294}}},
            {"TCG",{{2,126}}}
        }
    {}
};


BOOST_FIXTURE_TEST_SUITE( context_index_test, ContextFixture )

BOOST_AUTO_TEST_CASE(context_index_whole_genome)
{
    using namespace CLONES::Mutations;

    auto context_index = ContextIndex<>::build_index(FASTA_FILE);

    for (const auto& [context_test, positions_test]: test_positions) {
        std::set<CLONES::Mutations::GenomicPosition> positions;

        if (positions_test.size() != 0) {
            BOOST_CHECK_NO_THROW(positions = get_genomic_positions(context_index, context_test));
        } else {
            BOOST_CHECK_THROW(positions = get_genomic_positions(context_index, context_test), std::domain_error);
        }

        BOOST_CHECK_EQUAL(positions.size(), positions_test.size());

        auto it = positions.begin();
        auto it_test = positions_test.begin();
        for (; it != positions.end(); ++it, ++it_test) {
            BOOST_CHECK_EQUAL(*it, *it_test);
        }

    }
}

bool in_regions(const std::set<CLONES::Mutations::GenomicRegion>& genomic_regions,
               const CLONES::Mutations::GenomicPosition& genomic_position)
{
    for (const auto& genomic_region: genomic_regions) {
        if (genomic_region.contains(genomic_position)) {
            return true;
        }
    }

    return false;
}

BOOST_AUTO_TEST_CASE(context_index_regions)
{
    using namespace CLONES::Mutations;

    const std::set<GenomicRegion> regions{{{2,115}, 20}, {{1,5}, 73},
                                          {{2,247}, 11}};

    decltype(test_positions) in_context_positions;

    for (const auto& [context_test, positions_test]: test_positions) {
        for (const auto& position_test: positions_test) {
            if (!in_regions(regions, position_test)) {
                in_context_positions[context_test].insert(position_test);
            }
        }
    }

    auto context_index = ContextIndex<>::build_index(FASTA_FILE, regions);

    for (const auto& [context_test, positions_test]: in_context_positions) {
        std::set<CLONES::Mutations::GenomicPosition> positions;

        if (positions_test.size() != 0) {
            BOOST_CHECK_NO_THROW(positions = get_genomic_positions(context_index, context_test));
        } else {
            BOOST_CHECK_THROW(positions = get_genomic_positions(context_index, context_test), std::domain_error);
        }
        BOOST_CHECK_EQUAL(positions.size(), positions_test.size());

        auto it = positions.begin();
        auto it_test = positions_test.begin();
        for (; it != positions.end(); ++it, ++it_test) {
            BOOST_CHECK_EQUAL(*it, *it_test);
        }

    }
}

BOOST_AUTO_TEST_CASE(context_index_remove_insert)
{
    using namespace CLONES::Mutations;

    auto context_index = ContextIndex<>::build_index(FASTA_FILE);

    auto context = "CCT";

    BOOST_CHECK_EQUAL(context_index[context].size(), 8);

    std::vector<uint32_t> expected{8,14,20,26,38,67,88,372};

    for (size_t i=0; i<8; ++i) {
        BOOST_CHECK_EQUAL(context_index[context][i], expected[i]);
    }

    auto value = context_index.extract(context, 3);

    BOOST_CHECK_EQUAL(value, expected[3]);
    BOOST_CHECK_EQUAL(context_index[context].size(), 7);

    context_index.insert(context, value, 3);

    BOOST_CHECK_EQUAL(context_index[context].size(), 8);

    for (size_t i=0; i<8; ++i) {
        BOOST_CHECK_EQUAL(context_index[context][i], expected[i]);
    }
}

BOOST_AUTO_TEST_SUITE_END()
