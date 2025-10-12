/**
 * @file context_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing RACES::Mutations::ContextIndex class
 * @version 1.0
 * @date 2025-10-12
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE context_index

#include <boost/test/unit_test.hpp>

#include <filesystem>

#include "sbs_context_index.hpp"

BOOST_AUTO_TEST_CASE(context_index_creation)
{
    using namespace RACES::Mutations;

    BOOST_CHECK_NO_THROW(SBSContextIndex());

    std::mt19937_64 random_generator(0);

    const auto index_path = get_a_temporary_path("sbs_context_index_test",
                                                 std::filesystem::temp_directory_path());

    BOOST_CHECK_NO_THROW(SBSContextIndex<>::build(random_generator, index_path,
                                                  FASTA_FILE));

    std::filesystem::remove_all(index_path);

    std::set<GenomicRegion> regions{{{2,115}, 20},
                                    {{1,5}, 73},
                                    {{2,247}, 11}};

    BOOST_CHECK_NO_THROW(SBSContextIndex<>::build(random_generator, index_path,
                                                  FASTA_FILE, regions));

    BOOST_CHECK_THROW(SBSContextIndex<>::build(random_generator, "/TEST-ERROR",
                                               FASTA_FILE), std::runtime_error);
}

template<typename RANDOM_GENERATOR>
std::set<RACES::Mutations::GenomicPosition> get_genomic_positions(const RACES::Mutations::SBSContextIndex<RANDOM_GENERATOR>& context_index,
                                                                   const RACES::Mutations::SBSContext& mutational_context)
{
    std::set<RACES::Mutations::GenomicPosition> positions;

    for (const auto& g_pos: context_index[mutational_context]) {
        positions.insert(g_pos);
    }

    return positions;
}

struct ContextFixture
{
    using SBSContext = RACES::Mutations::SBSContext;
    using GenomicPosition = RACES::Mutations::GenomicPosition;

    std::map<SBSContext, std::set<GenomicPosition> > test_positions;

    ContextFixture():
        test_positions{
            {"ACT",{{1,76},{2,263},{3,5}}},
            {"GCG",{{1,30},{3,8}}},
            {"TCC",{{1,83},{2,295}}},
            {"TCT",{{1,61},{1,107},{2,163},{2,165}}},
            {"GCT",{{1,81},{2,127},{2,170},{2,293}}},
            {"TCG",{{2,125}}}
        }
    {}
};


BOOST_FIXTURE_TEST_SUITE( context_index_test, ContextFixture )

BOOST_AUTO_TEST_CASE(context_index_whole_genome)
{
    using namespace RACES::Mutations;

    std::mt19937_64 random_generator(0);

    const auto index_path = get_a_temporary_path("sbs_context_index_test",
                                                 std::filesystem::temp_directory_path());

    auto context_index = SBSContextIndex<>::build(random_generator, index_path,
                                                  FASTA_FILE);

    for (const auto& [context_test, positions_test]: test_positions) {
        std::set<RACES::Mutations::GenomicPosition> positions;

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

    std::filesystem::remove_all(index_path);
}

bool in_regions(const std::set<RACES::Mutations::GenomicRegion>& genomic_regions,
               const RACES::Mutations::GenomicPosition& genomic_position)
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
    using namespace RACES::Mutations;

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

    std::mt19937_64 random_generator(0);
    const auto index_path = get_a_temporary_path("sbs_context_index_test",
                                                 std::filesystem::temp_directory_path());

    auto context_index = SBSContextIndex<>::build(random_generator, index_path,
                                                  FASTA_FILE, regions);

    for (const auto& [context_test, positions_test]: in_context_positions) {
        std::set<RACES::Mutations::GenomicPosition> positions;

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

    std::filesystem::remove_all(index_path);
}

BOOST_AUTO_TEST_CASE(context_index_remove_insert)
{
    using namespace RACES::Mutations;

    std::mt19937_64 random_generator(0);
    const auto index_path = get_a_temporary_path("sbs_context_index_test",
                                                 std::filesystem::temp_directory_path());

    auto context_index = SBSContextIndex<>::build(random_generator, index_path,
                                                  FASTA_FILE);

    auto context = "CCT";

    auto& context_bucket = context_index[context];

    BOOST_CHECK_EQUAL(context_bucket.size(), 8);

    std::set<GenomicPosition> expected{{1, 7}, {1, 13}, {1, 19},
                                       {1, 25}, {1, 37}, {1, 66},
                                       {1, 87}, {2, 152}};
    for (size_t i=0; i<8; ++i) {
        BOOST_CHECK(expected.count(context_bucket[i])>0);
    }
    std::filesystem::remove_all(index_path);
}

BOOST_AUTO_TEST_SUITE_END()
