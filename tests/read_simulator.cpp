/**
 * @file read_simulator.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing read simulator statistics
 * @version 1.1
 * @date 2026-11-18
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
#define BOOST_TEST_MODULE read_simulator

#include <vector>
#include <string>
#include <memory>

#include <boost/test/unit_test.hpp>

#include "read_simulator.hpp"
#include "genome_fragment.hpp"

#include "error.hpp"

CLONES::Mutations::SID build_SNV(const CLONES::Mutations::ChrPosition position)
{
    using namespace CLONES::Mutations;

    return SID(1, position, 'A', 'T', Mutation::PASSENGER);
}

void add_SNVs(CLONES::Mutations::SequencingSimulations::ChrSampleStatistics& stats,
              const std::vector<CLONES::Mutations::ChrPosition>& positions)
{
    using namespace CLONES::Mutations::SequencingSimulations;

    for (const auto& position : positions) {
        const auto snv = build_SNV(position);

        stats.get_data().insert({snv, SIDData(snv, 1)});
    }
}

BOOST_AUTO_TEST_CASE(chr_statistics_merge)
{
    using namespace CLONES::Mutations;
    using namespace CLONES::Mutations::SequencingSimulations;

    // the SIDs of the two statistics interleave: some SIDs of the second
    // statistics precede, follow, or lay between those of the first one
    ChrSampleStatistics stats_a, stats_b;
    add_SNVs(stats_a, {100, 300});
    add_SNVs(stats_b, {50, 100, 200, 400});

    stats_a += stats_b;

    const auto& data = stats_a.get_data();
    BOOST_REQUIRE_EQUAL(data.size(), 5);

    std::vector<ChrPosition> positions;
    for (const auto& [mutation, SID_data] : data) {
        positions.push_back(mutation.position);
    }
    BOOST_CHECK(positions == std::vector<ChrPosition>({50, 100, 200, 300, 400}));

    // the occurrences of the common SID are summed
    BOOST_CHECK_EQUAL(data.at(build_SNV(100)).num_of_occurrences, 2);
    BOOST_CHECK_EQUAL(data.at(build_SNV(50)).num_of_occurrences, 1);
}

BOOST_AUTO_TEST_CASE(genome_fragment_at_chromosome_end)
{
    using namespace CLONES;
    using namespace CLONES::Mutations;

    // positions:                1234567890
    const std::string reference("ACGTACGTAC");

    // a deletion of the positions 6-8 anchored on position 5
    CLONES::map<GenomicPosition, std::shared_ptr<SID>> germline, somatic;
    somatic.emplace(GenomicPosition(1, 5),
                    std::make_shared<SID>(1, 5, "ACGT", "A", Mutation::PASSENGER));

    // a 6-base fragment beginning in position 3 would require the
    // reference base in position 11, which does not exist
    const GenomeFragment fragment(reference, germline, somatic, GenomicPosition(1, 3), 6);

    // the fragment is truncated at the reference end
    BOOST_CHECK_EQUAL(fragment.get_sequence(), "GTAAC");

    const auto covered = fragment.get_covered_reference_region();
    BOOST_CHECK_EQUAL(covered.begin(), 3);
    BOOST_CHECK_EQUAL(covered.end(), 10);
}
