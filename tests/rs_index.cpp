/**
 * @file rs_index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing CLONES::Mutations::RSIndex class
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
#define BOOST_TEST_MODULE rs_index

#include <cctype>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <iterator>

#include <boost/test/unit_test.hpp>

#include "rs_index.hpp"
#include "id_signature.hpp"
#include "fasta_chr_reader.hpp"

#include "error.hpp"

// the names of the 83 COSMIC indel types
std::vector<std::string> get_COSMIC_ID_type_names()
{
    std::vector<std::string> names;

    for (const std::string type : {"Del", "Ins"}) {
        for (const std::string base : {"C", "T"}) {
            for (size_t i=0; i<6; ++i) {
                names.push_back("1:" + type + ":" + base + ":" + std::to_string(i));
            }
        }

        for (size_t length=2; length<6; ++length) {
            for (size_t i=0; i<6; ++i) {
                names.push_back(std::to_string(length) + ":" + type + ":R:"
                                + std::to_string(i));
            }
        }
    }

    for (size_t length=2; length<6; ++length) {
        for (size_t i=1; i<length; ++i) {
            names.push_back(std::to_string(length) + ":Del:M:" + std::to_string(i));
        }
    }
    names.push_back("5:Del:M:5");

    return names;
}

BOOST_AUTO_TEST_CASE(indel_anchors_match_the_reference)
{
    using namespace CLONES::Mutations;
    using namespace CLONES::IO::FASTA;

    // the test FASTA file contains unknown (N) bases
    std::map<ChromosomeId, std::string> sequences;
    {
        ChromosomeData<Sequence> chr;
        Reader<ChromosomeData<Sequence>> reader(FASTA_FILE);
        while (reader.read(chr)) {
            sequences[chr.chr_id] = chr.nucleotides;
        }
    }

    // no progress bar is provided
    RSIndex index;
    BOOST_REQUIRE_NO_THROW(index = RSIndex::build_index(FASTA_FILE, 50, 1000000, 0));

    const auto names = get_COSMIC_ID_type_names();
    BOOST_REQUIRE_EQUAL(names.size(), 83);

    size_t num_of_repetitions{0};
    for (const auto& name : names) {
        const IDType type(name);

        while (index.count_available_for(type) > 0) {
            const auto repetition = index.select(type, true);
            ++num_of_repetitions;

            // indels are anchored on the base preceding the repetition
            const auto& g_position = repetition.g_position;

            BOOST_TEST_CONTEXT("type " << name << " at " << g_position) {
                BOOST_REQUIRE_GT(g_position.position, 1);

                // genomic positions are 1-based
                const auto& sequence = sequences.at(g_position.chr_id);
                const char ref_base = static_cast<char>(std::toupper(sequence[g_position.position-2]));
                const char anchor = static_cast<char>(std::toupper(repetition.prev_base));

                BOOST_CHECK_NE(anchor, 'N');
                BOOST_CHECK_EQUAL(anchor, ref_base);
            }
        }
    }

    BOOST_CHECK_GT(num_of_repetitions, 0);
}

// an indel-affected span: the type name, the chromosome, and the
// first and the last positions from the anchor to the repetition end
struct RepetitionSpan
{
    std::string type;
    CLONES::Mutations::ChromosomeId chr_id;
    CLONES::Mutations::ChrPosition first;
    CLONES::Mutations::ChrPosition last;

    bool operator<(const RepetitionSpan& other) const
    {
        return std::tie(type, chr_id, first, last)
                < std::tie(other.type, other.chr_id, other.first, other.last);
    }

    bool operator==(const RepetitionSpan& other) const
    {
        return std::tie(type, chr_id, first, last)
                == std::tie(other.type, other.chr_id, other.first, other.last);
    }
};

std::set<RepetitionSpan> collect_spans(CLONES::Mutations::RSIndex& index)
{
    using namespace CLONES::Mutations;

    std::set<RepetitionSpan> spans;
    for (const auto& name : get_COSMIC_ID_type_names()) {
        const IDType type(name);

        while (index.count_available_for(type) > 0) {
            const auto repetition = index.select(type, true);
            const auto& begin = repetition.g_position.position;
            const size_t length = std::max(static_cast<size_t>(repetition.num_of_repetitions)
                                                *repetition.unit.size(),
                                           static_cast<size_t>(1));

            spans.insert({name, repetition.g_position.chr_id, begin-1,
                          static_cast<ChrPosition>(begin+length-1)});
        }
    }

    return spans;
}

BOOST_AUTO_TEST_CASE(regions_to_avoid)
{
    using namespace CLONES::Mutations;

    auto all_index = RSIndex::build_index(FASTA_FILE, 50, 1000000, 0);
    const auto all_spans = collect_spans(all_index);
    BOOST_REQUIRE_GT(all_spans.size(), 10);

    // avoid three single bases and a longer region taken from repetitions
    CLONES::set<GenomicRegion> regions;
    auto span_it = all_spans.begin();
    for (size_t i=0; i<3; ++i, std::advance(span_it, all_spans.size()/4)) {
        regions.emplace(GenomicPosition(span_it->chr_id, span_it->last), 1);
    }
    regions.emplace(GenomicPosition(span_it->chr_id, span_it->first), 10);

    auto overlaps = [&regions](const RepetitionSpan& span) {
        for (const auto& region : regions) {
            if (region.get_chromosome_id() == span.chr_id
                    && region.begin() <= span.last && span.first <= region.end()) {
                return true;
            }
        }

        return false;
    };

    std::set<RepetitionSpan> expected;
    for (const auto& span : all_spans) {
        if (!overlaps(span)) {
            expected.insert(span);
        }
    }
    BOOST_REQUIRE_LT(expected.size(), all_spans.size());

    RSIndex index;
    BOOST_REQUIRE_NO_THROW(index = RSIndex::build_index(FASTA_FILE, regions, 50, 1000000, 0));

    // exactly the repetitions not affecting the regions to avoid are stored
    const auto spans = collect_spans(index);
    BOOST_CHECK(spans == expected);
}
