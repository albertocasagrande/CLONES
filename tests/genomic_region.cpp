/**
 * @file genomic_region.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Genomic region tests
 * @version 1.0
 * @date 2026-11-20
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
#define BOOST_TEST_MODULE genomic_region

#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "genomic_position.hpp"
#include "genomic_region.hpp"

using namespace CLONES::Mutations;


BOOST_AUTO_TEST_SUITE(genomic_region)

BOOST_AUTO_TEST_CASE(constructors)
{
    GenomicRegion empty;
    BOOST_CHECK_EQUAL(empty.size(), 0);

    GenomicRegion from_chr(2, 10);
    BOOST_CHECK_EQUAL(from_chr.get_chromosome_id(), 2);
    BOOST_CHECK_EQUAL(from_chr.begin(), 1);
    BOOST_CHECK_EQUAL(from_chr.end(), 10);
    BOOST_CHECK_EQUAL(from_chr.size(), 10);
    BOOST_CHECK_NO_THROW(GenomicRegion(3, 100, 110));
    BOOST_CHECK_NO_THROW(GenomicRegion(GenomicPosition(3, 100), 10));
    BOOST_CHECK_NO_THROW(GenomicRegion(3, 100, 100));

    GenomicRegion from_pos(GenomicPosition(3, 100), 5);
    BOOST_CHECK_EQUAL(from_pos.get_chromosome_id(), 3);
    BOOST_CHECK_EQUAL(from_pos.begin(), 100);
    BOOST_CHECK_EQUAL(from_pos.end(), 104);
    BOOST_CHECK_EQUAL(from_pos.size(), 5);
    BOOST_CHECK(from_pos.get_initial_position() == GenomicPosition(3, 100));
    BOOST_CHECK(from_pos.get_final_position() == GenomicPosition(3, 104));

    BOOST_CHECK_THROW(GenomicRegion(2, 0), std::domain_error);
    BOOST_CHECK_THROW(GenomicRegion(GenomicPosition(2, 10), 0), std::domain_error);
    BOOST_CHECK_THROW(GenomicRegion(2, 10, 9), std::domain_error);
}

BOOST_AUTO_TEST_CASE(equality)
{
    BOOST_CHECK(GenomicRegion(1, 10, 20) == GenomicRegion(1, 10, 20));
    BOOST_CHECK(GenomicRegion(1, 10, 20) != GenomicRegion(1, 10, 21));
    BOOST_CHECK(GenomicRegion(1, 10, 20) != GenomicRegion(1, 11, 20));
    BOOST_CHECK(GenomicRegion(1, 10, 20) != GenomicRegion(2, 10, 20));
}

BOOST_AUTO_TEST_CASE(follows_and_precedes)
{
    const auto A = GenomicRegion(1, 10, 20);

    BOOST_CHECK(GenomicRegion(1, 21, 30).follows(A));
    BOOST_CHECK(A.precedes(GenomicRegion(1, 21, 30)));

    BOOST_CHECK(!GenomicRegion(1, 22, 30).follows(A));
    BOOST_CHECK(!GenomicRegion(1, 20, 30).follows(A));
    BOOST_CHECK(!GenomicRegion(2, 21, 30).follows(A));
    BOOST_CHECK(!A.follows(GenomicRegion(1, 21, 30)));
}

BOOST_AUTO_TEST_CASE(contains_position)
{
    const auto A = GenomicRegion(1, 10, 20);

    BOOST_CHECK(A.contains(GenomicPosition(1, 10)));
    BOOST_CHECK(A.contains(GenomicPosition(1, 15)));
    BOOST_CHECK(A.contains(GenomicPosition(1, 20)));
    BOOST_CHECK(!A.contains(GenomicPosition(1, 9)));
    BOOST_CHECK(!A.contains(GenomicPosition(1, 21)));
    BOOST_CHECK(!A.contains(GenomicPosition(2, 15)));

    BOOST_CHECK(A.strictly_contains(GenomicPosition(1, 15)));
    BOOST_CHECK(!A.strictly_contains(GenomicPosition(1, 10)));
    BOOST_CHECK(!A.strictly_contains(GenomicPosition(1, 20)));

    BOOST_CHECK(!GenomicRegion().contains(GenomicPosition(0, 0)));
}

BOOST_AUTO_TEST_CASE(contains_region)
{
    const auto A = GenomicRegion(1, 10, 20);

    BOOST_CHECK(A.contains(A));
    BOOST_CHECK(A.contains(GenomicRegion(1, 12, 18)));
    BOOST_CHECK(A.contains(GenomicRegion(1, 10, 15)));
    BOOST_CHECK(!A.contains(GenomicRegion(1, 5, 15)));
    BOOST_CHECK(!A.contains(GenomicRegion(1, 15, 25)));
    BOOST_CHECK(!A.contains(GenomicRegion(2, 12, 18)));

    BOOST_CHECK(A.strictly_contains(GenomicRegion(1, 12, 18)));
    BOOST_CHECK(!A.strictly_contains(A));
    BOOST_CHECK(!A.strictly_contains(GenomicRegion(1, 10, 15)));
    BOOST_CHECK(!A.strictly_contains(GenomicRegion(1, 15, 20)));
}

BOOST_AUTO_TEST_CASE(ends_before_and_begins_after)
{
    const auto A = GenomicRegion(1, 10, 20);

    BOOST_CHECK(A.ends_before(GenomicPosition(1, 21)));
    BOOST_CHECK(!A.ends_before(GenomicPosition(1, 20)));
    BOOST_CHECK(!A.ends_before(GenomicPosition(2, 21)));

    BOOST_CHECK(A.begins_after(GenomicPosition(1, 9)));
    BOOST_CHECK(!A.begins_after(GenomicPosition(1, 10)));

    BOOST_CHECK(A.ends_before(GenomicRegion(1, 21, 30)));
    BOOST_CHECK(!A.ends_before(GenomicRegion(1, 20, 30)));
    BOOST_CHECK(GenomicRegion(1, 21, 30).begins_after(A));
    BOOST_CHECK(!GenomicRegion(1, 20, 30).begins_after(A));
}

BOOST_AUTO_TEST_CASE(overlaps)
{
    const auto A = GenomicRegion(1, 10, 20);

    BOOST_CHECK(A.overlaps(A));
    BOOST_CHECK(A.overlaps(GenomicRegion(1, 15, 30)));
    BOOST_CHECK(A.overlaps(GenomicRegion(1, 1, 15)));
    BOOST_CHECK(A.overlaps(GenomicRegion(1, 12, 18)));
    BOOST_CHECK(A.overlaps(GenomicRegion(1, 1, 30)));
    BOOST_CHECK(A.overlaps(GenomicRegion(1, 20, 30)));
    BOOST_CHECK(A.overlaps(GenomicRegion(1, 1, 10)));

    BOOST_CHECK(!A.overlaps(GenomicRegion(1, 21, 30)));
    BOOST_CHECK(!A.overlaps(GenomicRegion(1, 1, 9)));
    BOOST_CHECK(!A.overlaps(GenomicRegion(2, 10, 20)));

    std::vector<GenomicRegion> regions{GenomicRegion(1, 1, 5), GenomicRegion(1, 30, 40)};
    BOOST_CHECK(!A.overlaps(regions));
    regions.push_back(GenomicRegion(1, 18, 25));
    BOOST_CHECK(A.overlaps(regions));
}

BOOST_AUTO_TEST_CASE(split)
{
    auto A = GenomicRegion(1, 10, 20);

    const auto B = A.split(GenomicPosition(1, 15));
    BOOST_CHECK(A == GenomicRegion(1, 10, 14));
    BOOST_CHECK(B == GenomicRegion(1, 15, 20));

    auto C = GenomicRegion(1, 10, 20);
    BOOST_CHECK_THROW(C.split(GenomicPosition(1, 10)), std::domain_error);
    BOOST_CHECK_THROW(C.split(GenomicPosition(1, 21)), std::domain_error);
    BOOST_CHECK_THROW(C.split(GenomicPosition(2, 15)), std::domain_error);
}

BOOST_AUTO_TEST_CASE(join)
{
    auto A = GenomicRegion(1, 10, 14);
    auto B = GenomicRegion(1, 15, 20);

    A.join(B);
    BOOST_CHECK(A == GenomicRegion(1, 10, 20));
    BOOST_CHECK_EQUAL(B.size(), 0);

    // joining in the opposite order
    auto C = GenomicRegion(1, 15, 20);
    auto D = GenomicRegion(1, 10, 14);
    C.join(D);
    BOOST_CHECK(C == GenomicRegion(1, 10, 20));

    auto E = GenomicRegion(1, 10, 14);
    auto F = GenomicRegion(1, 16, 20);
    BOOST_CHECK_THROW(E.join(F), std::domain_error);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(intersection)

/**
 * @brief Check that the intersection of two regions is the expected one
 *
 * The test is performed in both argument orders.
 */
void check_intersection(const GenomicRegion& A, const GenomicRegion& B,
                        const GenomicRegion& expected)
{
    BOOST_CHECK_EQUAL(intersect(A, B), expected);
    BOOST_CHECK_EQUAL(intersect(B, A), expected);
}

/**
 * @brief Check that the intersection of two regions is empty
 *
 * The test is performed in both argument orders.
 */
void check_empty_intersection(const GenomicRegion& A, const GenomicRegion& B)
{
    BOOST_CHECK_EQUAL(intersect(A, B).size(), 0);
    BOOST_CHECK_EQUAL(intersect(B, A).size(), 0);
}

BOOST_AUTO_TEST_CASE(partial_overlap)
{
    check_intersection(GenomicRegion(1, 10, 20), GenomicRegion(1, 15, 30), GenomicRegion(1, 15, 20));
}

BOOST_AUTO_TEST_CASE(containment)
{
    check_intersection(GenomicRegion(1, 10, 20), GenomicRegion(1, 12, 18), GenomicRegion(1, 12, 18));
    check_intersection(GenomicRegion(1, 10, 20), GenomicRegion(1, 10, 15), GenomicRegion(1, 10, 15));
    check_intersection(GenomicRegion(1, 10, 20), GenomicRegion(1, 15, 20), GenomicRegion(1, 15, 20));
}

BOOST_AUTO_TEST_CASE(identical)
{
    check_intersection(GenomicRegion(1, 10, 20), GenomicRegion(1, 10, 20), GenomicRegion(1, 10, 20));
}

BOOST_AUTO_TEST_CASE(single_base_overlap)
{
    check_intersection(GenomicRegion(1, 10, 20), GenomicRegion(1, 20, 30), GenomicRegion(1, 20, 20));
    check_intersection(GenomicRegion(1, 1, 1), GenomicRegion(1, 1, 5), GenomicRegion(1, 1, 1));
}

BOOST_AUTO_TEST_CASE(no_overlap)
{
    check_empty_intersection(GenomicRegion(1, 10, 20), GenomicRegion(1, 21, 30));
    check_empty_intersection(GenomicRegion(1, 10, 20), GenomicRegion(1, 40, 50));
}

BOOST_AUTO_TEST_CASE(different_chromosomes)
{
    check_empty_intersection(GenomicRegion(1, 10, 20), GenomicRegion(2, 10, 20));
}

BOOST_AUTO_TEST_CASE(large_positions)
{
    const auto max_pos = std::numeric_limits<ChrPosition>::max();

    check_intersection(GenomicRegion(1, max_pos-1000, max_pos), GenomicRegion(1, max_pos-10, max_pos),
                       GenomicRegion(1, max_pos-10, max_pos));
}

BOOST_AUTO_TEST_CASE(consistent_with_overlaps)
{
    const auto A = GenomicRegion(1, 10, 20);

    for (ChrPosition begin=1; begin<30; ++begin) {
        for (ChrPosition end=begin; end<30; ++end) {
            const auto B = GenomicRegion(1, begin, end);
            const auto I = intersect(A, B);

            BOOST_CHECK_EQUAL(A.overlaps(B), I.size()!=0);

            if (I.size()!=0) {
                BOOST_CHECK(A.contains(I));
                BOOST_CHECK(B.contains(I));
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
