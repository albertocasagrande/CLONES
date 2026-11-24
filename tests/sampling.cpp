/**
 * @file sampling.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing sampling functions
 * @version 1.0
 * @date 2026-11-24
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
#define BOOST_TEST_MODULE sampling

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "sampling.hpp"

/**
 * @brief Check that a sample contains `k` distinct values in [0, n-1]
 *
 * @param sample is the sample to be checked
 * @param n is the number of values from which the sample was drawn
 * @param k is the number of sampled values
 */
void check_sample(const std::vector<size_t>& sample, const size_t n, const size_t k)
{
    BOOST_CHECK_EQUAL(sample.size(), k);

    const std::set<size_t> values(sample.begin(), sample.end());
    BOOST_CHECK_EQUAL(values.size(), k);

    if (!values.empty()) {
        BOOST_CHECK_LT(*values.rbegin(), n);
    }
}

BOOST_AUTO_TEST_CASE(sample_distinct_values)
{
    std::mt19937_64 generator(0);

    for (const auto& [n, k] : std::vector<std::pair<size_t, size_t>>{{1, 0}, {1, 1}, {5, 0},
                                                                     {5, 2}, {5, 5}, {100, 1},
                                                                     {100, 50}, {100, 99},
                                                                     {100, 100}, {1000, 10}}) {
        BOOST_TEST_CONTEXT("n " << n << ", k " << k) {
            for (size_t i=0; i<100; ++i) {
                check_sample(CLONES::sample_distinct(n, k, generator), n, k);
            }
        }
    }

    // empty samples
    BOOST_CHECK(CLONES::sample_distinct(0, 0, generator).empty());
    BOOST_CHECK(CLONES::sample_distinct(10, 0, generator).empty());
}

BOOST_AUTO_TEST_CASE(sample_distinct_all_values)
{
    std::mt19937_64 generator(0);

    // when k == n, the sample contains all the values in [0, n-1]
    for (const size_t n : {1, 2, 10, 1000}) {
        auto sample = CLONES::sample_distinct(n, n, generator);
        std::sort(sample.begin(), sample.end());

        std::vector<size_t> all_values(n);
        for (size_t i=0; i<n; ++i) {
            all_values[i] = i;
        }

        BOOST_TEST_CONTEXT("n " << n) {
            BOOST_CHECK(sample == all_values);
        }
    }
}

BOOST_AUTO_TEST_CASE(sample_distinct_large_n)
{
    std::mt19937_64 generator(0);

    // the cost does not depend on n: drawing few values from a huge set is
    // immediate
    for (const size_t n : {static_cast<size_t>(1e12), std::numeric_limits<size_t>::max()}) {
        BOOST_TEST_CONTEXT("n " << n) {
            for (size_t i=0; i<1000; ++i) {
                check_sample(CLONES::sample_distinct(n, 5, generator), n, 5);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(sample_distinct_invalid)
{
    std::mt19937_64 generator(0);

    BOOST_CHECK_THROW(CLONES::sample_distinct(5, 6, generator), std::invalid_argument);
    BOOST_CHECK_THROW(CLONES::sample_distinct(0, 1, generator), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(sample_distinct_reproducibility)
{
    // the same seed produces the same sample
    std::mt19937_64 generator_a(42), generator_b(42);

    for (size_t i=0; i<100; ++i) {
        BOOST_CHECK(CLONES::sample_distinct(1000, 20, generator_a)
                        == CLONES::sample_distinct(1000, 20, generator_b));
    }
}

BOOST_AUTO_TEST_CASE(sample_distinct_generators)
{
    // any uniform random bit generator can be used
    std::mt19937 mt19937(0);
    std::minstd_rand minstd_rand(0);
    std::ranlux24 ranlux24(0);

    check_sample(CLONES::sample_distinct(100, 10, mt19937), 100, 10);
    check_sample(CLONES::sample_distinct(100, 10, minstd_rand), 100, 10);
    check_sample(CLONES::sample_distinct(100, 10, ranlux24), 100, 10);
}

/**
 * @brief Compute the binomial coefficient
 *
 * @param n is the size of the set
 * @param k is the size of the subsets
 * @return the number of subsets of a set of size `n` having `k` elements
 */
size_t binomial(const size_t n, const size_t k)
{
    size_t result = 1;
    for (size_t i=1; i<=k; ++i) {
        result = result * (n - k + i) / i;
    }

    return result;
}

BOOST_AUTO_TEST_CASE(sample_distinct_uniform_subsets)
{
    // every k-subset of [0, n-1] has probability 1/binomial(n,k): Pearson's
    // chi-squared test with the critical values for a significance level of
    // 0.1% and binomial(n,k)-1 degrees of freedom
    struct Case {
        size_t n, k;
        double critical_value;
    };

    const std::vector<Case> cases{{5, 2, 27.877},    //  9 degrees of freedom
                                  {6, 3, 43.820},    // 19 degrees of freedom
                                  {7, 6, 22.458},    //  6 degrees of freedom
                                  {8, 1, 24.322}};   //  7 degrees of freedom

    const size_t num_of_samples = 200000;

    std::mt19937_64 generator(0);

    for (const auto& [n, k, critical_value] : cases) {
        std::map<std::vector<size_t>, size_t> counts;
        for (size_t i=0; i<num_of_samples; ++i) {
            auto sample = CLONES::sample_distinct(n, k, generator);
            std::sort(sample.begin(), sample.end());

            ++counts[sample];
        }

        const size_t num_of_subsets = binomial(n, k);
        const double expected = static_cast<double>(num_of_samples) / num_of_subsets;

        double chi_squared = 0;
        for (const auto& [subset, count] : counts) {
            chi_squared += (count - expected) * (count - expected) / expected;
        }

        BOOST_TEST_CONTEXT("n " << n << ", k " << k) {
            BOOST_CHECK_EQUAL(counts.size(), num_of_subsets);
            BOOST_CHECK_LT(chi_squared, critical_value);
        }
    }
}

BOOST_AUTO_TEST_CASE(sample_distinct_inclusion_probability)
{
    // every value is in the sample with probability k/n
    const size_t n = 50, k = 7, num_of_samples = 200000;

    std::mt19937_64 generator(0);

    std::vector<size_t> counts(n, 0);
    for (size_t i=0; i<num_of_samples; ++i) {
        for (const auto& value : CLONES::sample_distinct(n, k, generator)) {
            ++counts[value];
        }
    }

    const double p = static_cast<double>(k) / n;
    const double expected = num_of_samples * p;
    const double stddev = std::sqrt(num_of_samples * p * (1 - p));

    for (size_t value=0; value<n; ++value) {
        BOOST_TEST_CONTEXT("value " << value) {
            // within 5 standard deviations from the expected count
            BOOST_CHECK_SMALL(counts[value] - expected, 5 * stddev);
        }
    }
}
