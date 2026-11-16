/**
 * @file distributions.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing distribution functions and classes
 * @version 1.0
 * @date 2026-11-16
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
#define BOOST_TEST_MODULE distributions

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "distributions.hpp"

#include "error.hpp"

struct MeanStddev
{
    double mean;
    double stddev;
};

// aimed (mean, stddev) pairs for a cut equal to 150
const std::vector<MeanStddev> feasible_cases{{300, 50}, {350, 100}, {500, 30},
                                             {450, 10}, {200, 40}, {200, 45},
                                             {200, 49}};

const double cut = 150;

// exact mean and standard deviation of N(mu, sigma^2) left-truncated at cut
MeanStddev truncated_moments(const long double mu, const long double sigma,
                             const long double cut)
{
    const long double alpha = (cut - mu) / sigma;
    const long double tail = 0.5L * std::erfc(alpha / std::sqrt(2.0L));
    const long double pdf = std::exp(-alpha * alpha / 2)
                                / std::sqrt(2 * 3.14159265358979323846264338L);
    const long double lambda = pdf / tail;

    return {static_cast<double>(mu + sigma * lambda),
            static_cast<double>(sigma * std::sqrt(1 - lambda * (lambda - alpha)))};
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_moments)
{
    using namespace CLONES;

    for (const auto& [mean, stddev] : feasible_cases) {
        const auto normal = fit_left_truncated_normal(mean, stddev, cut);
        const auto moments = truncated_moments(normal.mean(), normal.stddev(), cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(moments.mean, mean, 0.1);
            BOOST_CHECK_CLOSE(moments.stddev, stddev, 0.1);
        }
    }
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_degenerate)
{
    using namespace CLONES;

    for (const double mean : {150.0, 300.0}) {
        const auto normal = fit_left_truncated_normal(mean, 0.0, cut);

        BOOST_CHECK_EQUAL(normal.mean(), mean);
        BOOST_CHECK_EQUAL(normal.stddev(), 0);
    }
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_infeasible)
{
    using namespace CLONES;

    // stddev not smaller than mean - cut
    BOOST_CHECK_THROW(fit_left_truncated_normal(200.0, 50.0, cut), std::invalid_argument);
    BOOST_CHECK_THROW(fit_left_truncated_normal(200.0, 60.0, cut), std::invalid_argument);

    // mean smaller than the cut, or equal to it with non-null stddev
    BOOST_CHECK_THROW(fit_left_truncated_normal(100.0, 10.0, cut), std::invalid_argument);
    BOOST_CHECK_THROW(fit_left_truncated_normal(150.0, 10.0, cut), std::invalid_argument);

    // negative stddev
    BOOST_CHECK_THROW(fit_left_truncated_normal(300.0, -10.0, cut), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(left_truncated_normal_creation)
{
    using namespace CLONES;

    BOOST_CHECK_NO_THROW(left_truncated_normal_distribution<double>());
    BOOST_CHECK_NO_THROW(left_truncated_normal_distribution<double>(0, 1, 2));
    BOOST_CHECK_NO_THROW(left_truncated_normal_distribution<double>(150, 0, 150));
    BOOST_CHECK_THROW(left_truncated_normal_distribution<double>(0, -1, 0),
                      std::invalid_argument);
    BOOST_CHECK_THROW(left_truncated_normal_distribution<double>(100, 0, 150),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(left_truncated_normal_sampling)
{
    using namespace CLONES;

    const size_t num_of_samples = 1000000;

    std::mt19937_64 generator(0);

    for (const auto& [mean, stddev] : feasible_cases) {
        const auto fitted = fit_left_truncated_normal(mean, stddev, cut);
        const left_truncated_normal_distribution<double> dist(fitted, cut);

        double sum{0}, sum_sq{0}, min_value{std::numeric_limits<double>::max()};
        for (size_t i=0; i<num_of_samples; ++i) {
            const double value = dist(generator);

            sum += value;
            sum_sq += value*value;
            min_value = std::min(min_value, value);
        }

        const double sample_mean = sum/num_of_samples;
        const double sample_stddev = std::sqrt(sum_sq/num_of_samples
                                               - sample_mean*sample_mean);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_GE(min_value, cut);
            BOOST_CHECK_CLOSE(sample_mean, mean, 0.5);
            BOOST_CHECK_CLOSE(sample_stddev, stddev, 1.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(left_truncated_normal_degenerate_sampling)
{
    using namespace CLONES;

    std::mt19937_64 generator(0);

    const left_truncated_normal_distribution<double> dist(300, 0, cut);
    for (size_t i=0; i<100; ++i) {
        BOOST_CHECK_EQUAL(dist(generator), 300);
    }
}
