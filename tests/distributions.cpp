/**
 * @file distributions.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing distribution functions and classes
 * @version 1.2
 * @date 2026-11-19
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

    using dist_type = left_truncated_normal_distribution<decltype(MeanStddev::mean)>;

    for (const auto& [mean, stddev] : feasible_cases) {
        const auto normal = dist_type::fit_untruncated(mean, stddev, cut);
        const auto moments = truncated_moments(normal.mean(), normal.stddev(), cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(moments.mean, mean, 0.1);
            BOOST_CHECK_CLOSE(moments.stddev, stddev, 0.1);
        }
    }
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_close_to_limit)
{
    using namespace CLONES;

    // stddev close to mean - cut: the standardized cut of the fitted
    // distribution is large
    const std::vector<MeanStddev> cases{{200, 49.9}, {200, 49.999}, {300, 149.9},
                                        {3000, 2849}, {1150, 1000 * (1 - 1e-8)},
                                        {1150, 1000 * (1 - 1e-12)}};

    for (const auto& [mean, stddev] : cases) {
        using dist_type = left_truncated_normal_distribution<decltype(MeanStddev::mean)>;

        const auto fitted = dist_type::fit_untruncated(mean, stddev, cut);
        const left_truncated_normal_distribution<double> dist(fitted, cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_GT(fitted.stddev(), 0);
            BOOST_CHECK_CLOSE(dist.mean(), mean, 1e-8);
            BOOST_CHECK_CLOSE(dist.stddev(), stddev, 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_degenerate)
{
    using namespace CLONES;

    for (const double mean : {150.0, 300.0}) {
        using dist_type = left_truncated_normal_distribution<double>;
        
        const auto normal = dist_type::fit_untruncated(mean, 0.0, cut);

        BOOST_CHECK_EQUAL(normal.mean(), mean);
        BOOST_CHECK_EQUAL(normal.stddev(), 0);
    }
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_infeasible)
{
    using namespace CLONES;

    using dist_type = left_truncated_normal_distribution<double>;

    // stddev not smaller than mean - cut
    BOOST_CHECK_THROW(dist_type::fit_untruncated(200.0, 50.0, cut), std::invalid_argument);
    BOOST_CHECK_THROW(dist_type::fit_untruncated(200.0, 60.0, cut), std::invalid_argument);

    // stddev too close to mean - cut to be fitted
    BOOST_CHECK_THROW(dist_type::fit_untruncated(200.0, std::nextafter(50.0, 0.0), cut),
                      std::invalid_argument);
    BOOST_CHECK_THROW(dist_type::fit_untruncated(1150.0, 1000 * (1 - 1e-15), cut),
                      std::invalid_argument);

    // mean smaller than the cut, or equal to it with non-null stddev
    BOOST_CHECK_THROW(dist_type::fit_untruncated(100.0, 10.0, cut), std::invalid_argument);
    BOOST_CHECK_THROW(dist_type::fit_untruncated(150.0, 10.0, cut), std::invalid_argument);

    // negative stddev
    BOOST_CHECK_THROW(dist_type::fit_untruncated(300.0, -10.0, cut), std::invalid_argument);
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
        using dist_type = left_truncated_normal_distribution<double>;

        const auto fitted = dist_type::fit_untruncated(mean, stddev, cut);
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

BOOST_AUTO_TEST_CASE(left_truncated_normal_mean)
{
    using namespace CLONES;

    for (const auto& [mean, stddev] : feasible_cases) {
        using dist_type = left_truncated_normal_distribution<double>;

        const auto fitted = dist_type::fit_untruncated(mean, stddev, cut);
        const left_truncated_normal_distribution<double> dist(fitted, cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(dist.mean(), mean, 0.1);
        }
    }

    // standardized cuts around the switch to the continued fraction (3) and
    // beyond the underflow of the exact formula (about 37)
    for (const double alpha : {-5.0, 0.0, 2.9, 3.1, 5.0, 29.0, 31.0, 50.0, 1e4}) {
        const left_truncated_normal_distribution<double> dist(0, 1, alpha);
        const double expected = (alpha > 31 ? alpha + 1/alpha - 2/(alpha*alpha*alpha)
                                 : truncated_moments(0, 1, alpha).mean);

        BOOST_TEST_CONTEXT("alpha " << alpha) {
            BOOST_CHECK_CLOSE(dist.mean(), expected, 1e-3);
            BOOST_CHECK_GT(dist.mean(), alpha);
        }
    }

    // degenerate distributions always return mu
    BOOST_CHECK_EQUAL(left_truncated_normal_distribution<double>(300, 0, cut).mean(), 300);
}

BOOST_AUTO_TEST_CASE(left_truncated_normal_stddev)
{
    using namespace CLONES;

    for (const auto& [mean, stddev] : feasible_cases) {
        using dist_type = left_truncated_normal_distribution<double>;

        const auto fitted = dist_type::fit_untruncated(mean, stddev, cut);
        const left_truncated_normal_distribution<double> dist(fitted, cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(dist.stddev(), stddev, 0.1);
        }
    }

    // standardized cuts around the switch to the continued fraction (3) and
    // beyond the underflow of the exact formula (about 37): for large alpha,
    // the standard deviation is about (1 - 3/alpha^2)/alpha
    for (const double alpha : {-5.0, 0.0, 2.9, 3.1, 5.0, 29.0, 31.0, 50.0, 1e4}) {
        const left_truncated_normal_distribution<double> dist(0, 1, alpha);
        const double expected = (alpha > 31 ? (1 - 3/(alpha*alpha))/alpha
                                 : truncated_moments(0, 1, alpha).stddev);

        BOOST_TEST_CONTEXT("alpha " << alpha) {
            BOOST_CHECK_CLOSE(dist.stddev(), expected, 1e-3);
            BOOST_CHECK_GT(dist.stddev(), 0);
            BOOST_CHECK_LT(dist.stddev(), 1);
        }
    }

    // the continued fraction agrees with the exact formula at the switch
    const double switch_alpha = 3;
    const left_truncated_normal_distribution<double> below(0, 1, std::nextafter(switch_alpha, 0.0)),
                                                     above(0, 1, switch_alpha);
    BOOST_CHECK_CLOSE(below.stddev(), above.stddev(), 1e-10);
    BOOST_CHECK_CLOSE(below.mean(), above.mean(), 1e-10);

    // degenerate distributions always return mu
    BOOST_CHECK_EQUAL(left_truncated_normal_distribution<double>(300, 0, cut).stddev(), 0);
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
