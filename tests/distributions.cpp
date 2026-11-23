/**
 * @file distributions.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing distribution functions and classes
 * @version 1.3
 * @date 2026-11-23
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
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>

#include "distributions.hpp"

#include "error.hpp"

struct MeanStddev
{
    double mean;
    double stddev;
};

using dist_type = CLONES::truncated_normal_distribution<double>;

const double inf = std::numeric_limits<double>::infinity();

// aimed (mean, stddev) pairs for a cut equal to 150
const std::vector<MeanStddev> feasible_cases{{300, 50}, {350, 100}, {500, 30},
                                             {450, 10}, {200, 40}, {200, 45},
                                             {200, 49}};

const double cut = 150;

// an upper cut for the two-sided cases
const double upper_cut = 600;

// aimed (mean, stddev) pairs for the interval [cut, upper_cut]
const std::vector<MeanStddev> two_sided_feasible_cases{{300, 50}, {300, 100}, {375, 129},
                                                       {200, 40}, {590, 9}, {500, 30},
                                                       {160, 9}};

// the floating-point type used to compute the exact moments
using exact_type = boost::multiprecision::cpp_bin_float_50;

/**
 * @brief Exact mean and standard deviation of N(mu, sigma^2) truncated to [a,b]
 *
 * The moments are computed from the closed-form formulas with 50 decimal
 * digits, so they do not suffer from the cancellation that affects the
 * double-precision formulas far in the tails. The interval is mirrored when
 * it lies in the left tail, so that the tail probabilities are computed by
 * erfc() rather than by subtracting values close to 1.
 */
MeanStddev truncated_moments(const double mu, const double sigma,
                             const double a, const double b)
{
    exact_type alpha = (exact_type(a) - mu) / sigma;
    exact_type beta = (exact_type(b) - mu) / sigma;
    int sign = 1;
    if (b <= mu) {
        std::swap(alpha, beta);
        alpha = -alpha;
        beta = -beta;
        sign = -1;
    }

    const exact_type pi = boost::math::constants::pi<exact_type>();

    auto pdf = [&pi](const exact_type& z) -> exact_type {
        return boost::multiprecision::isinf(z) ? exact_type(0)
                                               : exp(-z * z / 2) / sqrt(2 * pi);
    };
    auto upper_tail = [](const exact_type& z) -> exact_type {
        if (boost::multiprecision::isinf(z)) {
            return (z > 0 ? 0 : 1);
        }
        return boost::math::erfc(z / sqrt(exact_type(2))) / 2;
    };
    auto z_pdf = [&pdf](const exact_type& z) -> exact_type {
        return boost::multiprecision::isinf(z) ? exact_type(0) : exact_type(z * pdf(z));
    };

    const exact_type Z = upper_tail(alpha) - upper_tail(beta);
    const exact_type m = (pdf(alpha) - pdf(beta)) / Z;
    const exact_type variance = 1 + (z_pdf(alpha) - z_pdf(beta)) / Z - m * m;

    return {static_cast<double>(mu + sign * sigma * m),
            static_cast<double>(sigma * sqrt(variance))};
}

/**
 * @brief Exact CDF of N(mu, sigma^2) truncated to [a,b]
 */
long double truncated_cdf(const long double x, const long double mu,
                          const long double sigma, const long double a,
                          const long double b)
{
    auto upper_tail = [mu, sigma](const long double v) {
        return std::isinf(v) ? (v > 0 ? 0.0L : 1.0L)
                             : 0.5L * std::erfc((v - mu) / (sigma * std::sqrt(2.0L)));
    };

    return (upper_tail(a) - upper_tail(x)) / (upper_tail(a) - upper_tail(b));
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_moments)
{
    for (const auto& [mean, stddev] : feasible_cases) {
        const auto normal = dist_type::fit_untruncated(mean, stddev, cut);
        const auto moments = truncated_moments(normal.mean(), normal.stddev(), cut, inf);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(moments.mean, mean, 1e-8);
            BOOST_CHECK_CLOSE(moments.stddev, stddev, 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(fit_right_truncated_normal_moments)
{
    // mirror the left-truncated cases
    for (const auto& [mean, stddev] : feasible_cases) {
        const auto normal = dist_type::fit_untruncated(-mean, stddev, -inf, -cut);
        const auto left_normal = dist_type::fit_untruncated(mean, stddev, cut);
        const auto moments = truncated_moments(normal.mean(), normal.stddev(), -inf, -cut);

        BOOST_TEST_CONTEXT("aimed " << -mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(moments.mean, -mean, 1e-8);
            BOOST_CHECK_CLOSE(moments.stddev, stddev, 1e-8);
            BOOST_CHECK_CLOSE(normal.mean(), -left_normal.mean(), 1e-10);
            BOOST_CHECK_CLOSE(normal.stddev(), left_normal.stddev(), 1e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(fit_two_sided_truncated_normal_moments)
{
    for (const auto& [mean, stddev] : two_sided_feasible_cases) {
        const auto normal = dist_type::fit_untruncated(mean, stddev, cut, upper_cut);
        const auto moments = truncated_moments(normal.mean(), normal.stddev(),
                                               cut, upper_cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(moments.mean, mean, 1e-8);
            BOOST_CHECK_CLOSE(moments.stddev, stddev, 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_close_to_limit)
{
    // stddev close to mean - cut: the standardized cut of the fitted
    // distribution is large
    const std::vector<MeanStddev> cases{{200, 49.9}, {200, 49.999}, {300, 149.9},
                                        {3000, 2849}, {1150, 1000 * (1 - 1e-8)},
                                        {1150, 1000 * (1 - 1e-12)}};

    for (const auto& [mean, stddev] : cases) {
        const auto fitted = dist_type::fit_untruncated(mean, stddev, cut);
        const dist_type dist(fitted, cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_GT(fitted.stddev(), 0);
            BOOST_CHECK_CLOSE(dist.mean(), mean, 1e-8);
            BOOST_CHECK_CLOSE(dist.stddev(), stddev, 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(fit_two_sided_truncated_normal_close_to_limit)
{
    struct Case {
        double mean, stddev, a, b;
    };

    // stddev close to the one of the truncated exponential distribution
    // having the aimed mean, i.e., the limit for sigma -> +infinity
    const std::vector<Case> cases{{160, 9.99, cut, upper_cut},
                                  {375, 129.8, cut, upper_cut},
                                  {0.9, 0.0998, 0, 1},
                                  {1e6 + 0.9, 0.0998, 1e6, 1e6 + 1},
                                  {1150, 999.9, cut, 1e5}};

    for (const auto& [mean, stddev, a, b] : cases) {
        const auto fitted = dist_type::fit_untruncated(mean, stddev, a, b);
        const dist_type dist(fitted, a, b);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev << " in ["
                                    << a << "," << b << "]") {
            BOOST_CHECK_GT(fitted.stddev(), 0);
            BOOST_CHECK_CLOSE(dist.mean(), mean, 1e-8);
            BOOST_CHECK_CLOSE(dist.stddev(), stddev, 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(fit_truncated_normal_negligible_truncation)
{
    // no cut at all
    auto normal = dist_type::fit_untruncated(300.0, 50.0);
    BOOST_CHECK_EQUAL(normal.mean(), 300);
    BOOST_CHECK_EQUAL(normal.stddev(), 50);

    // cuts farther than 10 stddev from the mean
    normal = dist_type::fit_untruncated(300.0, 10.0, cut, upper_cut);
    BOOST_CHECK_EQUAL(normal.mean(), 300);
    BOOST_CHECK_EQUAL(normal.stddev(), 10);
}

BOOST_AUTO_TEST_CASE(fit_truncated_normal_degenerate)
{
    for (const double mean : {150.0, 300.0}) {
        const auto normal = dist_type::fit_untruncated(mean, 0.0, cut);

        BOOST_CHECK_EQUAL(normal.mean(), mean);
        BOOST_CHECK_EQUAL(normal.stddev(), 0);
    }

    for (const double mean : {cut, 300.0, upper_cut}) {
        const auto normal = dist_type::fit_untruncated(mean, 0.0, cut, upper_cut);

        BOOST_CHECK_EQUAL(normal.mean(), mean);
        BOOST_CHECK_EQUAL(normal.stddev(), 0);
    }
}

BOOST_AUTO_TEST_CASE(fit_left_truncated_normal_infeasible)
{
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

BOOST_AUTO_TEST_CASE(fit_two_sided_truncated_normal_infeasible)
{
    // stddev not smaller than mean - a or b - mean
    BOOST_CHECK_THROW(dist_type::fit_untruncated(160.0, 10.0, cut, upper_cut),
                      std::invalid_argument);
    BOOST_CHECK_THROW(dist_type::fit_untruncated(590.0, 10.0, cut, upper_cut),
                      std::invalid_argument);

    // stddev greater than the one of the truncated exponential distribution
    // with the aimed mean (about 0.09982) and than the one of the uniform
    // distribution (about 0.288675)
    BOOST_CHECK_THROW(dist_type::fit_untruncated(0.9, 0.0999, 0.0, 1.0),
                      std::invalid_argument);
    BOOST_CHECK_THROW(dist_type::fit_untruncated(0.5, 0.2887, 0.0, 1.0),
                      std::invalid_argument);

    // mean out of [a,b], or equal to a cut with non-null stddev
    BOOST_CHECK_THROW(dist_type::fit_untruncated(700.0, 10.0, cut, upper_cut),
                      std::invalid_argument);
    BOOST_CHECK_THROW(dist_type::fit_untruncated(upper_cut, 10.0, cut, upper_cut),
                      std::invalid_argument);
    BOOST_CHECK_THROW(dist_type::fit_untruncated(-100.0, 10.0, -inf, -cut),
                      std::invalid_argument);

    // the lower cut is greater than the upper cut
    BOOST_CHECK_THROW(dist_type::fit_untruncated(300.0, 10.0, upper_cut, cut),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(truncated_normal_creation)
{
    BOOST_CHECK_NO_THROW(dist_type());
    BOOST_CHECK_NO_THROW(dist_type(0, 1, 2));
    BOOST_CHECK_NO_THROW(dist_type(0, 1, -inf, 2));
    BOOST_CHECK_NO_THROW(dist_type(0, 1, 2, 3));
    BOOST_CHECK_NO_THROW(dist_type(0, 1, 2, 2));
    BOOST_CHECK_NO_THROW(dist_type(150, 0, 150));
    BOOST_CHECK_NO_THROW(dist_type(600, 0, 150, 600));
    BOOST_CHECK_NO_THROW(dist_type(std::normal_distribution<double>(0, 1)));
    BOOST_CHECK_NO_THROW(dist_type(std::normal_distribution<double>(0, 1), 2, 3));

    // negative sigma
    BOOST_CHECK_THROW(dist_type(0, -1, 0), std::invalid_argument);

    // null sigma and mu out of [a,b]
    BOOST_CHECK_THROW(dist_type(100, 0, 150), std::invalid_argument);
    BOOST_CHECK_THROW(dist_type(700, 0, 150, 600), std::invalid_argument);

    // a greater than b
    BOOST_CHECK_THROW(dist_type(0, 1, 3, 2), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(truncated_normal_getters)
{
    const dist_type dist(300, 50, cut, upper_cut);

    BOOST_CHECK_EQUAL(dist.get_mu(), 300);
    BOOST_CHECK_EQUAL(dist.get_sigma(), 50);
    BOOST_CHECK_EQUAL(dist.get_a(), cut);
    BOOST_CHECK_EQUAL(dist.get_b(), upper_cut);

    const dist_type unbounded(std::normal_distribution<double>(10, 2));
    BOOST_CHECK_EQUAL(unbounded.get_mu(), 10);
    BOOST_CHECK_EQUAL(unbounded.get_sigma(), 2);
}

/**
 * @brief Check that samples are in [a,b] and that their moments are the expected ones
 */
void check_sample_moments(const dist_type& dist, const size_t num_of_samples,
                          std::mt19937_64& generator)
{
    std::vector<double> samples(num_of_samples);
    for (auto& sample : samples) {
        sample = dist(generator);
    }

    const auto [min_value, max_value] = std::minmax_element(samples.begin(), samples.end());

    double sample_mean{0};
    for (const auto& sample : samples) {
        sample_mean += sample;
    }
    sample_mean /= num_of_samples;

    // two passes: E[X^2] - E[X]^2 cancels when the stddev is small with
    // respect to the mean
    double sample_variance{0};
    for (const auto& sample : samples) {
        sample_variance += (sample - sample_mean) * (sample - sample_mean);
    }
    const double sample_stddev = std::sqrt(sample_variance / num_of_samples);

    BOOST_CHECK_GE(*min_value, dist.get_a());
    BOOST_CHECK_LE(*max_value, dist.get_b());

    // the sample mean is within 5 standard errors from the mean
    BOOST_CHECK_SMALL(sample_mean - dist.mean(),
                      5 * dist.stddev() / std::sqrt(num_of_samples));
    BOOST_CHECK_CLOSE(sample_stddev, dist.stddev(), 1.0);
}

BOOST_AUTO_TEST_CASE(truncated_normal_sampling)
{
    const size_t num_of_samples = 1000000;

    std::mt19937_64 generator(0);

    for (const auto& [mean, stddev] : feasible_cases) {
        const dist_type dist(dist_type::fit_untruncated(mean, stddev, cut), cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            check_sample_moments(dist, num_of_samples, generator);
            BOOST_CHECK_CLOSE(dist.mean(), mean, 1e-8);
            BOOST_CHECK_CLOSE(dist.stddev(), stddev, 1e-8);
        }
    }

    for (const auto& [mean, stddev] : two_sided_feasible_cases) {
        const dist_type dist(dist_type::fit_untruncated(mean, stddev, cut, upper_cut),
                             cut, upper_cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev << " in ["
                                    << cut << "," << upper_cut << "]") {
            check_sample_moments(dist, num_of_samples, generator);
            BOOST_CHECK_CLOSE(dist.mean(), mean, 1e-8);
            BOOST_CHECK_CLOSE(dist.stddev(), stddev, 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(truncated_normal_sampling_distribution)
{
    struct Interval {
        double a, b;
    };

    // these intervals cover all the branches of the sampler: [alpha, beta]
    // containing 0 and wide or narrow, in the right tail and wide or narrow,
    // and in the left tail
    const std::vector<Interval> intervals{{-inf, inf}, {-1, 2}, {-0.5, 0.5}, {-3, 0.1},
                                          {0, inf}, {2, inf}, {1, 4}, {2, 2.3},
                                          {5, 5.05}, {-inf, -2}, {-2.3, -2}, {-inf, 0}};

    const size_t num_of_samples = 100000;

    // the critical value of the Kolmogorov-Smirnov statistic for a
    // significance level of 0.1%
    const double critical_value = 1.95 / std::sqrt(num_of_samples);

    std::mt19937_64 generator(0);

    for (const auto& [a, b] : intervals) {
        const double mu = 10, sigma = 3;
        const dist_type dist(mu, sigma, mu + sigma * a, mu + sigma * b);

        std::vector<double> samples(num_of_samples);
        for (auto& sample : samples) {
            sample = dist(generator);
        }
        std::sort(samples.begin(), samples.end());

        double ks_statistic = 0;
        for (size_t i=0; i<num_of_samples; ++i) {
            const auto cdf = truncated_cdf(samples[i], mu, sigma, dist.get_a(), dist.get_b());

            ks_statistic = std::max({ks_statistic,
                                     static_cast<double>(cdf - static_cast<long double>(i) / num_of_samples),
                                     static_cast<double>(static_cast<long double>(i + 1) / num_of_samples - cdf)});
        }

        BOOST_TEST_CONTEXT("standardized interval [" << a << "," << b << "]") {
            BOOST_CHECK_LT(ks_statistic, critical_value);
        }
    }
}

BOOST_AUTO_TEST_CASE(truncated_normal_far_tail_sampling)
{
    const size_t num_of_samples = 100000;

    std::mt19937_64 generator(0);

    for (const double alpha : {30.0, 1e4}) {
        const dist_type dist(0, 1, alpha);

        BOOST_TEST_CONTEXT("alpha " << alpha) {
            check_sample_moments(dist, num_of_samples, generator);
        }
    }
}

BOOST_AUTO_TEST_CASE(truncated_normal_mean)
{
    for (const auto& [mean, stddev] : feasible_cases) {
        const dist_type dist(dist_type::fit_untruncated(mean, stddev, cut), cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(dist.mean(), mean, 1e-8);
        }
    }

    // left-truncated standard normals, also beyond the underflow of the
    // double-precision formula (about 37)
    for (const double alpha : {-5.0, 0.0, 2.9, 3.1, 5.0, 29.0, 31.0, 50.0, 1e4}) {
        const dist_type dist(0, 1, alpha);

        BOOST_TEST_CONTEXT("alpha " << alpha) {
            // the mean may be close to 0: compare with respect to the stddev
            const auto expected = truncated_moments(0, 1, alpha, inf);
            BOOST_CHECK_SMALL(dist.mean() - expected.mean,
                              1e-12 * std::max(std::abs(expected.mean), expected.stddev));
            BOOST_CHECK_GT(dist.mean(), alpha);
        }
    }

    // two-sided truncations
    for (const auto& [a, b] : std::vector<std::pair<double, double>>{{-1, 2}, {-0.5, 0.5},
                                                                     {-3, 0.1}, {2, 2.3},
                                                                     {5, 5.05}, {-2.3, -2},
                                                                     {20, 25}, {-25, -20}}) {
        const dist_type dist(0, 1, a, b);

        BOOST_TEST_CONTEXT("[" << a << "," << b << "]") {
            const auto expected = truncated_moments(0, 1, a, b);
            BOOST_CHECK_SMALL(dist.mean() - expected.mean,
                              1e-12 * std::max(std::abs(expected.mean), expected.stddev));
            BOOST_CHECK_GT(dist.mean(), a);
            BOOST_CHECK_LT(dist.mean(), b);

            // mirroring the interval mirrors the mean
            BOOST_CHECK_CLOSE(dist_type(0, 1, -b, -a).mean(), -dist.mean(), 1e-10);
        }
    }

    // no truncation
    BOOST_CHECK_CLOSE(dist_type(300, 50).mean(), 300, 1e-10);
    BOOST_CHECK_CLOSE(dist_type(std::normal_distribution<double>(300, 50)).mean(), 300, 1e-10);

    // the cut is far from mu: mu + sigma*E[Z] would lose about 8 digits
    const dist_type far(-1e14, 1e10, cut);
    const auto far_expected = truncated_moments(-1e14, 1e10, cut, inf);
    BOOST_CHECK_CLOSE(far.mean(), far_expected.mean, 1e-10);
    BOOST_CHECK_CLOSE(far.stddev(), far_expected.stddev, 1e-10);

    // degenerate distributions always return mu
    BOOST_CHECK_EQUAL(dist_type(300, 0, cut).mean(), 300);
    BOOST_CHECK_EQUAL(dist_type(300, 0, cut, upper_cut).mean(), 300);
}

BOOST_AUTO_TEST_CASE(truncated_normal_stddev)
{
    for (const auto& [mean, stddev] : feasible_cases) {
        const dist_type dist(dist_type::fit_untruncated(mean, stddev, cut), cut);

        BOOST_TEST_CONTEXT("aimed " << mean << "±" << stddev) {
            BOOST_CHECK_CLOSE(dist.stddev(), stddev, 1e-8);
        }
    }

    // left-truncated standard normals, also beyond the underflow of the
    // double-precision formula (about 37)
    for (const double alpha : {-5.0, 0.0, 2.9, 3.1, 5.0, 29.0, 31.0, 50.0, 1e4}) {
        const dist_type dist(0, 1, alpha);

        BOOST_TEST_CONTEXT("alpha " << alpha) {
            BOOST_CHECK_CLOSE(dist.stddev(), truncated_moments(0, 1, alpha, inf).stddev, 1e-10);
            BOOST_CHECK_GT(dist.stddev(), 0);
            BOOST_CHECK_LT(dist.stddev(), 1);
        }
    }

    // two-sided truncations
    for (const auto& [a, b] : std::vector<std::pair<double, double>>{{-1, 2}, {-0.5, 0.5},
                                                                     {-3, 0.1}, {2, 2.3},
                                                                     {5, 5.05}, {-2.3, -2},
                                                                     {20, 25}, {-25, -20}}) {
        const dist_type dist(0, 1, a, b);

        BOOST_TEST_CONTEXT("[" << a << "," << b << "]") {
            BOOST_CHECK_CLOSE(dist.stddev(), truncated_moments(0, 1, a, b).stddev, 1e-10);

            // mirroring the interval does not change the standard deviation
            BOOST_CHECK_CLOSE(dist_type(0, 1, -b, -a).stddev(), dist.stddev(), 1e-10);
        }
    }

    // narrow intervals: the distribution is almost uniform
    const double width = 1e-3;
    BOOST_CHECK_CLOSE(dist_type(0, 1, -width/2, width/2).stddev(),
                      width / std::sqrt(12.0), 1e-4);

    // the moments are continuous where the computation switches from the
    // interval containing 0 to the right tail
    const dist_type below(0, 1, std::nextafter(0.0, -1.0)), above(0, 1, 0);
    BOOST_CHECK_CLOSE(below.stddev(), above.stddev(), 1e-10);
    BOOST_CHECK_CLOSE(below.mean(), above.mean(), 1e-10);

    // no truncation
    BOOST_CHECK_CLOSE(dist_type(300, 50).stddev(), 50, 1e-10);

    // degenerate distributions always return 0
    BOOST_CHECK_EQUAL(dist_type(300, 0, cut).stddev(), 0);
    BOOST_CHECK_EQUAL(dist_type(300, 0, cut, upper_cut).stddev(), 0);
}

BOOST_AUTO_TEST_CASE(truncated_normal_degenerate_sampling)
{
    std::mt19937_64 generator(0);

    for (const dist_type& dist : {dist_type(300, 0, cut), dist_type(300, 0, cut, upper_cut)}) {
        for (size_t i=0; i<100; ++i) {
            BOOST_CHECK_EQUAL(dist(generator), 300);
        }
    }
}
