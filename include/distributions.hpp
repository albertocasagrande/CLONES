/**
 * @file distributions.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines functions about probability distributions
 * @version 1.4
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

#ifndef __CLONES_DISTRIBUTIONS__
#define __CLONES_DISTRIBUTIONS__

#include <algorithm>
#include <array>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <random>
#include <numbers>
#include <limits>
#include <type_traits>

#include "utils.hpp"

#include "error.hpp"

namespace CLONES
{



/**
 * @brief Probability density function (PDF) of the standard normal distribution.
 *
 * Computes \f$ \varphi(z) = \frac{1}{\sqrt{2\pi}} e^{-z^2/2} \f$.
 *
 * @tparam T floating-point type.
 * @param z point at which the density is evaluated.
 * @return the standard normal density at @p z.
 */
template<typename T>
T norm_pdf(const T z)
{
    if (std::isinf(z)) {
        return 0.0;
    }

    return (std::exp(-0.5 * z * z) * std::numbers::inv_sqrtpi_v<T>) /
                std::numbers::sqrt2_v<T>;
}

/**
 * @brief Cumulative distribution function (CDF) of the standard normal distribution.
 *
 * Computes \f$ \Phi(z) = \frac{1}{2}\left(1 + \mathrm{erf}(z/\sqrt{2})\right) \f$
 * using `std::erf`.
 *
 * @tparam T floating-point type.
 * @param z point at which the CDF is evaluated.
 * @return the probability \f$ P(Z \le z) \f$ with \f$ Z \sim N(0,1) \f$.
 */
template<typename T>
inline T norm_cdf(const T z)
{

    if (z == -std::numeric_limits<double>::infinity()) {
        return 0.0;
    }

    if (z ==  std::numeric_limits<double>::infinity()) {
        return 1.0;
    }

    return 0.5 * (1.0 + std::erf(z / std::numbers::sqrt2_v<T>));
}

/**
 * @brief A truncated normal distribution
 *
 * C. P. Robert. "Simulation of truncated normal variables". Statistics
 * and Computing 5, 121–125 (1995).
 *
 * @tparam T is the floating-point type
 */
template<typename T = double>
  requires std::is_floating_point_v<T>
class truncated_normal_distribution
{
    T mu;       //!< The mean of the non-truncated normal distribution
    T sigma;    //!< The standard deviation of the non-truncated normal distribution
    T a;        //!< The lower cut value
    T b;        //!< The upper cut value
    T alpha;    //!< The canonized lower cut value in \f$ N(0, 1) \f$
    T beta;     //!< The canonized upper cut value in \f$ N(0, 1) \f$

    /**
     * @brief Sample \f$N(0, 1)\f$ truncated to \f$[\alpha, \beta]\f$
     * 
     * This method samples \f$ X \sim N(0, 1) \f$ conditioned on
     * \f$X \in [\alpha, \beta]\f$.
     *
     * C. P. Robert. "Simulation of truncated normal variables". Statistics
     * and Computing 5, 121–125 (1995).
     *
     * @tparam RANDOM_GENERATOR is the random generator type
     * @param random_gen is a random generator
     * @param alpha is the lower bound cut
     * @param beta is the upper bound cut
     * @return a sample from the normal distribution \f$N(0,1)\f$
     *    greater than or equal to @p alpha and smaller than or
     *    equal to @p beta.
     */
    template<typename RANDOM_GENERATOR>
    static T canonical_sampler(RANDOM_GENERATOR& random_gen,
                               const T alpha, const T beta)
    {
        constexpr T sqrt2pi = std::numbers::sqrt2_v<T> /
                                    std::numbers::inv_sqrtpi_v<T>;

        // if beta <= 0 
        if (beta <= 0) {                      
            // mirror the sampling
            return -canonical_sampler(random_gen, -beta, -alpha);
        }

        std::uniform_real_distribution<T> u_sampler_01(0, 1);

        // if 0 \in (alpha, beta)
        if (alpha < 0) {
            T z;

            if (beta - alpha >= sqrt2pi) {
                // [alpha, beta] is larger than sqrt(2*pi) and the
                // probability of getting a sample in [alpha, beta]
                // by sampling N(0,1) is almost 1/2. As a matter
                // of the fact, Phi(sqrt(2*pi))-Phi(0) \approx 0.4939.
                // Hence, the expected number of draws of N(0,1) to
                // satisfy X in [alpha, beta] is 2.
                std::normal_distribution<T> n_sampler_01(0, 1);

                do {
                    z = n_sampler_01(random_gen);
                } while (z < alpha || z > beta);
            } else {
                // [alpha, beta] is narrow so we can sample z from [alpha, beta]
                // and return it when x \sim [0, 1] <= e^{-(z^2)/2}
                // (see Robert (1995)). The distribution will be the same
                // as N(0, 1) truncated in [alpha, beta].
                std::uniform_real_distribution<T> u_sampler_ab(alpha, beta);

                do {
                    z = u_sampler_ab(random_gen);
                } while (u_sampler_01(random_gen) > std::exp(-z * z / 2));
            }

            return z;
        }
        
        // if 0 <= alpha < beta: right tail

        // using a shifted exponential lambda*e^{-lambda(z-alpha)}
        // determine the proper lambda 
        const T alpha_square = alpha * alpha;
        const T s = std::sqrt(alpha_square + 4);
        const T lambda = (alpha + s) / 2;

        // determine the beta's cutoff to establish whether
        // the interval [alpha, beta] is wide or narrow
        const T cutoff = alpha + 2 * std::sqrt(std::numbers::e_v<T>) / (alpha + s)
              * std::exp((alpha_square - alpha * s) / 4);

        T z;
        if (beta > cutoff) {
            // [alpha, beta] is wide. Sample from alpha + E(lambda)
            // and proceed as above (see Robert (1995)).
            std::exponential_distribution<T> exp_sampler(lambda);

            T d;
            do {
                do {
                    z = alpha + exp_sampler(random_gen);
                } while (z > beta);
                d = z - lambda;
            } while (u_sampler_01(random_gen) > std::exp(-d * d / 2));
        } else {
            // [alpha, beta] is narrow. Sample from U(alpha, beta)
            // and proceed as above (see Robert (1995)).
            std::uniform_real_distribution<T> u_sampler_ab(alpha, beta);
            do {
                z = u_sampler_ab(random_gen);
            } while (u_sampler_01(random_gen) > std::exp((alpha_square - z * z) / 2));
        }
        return z;
    }

    /**
     * @brief Canonize a value
     * 
     * This method computes the value \f$ v\f$ whose
     * probability in \f$ N(0, 1)\f$ equals the probability
     * of @p value in \f$ N(\mu, \sigma^2)\f$.
     * 
     * @param value is the value to be canonized
     * @return the value \f$ v\f$ whose probability in
     *   \f$ N(0, 1)\f$ equals the probability of @p value
     *   in \f$ N(\mu, \sigma^2)\f$.
     */
    inline T canonize(const T value) const
    {
        return (value - get_mu()) / get_sigma();
    }

    /**
     * @brief Lift a value
     * 
     * This method computes the value \f$ v\f$ whose
     * probability in \f$ N(\mu, \sigma^2)\f$ equals the 
     * probability of @p value in \f$ N(0, 1)\f$.
     * 
     * @param value is the value to be canonized
     * @return the value \f$ v\f$ whose probability in 
     *   \f$ N(\mu, \sigma^2)\f$ equals the probability
     *   of @p value in \f$ N(0, 1)\f$.
     */
    inline T lift(const T value) const
    {
        return get_mu() + value * get_sigma();
    }

    /**
     * @brief Mean excess and variance of the standard normal distribution
     *        truncated to \f$[\alpha, \alpha + w]\f$
     *
     * Let \f$ Z \sim N(0,1) \f$ conditioned on
     * \f$ Z \in [\alpha, \alpha + w] \f$. This method computes
     * \f$ E[Z] - \alpha \f$ and \f$ Var[Z] \f$ by Gauss-Legendre
     * quadrature of the density about the point \f$ z_0 \f$ of the
     * interval closest to 0, where the mass concentrates. The integrand is
     * normalized by \f$ \varphi(z_0) \f$, so it does not underflow in the
     * tails, and the quadrature nodes are expressed as distances from
     * \f$ z_0 \f$, so the moments do not suffer from cancellation even when
     * \f$ |\alpha| \f$ is large or \f$ w \f$ is small.
     *
     * @param alpha is the lower cut. It must be finite.
     * @param width is the width \f$ w \f$ of the truncation interval. It
     *      must be positive and it may be \f$ +\infty \f$.
     * @return the array \f$ [E[Z] - \alpha, Var[Z]] \f$.
     */
    static std::array<T, 2> std_truncated_excess_and_variance(const T alpha, const T width)
    {
        // 8-point Gauss-Legendre nodes and weights on [-1, 1]
        constexpr std::array<T, 4> gl_nodes{0.1834346424956498, 0.5255324099163290,
                                            0.7966664774136267, 0.9602898564975363};
        constexpr std::array<T, 4> gl_weights{0.3626837833783620, 0.3137066458778873,
                                              0.2223810344533745, 0.1012285362903763};
        constexpr size_t num_of_panels = 32;
        constexpr size_t num_of_points = num_of_panels * 2 * gl_nodes.size();

        // z0 is the point of [alpha, alpha + width] closest to 0, and
        // [d_low, d_high] is the interval as distances from z0
        T z0, d_low, d_high;
        if (alpha >= 0) {
            z0 = alpha;
            d_low = 0;
            d_high = width;
        } else if (alpha + width <= 0) {
            z0 = alpha + width;
            d_low = -width;
            d_high = 0;
        } else {
            z0 = 0;
            d_low = alpha;
            d_high = alpha + width;
        }

        // from z0, the density decays at least as fast as exp(-|z0| d) and
        // exp(-d^2/2), so its mass farther than half_width is negligible
        const T half_width = std::min(T(12), 40 / std::max(T(1), std::abs(z0)));
        d_low = std::max(d_low, -half_width);
        d_high = std::min(d_high, half_width);

        const T panel_size = (d_high - d_low) / num_of_panels;

        std::array<T, num_of_points> d, w;
        size_t i = 0;
        for (size_t panel = 0; panel < num_of_panels; ++panel) {
            const T center = d_low + (panel + T(0.5)) * panel_size;
            for (size_t j = 0; j < gl_nodes.size(); ++j) {
                for (const T sign : {T(-1), T(1)}) {
                    // the distance from z0 and phi(z0 + d)/phi(z0)
                    d[i] = center + sign * gl_nodes[j] * panel_size / 2;
                    w[i] = gl_weights[j] * std::exp(-d[i] * (d[i] + 2 * z0) / 2);
                    ++i;
                }
            }
        }

        T mass = 0, d_mean = 0;
        for (i = 0; i < num_of_points; ++i) {
            mass += w[i];
            d_mean += w[i] * d[i];
        }
        d_mean /= mass;

        T variance = 0;
        for (i = 0; i < num_of_points; ++i) {
            variance += w[i] * (d[i] - d_mean) * (d[i] - d_mean);
        }

        // E[Z] - alpha = (z0 - alpha) + E[Z - z0]
        const T z0_excess = (alpha >= 0 ? 0 : (alpha + width <= 0 ? width : -alpha));

        return {z0_excess + d_mean, variance / mass};
    }

    /**
     * @brief Fits a left-truncated normal distribution
     *
     * See `fit_untruncated()`.
     *
     * @param mean is desired mean of the truncated distribution.
     * @param stddev is desired standard deviation of the truncated distribution.
     *      It must be positive.
     * @param a is left truncation value. It must be finite and smaller than
     *      @p mean.
     * @param max_alpha is the largest standardized cut that can be fitted.
     * @return the (non-truncated) normal distribution with parameters
     *         \f$ \mu \f$ and \f$ \sigma \f$.
     */
    static std::normal_distribution<T>
    fit_left_truncated(const T mean, const T stddev, const T a, const T max_alpha)
    {
        // the ratio stddev/(mean - a) of a left-truncated distribution
        // increases from 0 (alpha -> -inf) to 1 (alpha -> +inf, the
        // distribution looks exponential), and 1 - ratio ~ 1/alpha^2
        const T target_ratio = stddev / (mean - a);

        if (target_ratio >= 1) {
            throw CLONES::Error<std::invalid_argument>("The aimed standard deviation must be "
                                                       "smaller than the difference between "
                                                       "the aimed mean and the cut.");
        }

        auto eq = [target_ratio](const T alpha) -> T {
            const auto [excess, variance]
                = std_truncated_excess_and_variance(alpha, std::numeric_limits<T>::infinity());

            return std::sqrt(variance) / excess / target_ratio - 1;
        };

        constexpr T alpha_min = -10;

        // eq increases: if it is non-negative at alpha_min, the cut is more
        // than 10 sigma below mu and the truncation is numerically irrelevant
        if (eq(alpha_min) >= 0) {
            return std::normal_distribution<T>(mean, stddev);
        }

        // since 1 - ratio ~ 1/alpha^2, the root is about 1/sqrt(1 - target_ratio)
        T alpha_max = 8;
        while (eq(alpha_max) < 0) {
            alpha_max *= 2;

            if (alpha_max > max_alpha) {
                throw CLONES::Error<std::invalid_argument>("The aimed standard deviation is "
                                                           "too close to the difference "
                                                           "between the aimed mean and the cut.");
            }
        }

        const T alpha = brent_root<decltype(eq), T, T>(eq, alpha_min, alpha_max);

        const auto excess = std_truncated_excess_and_variance(alpha,
                                                              std::numeric_limits<T>::infinity())[0];
        const T sigma = (mean - a) / excess;

        return std::normal_distribution<T>(a - alpha * sigma, sigma);
    }

    /**
     * @brief Finds the standardized lower cut that gives a mean excess
     *
     * @param width is the standardized width \f$ \beta - \alpha \f$ of the
     *      truncation interval. It must be positive and finite.
     * @param excess is the aimed mean excess. It must be in
     *      \f$ (0, width) \f$.
     * @return the standardized lower cut \f$ \alpha \f$ such that the
     *      standard normal distribution truncated to
     *      \f$[\alpha, \alpha + width]\f$ has mean \f$ \alpha + excess \f$.
     */
    static T find_alpha(const T width, const T excess)
    {
        // the mean excess decreases from width (alpha -> -inf) to 0
        // (alpha -> +inf), and it is width/2 when alpha = -width/2
        auto eq = [width, excess](const T alpha) -> T {
            return std_truncated_excess_and_variance(alpha, width)[0] / excess - 1;
        };

        const T alpha_0 = -width / 2;
        const T direction = (eq(alpha_0) > 0 ? 1 : -1);

        T alpha_1 = alpha_0, step = 1;
        do {
            alpha_1 = alpha_0 + direction * step;
            step *= 2;

            if (!std::isfinite(step)) {
                throw CLONES::Error<std::runtime_error>("The root search failed.");
            }
        } while (direction * eq(alpha_1) > 0);

        return brent_root<decltype(eq), T, T>(eq, std::min(alpha_0, alpha_1),
                                                      std::max(alpha_0, alpha_1));
    }

    /**
     * @brief Fits a normal distribution truncated to a finite interval
     *
     * See `fit_untruncated()`.
     *
     * @param mean is desired mean of the truncated distribution.
     * @param stddev is desired standard deviation of the truncated distribution.
     *      It must be positive.
     * @param a is left truncation value. It must be finite and smaller than
     *      @p mean.
     * @param b is right truncation value. It must be finite and greater than
     *      @p mean.
     * @param max_sigma_ratio is the largest ratio between \f$ \sigma \f$ and
     *      @p stddev that can be fitted.
     * @return the (non-truncated) normal distribution with parameters
     *         \f$ \mu \f$ and \f$ \sigma \f$.
     */
    static std::normal_distribution<T>
    fit_truncated(const T mean, const T stddev, const T a, const T b,
                  const T max_sigma_ratio)
    {
        // truncating to [a,b] gives a smaller standard deviation than
        // truncating to [a, +inf) or (-inf, b]
        if (stddev >= mean - a || stddev >= b - mean) {
            throw CLONES::Error<std::invalid_argument>("The aimed standard deviation must be "
                                                       "smaller than the differences between "
                                                       "the aimed mean and the cuts.");
        }

        const T interval_size = b - a;

        // The fitted distribution is parameterized by the logarithm of the
        // standardized width w = (b - a)/sigma. For every w, alpha(w) is
        // chosen to get the aimed mean.
        auto get_alpha = [mean, a, interval_size](const T width) {
            return find_alpha(width, (mean - a) * width / interval_size);
        };

        // the standard deviation increases in sigma, i.e., decreases in w
        auto eq = [&get_alpha, stddev, interval_size](const T log_width) -> T {
            const T width = std::exp(log_width);
            const T alpha = get_alpha(width);
            const T variance = std_truncated_excess_and_variance(alpha, width)[1];

            return interval_size * std::sqrt(variance) / (width * stddev) - 1;
        };

        // truncating reduces the standard deviation, so sigma >= stddev
        const T log_width_max = std::log(interval_size / stddev);

        if (eq(log_width_max) >= 0) {
            return std::normal_distribution<T>(a - get_alpha(std::exp(log_width_max)) * stddev,
                                               stddev);
        }

        T log_width_min = log_width_max;
        {
            const T log_width_limit = log_width_max - std::log(max_sigma_ratio);
            do {
                log_width_min -= std::log(T(2));

                if (log_width_min < log_width_limit) {
                    throw CLONES::Error<std::invalid_argument>("The aimed standard deviation "
                                                               "is too large for the "
                                                               "truncation interval.");
                }
            } while (eq(log_width_min) < 0);
        }

        const T log_width = brent_root<decltype(eq), T, T>(eq, log_width_min,
                                                                   log_width_max);

        const T width = std::exp(log_width);
        const T sigma = interval_size / width;

        return std::normal_distribution<T>(a - get_alpha(width) * sigma, sigma);
    }

    /**
     * @brief Get the moments of the truncated distribution about a reference point
     *
     * Let \f$ Z \sim N(0,1) \f$ conditioned on \f$ Z \in [\alpha, \beta] \f$.
     * This method selects a reference point \f$ x_r \f$ among @p a, @p b,
     * and \f$ \mu \f$ near which the mass of the truncated distribution
     * concentrates, and it computes the standardized shift
     * \f$ s \f$ of the mean from \f$ x_r \f$ and \f$ Var[Z] \f$, so that
     * the mean is \f$ x_r + \sigma s \f$ and the standard deviation is
     * \f$ \sigma\sqrt{Var[Z]} \f$. Since \f$ \sigma s \f$ is small with
     * respect to \f$ x_r \f$, the mean does not suffer from the cancellation
     * of \f$ \mu + \sigma E[Z] \f$, which occurs when the cuts are far from
     * \f$ \mu \f$.
     *
     * The distribution must have a positive standard deviation.
     *
     * @return the array \f$ [x_r, s, Var[Z]] \f$.
     */
    std::array<T, 3> get_moments_about_reference() const
    {
        // the interval is in the right tail: the mass is close to a
        if (alpha >= 0) {
            const auto [excess, variance]
                = std_truncated_excess_and_variance(alpha, (b - a) / sigma);

            return {a, excess, variance};
        }

        // the interval is in the left tail: the mass is close to b
        if (beta <= 0) {
            const auto [excess, variance]
                = std_truncated_excess_and_variance(-beta, (b - a) / sigma);

            return {b, -excess, variance};
        }

        // 0 is in [alpha, beta]: the mass is close to mu and it is negligible
        // farther than 40 standard deviations, so the cuts can be clipped
        // there without affecting the moments. This also avoids the
        // cancellation between E[Z] - alpha and alpha when alpha is very
        // negative.
        constexpr T max_distance = 40;
        const T clipped_alpha = std::max(alpha, -max_distance);
        const T clipped_beta = std::min(beta, max_distance);

        const auto [excess, variance]
            = std_truncated_excess_and_variance(clipped_alpha, clipped_beta - clipped_alpha);

        return {mu, clipped_alpha + excess, variance};
    }

public:
    using result_type = T;

    /**
     * @brief A constructor
     *
     * @param[in] mu is the mean of the non-truncated normal distribution
     * @param[in] sigma is the standard deviation of the non-truncated normal
     *      distribution
     * @param[in] a is the lower truncation value
     * @param[in] b is the upper truncation value
     * @throw std::invalid_argument if @p sigma is negative, if @p a is
     *      greater than @p b, or if @p sigma is null and @p mu is lower
     *      than @p a or greater than @p b
     */
    explicit truncated_normal_distribution(const T mu = 0, const T sigma = 0,
                                           const T a = -std::numeric_limits<T>::infinity(),
                                           const T b = std::numeric_limits<T>::infinity()):
        mu{mu}, sigma{sigma}, a{a}, b{b}
    {
        if (sigma < 0) {
            throw CLONES::Error<std::invalid_argument>("The standard deviation "
                                                       "must be non-negative.");
        }

        if (a > b) {
            throw CLONES::Error<std::invalid_argument>("The lower cut must not be "
                                                       "greater than the upper cut.");
        }

        if (sigma == 0 && (mu < a || mu > b) ) {
            throw CLONES::Error<std::invalid_argument>("When the standard deviation "
                                                       "is null, the mean must "
                                                       "be in the interval [a,b].");
        }

        if (sigma != 0) {
            alpha = canonize(a);
            beta = canonize(b);
        } else {
            alpha = a - mu;
            beta = b - mu;
        }
    }

    /**
     * @brief Build the truncated normal distribution of a normal distribution
     *
     * @param[in] normal_distribution is a normal distribution
     * @param[in] a is the lower truncation value
     * @param[in] b is the upper truncation value
     */
    explicit truncated_normal_distribution(const std::normal_distribution<T>& normal_distribution,
                                           const T a = std::numeric_limits<T>::lowest(),
                                           const T b = std::numeric_limits<T>::max()):
        truncated_normal_distribution{normal_distribution.mean(),
                                      normal_distribution.stddev(), a, b}
    {}

    /**
     * @brief Sample a value
     *
     * @tparam RANDOM_GENERATOR is the type of the random generator
     * @param[in, out] generator is a random generator
     * @return a value sampled according to the distribution. It is
     *      always greater than or equal to the cut
     */
    template<typename RANDOM_GENERATOR>
    T operator()(RANDOM_GENERATOR& generator) const
    {
        if (sigma == 0) {
            return mu;
        }

        return lift(canonical_sampler(generator, alpha, beta));
    }

    /**
     * @brief Get the mean of the non-truncated normal distribution
     *
     * @return a constant reference to the mean of the non-truncated
     *    normal distribution
     */
    inline const T& get_mu() const
    {
        return mu;
    }

    /**
     * @brief Get the standard deviation of the non-truncated normal distribution
     *
     * @return a constant reference to the standard deviation of the non-truncated
     *    normal distribution
     */
    inline const T& get_sigma() const
    {
        return sigma;
    }

    /**
     * @brief Get the lower cut value
     *
     * @return a constant reference to the lower cut value
     */
    inline const T& get_a() const
    {
        return a;
    }

    /**
     * @brief Get the upper cut value
     *
     * @return a constant reference to the upper cut value
     */
    inline const T& get_b() const
    {
        return b;
    }

    /**
     * @brief Get the mean of the truncated distribution
     *
     * @return the mean of the values sampled by the distribution
     */
    T mean() const
    {
        if (sigma == 0) {
            return mu;
        }

        const auto [reference, shift, variance] = get_moments_about_reference();

        return reference + sigma * shift;
    }

    /**
     * @brief Get the standard deviation of the truncated distribution
     *
     * @return the standard deviation of the values sampled by the distribution
     */
    T stddev() const
    {
        if (sigma == 0) {
            return 0;
        }

        return sigma * std::sqrt(get_moments_about_reference()[2]);
    }

    /**
     * @brief Finds the normal distribution that, truncated to \f$[a,b]\f$,
     *        has the given mean and standard deviation.
     *
     * Looks for the parameters \f$ \mu \f$ and \f$ \sigma \f$ of a normal
     * distribution such that the distribution obtained by truncating it in
     * \f$[a,b]\f$ has mean @p mean and standard deviation @p stddev.
     *
     * Let \f$ \alpha = (a - \mu)/\sigma \f$ and \f$ \beta = (b - \mu)/\sigma \f$
     * be the standardized cuts, and let \f$ \delta(\alpha, \beta) \f$ and
     * \f$ v(\alpha, \beta) \f$ be the mean excess over \f$ \alpha \f$ and the
     * variance of the standard normal distribution truncated to
     * \f$[\alpha, \beta]\f$. Then, the truncated distribution has mean
     * \f$ a + \sigma\delta(\alpha, \beta) \f$ and standard deviation
     * \f$ \sigma\sqrt{v(\alpha, \beta)} \f$.
     *
     * - When @p b is \f$ +\infty \f$, the ratio
     *   \f$ \sqrt{v(\alpha, \infty)}/\delta(\alpha, \infty) \f$ only depends
     *   on \f$ \alpha \f$, and it increases from 0 to 1. So, \f$ \alpha \f$ is
     *   the root of \f$ \sqrt{v}/\delta = stddev/(mean - a) \f$, and
     *   \f$ \sigma = (mean - a)/\delta(\alpha, \infty) \f$. When @p a is
     *   \f$ -\infty \f$, the problem is mirrored.
     * - When both @p a and @p b are finite, for any standardized width
     *   \f$ w = (b - a)/\sigma \f$, the mean decreases in \f$ \alpha \f$, so
     *   \f$ \alpha(w) \f$ is found by a root search. Then, \f$ w \f$ is the
     *   root of \f$ \sigma\sqrt{v(\alpha(w), \alpha(w) + w)} = stddev \f$,
     *   whose left-hand side increases in \f$ \sigma \f$.
     *
     * When the cuts are more than 10 @p stddev away from @p mean, the
     * truncation is numerically irrelevant and \f$ N(mean, stddev^2) \f$ is
     * returned.
     *
     * @param mean is desired mean of the truncated distribution.
     * @param stddev is desired standard deviation of the truncated distribution.
     * @param a is left truncation value.
     * @param b is right truncation value.
     * @return the (non-truncated) normal distribution with parameters
     *         \f$ \mu \f$ and \f$ \sigma \f$.
     * @throw std::invalid_argument if @p stddev is negative, if @p a is greater
     *        than @p b, if @p mean is not in \f$[a,b]\f$, or it is equal to
     *        @p a or @p b while @p stddev is positive, or if no normal
     *        distribution truncated to \f$[a,b]\f$ has standard deviation
     *        @p stddev and mean @p mean, or the fitted one is too far from
     *        the cuts to be computed.
     * @throw std::runtime_error if the root search fails.
     */
    static std::normal_distribution<T>
    fit_untruncated(const T mean, const T stddev,
                    const T a = -std::numeric_limits<T>::infinity(),
                    const T b = std::numeric_limits<T>::infinity())
    {
        if (stddev < 0) {
            throw CLONES::Error<std::invalid_argument>("The standard deviation "
                                                       "must be non-negative.");
        }

        if (a > b) {
            throw CLONES::Error<std::invalid_argument>("The lower cut must not be "
                                                       "greater than the upper cut.");
        }

        // a distribution truncated to [a,b] only takes values in [a,b], so
        // its mean is in [a,b] too, and it is a or b only if stddev is null
        if (mean < a || mean > b || ((mean == a || mean == b) && stddev > 0)) {
            throw CLONES::Error<std::invalid_argument>("The aimed mean must be "
                                                       "in the interval (a,b) or, "
                                                       "when the standard deviation "
                                                       "is null, in [a,b].");
        }

        if (stddev == 0) {
            return std::normal_distribution<T>(mean, stddev);
        }

        // the truncation is numerically irrelevant
        constexpr T negligible_distance = 10;
        if (mean - a >= negligible_distance * stddev
                && b - mean >= negligible_distance * stddev) {
            return std::normal_distribution<T>(mean, stddev);
        }

        // reduce to the case in which a is finite by mirroring
        if (std::isinf(a)) {
            const auto mirrored = fit_untruncated(-mean, stddev, -b, -a);

            return std::normal_distribution<T>(-mirrored.mean(), mirrored.stddev());
        }

        // 1e7 is the largest ratio between sigma and stddev that the
        // root searches below can handle
        constexpr T max_sigma_ratio = 1e7;

        if (std::isinf(b)) {
            return fit_left_truncated(mean, stddev, a, max_sigma_ratio);
        }

        return fit_truncated(mean, stddev, a, b, max_sigma_ratio);
    }
};

} // CLONES

#endif // __CLONES_DISTRIBUTIONS__
