/**
 * @file distributions.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines functions about probability distributions
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

#ifndef __CLONES_DISTRIBUTIONS__
#define __CLONES_DISTRIBUTIONS__

#include <cmath>
#include <stdexcept>
#include <random>
#include <numbers>

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
    return std::exp(-0.5 * z * z) / std::sqrt(2.0 * std::numbers::pi);
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
T norm_cdf(const T z)
{
    return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
}

/**
 * @brief Finds the normal distribution that, left-truncated at @p cut, has
 *        the given mean and standard deviation.
 *
 * Looks for the parameters \f$ \mu \f$ and \f$ \sigma \f$ of a normal
 * distribution such that the distribution obtained by left-truncating it at
 * @p cut has mean @p mean and standard deviation @p stddev.
 *
 * Letting \f$ \alpha = (cut - \mu)/\sigma \f$ and
 * \f$ \lambda(\alpha) = \varphi(\alpha) / (1 - \Phi(\alpha)) \f$
 * (the inverse Mills ratio), the function numerically solves for
 * \f$ \alpha \f$ the equation
 * \f[ \frac{\sqrt{1 - \lambda(\alpha)(\lambda(\alpha) - \alpha)}}
 *        {\lambda(\alpha) - \alpha} = \frac{aimed\_stddev}{aimed\_mean - cut} \f]
 * in the interval [-10, 8], and then derives
 * \f$ \sigma = (aimed\_mean - cut)/(\lambda(\alpha) - \alpha) \f$ and
 * \f$ \mu = cut - \alpha\sigma \f$.
 *
 * @tparam T floating-point type.
 * @param mean is desired mean of the truncated distribution.
 * @param stddev is desired standard deviation of the truncated distribution.
 * @param cut is left truncation value.
 * @return the (non-truncated) normal distribution with parameters
 *         \f$ \mu \f$ and \f$ \sigma \f$.
 * @throw std::invalid_argument if @p mean is not greater than @p cut.
 * @throw std::runtime_error if the root search fails.
 */
template<typename T>
std::normal_distribution<T>
fit_left_truncated_normal(const T mean, const T stddev, const T cut)
{
    if (stddev < 0) {
        throw CLONES::Error<std::invalid_argument>("The standard deviation "
                                                   "must be non-negative.");
    }

    // A distribution left-truncated at cut only takes values greater than
    // cut, so its mean must be greater than cut too
    if (mean < cut || (mean == cut && stddev > 0)) {
        throw CLONES::Error<std::invalid_argument>("The aimed mean must be "
                                                   "greater than the cut or, "
                                                   "when the standard deviation "
                                                   "is null, equal it.");
    }

    // If the standard deviation is 0 and mean >= cut, then
    // the aimed distribution is N(mean, 0)
    if (stddev == 0) {
        return std::normal_distribution<T>(mean, stddev);
    }

    // For X ~ N(mu, sigma^2) truncated to X > cut, with
    // alpha = (cut - mu) / sigma and lambda = pdf(alpha) / (1 - cdf(alpha)):
    //
    //   E[X]   = mu + sigma * lambda      =>  E[X] - cut = sigma * (lambda - alpha)
    //   Var[X] = sigma^2 * (1 - lambda * (lambda - alpha))
    //
    // Hence, the ratio stddev / (mean - cut) depends on alpha alone:
    //
    //   ratio(alpha) = sqrt(1 - lambda * (lambda - alpha)) / (lambda - alpha)
    //
    // ratio(alpha) increases from 0 (alpha -> -inf, almost no truncation) to
    // 1 (alpha -> +inf, the tail looks exponential). So, alpha can be found
    // as the only root of ratio(alpha) - target_ratio, provided that
    // 0 < target_ratio < 1.
    const T target_ratio = stddev / (mean - cut);

    // ratio(alpha) < 1 for every alpha, so no root exists
    if (target_ratio >= 1) {
        throw CLONES::Error<std::invalid_argument>("The aimed standard deviation must be "
                                                   "smaller than the difference between "
                                                   "the aimed mean and the cut.");
    }

    // eq(alpha) = ratio(alpha) - target_ratio. The guards below deal with
    // the loss of precision of 1 - norm_cdf(alpha) for large alpha.
    auto eq = [target_ratio](const T& alpha) {
        // outside of the search interval: never evaluated by the root finder
        if (alpha > 8) {
            return -target_ratio;
        }

        // the tail probability is negligible: return a positive value, since
        // ratio(alpha) is close to 1 > target_ratio there
        const auto denom_cdf = 1 - norm_cdf(alpha);
        if (denom_cdf <= 1e-15) {
            return static_cast<T>(1);
        }

        // lambda is the inverse Mills ratio; diff = (E[X] - cut) / sigma,
        // which is positive in exact arithmetic
        const auto lambda_a = norm_pdf(alpha) / denom_cdf;
        const auto diff = lambda_a - alpha;

        // avoid dividing by (almost) zero
        if (std::abs(diff) < 1e-12) {
            return static_cast<T>(1e12);
        }

        // term = Var[X] / sigma^2, which lies in (0, 1) in exact arithmetic;
        // it can be negative only because of rounding errors
        const auto term = 1 - lambda_a * diff;
        if (term < 0) {
            return -target_ratio;
        }

        return (std::sqrt(term) / diff) - target_ratio;
    };

    constexpr T alpha_min = -10;

    // eq is increasing: if it is non-negative at alpha_min, the root lies below
    // alpha_min, i.e., the cut is more than 10 sigma below mu and truncation is
    // numerically irrelevant
    if (eq(alpha_min)>=0) {
        return std::normal_distribution<T>(mean, stddev);
    }

    // Numerically search for alpha in the safe interval [-10.0, 8.0]
    const auto alpha = brent_root(eq, alpha_min, 8.0);

    // Compute the final parameters: sigma from E[X] - cut = sigma * (lambda - alpha)
    // and mu from the definition of alpha
    const auto cdf = norm_cdf(alpha);
    const auto pdf = norm_pdf(alpha);
    const auto lambda_a = pdf / (1 - cdf);

    const auto sigma = (mean - cut) / (lambda_a - alpha);
    const auto mu = cut - alpha * sigma;

    return std::normal_distribution(mu, sigma);
}

/**
 * @brief A left-truncated normal distribution
 *
 * This class samples \f$ X \sim N(\mu, \sigma^2) \f$ conditioned on
 * \f$ X \geq cut \f$. Let \f$ \alpha = (cut - \mu)/\sigma \f$ be the
 * standardized cut. When \f$ \alpha < 0 \f$, samples are produced by
 * rejecting the standard normal samples below \f$ \alpha \f$ (acceptance
 * probability \f$ \geq 1/2 \f$). Otherwise, the exponential-proposal
 * rejection method of Robert (1995) is used (acceptance probability
 * \f$ \gtrsim 0.76 \f$). Hence, sampling is efficient for any cut.
 *
 * C. P. Robert. "Simulation of truncated normal variables". Statistics
 * and Computing 5, 121–125 (1995).
 *
 * @tparam T is the floating-point type
 */
template<typename T = double>
class left_truncated_normal_distribution
{
    T mu;       //!< The mean of the non-truncated normal distribution
    T sigma;    //!< The standard deviation of the non-truncated normal distribution
    T cut;      //!< The truncation value

    T alpha;    //!< The standardized cut, i.e., \f$ (cut - \mu)/\sigma \f$
    T lambda;   //!< The optimal rate of the exponential proposal when \f$ \alpha \geq 0 \f$

public:
    using result_type = T;

    /**
     * @brief A constructor
     *
     * @param[in] mu is the mean of the non-truncated normal distribution
     * @param[in] sigma is the standard deviation of the non-truncated normal
     *      distribution
     * @param[in] cut is the truncation value
     * @throw std::invalid_argument if @p sigma is negative or if @p sigma
     *      is null and @p mu is smaller than @p cut
     */
    explicit left_truncated_normal_distribution(const T mu = 0, const T sigma = 0,
                                                const T cut = 0):
        mu{mu}, sigma{sigma}, cut{cut},
        alpha{sigma > 0 ? (cut - mu) / sigma : 0},
        lambda{(alpha + std::sqrt(alpha * alpha + 4)) / 2}
    {
        if (sigma < 0) {
            throw CLONES::Error<std::invalid_argument>("The standard deviation "
                                                       "must be non-negative.");
        }

        if (sigma == 0 && mu < cut) {
            throw CLONES::Error<std::invalid_argument>("When the standard deviation "
                                                       "is null, the mean must not "
                                                       "be smaller than the cut.");
        }
    }

    /**
     * @brief Build the left-truncated normal distribution of a normal distribution
     *
     * @param[in] normal_distribution is a normal distribution
     * @param[in] cut is the truncation value
     */
    left_truncated_normal_distribution(const std::normal_distribution<T>& normal_distribution,
                                       const T cut):
        left_truncated_normal_distribution{normal_distribution.mean(),
                                           normal_distribution.stddev(), cut}
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

        // naive rejection: the acceptance probability is at least 1/2
        if (alpha < 0) {
            std::normal_distribution<T> std_normal;
            T z;
            do {
                z = std_normal(generator);
            } while (z < alpha);

            // max() guards against rounding errors in mu + sigma * z
            return std::max(cut, mu + sigma * z);
        }

        // exponential proposal z = alpha + E, E ~ Exp(lambda), accepted
        // with probability exp(-(z - lambda)^2 / 2)
        std::exponential_distribution<T> exp_dist(lambda);
        std::uniform_real_distribution<T> uniform_dist(0, 1);
        while (true) {
            const T z = alpha + exp_dist(generator);
            const T d = z - lambda;
            if (uniform_dist(generator) <= std::exp(-d * d / 2)) {
                return std::max(cut, mu + sigma * z);
            }
        }
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
     * @brief Get the truncation value
     *
     * @return a constant reference to the truncation value
     */
    inline const T& get_cut() const
    {
        return cut;
    }

    /**
     * @brief Get the mean of the left-truncated distribution
     *
     * The mean is \f$ \mu + \sigma\lambda(\alpha) \f$, where
     * \f$ \lambda(\alpha) = \varphi(\alpha) / (1 - \Phi(\alpha)) \f$ is the
     * inverse Mills ratio of the standardized cut \f$ \alpha \f$.
     *
     * @return the mean of the values sampled by the distribution
     */
    T mean() const
    {
        if (sigma == 0) {
            return mu;
        }

        return mu + sigma * inverse_mills_ratio();
    }

    /**
     * @brief Get the standard deviation of the left-truncated distribution
     *
     * The standard deviation is
     * \f$ \sigma\sqrt{1 - \lambda(\alpha)(\lambda(\alpha) - \alpha)} \f$, where
     * \f$ \lambda(\alpha) = \varphi(\alpha) / (1 - \Phi(\alpha)) \f$ is the
     * inverse Mills ratio of the standardized cut \f$ \alpha \f$.
     *
     * @return the standard deviation of the values sampled by the distribution
     */
    T stddev() const
    {
        if (sigma == 0) {
            return 0;
        }

        // for large alpha, 1 - lambda * (lambda - alpha) suffers from
        // cancellation: use its asymptotic expansion
        // 1/alpha^2 - 6/alpha^4 + 50/alpha^6
        if (alpha > asymptotic_alpha) {
            const T u2 = 1 / (alpha * alpha);

            return sigma * std::sqrt(u2 * (1 - u2 * (6 - 50 * u2)));
        }

        const T lambda_a = inverse_mills_ratio();

        // the variance ratio lies in (0, 1): max() guards against rounding errors
        return sigma * std::sqrt(std::max(T(0), 1 - lambda_a * (lambda_a - alpha)));
    }

private:
    //! The standardized cut above which asymptotic expansions are used
    static constexpr T asymptotic_alpha = 30;

    /**
     * @brief Compute the inverse Mills ratio of the standardized cut
     *
     * @return \f$ \lambda(\alpha) = \varphi(\alpha) / (1 - \Phi(\alpha)) \f$
     */
    T inverse_mills_ratio() const
    {
        // for large alpha, phi(alpha) and 1 - Phi(alpha) underflow: use the
        // asymptotic expansion alpha + 1/alpha - 2/alpha^3
        if (alpha > asymptotic_alpha) {
            return alpha + 1 / alpha - 2 / (alpha * alpha * alpha);
        }

        // erfc avoids the cancellation of 1 - norm_cdf(alpha) for large alpha
        const T tail = std::erfc(alpha / std::sqrt(T(2))) / 2;

        return norm_pdf(alpha) / tail;
    }
};

} // CLONES

#endif // __CLONES_DISTRIBUTIONS__
