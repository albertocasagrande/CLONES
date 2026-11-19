/**
 * @file distributions.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines functions about probability distributions
 * @version 1.3
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

#include <algorithm>
#include <array>
#include <cstddef>
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

    //! The standardized cut from which the continued fractions are used
    static constexpr T continued_fraction_alpha = 3;
    
    /**
     * @brief Evaluates the continued fraction of the left-truncated standard
     *        normal distribution.
     *
     * Computes the first three terms of the sequence
     * \f$ T_n = \alpha + (n+1)/T_{n+1} \f$, whose limit satisfies
     * \f$ \lambda(\alpha) - \alpha = 1/T_1 \f$, where \f$ \lambda(\alpha) \f$
     * is the inverse Mills ratio. The continued fraction is truncated after
     * 64 terms, which are enough to reach the double precision for
     * \f$ \alpha \geq 3 \f$.
     *
     * @param alpha is the standardized cut. It must be positive.
     * @return the array \f$ [T_1, T_2, T_3] \f$.
     */
    static std::array<T, 3> left_truncated_std_normal_cf(const T alpha)
    {
        constexpr size_t num_of_terms = 64;
    
        std::array<T, 3> terms;
        T t = alpha;
        for (size_t n = num_of_terms; n > 0; --n) {
            t = alpha + static_cast<T>(n + 1) / t;
            if (n <= terms.size()) {
                terms[n - 1] = t;
            }
        }
    
        return terms;
    }
    
    /**
     * @brief Mean excess of the standard normal distribution left-truncated at
     *        @p alpha.
     *
     * Computes \f$ E[Z - \alpha \mid Z \geq \alpha] = \lambda(\alpha) - \alpha \f$,
     * where \f$ Z \sim N(0,1) \f$ and
     * \f$ \lambda(\alpha) = \varphi(\alpha) / (1 - \Phi(\alpha)) \f$ is the
     * inverse Mills ratio. For \f$ \alpha \geq 3 \f$, where computing
     * \f$ \lambda(\alpha) - \alpha \f$ suffers from cancellation and
     * \f$ \varphi(\alpha) \f$ eventually underflows, the continued fraction
     * \f$ 1/(\alpha + 2/(\alpha + 3/(\alpha + \ldots))) \f$ is used.
     *
     * @param alpha is the standardized cut.
     * @return the mean excess \f$ \lambda(\alpha) - \alpha \f$, which is positive.
     */
    static T left_truncated_std_normal_excess(const T alpha)
    {
        if (alpha >= continued_fraction_alpha) {
            return 1 / left_truncated_std_normal_cf(alpha)[0];
        }
    
        // erfc avoids the cancellation of 1 - norm_cdf(alpha)
        const T tail = std::erfc(alpha / std::sqrt(T(2))) / 2;
    
        return norm_pdf(alpha) / tail - alpha;
    }
    
    /**
     * @brief Variance of the standard normal distribution left-truncated at
     *        @p alpha.
     *
     * Computes \f$ Var[Z \mid Z \geq \alpha] =
     * 1 - \lambda(\alpha)(\lambda(\alpha) - \alpha) \f$, where
     * \f$ Z \sim N(0,1) \f$ and \f$ \lambda(\alpha) \f$ is the inverse Mills
     * ratio. For \f$ \alpha \geq 3 \f$, where this expression suffers from
     * cancellation, the variance is computed as
     * \f$ (\alpha + 4/T_2 - 3/T_3)/(T_2 T_1^2) \f$, where
     * \f$ T_n = \alpha + (n+1)/T_{n+1} \f$ are the terms of the continued
     * fraction of \f$ \lambda(\alpha) - \alpha = 1/T_1 \f$.
     *
     * @param alpha is the standardized cut.
     * @return the variance, which lies in \f$ (0, 1) \f$.
     */
    static T left_truncated_std_normal_variance(const T alpha)
    {
        if (alpha >= continued_fraction_alpha) {
            const auto [t1, t2, t3] = left_truncated_std_normal_cf(alpha);
    
            return (alpha + 4 / t2 - 3 / t3) / (t2 * t1 * t1);
        }
    
        const T excess = left_truncated_std_normal_excess(alpha);
    
        // max() guards against rounding errors
        return std::max(T(0), 1 - (alpha + excess) * excess);
    }

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

        return cut + sigma * left_truncated_std_normal_excess(alpha);
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

        return sigma * std::sqrt(left_truncated_std_normal_variance(alpha));
    }

    /**
     * @brief Finds the normal distribution that, left-truncated at @p cut, has
     *        the given mean and standard deviation.
     *
     * Looks for the parameters \f$ \mu \f$ and \f$ \sigma \f$ of a normal
     * distribution such that the distribution obtained by left-truncating it at
     * @p cut has mean @p mean and standard deviation @p stddev.
     *
     * Letting \f$ \alpha = (cut - \mu)/\sigma \f$, \f$ \delta(\alpha) \f$ be
     * the mean excess (see left_truncated_std_normal_excess()) and
     * \f$ v(\alpha) \f$ be the variance (see left_truncated_std_normal_variance())
     * of the standard normal distribution left-truncated at \f$ \alpha \f$, the
     * function numerically solves for \f$ \alpha \f$ the equation
     * \f[ \frac{\sqrt{v(\alpha)}}{\delta(\alpha)} = \frac{stddev}{mean - cut} \f]
     * and then derives \f$ \sigma = (mean - cut)/\delta(\alpha) \f$ and
     * \f$ \mu = cut - \alpha\sigma \f$. When the solution is smaller than
     * \f$ -10 \f$, the truncation is numerically irrelevant and
     * \f$ N(mean, stddev^2) \f$ is returned.
     *
     * @param mean is desired mean of the truncated distribution.
     * @param stddev is desired standard deviation of the truncated distribution.
     * @param cut is left truncation value.
     * @return the (non-truncated) normal distribution with parameters
     *         \f$ \mu \f$ and \f$ \sigma \f$.
     * @throw std::invalid_argument if @p stddev is negative, if @p mean is
     *        smaller than @p cut, or equal to it while @p stddev is positive,
     *        or if @p stddev is not smaller than @p mean - @p cut, or so
     *        close to it that the distribution cannot be fitted.
     * @throw std::runtime_error if the root search fails.
     */
    static std::normal_distribution<T>
    fit_untruncated(const T mean, const T stddev, const T cut)
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
        // alpha = (cut - mu) / sigma:
        //
        //   E[X] - cut = sigma * excess(alpha)
        //   Var[X]     = sigma^2 * variance(alpha)
        //
        // Hence, the ratio stddev / (mean - cut) depends on alpha alone:
        //
        //   ratio(alpha) = sqrt(variance(alpha)) / excess(alpha)
        //
        // ratio(alpha) increases from 0 (alpha -> -inf, almost no truncation) to
        // 1 (alpha -> +inf, the tail looks exponential), and
        // 1 - ratio(alpha) ~ 1/alpha^2 for large alpha. So, alpha can be found
        // as the only root of ratio(alpha) - target_ratio, provided that
        // 0 < target_ratio < 1.
        const T target_ratio = stddev / (mean - cut);
    
        // ratio(alpha) < 1 for every alpha, so no root exists
        if (target_ratio >= 1) {
            throw CLONES::Error<std::invalid_argument>("The aimed standard deviation must be "
                                                       "smaller than the difference between "
                                                       "the aimed mean and the cut.");
        }
    
        auto eq = [target_ratio](const T& alpha) {
            return std::sqrt(left_truncated_std_normal_variance(alpha))
                    / left_truncated_std_normal_excess(alpha) - target_ratio;
        };
    
        constexpr T alpha_min = -10;
    
        // eq is increasing: if it is non-negative at alpha_min, the root lies below
        // alpha_min, i.e., the cut is more than 10 sigma below mu and truncation is
        // numerically irrelevant
        if (eq(alpha_min)>=0) {
            return std::normal_distribution<T>(mean, stddev);
        }
    
        // look for an upper bound of the root: since 1 - ratio(alpha) ~ 1/alpha^2,
        // the root is about 1/sqrt(1 - target_ratio)
        constexpr T alpha_limit = 1e7;
        T alpha_max = 8;
        while (eq(alpha_max) < 0) {
            alpha_max *= 2;
    
            // ratio(alpha) cannot be distinguished from 1 any more
            if (alpha_max > alpha_limit) {
                throw CLONES::Error<std::invalid_argument>("The aimed standard deviation is "
                                                           "too close to the difference "
                                                           "between the aimed mean and the cut.");
            }
        }
    
        const auto alpha = brent_root(eq, alpha_min, alpha_max);
    
        // Compute the final parameters: sigma from E[X] - cut = sigma * excess(alpha)
        // and mu from the definition of alpha
        const auto sigma = (mean - cut) / left_truncated_std_normal_excess(alpha);
        const auto mu = cut - alpha * sigma;
    
        return std::normal_distribution(mu, sigma);
    }
};

} // CLONES

#endif // __CLONES_DISTRIBUTIONS__
