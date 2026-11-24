/**
 * @file sampling.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines sampling functions
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

#ifndef __CLONES_SAMPLING__
#define __CLONES_SAMPLING__

#include <cstddef>
#include <vector>
#include <random>
#include <unordered_set>

#include "error.hpp"

namespace CLONES
{

/**
 * @brief Sample k distinct naturals from [0, n-1]
 *
 * This function uniformly samples a k-subset of \f$\{0, \ldots, n-1\}\f$,
 * i.e., every subset of \f$\{0, \ldots, n-1\}\f$ having \f$k\f$ elements
 * is returned with probability \f$1/\binom{n}{k}\f$.
 *
 * It implements Floyd's algorithm [1]: for every \f$j\f$ from \f$n-k\f$ to
 * \f$n-1\f$, a value \f$t\f$ is uniformly drawn from \f$\{0, \ldots, j\}\f$;
 * \f$t\f$ is added to the sample if it has not been chosen yet, otherwise
 * \f$j\f$ is added. Since \f$j\f$ cannot have been chosen in the previous
 * steps, exactly one new value is added at each step, and, by induction on
 * \f$j\f$, the values chosen up to step \f$j\f$ form a uniformly distributed
 * subset of \f$\{0, \ldots, j\}\f$.
 *
 * Differently from selection sampling (e.g., Algorithm S in [2],
 * Section 3.4.2, and `std::sample()`), whose cost is linear in \f$n\f$, the
 * algorithm draws exactly \f$k\f$ random values. Hence, it takes
 * \f$O(k)\f$ expected time and \f$O(k)\f$ memory, whatever \f$n\f$ is.
 *
 * The set of the sampled values is uniformly distributed, but their order
 * in the returned vector is not: for instance, the value \f$n-1\f$ can only
 * occupy the last position. If a random order or a sorted vector is
 * required, the result must be shuffled or sorted.
 *
 * [1] J. Bentley and B. Floyd. "Programming pearls: A sample of brilliance".
 *     Communications of the ACM 30(9), 754–757 (1987).
 *     https://doi.org/10.1145/30401.315746
 *
 * [2] D. E. Knuth. "The Art of Computer Programming, Volume 2:
 *     Seminumerical Algorithms". 3rd edition, Addison-Wesley (1997).
 *
 * @tparam RANDOM_GENERATOR is the random number generator type. It must
 *      satisfy `std::uniform_random_bit_generator`
 * @param[in] n is the number of values from which the sample is drawn,
 *      i.e., the values are drawn from \f$\{0, \ldots, n-1\}\f$
 * @param[in] k is the number of values to be sampled. It must not be
 *      greater than @p n
 * @param[in,out] generator is the random number generator
 * @return a vector containing @p k distinct values in
 *      \f$\{0, \ldots, n-1\}\f$. The set of these values is uniformly
 *      distributed among the subsets of \f$\{0, \ldots, n-1\}\f$ having
 *      @p k elements, while their order is not uniformly distributed
 * @throw std::invalid_argument if @p k is greater than @p n
 */
template<typename RANDOM_GENERATOR>
    requires std::uniform_random_bit_generator<RANDOM_GENERATOR>
std::vector<size_t> sample_distinct(const size_t n, const size_t k, RANDOM_GENERATOR& generator)
{
    if (k>n) {
        throw CLONES::Error<std::invalid_argument>("k must be smaller or equal to n");
    }

    std::unordered_set<size_t> chosen;
    std::vector<size_t> sample;
    sample.reserve(k);

    for (size_t j = n - k; j < n; ++j) {
        const size_t t = std::uniform_int_distribution<size_t>(0, j)(generator);

        if (chosen.insert(t).second) {
            // if t has never been chosen, push it into sample

            sample.push_back(t);
        } else {
            // if t was already chosen, j cannot have been
            chosen.insert(j);

            sample.push_back(j);
        }
    }

    return sample;
}

}

#endif // __CLONES_SAMPLING__