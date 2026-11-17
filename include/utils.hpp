/**
 * @file utils.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines utility functions
 * @version 1.8
 * @date 2026-11-17
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

#ifndef __CLONES_UTILS__
#define __CLONES_UTILS__

#include <string>
#include <filesystem>
#include <map>
#include <set>
#include <algorithm>
#include <ostream>
#include <limits>

#include <cmath>
#include <stdexcept>
#include <concepts>
#include <cstddef>
#include <random>
#include <type_traits>

#include<regex>
#include<clocale>

#include "error.hpp"

/**
 * @brief Build a std::regex
 *
 * This method implements a workaround for creating `std::regex`
 * objects in unsupported Window locales.
 *
 * @param regex_str is the regular expression used to build
 *          the std::regex object.
 */
inline std::regex build_regex(const char* regex_str)
{
#if defined(__WIN32__) && defined(__GNUG__) && __GNUC__ <= 14
    // saving current locale
    const std::string orig_locale = std::setlocale(LC_ALL, nullptr);

    // setting a temporary locale
    std::setlocale(LC_ALL, "C");

    std::regex re(regex_str);

    // resetting original locale
    std::setlocale(LC_ALL, orig_locale.c_str());
#else
    std::regex re(regex_str);
#endif

    return re;
}

/**
 * @brief Build a std::regex
 *
 * This method implements a workaround for creating `std::regex`
 * objects in unsupported Window locales.
 *
 * @param regex_str is the regular expression used to build
 *          the std::regex object.
 */
inline std::regex build_regex(const std::string& regex_str)
{
    return build_regex(regex_str.c_str());
}

/**
 * @brief Convert a path to a string
 *
 * @param fs_path is the path to be converted
 * @return a string representation of the path
 */
inline std::string to_string(const std::filesystem::path& fs_path)
{
#if defined(__WIN32__) || defined(__WIN64__)
    std::wstring w_path(fs_path);

    return std::string(w_path.begin(), w_path.end());
#else
    return fs_path;
#endif
}

/**
 * @brief Compute the ceiling of a division
 *
 * @tparam TYPE is the type of the operands
 * @param x is the dividend
 * @param y is the divisor
 * @return the ceiling of the division `x/y`
 */
template<typename TYPE>
inline static TYPE ceil_div(const TYPE& x, const TYPE& y)
{
    return 1 + ((x - 1) / y);
}

/**
 * @brief Get the a temporary path
 *
 * The basename of the returned path is @p prefix followed by a random
 * 64-bit hexadecimal suffix, so that concurrent processes and threads
 * get different paths with overwhelming probability. However, the
 * returned path is not reserved: use `TemporaryDirectory` when an
 * exclusive location is required.
 *
 * @param[in] prefix is the prefix of the path basename
 * @param[in] parent_dir is the parent directory of the resulting path
 * @return a path that does not exists, whose parent directory
 *      is @p parent_dir, and whose basename prefix is @p prefix
 */
std::filesystem::path
get_a_temporary_path(const std::string& prefix="CLONES_tmp_",
                     const std::filesystem::path parent_dir=std::filesystem::temp_directory_path());

/**
 * @brief A temporary directory removed on destruction
 *
 * The directory is atomically created by the constructor, so it is
 * exclusively owned by the object even when several processes or
 * threads create temporary directories in the same parent directory
 * at the same time. The destructor removes the directory together
 * with its content.
 */
class TemporaryDirectory
{
    std::filesystem::path dir_path;     //!< The path of the directory
public:
    /**
     * @brief Create a temporary directory
     *
     * @param[in] prefix is the prefix of the directory basename
     * @param[in] parent_dir is the parent directory of the temporary directory
     * @throw CLONES::Error<std::domain_error> if @p parent_dir is not an
     *      existing directory
     * @throw CLONES::Error<std::runtime_error> if the directory cannot be
     *      created
     */
    explicit TemporaryDirectory(const std::string& prefix="CLONES_tmp_dir_",
                                const std::filesystem::path& parent_dir=std::filesystem::temp_directory_path());

    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

    /**
     * @brief Get the path of the temporary directory
     *
     * @return the path of the temporary directory
     */
    inline const std::filesystem::path& path() const
    {
        return dir_path;
    }

    /**
     * @brief Destroy the temporary directory and its content
     */
    ~TemporaryDirectory();
};

/**
 * @brief Finds a root of a continuous function using Brent's method.
 *
 * Searches for a zero of @p func in the interval [@p low, @p high]. The
 * function values at the endpoints must have opposite signs (or be zero).
 * At each step the method tries inverse quadratic interpolation or the
 * secant method, and falls back to bisection whenever the interpolated
 * point is not acceptable. Hence, it keeps the robustness of the
 * bisection algorithm while usually converging superlinearly.
 *
 * The iteration stops when \f$ |f(b)| < tol \f$, when the bracketing
 * interval is narrower than about @p tol, or after @p max_iter iterations.
 *
 * @tparam FUNC type of the callable object.
 * @tparam FUNC_DOMAIN type of the function domain.
 * @tparam FUNC_CODOMAIN type of the function codomain (default `double`).
 * @param func function whose root is sought.
 * @param low lower endpoint of the interval.
 * @param high upper endpoint of the interval.
 * @param max_iter maximum number of iterations (default 100).
 * @param tol tolerance on the function value and on the interval
 *            width (default 1e-12).
 * @return an approximation of the root in [@p low, @p high].
 * @throw std::runtime_error if \f$ f(low) \cdot f(high) > 0 \f$, i.e. if a
 *        root is not guaranteed to lie in the interval.
 */
template <typename FUNC, typename FUNC_DOMAIN, typename FUNC_CODOMAIN=double>
requires std::invocable<FUNC, FUNC_DOMAIN> &&
         std::same_as<std::invoke_result_t<FUNC, FUNC_DOMAIN>, FUNC_CODOMAIN>
FUNC_DOMAIN brent_root(FUNC func, FUNC_DOMAIN low, FUNC_DOMAIN high,
                       const size_t max_iter = 100,
                       const FUNC_CODOMAIN tol = 1e-12)
{
    // b is the current estimate, a the previous one, and [b, c] (or [c, b])
    // always brackets the root
    auto a = low, b = high;
    auto f_a = func(a), f_b = func(b);

    if (f_a * f_b > 0) {
        throw std::runtime_error("The root is not in the interval ["
                                 + std::to_string(low) + ","
                                 + std::to_string(high) + "]");
    }

    auto c = a;
    auto f_c = f_a;
    auto step = b - a, prev_step = step;

    for (size_t i = 0; i < max_iter; ++i) {
        if (f_b * f_c > 0) {
            c = a;
            f_c = f_a;
            step = prev_step = b - a;
        }

        // ensure that b is the best estimate so far
        if (std::abs(f_c) < std::abs(f_b)) {
            a = b;
            b = c;
            c = a;
            f_a = f_b;
            f_b = f_c;
            f_c = f_a;
        }

        const FUNC_DOMAIN tol_b = 2 * std::numeric_limits<FUNC_DOMAIN>::epsilon()
                                  * std::abs(b) + tol / 2;
        const auto half_width = (c - b) / 2;

        if (std::abs(half_width) <= tol_b || std::abs(f_b) < tol) {
            return b;
        }

        if (std::abs(prev_step) >= tol_b && std::abs(f_a) > std::abs(f_b)) {
            // try interpolation: the new step is p/q
            const auto s = f_b / f_a;
            FUNC_DOMAIN p, q;
            if (a == c) {
                // secant method
                p = 2 * half_width * s;
                q = 1 - s;
            } else {
                // inverse quadratic interpolation
                const auto r_a = f_a / f_c;
                const auto r_b = f_b / f_c;
                p = s * (2 * half_width * r_a * (r_a - r_b) - (b - a) * (r_b - 1));
                q = (r_a - 1) * (r_b - 1) * (s - 1);
            }

            if (p > 0) {
                q = -q;
            } else {
                p = -p;
            }

            // accept the interpolation only if it falls inside the bracket
            // and the steps are shrinking fast enough
            if (2 * p < std::min(3 * half_width * q - std::abs(tol_b * q),
                                 std::abs(prev_step * q))) {
                prev_step = step;
                step = p / q;
            } else {
                step = prev_step = half_width;
            }
        } else {
            // fall back to bisection
            step = prev_step = half_width;
        }

        a = b;
        f_a = f_b;
        if (std::abs(step) > tol_b) {
            b += step;
        } else {
            b += (half_width > 0 ? tol_b : -tol_b);
        }
        f_b = func(b);
    }
    return b;
}

namespace std
{

/**
 * @brief Test whether a string is a suffix of another string
 *
 * @param suffix is the suffix to be tested
 * @param text is the test whose suffix is tested
 * @return `true` if and only if `suffix` is a suffix of `text`
 */
inline bool is_suffix_of(const std::string& suffix, const std::string& text) {
    return ((text.size()>=suffix.size())
            && (text.compare(text.size()-suffix.size(), suffix.size(), suffix)));
}

/**
 * @brief Get the union of two sets
 *
 * @tparam T is the type of the set values
 * @param A is a set of T values
 * @param B is a set of T values
 * @return the set \f$A \cup B\f$
 */
template<typename T>
std::set<T> get_union(const std::set<T>& A, const std::set<T>& B)
{
    std::set<T> C;

    std::set_union(A.begin(), A.end(), B.begin(), B.end(),
                   std::inserter(C, C.end()));

    return C;
}

/**
 * @brief Test whether a value is the maximum of its type
 *
 * @tparam VALUE_TYPE is the type of the value to be tested
 * @param value is the value to be tested
 * @return `true` if and only if `value` is the maximum
 *   among all values representable by `value` type
 */
template<typename VALUE_TYPE>
inline bool is_max(const VALUE_TYPE& value)
{
    return std::numeric_limits<VALUE_TYPE>::max() == value;
}

/**
 * @brief Stream a map
 *
 * @tparam KEYS is the type of the map keys
 * @tparam VALUES is the type of the map values
 * @param os is the output stream
 * @param m is the map to be streamed
 * @return a reference to the updated stream
 */
template<typename KEYS, typename VALUES>
std::ostream& operator<<(std::ostream& os, const std::map<KEYS, VALUES>& m)
{
  os << "{";
  std::string sep = "";
  for (const auto& [key, value]: m) {
    os << sep << "\"" << key << "\": " << value;
    sep = ", ";
  }
  os << "}";

  return os;
}

}   // std

#endif // __CLONES_UTILS__
