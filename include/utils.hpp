/**
 * @file utils.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines utility functions
 * @version 1.4
 * @date 2025-10-11
 *
 * @copyright Copyright (c) 2023-2025
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

#ifndef __RACES_UTILS__
#define __RACES_UTILS__

#include <string>
#include <filesystem>
#include <map>
#include <set>
#include <algorithm>
#include <ostream>
#include <concepts>
#include <type_traits>

#include<regex>
#include<clocale>

namespace RACES
{

/**
 * @brief A structure defined to deal exclusively with non-empty packs
 */
template<typename...> struct last_arg_type; // primary has no ::type

/**
 * @brief A specialization of `last_arg_type` for packs with one type
 *
 * @tparam T is the only type of the pack
 */
template<typename T>
struct last_arg_type<T> { using type = T; };

/**
 * @brief A specialization of `last_arg_type` for packs with at least one type
 *
 * @tparam T is the first type of the pack
 * @tparam REST are the remaining pack
 */
template<typename T, typename... REST>
struct last_arg_type<T, REST...> : last_arg_type<REST...> { };

/**
 * @brief Concept that is true only when the pack is non-empty and the last arg is not T
 *
 * @tparam T is the non-aimed type of the last parameter
 * @tparam ARGS is the type of the parameters
 */
template<typename T, typename... ARGS>
concept last_is_not = requires { typename last_arg_type<ARGS...>::type; } &&
                        !std::same_as<std::remove_cvref_t<typename last_arg_type<ARGS...>::type>, T>;

} // RACES

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
 * @param[in] prefix is the prefix of the path basename
 * @param[in] parent_dir is the parent directory of the resulting path
 * @return a path that does not exists, whose parent directory
 *      is `parent_dir`, and whose basename prefix is `prefix`
 */
std::filesystem::path
get_a_temporary_path(const std::string& prefix="RACES_tmp_",
                     const std::filesystem::path parent_dir=std::filesystem::temp_directory_path());

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

#endif // __RACES_UTILS__
