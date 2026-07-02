/**
 * @file error.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines CLONES errors
 * @version 1.3
 * @date 2026-07-02
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

#ifndef __CLONES_ERROR__
#define __CLONES_ERROR__

#include <cstdint>

#include <string>
#include <sstream>
#include <concepts>
#include <exception>
#include <source_location>

namespace CLONES
{

/**
 * @brief An exception wrapper storing the location of object creation
 *
 * @tparam EXCEPTION_TYPE is the base type for the exception
 */
template<typename EXCEPTION_TYPE>
  requires std::derived_from<EXCEPTION_TYPE, std::exception>
class Error : public EXCEPTION_TYPE
{
    std::source_location location;  //!< The object creation location
public:

    /**
     * @brief A constructor wrapping an already existing exception
     *
     * @param ex is the exception to be wrapped
     * @param location is the location associated to the error
     */
    explicit Error(const EXCEPTION_TYPE& ex,
                   std::source_location location = std::source_location::current()):
        EXCEPTION_TYPE(ex.what()), location(location)
    {}

    /**
     * @brief A constructor
     *
     * @param msg is the error message
     * @param location is the location associated to the error
     */
    explicit Error(const std::string& msg,
                   std::source_location location = std::source_location::current()):
        EXCEPTION_TYPE(msg), location(location)
    {}

    /**
     * @brief A constructor
     *
     * @param msg is the error message
     * @param location is the location associated to the error
     */
    explicit Error(const char* msg,
                   std::source_location location = std::source_location::current()):
        EXCEPTION_TYPE(msg), location(location)
    {}

    /**
     * @brief Get the error message
     *
     * @return The error message
     */
    inline const char* what() const noexcept override
    {
        return EXCEPTION_TYPE::what();
    }

    /**
     * @brief Get the debug message
     *
     * @return The debug message
     */
    inline std::string debug_msg() const noexcept
    {
        std::ostringstream oss;

        oss << function_name() << " (" << file_name() << ":"
            << line() << "): " << what();

        return oss.str();
    }

    /**
     * @brief Get the name of the function in which the error was created
     *
     * @return The name of the function in which the error was created
     */
    inline const char* function_name() const noexcept
    {
        return location.function_name();
    }

    /**
     * @brief Get the name of the file in which the error was created
     *
     * @return The name of the file in which the error was created
     */
    inline const char* file_name() const noexcept
    {
        return location.file_name();
    }

    /**
     * @brief Get the line number in which the error was created
     *
     * @return The line number in which the error was created
     */
    inline constexpr std::uint_least32_t line() const noexcept
    {
        return location.line();
    }
};

} // CLONES

#endif // __CLONES_ERROR__