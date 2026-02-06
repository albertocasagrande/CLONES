/**
 * @file warning.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines warning functions
 * @version 1.1
 * @date 2026-02-06
 *
 * @copyright Copyright (c) 2025
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

#ifndef __CLONES_WARNING__
#define __CLONES_WARNING__

#include <string>
#include <functional>

namespace CLONES
{

    /**
     * @brief The warning types
     */
    enum class WarningType
    {
        NO_MUT_FOR_CONTEXT,    //!< No mutations are available for a context
        NO_MUT_FOR_RPATTERN    //!< No mutations are available for a repetition pattern
    };

    /**
     * @brief The warning function type
     */
    using WarningFunction = std::function<void(const WarningType, const std::string)>;

    /**
     * @brief The default warning function
     * 
     * @param type is the warning type
     * @param message is the warning message
     */
    void warning(const WarningType type, const std::string message);
}

#endif // __CLONES_WARNING__