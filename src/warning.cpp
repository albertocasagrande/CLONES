/**
 * @file warning.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement the default warning function
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

#include <iostream>

#include "warning.hpp"

namespace CLONES
{

void warning(const WarningType type, const std::string message)
{
    std::cerr << "Warning: " << message;
    switch(type) {
        case WarningType::NO_MUT_FOR_CONTEXT:
            std::cerr << " Please try decreasing `sampling_rate` "
                      << "when building ContextIndex." << std::endl;
            break;
        case WarningType::NO_MUT_FOR_RPATTERN:
            std::cerr << " Please try increasing `max_stored_repetitions` "
                      << "when building RSIndex." << std::endl;
            break;
    }
}

}
