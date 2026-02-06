/**
 * @file utils.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement utility functions
 * @version 1.1
 * @date 2026-02-06
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


#include <string>
#include <filesystem>
#include <sstream>

#include "utils.hpp"

std::filesystem::path
get_a_temporary_path(const std::string& prefix,
                     const std::filesystem::path parent_dir)
{
    if (!std::filesystem::exists(parent_dir)) {
        std::ostringstream oss;

        oss << "get_a_temporary_path: the directory \""
            << to_string(parent_dir)
            << "\" does not exist.";

        throw std::domain_error(oss.str());
    }

    if (!std::filesystem::is_directory(parent_dir)) {
        std::ostringstream oss;

        oss << "get_a_temporary_path: the directory \""
            << to_string(parent_dir)
            << "\" is not a directory.";

        throw std::domain_error(oss.str());
    }

    size_t counter{0};
    std::filesystem::path tmp_path;
    do {
        tmp_path = parent_dir / (prefix + std::to_string(counter));
        ++counter;
    } while (std::filesystem::exists(tmp_path));

    return tmp_path;
}