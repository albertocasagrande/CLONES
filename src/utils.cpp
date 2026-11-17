/**
 * @file utils.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement utility functions
 * @version 1.3
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


#include <string>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <cstdint>
#include <system_error>

#include "utils.hpp"

#include "error.hpp"


namespace
{

/**
 * @brief Get a random 64-bit hexadecimal string
 *
 * Each thread has its own generator, seeded by both `std::random_device`
 * and the current time, so that different processes and threads produce
 * different sequences.
 *
 * @return a random 16-character hexadecimal string
 */
std::string get_random_suffix()
{
    thread_local std::mt19937_64 generator = [] {
        std::random_device device;

        const auto now = static_cast<uint64_t>(std::chrono::high_resolution_clock::now()
                                                    .time_since_epoch().count());

        std::seed_seq seed{device(), device(), device(), device(),
                           static_cast<uint32_t>(now), static_cast<uint32_t>(now >> 32)};

        return std::mt19937_64(seed);
    }();

    std::ostringstream oss;
    oss << std::hex << std::setw(16) << std::setfill('0') << generator();

    return oss.str();
}

void validate_parent_directory(const std::filesystem::path& parent_dir)
{
    if (!std::filesystem::exists(parent_dir)) {
        throw CLONES::Error<std::domain_error>("The directory \""
                                               + to_string(parent_dir)
                                               + "\" does not exist.");
    }

    if (!std::filesystem::is_directory(parent_dir)) {
        throw CLONES::Error<std::domain_error>("The directory \""
                                               + to_string(parent_dir)
                                               + "\" is not a directory.");
    }
}

}   // namespace

std::filesystem::path
get_a_temporary_path(const std::string& prefix,
                     const std::filesystem::path parent_dir)
{
    validate_parent_directory(parent_dir);

    std::filesystem::path tmp_path;
    do {
        tmp_path = parent_dir / (prefix + get_random_suffix());
    } while (std::filesystem::exists(tmp_path));

    return tmp_path;
}

TemporaryDirectory::TemporaryDirectory(const std::string& prefix,
                                       const std::filesystem::path& parent_dir)
{
    validate_parent_directory(parent_dir);

    // `create_directory()` returns `true` only for the caller that
    // actually created the directory: this makes the reservation atomic
    constexpr size_t max_attempts = 100;
    for (size_t i=0; i<max_attempts; ++i) {
        const auto candidate = parent_dir / (prefix + get_random_suffix());

        std::error_code error;
        if (std::filesystem::create_directory(candidate, error)) {
            dir_path = candidate;

            return;
        }

        // the candidate already exists: try another name
        if (error && error != std::errc::file_exists) {
            throw CLONES::Error<std::runtime_error>("Cannot create the directory \""
                                                    + to_string(candidate) + "\": "
                                                    + error.message());
        }
    }

    throw CLONES::Error<std::runtime_error>("Cannot create a temporary directory in \""
                                            + to_string(parent_dir) + "\".");
}

TemporaryDirectory::~TemporaryDirectory()
{
    std::error_code error;

    std::filesystem::remove_all(dir_path, error);
}