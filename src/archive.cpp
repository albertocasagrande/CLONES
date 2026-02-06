/**
 * @file archive.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements some archive classes and their methods
 * @version 1.3
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

#include "archive.hpp"

namespace CLONES
{

namespace Archive {


WrongFileFormatDescr::WrongFileFormatDescr(const std::string& expected_descriptor,
                                           const std::string& read_descriptor,
                                           const std::filesystem::path& filepath):
    expected_descr(expected_descriptor), read_descr(read_descriptor),
    filepath(std::filesystem::absolute(filepath))
{}

WrongFileFormatVersion::WrongFileFormatVersion(const uint8_t& expected_version,
                                               const uint8_t& read_version,
                                               const std::filesystem::path& filepath):
    expected_version(expected_version), read_version(read_version),
    filepath(std::filesystem::absolute(filepath))
{}


namespace Basic {

Basic::Basic():
    fs()
{}

Basic::Basic(std::filesystem::path position, std::ios_base::openmode mode):
    fs(position, mode), filepath(position)
{}

std::streampos Basic::size()
{
    // store current position
    const auto pos = fs.tellg();

    // seek the end
    fs.seekg(0, std::ios::end);

    // read the number of bytes from the beginning
    const auto file_bytes = fs.tellg();

    // move to the original position
    fs.seekg(pos);

    return file_bytes;
}

Basic::~Basic()
{
    if (is_open()) {
        close();
    }
}

Out::Out():
    Basic()
{}

Out::Out(std::filesystem::path position):
    Basic(position, std::fstream::out)
{}

Out::Out(std::filesystem::path position, std::ios_base::openmode mode):
    Basic(position, mode | std::fstream::out)
{}

In::In(std::filesystem::path position):
    Basic(position, std::fstream::in)
{}

In::In(std::filesystem::path position, std::ios_base::openmode mode):
    Basic(position, mode | std::fstream::in)
{}

ProgressViewer::ProgressViewer():
    progress_bar(nullptr), total_steps{0}, next_percentage{0}, performed_steps{0}
{}

void ProgressViewer::initialize(CLONES::UI::ProgressBar* progress_bar, const size_t total_steps)
{
    this->progress_bar = progress_bar;
    this->total_steps = total_steps;
    performed_steps = 0;
    next_percentage = total_steps/100;
}

void ProgressViewer::advance(const size_t& steps)
{
    if (progress_bar != nullptr) {
        performed_steps += steps;

        if (performed_steps > next_percentage) {
            progress_bar->set_progress(progress_bar->get_progress()+1);
            next_percentage = (progress_bar->get_progress()+1)*total_steps/100;
        }
    }
}

void ProgressViewer::reset()
{
    progress_bar = nullptr;
}

}  // Basic

namespace Binary {

Out::Out():
    Archive::Basic::Out(), ProgressViewer()
{}

Out::Out(std::filesystem::path position):
    Archive::Basic::Out(position, std::fstream::binary), ProgressViewer()
{}

Out::Out(std::filesystem::path position, std::ios_base::openmode mode):
    Archive::Basic::Out(position, mode | std::fstream::binary), ProgressViewer()
{}

ByteCounter::ByteCounter(CLONES::UI::ProgressBar* progress_bar):
    Out(), bytes(0), progress_bar(progress_bar)
{}

In::In(std::filesystem::path position):
    Archive::Basic::In(position, std::fstream::binary), ProgressViewer()
{}

In::In(std::filesystem::path position, std::ios_base::openmode mode):
    Archive::Basic::In(position, std::fstream::binary | mode), ProgressViewer()
{}

}

}

}
