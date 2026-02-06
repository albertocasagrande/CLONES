/**
 * @file sequencer.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements sequencer models
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

#include <cmath>

#include "sequencer.hpp"

namespace CLONES
{

namespace Sequencers
{

char SangerQualityCodec::encode(const double error_prob)
{
    if ((error_prob > 1) || (error_prob < 0)) {
        throw std::domain_error("The error probability must be a value in [0, 1]: got "
                                + std::to_string(error_prob));
    }

    // Simulating quality scores for wrongly sequenced bases
    if (error_prob == 0) {
        return max_quality_code();
    }

    const int qual = - 10 * std::log10(error_prob);

    if (qual > max_quality_value()) {
        return max_quality_code();
    }

    if (qual < min_quality_value()) {
        return min_quality_code();
    }

    return static_cast<char>(qual) + base_quality_code();
}


double SangerQualityCodec::decode(const char quality_code)
{
    const double quality_value = quality_code-base_quality_code();

    if ((quality_code > max_quality_code()) || (quality_code < min_quality_code())) {
        throw std::domain_error("The quality code must be a character in ['"
                                + std::string(1, min_quality_code())
                                + "' (ASCII-"+ std::to_string(min_quality_value())+"), "
                                + std::string(1, max_quality_code())
                                + "' (ASCII-"+ std::to_string(max_quality_value())+")]"
                                + ": got " + std::string(1, quality_code)
                                + "' (ASCII-"+ std::to_string(quality_value) +")");
    }

    return std::pow(static_cast<double>(10), -quality_value/10);
}

}   // Sequencers

}   // CLONES
