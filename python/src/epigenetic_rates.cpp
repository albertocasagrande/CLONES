/**
 * @file epigenetic_rates.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement the Python wrapper class and functions for `CLONES::EpigeneticRates`
 * @version 1.0
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

#include "epigenetic_rates.hpp"

std::shared_ptr<CLONES::Mutants::EpigeneticRates>
EpigeneticRatesWrapper::create(boost::python::list const& rates_list)
{
    using namespace CLONES::Mutants;

    namespace bp = boost::python;

    const auto size = bp::len(rates_list);
    if (size==1) {
        return std::make_shared<EpigeneticRates>(bp::extract<double>(rates_list[0]));
    }
    if (size==2) {
        return std::make_shared<EpigeneticRates>(bp::extract<double>(rates_list[0]),
                                                 bp::extract<double>(rates_list[1]));
    }

    throw std::domain_error("EpigeneticRates cannot be built by using " +
                            std::to_string(size) + " doubles");
}
