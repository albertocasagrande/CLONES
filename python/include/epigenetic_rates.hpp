/**
 * @file epigenetic_rates.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define the Python wrapper class and functions for `CLONES::EpigeneticRates`
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

#ifndef __CLONES_PYTHON_EPIGENETIC_RATES__
#define __CLONES_PYTHON_EPIGENETIC_RATES__

#include <memory>

#include <boost/python.hpp>

#include "mutant_properties.hpp"


using namespace boost::python;

/**
 * @brief A wrapper structure for `CLONES::EpigeneticRates`
 */
struct EpigeneticRatesWrapper
{
    /**
     * @brief Create a new `CLONES::EpigeneticRates` object
     *
     * This method creates a new `CLONES::EpigeneticRates` object
     * by using as parameters the values stored in Python list.
     * If the list contains 1 or 2 real values, then this
     * method creates the object. Otherwise, it throws a
     * domain exception.
     *
     * @param rates_list is a Python list of object
     * @return a shared pointer to the newly created
     *        `CLONES::EpigeneticRates` object
     */
    static std::shared_ptr<CLONES::Mutants::EpigeneticRates>
    create(boost::python::list const& rates_list);

    /**
     * @brief Set the methylation rate of a `CLONES::EpigeneticRates` object
     *
     * @param epigenetic_rates is the `CLONES::EpigeneticRates` object
     * @param value is the methylation rate value to be set
     */
    inline static
    void set_methylation_rate(CLONES::Mutants::EpigeneticRates *epigenetic_rates, const double& value)
    {
        epigenetic_rates->set_methylation_rate(value);
    }

    /**
     * @brief Set the demethylation rate of a `CLONES::EpigeneticRates` object
     *
     * @param epigenetic_rates is the `CLONES::EpigeneticRates` object
     * @param value is the demethylation rate value to be set
     */
    inline static
    void set_demethylation_rate(CLONES::Mutants::EpigeneticRates *epigenetic_rates, const double& value)
    {
        epigenetic_rates->set_demethylation_rate(value);
    }
};

#endif // __CLONES_PYTHON_EPIGENETIC_RATES__
