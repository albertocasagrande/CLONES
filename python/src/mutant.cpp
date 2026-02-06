/**
 * @file mutant.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements the Python wrapper class and functions for `MutantProperties`
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

#include "cell_event.hpp"

#include "mutant.hpp"
#include "epigenetic_rates.hpp"

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

std::shared_ptr<CLONES::Mutants::MutantProperties>
CloneWrapper::create(std::string const& name, boost::python::list const& epigenetic_rates)
{
    using namespace CLONES::Mutants;

    namespace bp = boost::python;

    std::vector<EpigeneticRates> c_epigenetic_rates;

    for (auto i=0; i<bp::len(epigenetic_rates); ++i) {
        bp::extract<bp::list> o_extractor(epigenetic_rates[i]);

        if (o_extractor.check()) {
            c_epigenetic_rates.push_back(*(EpigeneticRatesWrapper::create(o_extractor)));
        } else {
            bp::extract<EpigeneticRates> o_extractor(epigenetic_rates[i]);

            if (o_extractor.check()) {
                c_epigenetic_rates.push_back(o_extractor);
            } else {
                throw std::domain_error("The list must contain EpigeneticRates");
            }
        }
    }

    return std::make_shared<MutantProperties>(name, c_epigenetic_rates);
}

void CloneWrapper::set_rates(CLONES::Mutants::MutantProperties *mutant,
                                const std::string& methylation_signature,
                                boost::python::dict const& rates)
{
    using namespace CLONES::Mutants;

    namespace bp = boost::python;

    std::map<CellEventType, double> c_rates;

    auto rates_list = bp::list(rates.items());
    for (auto i=0; i<len(rates_list); ++i) {
        const unsigned int event_type_num = bp::extract<unsigned int>(rates_list[i][0]);
        const CellEventType event_type_id = static_cast<CellEventType>(event_type_num);
        c_rates[event_type_id] = bp::extract<double>(rates_list[i][1]);
    }

    (*mutant)[methylation_signature].set_rates(c_rates);
}
