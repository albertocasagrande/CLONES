/**
 * @file simulation_wrapper.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements the Python wrapper class and functions for `Simulation`
 * @version 1.0
 * @date 2026-02-06
 *
 * @copyright Copyright (c) 2023-2026-2024
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

#include <shared_mutex>

#include "progress_bar.hpp"
#include "tissue_plotter.hpp"
#include "simulation_wrapper.hpp"
#include "ending_conditions.hpp"

#ifdef WITH_SDL2
#include "SDL_plot.hpp"
#endif

SimulationWrapper::_SimulationWrapper::_SimulationWrapper(int random_seed):
    simulation(random_seed)
{}

SimulationWrapper::SimulationWrapper(int random_seed):
    obj_ptr(std::make_shared<SimulationWrapper::_SimulationWrapper>(random_seed))
{}

void SimulationWrapper::schedule_mutation(const CLONES::Mutants::MutantProperties& src,
                                          const CLONES::Mutants::MutantProperties& dst,
                                          const CLONES::Time time)
{
    std::unique_lock lock(obj_ptr->s_mutex);

    obj_ptr->simulation.schedule_mutation(src, dst, time);
}

struct PythonEndTest : public CLONES::Mutants::Evolutions::TimeTest
{
    /**
     * @brief The empty constructor
     */
    explicit PythonEndTest(const CLONES::Time& time):
        CLONES::Mutants::Evolutions::TimeTest(time)
    {}

    /**
     * @brief Establish whether the simulation must end
     *
     * @param simulation is the considered simulation
     * @return `true` if and only if a signal has been sent to the Python process
     */
    inline bool operator()(const CLONES::Mutants::Evolutions::Simulation& simulation)
    {
        return CLONES::Mutants::Evolutions::TimeTest::operator()(simulation) || PyErr_CheckSignals() == -1;
    }
};

void SimulationWrapper::run_up_to(const CLONES::Time& final_time, const bool quiet,
                                  const bool plot)
{
    using namespace CLONES::UI;

    std::unique_lock lock(obj_ptr->s_mutex);

#ifdef WITH_SDL2

    TissuePlotter<SDLWindow>* plotter{nullptr};

    if (plot) {
        plotter = new TissuePlotter<SDLWindow>(obj_ptr->simulation.tissue());
        plotter->set_frames_per_second(1);
    }
#else
    (void)plot;

    TissuePlotter<Plot2DWindow>* plotter{nullptr};
#endif

    ProgressBar* bar{nullptr};

    PythonEndTest done(final_time);
    if (quiet) {
        obj_ptr->simulation.run(done, plotter, bar);
    } else {
        bar = new ProgressBar(std::cout);
        obj_ptr->simulation.run(done, plotter, bar);

        delete bar;
    }

    if (plotter != nullptr) {
#ifdef WITH_SDL2
        while (plotter->waiting_end()) {
            plotter->plot(obj_ptr->simulation.get_statistics());
        }
#endif

        delete plotter;
    }
}

const CLONES::Time& SimulationWrapper::get_time() const
{
    std::shared_lock lock(obj_ptr->s_mutex);

    return obj_ptr->simulation.get_time();
}

void SimulationWrapper::add_mutant(const CLONES::Mutants::MutantProperties& mutant)
{
    std::unique_lock lock(obj_ptr->s_mutex);

    obj_ptr->simulation.add_mutant(mutant);
}

CLONES::Mutants::Evolutions::PositionInTissue
from_Python_list_to_position(boost::python::list const& position, const uint8_t num_of_dimensions)
{
    namespace bp = boost::python;
    using namespace CLONES::Mutants::Evolutions;

    std::vector<AxisPosition> c_position;
    try {
        c_position = std::vector<AxisPosition>{
            bp::stl_input_iterator<AxisPosition>(position),
            bp::stl_input_iterator<AxisPosition>()
        };
    } catch (...) {
        throw std::domain_error("Expected a list of integer numbers");
    }

    if (c_position.size() != num_of_dimensions) {
        auto dim_text = std::to_string(num_of_dimensions);
        throw std::domain_error("Expected a " + dim_text + "D position, "
                                + "i.e., a list having " + dim_text+ " elements");
    }

    PositionInTissue pos;

    if (num_of_dimensions==2) {
        return PositionInTissue(c_position[0], c_position[1]);
    }

    return PositionInTissue(c_position[0], c_position[1], c_position[2]);
}


void SimulationWrapper::place_cell(const CLONES::Mutants::MutantProperties& mutant,
                                   const std::string& methylation_signature,
                                   boost::python::list const& position)
{
    const auto num_of_dimensions = obj_ptr->simulation.tissue().num_of_dimensions();
    auto c_position = from_Python_list_to_position(position, num_of_dimensions);

    {
        std::unique_lock lock(obj_ptr->s_mutex);
        obj_ptr->simulation.place_cell(mutant[methylation_signature], c_position);
    }

}

void SimulationWrapper::set_tissue(const std::string& name, boost::python::list const& sizes_list)
{
    namespace bp = boost::python;
    using namespace CLONES::Mutants::Evolutions;

    std::vector<AxisSize> c_sizes;
    try {
        c_sizes = std::vector<AxisSize>{
            bp::stl_input_iterator<AxisSize>(sizes_list),
            bp::stl_input_iterator<AxisSize>()
        };
    } catch (...) {
        throw std::domain_error("Expected a list of sizes");
    }

    {
        std::unique_lock lock(obj_ptr->s_mutex);
        obj_ptr->simulation.set_tissue(name, c_sizes);
    }
}
