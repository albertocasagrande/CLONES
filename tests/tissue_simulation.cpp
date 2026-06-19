/**
 * @file tissue_simulation.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Some archive tests
 * @version 1.0
 * @date 2026-06-19
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tissue_simulation

#include <sstream>
#include <filesystem>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "tissue_simulation.hpp"
#include "ending_conditions.hpp"

#include "test_common.hpp"

struct ArchiveFixture {
    CLONES::Mutants::Evolutions::TissueSimulation tissue_simulation;

    ArchiveFixture():
        tissue_simulation()
    {}

    ~ArchiveFixture()
    {
        std::filesystem::remove_all(tissue_simulation.get_logger().get_directory());
    }
};

BOOST_AUTO_TEST_CASE(build_simulation)
{
    using namespace CLONES::Mutants::Evolutions;

    BOOST_CHECK_NO_THROW(TissueSimulation tissue_simulation);
}

BOOST_FIXTURE_TEST_SUITE( simulatedData, ArchiveFixture )

BOOST_AUTO_TEST_CASE(set_tissue)
{
    BOOST_CHECK_NO_THROW(tissue_simulation.set_tissue("Liver", {1000,1000}));
}

BOOST_AUTO_TEST_CASE(add_mutant)
{
    using namespace CLONES;

    BOOST_CHECK_NO_THROW(tissue_simulation.add_mutant("A"));
    BOOST_CHECK_NO_THROW(tissue_simulation.add_mutant("B"));

    BOOST_CHECK_THROW(tissue_simulation.add_mutant("A"), Error<std::runtime_error>);
}

BOOST_AUTO_TEST_CASE(add_epigenetic_state)
{
    using namespace CLONES;

    BOOST_CHECK_NO_THROW(tissue_simulation.add_epigenetic_states({"E1", "E2", "E3"}));

    BOOST_CHECK_THROW(tissue_simulation.add_epigenetic_state("E1"),
                      Error<std::runtime_error>);

    BOOST_CHECK_THROW(tissue_simulation.add_epigenetic_states({"E5", "E1"}),
                      Error<std::runtime_error>);
}

BOOST_AUTO_TEST_CASE(add_mutant_epigenetic_state)
{
    BOOST_CHECK_NO_THROW(tissue_simulation.add_mutant("A"));
    BOOST_CHECK_NO_THROW(tissue_simulation.add_mutant("B"));
    BOOST_CHECK_NO_THROW(tissue_simulation.add_epigenetic_states({"E1", "E2", "E3"}));
}

void set_rate(CLONES::Mutants::Evolutions::TissueSimulation& tissue_simulation)
{
    using namespace CLONES::Mutants;
    using namespace CLONES::Mutants::Evolutions;

    BOOST_CHECK_NO_THROW(tissue_simulation.set_tissue("Liver", {1000,1000}));
    BOOST_CHECK_NO_THROW(tissue_simulation.add_mutant("A"));
    BOOST_CHECK_NO_THROW(tissue_simulation.add_mutant("B"));

    BOOST_CHECK_NO_THROW(tissue_simulation.add_epigenetic_states({"E1", "E2", "E3"}));

    Tissue::mutant_view A, B;

    BOOST_CHECK_NO_THROW(A = tissue_simulation.tissue().get_mutant_view("A"));

    Species& A_E1{A["E1"]}, A_E2{A["E2"]}, A_E3{A["E3"]};
    BOOST_CHECK_NO_THROW(A_E1.set_rate(CellEventType::DEATH, 0.1));
    BOOST_CHECK_NO_THROW(A_E1.set_rate(CellEventType::DUPLICATION, 0.3));
    BOOST_CHECK_NO_THROW(A_E1.set_rate(CellEventType::DUP_AND_EPI_SWITCH, A_E2, 0.01));
    BOOST_CHECK_NO_THROW(A_E1.set_rate(CellEventType::DUP_AND_EPI_SWITCH, A_E3, 0.01));

    BOOST_CHECK_NO_THROW(A_E2.set_rate(CellEventType::DEATH, 0.1));
    BOOST_CHECK_NO_THROW(A_E2.set_rate(CellEventType::DUPLICATION, 0.45));
    BOOST_CHECK_NO_THROW(A_E2.set_rate(CellEventType::DUP_AND_EPI_SWITCH, A_E1, 0.01));
    BOOST_CHECK_NO_THROW(A_E2.set_rate(CellEventType::DUP_AND_EPI_SWITCH, A_E3, 0.02));

    BOOST_CHECK_NO_THROW(A_E3.set_rate(CellEventType::DEATH, 0.08));
    BOOST_CHECK_NO_THROW(A_E3.set_rate(CellEventType::DUPLICATION, 0.3));
    BOOST_CHECK_NO_THROW(A_E3.set_rate(CellEventType::DUP_AND_EPI_SWITCH, A_E2, 0.02));
    BOOST_CHECK_NO_THROW(A_E3.set_rate(CellEventType::DUP_AND_EPI_SWITCH, A_E3, 0.01));

    BOOST_CHECK_NO_THROW(B = tissue_simulation.tissue().get_mutant_view("B"));

    Species& B_E1{B["E1"]}, B_E2{B["E2"]}, B_E3{B["E3"]};
    BOOST_CHECK_NO_THROW(B_E1.set_rate(CellEventType::DEATH, 0.1));
    BOOST_CHECK_NO_THROW(B_E1.set_rate(CellEventType::DUPLICATION, 0.2));
    BOOST_CHECK_NO_THROW(B_E1.set_rate(CellEventType::DUP_AND_EPI_SWITCH, B_E2, 0.01));
    BOOST_CHECK_NO_THROW(B_E1.set_rate(CellEventType::DUP_AND_EPI_SWITCH, B_E3, 0.01));

    BOOST_CHECK_NO_THROW(B_E2.set_rate(CellEventType::DEATH, 0.1));
    BOOST_CHECK_NO_THROW(B_E2.set_rate(CellEventType::DUPLICATION, 0.2));
    BOOST_CHECK_NO_THROW(B_E2.set_rate(CellEventType::DUP_AND_EPI_SWITCH, B_E1, 0.01));
    BOOST_CHECK_NO_THROW(B_E2.set_rate(CellEventType::DUP_AND_EPI_SWITCH, B_E3, 0.02));

    BOOST_CHECK_NO_THROW(B_E3.set_rate(CellEventType::DEATH, 0.09));
    BOOST_CHECK_NO_THROW(B_E3.set_rate(CellEventType::DUPLICATION, 0.18));
    BOOST_CHECK_NO_THROW(B_E3.set_rate(CellEventType::DUP_AND_EPI_SWITCH, B_E2, 0.02));
    BOOST_CHECK_NO_THROW(B_E3.set_rate(CellEventType::DUP_AND_EPI_SWITCH, B_E3, 0.01));
}

BOOST_AUTO_TEST_CASE(set_rates)
{
    set_rate(tissue_simulation);
}

BOOST_AUTO_TEST_CASE(place_cell)
{
    tissue_simulation.set_tissue("Liver", {1000,1000})
                     .add_mutant("A")
                     .add_mutant("B")
                     .add_epigenetic_states({"E1", "E2", "E3"});

    auto B = tissue_simulation.tissue().get_mutant_view("B");

    BOOST_CHECK_NO_THROW(tissue_simulation.place_cell(B["E1"], {250, 500}));
}

BOOST_AUTO_TEST_CASE(schedule_mutation)
{
    tissue_simulation.set_tissue("Liver", {1000,1000})
                     .add_mutant("A")
                     .add_mutant("B")
                     .add_epigenetic_states({"E1", "E2", "E3"});

    BOOST_CHECK_NO_THROW(tissue_simulation.schedule_mutation("B", "A", 50));
}

BOOST_AUTO_TEST_CASE(simulate)
{
    using namespace CLONES::Mutants;

    tissue_simulation.set_tissue("Liver", {1000,1000})
                    .add_mutant("A")
                    .add_mutant("B")
                    .add_epigenetic_states({"E1", "E2", "E3"});

    auto A = tissue_simulation.tissue().get_mutant_view("A");

    A["E1"].set_rate(CellEventType::DEATH, 0.1)
           .set_rate(CellEventType::DUPLICATION, 0.3)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, A["E2"], 0.01)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, A["E3"], 0.01);
    A["E2"].set_rate(CellEventType::DEATH, 0.1)
           .set_rate(CellEventType::DUPLICATION, 0.45)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, A["E1"], 0.01)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, A["E3"], 0.02);
    A["E3"].set_rate(CellEventType::DEATH, 0.08)
           .set_rate(CellEventType::DUPLICATION, 0.3)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, A["E2"], 0.02)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, A["E3"], 0.01);

    auto B = tissue_simulation.tissue().get_mutant_view("B");

    B["E1"].set_rate(CellEventType::DEATH, 0.1)
           .set_rate(CellEventType::DUPLICATION, 0.2)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, B["E2"], 0.01)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, B["E3"], 0.01);
    B["E2"].set_rate(CellEventType::DEATH, 0.1)
           .set_rate(CellEventType::DUPLICATION, 0.2)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, B["E1"], 0.01)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, B["E3"], 0.02);
    B["E3"].set_rate(CellEventType::DEATH, 0.09)
           .set_rate(CellEventType::DUPLICATION, 0.18)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, B["E2"], 0.02)
           .set_rate(CellEventType::DUP_AND_EPI_SWITCH, B["E3"], 0.01);

    tissue_simulation.place_cell(B["E1"], {250, 500})
                     .schedule_mutation("B", "A", 50);

    tissue_simulation.death_activation_level = 100;
    tissue_simulation.storage_enabled = false;

    CLONES::Mutants::Evolutions::TimeTest done(70);

    BOOST_CHECK_NO_THROW(tissue_simulation.run(done));
}

BOOST_AUTO_TEST_SUITE_END()
