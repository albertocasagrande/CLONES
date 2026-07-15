/**
 * @file archive.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Some archive tests
 * @version 1.7
 * @date 2026-07-14
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
#define BOOST_TEST_MODULE archive

#include <sstream>
#include <filesystem>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "archive.hpp"
#include "tissue_simulation.hpp"
#include "ending_conditions.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include "error.hpp"

#include "test_common.hpp"

struct ArchiveFixture {
    double time_horizon;

    CLONES::Mutants::Evolutions::TissueSimulation tissue_simulation;

    ArchiveFixture():
        time_horizon(70), tissue_simulation()
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
                         .schedule_mutation(B, A, 50);

        tissue_simulation.death_activation_level = 100;
        tissue_simulation.storage_enabled = false;

        CLONES::Mutants::Evolutions::TimeTest done(time_horizon);

        tissue_simulation.run(done);
    }

    ~ArchiveFixture()
    {
        std::filesystem::remove_all(tissue_simulation.get_logger().get_directory());
    }
};

template<typename T>
void basic_type_test(const std::vector<T>& to_save)
{
    auto filename = get_a_temporary_path();
    {
        CLONES::Archive::Binary::Out o_archive(filename);

        for (const auto& value : to_save) {
            o_archive & value;
        }
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        for (const auto& value : to_save) {
            T read_value;

            i_archive & read_value;

            BOOST_CHECK(read_value==value);
        }
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(binary_size_t)
{
    basic_type_test<size_t>({3, 4, 0, 5});
}

BOOST_AUTO_TEST_CASE(binary_double)
{
    basic_type_test<double>({-3.3, 1/4, -0, 2/3});
}

BOOST_AUTO_TEST_CASE(binary_char)
{
    basic_type_test<char>({'a', '.', '\n', '\0'});
}

BOOST_AUTO_TEST_CASE(binary_string)
{
    basic_type_test<std::string>({"string", "", "\n", "\"ci\0ao\""});
}

BOOST_AUTO_TEST_CASE(binary_vector)
{
    basic_type_test<std::vector<size_t>>({{},{3, 4}, {}, {0}, {5}});
    basic_type_test<std::vector<double>>({{-3.3}, {}, {1/4, -0},{}, {2/3}});
    basic_type_test<std::vector<char>>({{},{},{'a', '.'}, {}, {'\n', '\0'},{}});
    basic_type_test<std::vector<std::string>>({{},{"string", ""}, {}, {"\n"}, {"\"ci\0ao\""}});
}

struct TestData
{
    size_t test_member[1000];

    TestData()
    {
        for (size_t i=0; i<1000; ++i) {
            test_member[i] = i;
        }
    }

    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<CLONES::Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        for (size_t i=0; i<1000; ++i) {
            archive & test_member[i];
        }
    }

    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<CLONES::Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static TestData load(ARCHIVE& archive)
    {
        TestData test_data;

        for (size_t i=0; i<1000; ++i) {
            archive & test_data.test_member[i];
        }

        return test_data;
    }

    bool operator==(const TestData& a) const
    {
        for (size_t i=0; i<1000; ++i) {
            if (test_member[i]!=a.test_member[i]) {
                return false;
            }
        }

        return true;
    }

    inline bool operator!=(const TestData& a) const
    {
        return !(*this==a);
    }
};

BOOST_AUTO_TEST_CASE(binary_dynamic_memory_object)
{
    std::filesystem::path filename = get_a_temporary_path();

    std::shared_ptr<int> a=std::make_shared<int>(3);
    std::vector<std::shared_ptr<int>> v{a, std::make_shared<int>(4)};
    std::shared_ptr<std::vector<int>> sv=std::make_shared<std::vector<int>>(10000);
    std::shared_ptr<TestData> a_struct=std::make_shared<TestData>();
    std::vector<std::shared_ptr<TestData>> v_struct{a_struct, std::make_shared<TestData>()};
    std::shared_ptr<std::vector<TestData>> sv_struct=std::make_shared<std::vector<TestData>>(10000);
    {
        CLONES::Archive::Binary::Out o_archive(filename);

        o_archive & a & a & a & a & a
                  & v & v & sv & sv
                  & a_struct & a_struct & a_struct & a_struct & a_struct
                  & v_struct & v_struct & sv_struct & sv_struct;
    }

    {
        std::shared_ptr<int> loaded_a1, loaded_a2, loaded_a3, loaded_a4, loaded_a5;
        std::vector<std::shared_ptr<int>> loaded_v1, loaded_v2;
        std::shared_ptr<std::vector<int>> loaded_sv1, loaded_sv2;
        std::shared_ptr<TestData> loaded_a_struct1, loaded_a_struct2, loaded_a_struct3,
                                  loaded_a_struct4, loaded_a_struct5;
        std::vector<std::shared_ptr<TestData>> loaded_v_struct1, loaded_v_struct2;
        std::shared_ptr<std::vector<TestData>> loaded_sv_struct1, loaded_sv_struct2;

        CLONES::Archive::Binary::In i_archive(filename);

        i_archive & loaded_a1 & loaded_a2 & loaded_a3
                  & loaded_a4 & loaded_a5 & loaded_v1
                  & loaded_v2 & loaded_sv1 & loaded_sv2
                  & loaded_a_struct1 & loaded_a_struct2
                  & loaded_a_struct3 & loaded_a_struct4
                  & loaded_a_struct5 & loaded_v_struct1
                  & loaded_v_struct2 & loaded_sv_struct1
                  & loaded_sv_struct2;

        BOOST_CHECK(*loaded_a1==*a);
        BOOST_CHECK(loaded_a2.get()==loaded_a1.get());
        BOOST_CHECK(loaded_a3.get()==loaded_a1.get());
        BOOST_CHECK(loaded_a4.get()==loaded_a1.get());
        BOOST_CHECK(loaded_a5.get()==loaded_a1.get());
        BOOST_CHECK(loaded_v1.size()==v.size());
        BOOST_CHECK(loaded_v1[0].get()==loaded_a1.get());
        BOOST_CHECK(loaded_v2[0].get()==loaded_a1.get());
        BOOST_CHECK(loaded_sv1->size()==sv->size());
        BOOST_CHECK(loaded_sv1.get()==loaded_sv2.get());
        BOOST_CHECK(*loaded_a_struct1==*a_struct);
        BOOST_CHECK(loaded_a_struct2.get()==loaded_a_struct1.get());
        BOOST_CHECK(loaded_a_struct3.get()==loaded_a_struct1.get());
        BOOST_CHECK(loaded_a_struct4.get()==loaded_a_struct1.get());
        BOOST_CHECK(loaded_a_struct5.get()==loaded_a_struct1.get());
        BOOST_CHECK(loaded_v_struct1[0].get()==loaded_a_struct1.get());
        BOOST_CHECK(loaded_v_struct2[0].get()==loaded_a_struct1.get());
        BOOST_CHECK(loaded_sv_struct1->size()==sv_struct->size());
        BOOST_CHECK(loaded_sv_struct2.get()==loaded_sv_struct1.get());
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(binary_list)
{
    basic_type_test<std::list<size_t>>({{},{3, 4}, {}, {0}, {5}});
    basic_type_test<std::list<double>>({{-3.3}, {}, {1/4, -0},{}, {2/3}});
    basic_type_test<std::list<char>>({{},{},{'a', '.'}, {}, {'\n', '\0'},{}});
    basic_type_test<std::list<std::string>>({{},{"string", ""}, {}, {"\n"}, {"\"ci\0ao\""}});
}

BOOST_AUTO_TEST_CASE(binary_map)
{
    std::map<std::string, std::vector<std::vector<double>>> to_save{{"ciao", {{}, {-0.3}, {5,6}}},
                                                                    {"", {}},
                                                                    {"ult\0imo", {{-3/7},{}}}};
    auto filename = get_a_temporary_path();
    {
        CLONES::Archive::Binary::Out o_archive(filename);

        o_archive & to_save;
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        std::map<std::string, std::vector<std::vector<double>>> read_value;

        i_archive & read_value;

        BOOST_CHECK(read_value==to_save);
    }

    std::filesystem::remove(filename);
}


BOOST_AUTO_TEST_CASE(binary_timed_mutation)
{
    using namespace CLONES::Mutants::Evolutions;

    std::vector<TimedEvent> to_save{{5,TissueSimulationEventWrapper(Mutation(0,1))},
                                    {3.5,TissueSimulationEventWrapper(Mutation(1,7))},
                                    {8.1,TissueSimulationEventWrapper(Mutation(2,1))}};

    auto filename = get_a_temporary_path();
    {
        CLONES::Archive::Binary::Out o_archive(filename);

        for (const auto& value : to_save) {
            o_archive & value;
        }
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        for (auto& value : to_save) {
            TimedEvent read_value = TimedEvent::load(i_archive);

            BOOST_CHECK(read_value==value);
        }
    }

    std::filesystem::remove(filename);
}


BOOST_AUTO_TEST_CASE(binary_timed_mutation_queue)
{
    using namespace CLONES::Mutants::Evolutions;

    std::vector<TimedEvent> to_save{{5,TissueSimulationEventWrapper(Mutation(0,1))},
                                    {3.5,TissueSimulationEventWrapper(Mutation(1,7))},
                                    {8.1,TissueSimulationEventWrapper(Mutation(2,1))}};

    using PriorityQueue = std::priority_queue<TimedEvent,
                                              std::vector<TimedEvent>,
                                              std::greater<TimedEvent>>;

    PriorityQueue queue(to_save.begin(), to_save.end());

    auto filename = get_a_temporary_path();
    {
        CLONES::Archive::Binary::Out o_archive(filename);

        o_archive & queue;
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        PriorityQueue i_queue;

        i_archive & i_queue;

        BOOST_CHECK(i_queue.size()==queue.size());

        while (!i_queue.empty() && !queue.empty()) {

            BOOST_CHECK(i_queue.top()==queue.top());
            i_queue.pop();
            queue.pop();
        }

        BOOST_CHECK(i_queue.size()==queue.size());
    }

    std::filesystem::remove(filename);
}

BOOST_FIXTURE_TEST_SUITE( simulatedData, ArchiveFixture )

BOOST_AUTO_TEST_CASE(binary_tissue)
{
    using namespace CLONES::Mutants::Evolutions;

    auto filename = get_a_temporary_path();

    {
        CLONES::Archive::Binary::Out o_archive(filename);

        o_archive & tissue_simulation.tissue();
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        Tissue in_tissue = Tissue::load(i_archive);

        BOOST_CHECK(tissue_simulation.tissue()==in_tissue);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(tissue_simulation_statistics)
{
    using namespace CLONES::Mutants::Evolutions;

    auto filename = get_a_temporary_path();

    {
        CLONES::Archive::Binary::Out o_archive(filename);

        o_archive & tissue_simulation.get_statistics();
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        TissueStatistics statistics = TissueStatistics::load(i_archive);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(tissue_simulation_tissue)
{
    using namespace CLONES::Mutants::Evolutions;

    auto filename = get_a_temporary_path();

    {
        CLONES::Archive::Binary::Out o_archive(filename);

        o_archive & tissue_simulation;
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        TissueSimulation in_simulation = TissueSimulation::load(i_archive);

        BOOST_CHECK(tissue_simulation==in_simulation);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_SUITE_END()
