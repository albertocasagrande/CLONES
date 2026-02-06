/**
 * @file archive.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Some archive tests
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE archive

#include <sstream>
#include <filesystem>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "archive.hpp"
#include "simulation.hpp"
#include "ending_conditions.hpp"
#include "logger.hpp"


struct ArchiveFixture {
    long double time_horizon;

    CLONES::Mutants::Evolutions::Simulation simulation;

    ArchiveFixture():
        time_horizon(70), simulation()
    {
        using namespace CLONES::Mutants;

        MutantProperties A("A",{{0.01,0.01}});
        A["-"].set_rates({{CellEventType::DEATH, 0.1},
                          {CellEventType::DUPLICATION, 0.3}});
        A["+"].set_rates({{CellEventType::DEATH, 0.1},
                          {CellEventType::DUPLICATION, 0.45}});

        MutantProperties B("B",{{0.01,0.01}});
        B["-"].set_rates({{CellEventType::DEATH, 0.1},
                          {CellEventType::DUPLICATION, 0.2}});
        B["+"].set_rates({{CellEventType::DEATH, 0.01},
                          {CellEventType::DUPLICATION, 0.02}});

        simulation.set_tissue("Liver", {1000,1000})
                  .add_mutant(A)
                  .add_mutant(B)
                  .place_cell(B["-"], {250, 500})
                  .schedule_mutation(B,A,50);

        simulation.death_activation_level = 100;
        simulation.storage_enabled = false;

        CLONES::Mutants::Evolutions::TimeTest done(time_horizon);

        simulation.run(done);
    }

    ~ArchiveFixture()
    {}
};

std::filesystem::path get_a_temporary_path()
{
    namespace fs = std::filesystem;

    fs::path tmp_dir = fs::temp_directory_path();

    std::string basename{"file_"};
    std::string name = basename + "0";
    size_t counter = 0;

    while (fs::exists(tmp_dir / name)) {
        std::ostringstream oss;

        oss << basename << ++counter;

        name = oss.str();
    }

    return tmp_dir / name;
}

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

template<typename T>
bool operator==(const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size()!=b.size()) {
        return false;
    }

    auto a_it=a.begin(), b_it=b.begin();
    for (; a_it!=a.end(); ++a_it, ++b_it) {
        if (*a_it!=*b_it) {
            return false;
        }
    }

    return true;
}

template<typename KEY, typename T>
bool operator==(const std::map<KEY,T>& a, const std::map<KEY,T>& b)
{
    if (a.size()!=b.size()) {
        return false;
    }

    auto a_it=a.begin(), b_it=b.begin();
    for (; a_it!=a.end(); ++a_it, ++b_it) {
        if (a_it->first!=b_it->first || a_it->second!=b_it->second) {
            return false;
        }
    }

    return true;
}

bool operator==(const CLONES::Mutants::Cell& a, const CLONES::Mutants::Cell& b)
{
    return (a.get_id()==b.get_id() &&
            a.get_parent_id()==b.get_parent_id() &&
            a.get_species_id()==b.get_species_id());
}

bool operator==(const CLONES::Mutants::Evolutions::CellInTissue& a, const CLONES::Mutants::Evolutions::CellInTissue& b)
{
    using namespace  CLONES::Mutants;

    return (static_cast<const Cell&>(a)==static_cast<const Cell&>(b) &&
            static_cast<const Evolutions::PositionInTissue&>(a)==static_cast<const Evolutions::PositionInTissue&>(b));
}

bool operator==(const CLONES::Mutants::SpeciesProperties& a, const CLONES::Mutants::SpeciesProperties& b)
{
    return (a.get_name()==b.get_name() &&
            a.get_id()==b.get_id() &&
            a.get_mutant_id()==b.get_mutant_id() &&
            a.get_rates()==b.get_rates() &&
            a.get_epigenetic_switch_rates()==b.get_epigenetic_switch_rates());
}

inline bool operator!=(const CLONES::Mutants::SpeciesProperties& a, const CLONES::Mutants::SpeciesProperties& b)
{
    return !(a==b);
}

bool operator==(const CLONES::Mutants::Evolutions::Species& a, const CLONES::Mutants::Evolutions::Species& b)
{
    using namespace CLONES::Mutants;

    if (static_cast<const SpeciesProperties&>(a)!=static_cast<const SpeciesProperties&>(b)) {
        return false;
    }

    auto cell_a_it = a.begin(), cell_b_it = b.begin();

    for (; cell_a_it != a.end(); ++cell_a_it, ++cell_b_it) {
        if (*cell_a_it != *cell_b_it) {
            return false;
        }
    }

    return true;
}

bool operator==(const CLONES::Mutants::Evolutions::Tissue& a, const CLONES::Mutants::Evolutions::Tissue& b)
{
    using namespace CLONES::Mutants;

    if (a.get_name()!=b.get_name()) {
        return false;
    }

    std::set<MutantId> mutant_ids;

    auto a_it = a.begin(), b_it = b.begin();

    for (; a_it != a.end(); ++a_it, ++b_it) {

        if (*a_it != *b_it) {
            return false;
        }
        mutant_ids.insert(a_it->get_mutant_id());
    }

    for (const auto& mutant_id : mutant_ids) {
        auto a_view = a.get_mutant_species(mutant_id);
        auto b_view = b.get_mutant_species(mutant_id);

        auto a_it = a_view.begin(), b_it = b_view.begin();

        for (; a_it != a_view.end(); ++a_it, ++b_it) {

            if (a_it->get_id() != b_it->get_id()) {
                return false;
            }
        }
    }

    return true;
}

bool operator==(const CLONES::Mutants::Evolutions::Simulation& a, const CLONES::Mutants::Evolutions::Simulation& b)
{
    return a.get_time()==b.get_time() &&
           a.tissue()==b.tissue() &&
           a.death_activation_level==b.death_activation_level;
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

    std::vector<TimedEvent> to_save{{5,SimulationEventWrapper(Mutation(0,1))},
                                    {3.5,SimulationEventWrapper(Mutation(1,7))},
                                    {8.1,SimulationEventWrapper(Mutation(2,1))}};

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

    std::vector<TimedEvent> to_save{{5,SimulationEventWrapper(Mutation(0,1))},
                                    {3.5,SimulationEventWrapper(Mutation(1,7))},
                                    {8.1,SimulationEventWrapper(Mutation(2,1))}};

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

BOOST_AUTO_TEST_CASE(binary_epigenetic_mutant)
{
    using namespace CLONES::Mutants;

    MutantProperties to_save("A",{{0.01,0.01},{0.01,0.01}});
    to_save["--"].set_rates({{CellEventType::DEATH, 0.1},
                             {CellEventType::DUPLICATION, 0.3}});
    to_save["+-"].set_rates({{CellEventType::DEATH, 0.1},
                             {CellEventType::DUPLICATION, 0.45}});
    to_save["-+"].set_rates({{CellEventType::DEATH, 0.1},
                             {CellEventType::DUPLICATION, 0.2}});
    to_save["+-"].set_rates({{CellEventType::DEATH, 0.01},
                             {CellEventType::DUPLICATION, 0.02}});

    auto filename = get_a_temporary_path();
    {
        CLONES::Archive::Binary::Out o_archive(filename);

        for (const auto& species : to_save.get_species()) {
            o_archive & species;
        }
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        for (const auto& species : to_save.get_species()) {
            SpeciesProperties in_species =  SpeciesProperties::load(i_archive);

            BOOST_CHECK(species==in_species);
        }
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

        o_archive & simulation.tissue();
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        Tissue in_tissue = Tissue::load(i_archive);

        BOOST_CHECK(simulation.tissue()==in_tissue);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(simulation_statistics)
{
    using namespace CLONES::Mutants::Evolutions;

    auto filename = get_a_temporary_path();

    {
        CLONES::Archive::Binary::Out o_archive(filename);

        o_archive & simulation.get_statistics();
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        TissueStatistics statistics = TissueStatistics::load(i_archive);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(simulation_tissue)
{
    using namespace CLONES::Mutants::Evolutions;

    auto filename = get_a_temporary_path();

    {
        CLONES::Archive::Binary::Out o_archive(filename);

        o_archive & simulation;
    }

    {
        CLONES::Archive::Binary::In i_archive(filename);

        Simulation in_simulation = Simulation::load(i_archive);

        BOOST_CHECK(simulation==in_simulation);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_SUITE_END()
