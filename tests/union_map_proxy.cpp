/**
 * @file union_map_proxy.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing CLONES::union_map_proxy class
 * @version 1.0
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
#define BOOST_TEST_MODULE union_map_proxy

#include <boost/test/unit_test.hpp>

#include <string>

#include "union_map_proxy.hpp"

struct union_fixture
{
    using map_type = std::map<int, std::string>;
    using result_list_type =  std::list<std::pair<int, std::string>>;

    struct test_type
    {
        map_type map1;
        map_type map2;

        result_list_type result;

        test_type(map_type&& map1, map_type&& map2, result_list_type&& result):
            map1{std::move(map1)}, map2{std::move(map2)}, result{std::move(result)}
        {}
    };

    std::list<test_type> test_union_map_proxys;

    union_fixture():
        test_union_map_proxys{{{{1, "a"}, {2, "b"}, {3, "c"}, {4, "d"}},
                         {{1, "A"}, {5, "E"}},
                         {{1, "a"}, {1, "A"}, {2, "b"}, {3, "c"}, {4, "d"}, {5, "E"}}},
                        {{{1, "a"}, {2, "b"}, {3, "c"}, {5, "d"}},
                         {{1, "A"}, {5, "E"}},
                         {{1, "a"}, {1, "A"}, {2, "b"}, {3, "c"}, {5, "d"}, {5, "E"}}},
                        {{{1, "a"}, {2, "b"}, {3, "c"}, {5, "d"}},
                         {{1, "A"}, {4, "E"}},
                         {{1, "a"}, {1, "A"}, {2, "b"}, {3, "c"}, {4, "E"}, {5, "d"}}},
                        {{{1, "A"}, {5, "E"}},
                         {{1, "a"}, {2, "b"}, {3, "c"}, {4, "d"}},
                         {{1, "A"}, {1, "a"}, {2, "b"}, {3, "c"}, {4, "d"}, {5, "E"}}},
                        {{{1, "A"}, {5, "E"}},
                         {{1, "a"}, {2, "b"}, {3, "c"}, {5, "d"}},
                         {{1, "A"}, {1, "a"}, {2, "b"}, {3, "c"}, {5, "E"}, {5, "d"}}},
                        {{{1, "A"}, {4, "E"}},
                         {{1, "a"}, {2, "b"}, {3, "c"}, {5, "d"}},
                         {{1, "A"}, {1, "a"}, {2, "b"}, {3, "c"}, {4, "E"}, {5, "d"}}}}
    {}
};


BOOST_FIXTURE_TEST_SUITE( union_test, union_fixture )

BOOST_AUTO_TEST_CASE(union_map_proxy_creation)
{
    using namespace CLONES;

    using union_map_proxy_type = union_map_proxy<int, std::string>;

    for (const auto& test : test_union_map_proxys) {
        BOOST_CHECK_NO_THROW(union_map_proxy_type(test.map1, test.map2));
    }
}

BOOST_AUTO_TEST_CASE(union_map_proxy_getting_begin)
{
    using namespace CLONES;

    for (const auto& test : test_union_map_proxys) {
        union_map_proxy<int, std::string> obj{test.map1, test.map2};

        BOOST_CHECK_NO_THROW(obj.begin());
    }
}

BOOST_AUTO_TEST_CASE(union_map_proxy_getting_end)
{
    using namespace CLONES;

    for (const auto& test : test_union_map_proxys) {
        union_map_proxy<int, std::string> obj{test.map1, test.map2};

        BOOST_CHECK_NO_THROW(obj.end());
    }
}

BOOST_AUTO_TEST_CASE(union_map_proxy_forward_iteration)
{
    using namespace CLONES;

    for (const auto& test : test_union_map_proxys) {
        union_map_proxy<int, std::string> union_obj{test.map1, test.map2};

        auto result_it = test.result.begin();
        auto union_it = union_obj.begin();

        for (;result_it != test.result.end(); ++union_it, ++result_it) {
            BOOST_CHECK_EQUAL(result_it->first, union_it->first);
            BOOST_CHECK_EQUAL(result_it->second, union_it->second);
        }

        BOOST_CHECK_EQUAL(union_it==union_obj.end(), true);
    }
}

BOOST_AUTO_TEST_CASE(union_map_proxy_backward_iteration)
{
    using namespace CLONES;

    for (const auto& test : test_union_map_proxys) {
        union_map_proxy<int, std::string> union_obj{test.map1, test.map2};

        auto result_it = test.result.end();
        auto union_it = union_obj.end();

        do {
            --result_it;
            --union_it;

            BOOST_CHECK_EQUAL(result_it->first, union_it->first);
            BOOST_CHECK_EQUAL(result_it->second, union_it->second);
        } while (result_it != test.result.begin());

        BOOST_CHECK_EQUAL(union_it==union_obj.begin(), true);
    }
}

BOOST_AUTO_TEST_CASE(union_map_proxy_find_iteration)
{
    using namespace CLONES;

    for (const auto& test : test_union_map_proxys) {
        union_map_proxy<int, std::string> union_obj{test.map1, test.map2};

        if (test.result.size() > 0) {
            int max_key{test.result.front().first};

            union_map_proxy<int, std::string>::const_iterator union_it;

            auto result_it = test.result.begin();
            while (result_it != test.result.end()) {

                const auto key = result_it->first;
                max_key = std::max(max_key, key);

                BOOST_CHECK_NO_THROW(union_it = union_obj.find(key));

                BOOST_CHECK_EQUAL(key, union_it->first);
                BOOST_CHECK_EQUAL(result_it->second, union_it->second);

                bool repeat{true};
                do {
                    ++result_it;
                    if (result_it != test.result.end()
                             && result_it->first == key) {

                        ++union_it;
                        BOOST_CHECK(union_it != union_obj.end());
                        BOOST_CHECK_EQUAL(key, union_it->first);
                        BOOST_CHECK_EQUAL(result_it->second, union_it->second);
                    } else {
                        repeat = false;
                    }
                } while (repeat);
            }

            BOOST_CHECK_NO_THROW(union_it = union_obj.find(max_key+1));
            BOOST_CHECK_EQUAL(union_it == union_obj.end(), true);
        }
    }
}

BOOST_AUTO_TEST_CASE(union_map_proxy_lower_bound_iteration)
{
    using namespace CLONES;

    for (const auto& test : test_union_map_proxys) {
        union_map_proxy<int, std::string> union_obj{test.map1, test.map2};

        if (test.result.size() > 0) {
            int min_key{test.result.front().first};
            int max_key{min_key};

            for (const auto& [key, value] : test.result) {
                min_key = std::min(min_key, key);
                max_key = std::max(max_key, key);
            }

            auto result_it = test.result.begin();
            for (int key=min_key-1; key<max_key+2; ++key) {
                while (result_it != test.result.end()
                        && result_it->first < key) {
                    ++result_it;
                }
                union_map_proxy<int, std::string>::const_iterator union_it;
                BOOST_CHECK_NO_THROW(union_it = union_obj.lower_bound(key));

                BOOST_CHECK_EQUAL(result_it->first, union_it->first);
                BOOST_CHECK_EQUAL(result_it->second, union_it->second);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(union_map_proxy_upper_bound_iteration)
{
    using namespace CLONES;

    for (const auto& test : test_union_map_proxys) {
        union_map_proxy<int, std::string> union_obj{test.map1, test.map2};

        if (test.result.size() > 0) {
            int min_key{test.result.front().first};
            int max_key{min_key};

            for (const auto& [key, value] : test.result) {
                min_key = std::min(min_key, key);
                max_key = std::max(max_key, key);
            }

            auto result_it = test.result.begin();
            for (int key=min_key-1; key<max_key+2; ++key) {

                while (result_it != test.result.end()
                        && result_it->first <= key) {
                    ++result_it;
                }
                union_map_proxy<int, std::string>::const_iterator union_it;
                BOOST_CHECK_NO_THROW(union_it = union_obj.upper_bound(key));

                BOOST_CHECK_EQUAL(result_it->first, union_it->first);
                BOOST_CHECK_EQUAL(result_it->second, union_it->second);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
