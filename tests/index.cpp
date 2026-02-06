/**
 * @file index.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Index tests
 * @version 1.2
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
#define BOOST_TEST_MODULE index

#include <set>
#include <map>
#include <random>
#include <type_traits>
#include <functional>
#include <utility>
#include <random>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "index.hpp"
#include "genomic_position.hpp"

#include "utils.hpp"

#define DATASET_SIZE 10000
#define WRITE_CACHE_SIZE 7000
#define READ_CACHE_SIZE 9000
#define NUM_OF_CHOICES 800

using namespace CLONES::Mutations;

typedef boost::mpl::list<
    std::pair<GenomicPosition, GenomicPosition>,
    std::pair<size_t, GenomicPosition>,
    std::pair<GenomicPosition, size_t>,
    std::pair<size_t, size_t>
> test_types;

template<typename TYPE>
inline TYPE create_data(const size_t& i)
{
    return i;
}

template<>
inline CLONES::Mutations::GenomicPosition create_data(const size_t& i)
{
    return CLONES::Mutations::GenomicPosition(i%22, static_cast<uint32_t>(i));
}

template<class PAIR>
struct IndexFixture
{
    using key_type  = typename PAIR::first_type;
    using value_type = typename PAIR::second_type;

    using IndexBuilderType = CLONES::Archive::IndexBuilder<key_type, value_type>;

    std::filesystem::path index_path;
    std::map<key_type, std::set<value_type>> dataset;

    IndexFixture():
        index_path(get_a_temporary_path())
    {
        IndexBuilderType builder(index_path, WRITE_CACHE_SIZE);
        for (size_t i=0; i<DATASET_SIZE; ++i) {

            key_type key = create_data<key_type>(i%10);
            value_type value = create_data<value_type>(i);
            BOOST_CHECK_NO_THROW(builder.insert(key, value));

            auto found = dataset.find(key);
            if (found != dataset.end()) {
                found->second.insert(std::move(value));
            } else {
                dataset.emplace(std::move(key),
                                std::set<value_type>({std::move(value)}));
            }
        }

        std::mt19937_64 random_generator(0);

        builder.shuffle(random_generator);
    }

    ~IndexFixture()
    {
        std::filesystem::remove_all(index_path);
    }
};


template<class KEY, class VALUE>
void create_index()
{
    using namespace CLONES::Archive;
    namespace fs = std::filesystem;

    fs::path index_path = get_a_temporary_path();

    BOOST_CHECK(!fs::exists(index_path));

    IndexBuilder<KEY, VALUE> builder(index_path,
                                     WRITE_CACHE_SIZE);

    BOOST_CHECK(fs::exists(index_path));

    fs::remove(index_path);

    using IndexBuilderType = IndexBuilder<KEY,VALUE>;

    BOOST_CHECK_THROW(IndexBuilderType("/"), std::domain_error);

    BOOST_CHECK_THROW(IndexBuilderType("/Pippo"), std::filesystem::filesystem_error);

    BOOST_CHECK_THROW(IndexBuilderType(get_a_temporary_path(), 0),
                      std::domain_error);
}

template<class KEY, class VALUE>
void test_random_function(const std::filesystem::path& index_path,
                          const std::map<KEY, std::set<VALUE>>& dataset,
                          const std::function<void(std::mt19937_64&, const std::filesystem::path&, const std::map<KEY, std::set<VALUE>>&)>& test_function)
{
    {
        std::mt19937_64 random_generator(0);

        test_function(random_generator, index_path, dataset);
    }

    {
        std::mt19937_64 random_generator(2);

        test_function(random_generator, index_path, dataset);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(create_bucket_T, T, test_types)
{
    using namespace CLONES::Archive;
    namespace fs = std::filesystem;

    fs::path index_path = get_a_temporary_path();

    BOOST_CHECK(!fs::exists(index_path));

    IndexBuilder<typename T::first_type,
                 typename T::second_type> builder(index_path,
                                                  WRITE_CACHE_SIZE);

    BOOST_CHECK(fs::exists(index_path));

    fs::remove(index_path);

    using IndexBuilderType = IndexBuilder<typename T::first_type,
                                          typename T::second_type>;

    BOOST_CHECK_THROW(IndexBuilderType("/"), std::domain_error);

    BOOST_CHECK_THROW(IndexBuilderType("/Pippo"),
                      std::filesystem::filesystem_error);

    BOOST_CHECK_THROW(IndexBuilderType(get_a_temporary_path(), 0),
                      std::domain_error);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(load_index_T, T, test_types, IndexFixture<T>)
{
    using namespace CLONES::Archive;

    using IndexType = IndexReader<typename T::first_type,
                                  typename T::second_type,
                                  std::mt19937_64>;

    BOOST_CHECK_NO_THROW(IndexType{this->index_path});
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(index_key_bucket_access_T, T, test_types, IndexFixture<T>)
{
    using namespace CLONES::Archive;

    using IndexType = IndexReader<typename T::first_type,
                                  typename T::second_type,
                                  std::mt19937_64>;
    IndexType index(this->index_path, READ_CACHE_SIZE);

    BOOST_CHECK(this->dataset.size()==index.num_of_keys());

    BOOST_CHECK_NO_THROW(index.get_keys());
    for (const auto& [key, values]: this->dataset) {
        BOOST_CHECK(index[key].size()==values.size());
    }
}


template<class KEY, class VALUE>
void index_extract(std::mt19937_64& random_generator,
                   const std::filesystem::path& index_path,
                   const std::map<KEY, std::set<VALUE>>& dataset)
{
    using namespace CLONES::Archive;
    std::map<KEY, std::set<VALUE>> local_dataset(dataset);

    IndexReader<KEY, VALUE, std::mt19937_64> index(index_path, READ_CACHE_SIZE);

    for (auto& [key, values]: local_dataset) {
        size_t num_of_values = values.size();
        for (size_t i=0;i<num_of_values; ++i) {
            VALUE value;
            BOOST_CHECK_NO_THROW(value = index.extract(random_generator, key));

            auto found = values.find(value);

            BOOST_CHECK(found != values.end());

            if (found != values.end()) {
                values.erase(found);
            }
        }

        BOOST_CHECK(values.size() == 0);
        BOOST_CHECK(index.extractable_for(key) == 0);
        BOOST_CHECK(index.num_of_values(key) == dataset.at(key).size());
    }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(index_extract_T, T, test_types, IndexFixture<T>)
{
    using key_type = typename T::first_type;
    using value_type = typename T::second_type;

    test_random_function<key_type, value_type>(this->index_path, this->dataset, index_extract<key_type, value_type>);
}

template<class KEY, class VALUE>
void index_choose(std::mt19937_64& random_generator,
                  const std::filesystem::path& index_path,
                  const std::map<KEY, std::set<VALUE>>& dataset)
{
    using namespace CLONES::Archive;

    IndexReader<KEY, VALUE, std::mt19937_64> index(index_path, READ_CACHE_SIZE);

    for (const auto& [key, values]: dataset) {
        for (size_t i=0;i<NUM_OF_CHOICES; ++i) {
            VALUE value;
            BOOST_CHECK_NO_THROW(value = index.choose(random_generator, key));

            auto found = values.find(value);

            BOOST_CHECK(found != values.end());
        }

        BOOST_CHECK(index.extractable_for(key) == values.size());
        BOOST_CHECK(index.num_of_values(key) == values.size());
    }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(index_choose_T, T, test_types, IndexFixture<T>)
{
    using key_type = typename T::first_type;
    using value_type = typename T::second_type;

    test_random_function<key_type, value_type>(this->index_path, this->dataset, index_choose<key_type, value_type>);
}

template<class KEY, class VALUE>
void index_extract_class(std::mt19937_64& random_generator,
                         const std::filesystem::path& index_path,
                         const std::map<KEY, std::set<VALUE>>& dataset)
{
    using namespace CLONES::Archive;
    std::map<KEY, std::set<VALUE>> local_dataset(dataset);

    IndexReader<KEY, VALUE, std::mt19937_64> index(index_path, READ_CACHE_SIZE);

    std::set<KEY> key_done;

    for (const auto key : index.get_keys()) {
        if (!key_done.contains(key)) {
            auto class_keys = CLONES::partition<KEY>::get_class_of(key);
            key_done.insert(class_keys.begin(), class_keys.end());

            std::set<KEY> class_set(class_keys.begin(), class_keys.end());

            size_t num_of_values;

            BOOST_CHECK_NO_THROW(num_of_values = index.num_of_class_values(key));

            for (size_t i=0; i<num_of_values; ++i) {
                std::pair<KEY, VALUE> extracted;

                BOOST_CHECK_NO_THROW(extracted = index.extract_from_class(random_generator, key));

                auto found = class_set.find(extracted.first);

                BOOST_CHECK(found != class_set.end());

                if (found != class_set.end()) {
                    auto found_dataset = local_dataset.find(extracted.first);

                    BOOST_CHECK(found_dataset != local_dataset.end());

                    if (found_dataset != local_dataset.end()) {
                        std::set<VALUE>& bucket_set = found_dataset->second;

                        auto set_found = bucket_set.find(extracted.second);
                        BOOST_CHECK(set_found != bucket_set.end());

                        if (set_found != bucket_set.end()) {
                            bucket_set.erase(set_found);
                        }
                    }
                }
            }

            size_t total_dataset_class{0};
            for (const auto& class_key : class_set) {
                auto found_dataset = local_dataset.find(class_key);

                if (found_dataset != local_dataset.end()) {
                    total_dataset_class += found_dataset->second.size();
                }
            }

            BOOST_CHECK(index.extractable_from_class(key) == 0);
            BOOST_CHECK(total_dataset_class == 0);
        }

        return;
    }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(index_extract_class_T, T, test_types, IndexFixture<T>)
{
    using key_type = typename T::first_type;
    using value_type = typename T::second_type;

    test_random_function<key_type, value_type>(this->index_path, this->dataset, index_extract_class<key_type, value_type>);
}

template<class KEY, class VALUE>
void index_choose_class(std::mt19937_64& random_generator,
                        const std::filesystem::path& index_path,
                        const std::map<KEY, std::set<VALUE>>& dataset)
{
    using namespace CLONES::Archive;
    IndexReader<KEY, VALUE, std::mt19937_64> index(index_path, READ_CACHE_SIZE);

    std::set<KEY> key_done;

    size_t total_available{0};
    for (const auto key : index.get_keys()) {
        if (!key_done.contains(key)) {
            auto class_keys = CLONES::partition<KEY>::get_class_of(key);
            key_done.insert(class_keys.begin(), class_keys.end());

            std::set<KEY> class_set(class_keys.begin(), class_keys.end());

            for (size_t i=0;i<NUM_OF_CHOICES; ++i) {
                std::pair<KEY, VALUE> choosen;

                BOOST_CHECK_NO_THROW(choosen = index.choose_from_class(random_generator, key));

                auto found = class_set.find(choosen.first);

                BOOST_CHECK(found != class_set.end());

                if (found != class_set.end()) {
                    auto found_dataset = dataset.find(choosen.first);

                    BOOST_CHECK(found_dataset != dataset.end());

                    if (found_dataset != dataset.end()) {
                        const std::set<VALUE>& bucket_set = found_dataset->second;

                        const auto set_found = bucket_set.find(choosen.second);
                        BOOST_CHECK(set_found != bucket_set.end());
                    }
                }
            }

            size_t total_dataset_class{0};
            for (const auto& class_key : class_set) {
                auto found_dataset = dataset.find(class_key);

                if (found_dataset != dataset.end()) {
                    total_dataset_class += found_dataset->second.size();
                }
            }

            total_available += index.extractable_from_class(key);
            BOOST_CHECK(index.extractable_from_class(key) == index.num_of_class_values(key));
            BOOST_CHECK(total_dataset_class == index.num_of_class_values(key));
        }
    }

    BOOST_CHECK(total_available==DATASET_SIZE);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(index_choose_class_T, T, test_types, IndexFixture<T>)
{
    using key_type = typename T::first_type;
    using value_type = typename T::second_type;

    test_random_function<key_type, value_type>(this->index_path, this->dataset, index_choose_class<key_type, value_type>);
}
