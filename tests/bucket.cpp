/**
 * @file bucket.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Bucket tests
 * @version 1.3
 * @date 2025-11-27
 *
 * @copyright Copyright (c) 2023-2025
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
#define BOOST_TEST_MODULE bucket

#include <stdint.h>
#include <type_traits>
#include <functional>
#include <random>
#include <set>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "utils.hpp"
#include "genomic_position.hpp"

#include "bucket.hpp"

#define DEFAULT_DATASET_SIZE 10000
#define DEFAULT_WRITE_CACHE_SIZE 700
#define DEFAULT_READ_CACHE_SIZE 900

typedef boost::mpl::list<size_t, RACES::Mutations::GenomicPosition> test_types;

template<typename TYPE>
inline TYPE create_data(const size_t& i)
{
    return i;
}

template<>
inline RACES::Mutations::GenomicPosition create_data(const size_t& i)
{
    return RACES::Mutations::GenomicPosition(i%22, static_cast<uint32_t>(i));
}

template<typename TYPE>
struct BucketFixture
{
    std::filesystem::path bucket_filepath;
    std::set<TYPE> dataset;

    BucketFixture():
        bucket_filepath{get_a_temporary_path()}
    {
        RACES::Archive::BucketWriter<TYPE> bucket_writer{bucket_filepath, DEFAULT_WRITE_CACHE_SIZE};

        for (size_t i=0; i<DEFAULT_DATASET_SIZE; ++i) {
            TYPE data = create_data<TYPE>(i);
            BOOST_CHECK_NO_THROW(bucket_writer.push_back(data));
            dataset.emplace(std::move(data));
        }

        bucket_writer.flush();
    }

    ~BucketFixture()
    {
        std::filesystem::remove(bucket_filepath);
    }
};

template<typename TYPE, typename RANDOM_GENERATOR>
TYPE
test_random_tour_on(const RACES::Archive::BucketRandomTour<TYPE, RANDOM_GENERATOR>& tour,
                    const std::set<TYPE>& dataset)
{
    std::set<TYPE> local_dataset(dataset.begin(), dataset.end());

    TYPE last_value;
    for (TYPE value : tour) {
        auto found = local_dataset.find(value);
        BOOST_CHECK(found != local_dataset.end());

        last_value = value;
        if (found != local_dataset.end()) {
            local_dataset.erase(found);
        }
    }

    BOOST_CHECK(local_dataset.empty());

    return last_value;
}

template<typename TYPE>
void shuffle_bucket(const std::filesystem::path& bucket_filepath,
                    const std::set<TYPE>& dataset,
                    const size_t cache_size)
{
    using namespace RACES::Archive;
    namespace fs = std::filesystem;

    BucketWriter<TYPE> bucket_writer{bucket_filepath, cache_size};

    std::mt19937_64 gen(0);

    if (cache_size < 2*sizeof(TYPE)) {
        throw std::domain_error("Minimum cache size (i.e., " + std::to_string(sizeof(TYPE))
                                + ") not respected. Parameter \"cache_size\" set to "
                                + std::to_string(cache_size) + ".");
    }

    {
        // minimum buffer size not respected
        BOOST_CHECK_THROW(bucket_writer.shuffle(gen, sizeof(TYPE),
                                                fs::temp_directory_path()), std::domain_error);
    }

    {
        BOOST_CHECK_NO_THROW(bucket_writer.shuffle(gen, cache_size,
                                                   fs::temp_directory_path()));
    }

    {
        BucketReader<TYPE> bucket_reader(bucket_filepath, DEFAULT_READ_CACHE_SIZE);
        std::set<TYPE> local_dataset(dataset.begin(), dataset.end());

        BOOST_CHECK(bucket_reader.size()==dataset.size());

        for (const TYPE& data: bucket_reader) {
            auto found = local_dataset.find(data);
            BOOST_CHECK(found != dataset.end());

            if (found != local_dataset.end()) {
                local_dataset.erase(found);
            }
        }
        BOOST_CHECK(local_dataset.empty());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(create_bucket_T, T, test_types)
{
    using namespace RACES::Archive;
    namespace fs = std::filesystem;

    fs::path bucket_filepath = get_a_temporary_path();

    BOOST_CHECK(!fs::exists(bucket_filepath));

    BucketWriter<T> bucket_writer(bucket_filepath);

    BOOST_CHECK(fs::exists(bucket_filepath));
    BOOST_CHECK(bucket_writer.get_path()==bucket_filepath);

    bucket_writer.flush();

    fs::remove(bucket_writer.get_path());

    BOOST_CHECK_THROW(BucketWriter<T>("/"), std::domain_error);

    BOOST_CHECK_THROW(BucketWriter<T>(get_a_temporary_path(), 0),
                      std::domain_error);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(load_bucket_T, T, test_types, BucketFixture<T>)
{
    using namespace RACES::Archive;

    BucketReader<T> bucket_reader(this->bucket_filepath);
    BOOST_CHECK(bucket_reader.size()==DEFAULT_DATASET_SIZE);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(sequential_bucket_T, T, test_types, BucketFixture<T>)
{
    using namespace RACES::Archive;

    BucketReader<T> bucket_reader(this->bucket_filepath, DEFAULT_READ_CACHE_SIZE);

    BOOST_CHECK(bucket_reader.size()==DEFAULT_DATASET_SIZE);

    size_t i=0;
    for (const auto& value: bucket_reader) {
        BOOST_CHECK(create_data<T>(i)==value);
        ++i;
    }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(random_io_bucket_T, T, test_types, BucketFixture<T>)
{
    using namespace RACES::Archive;

    BucketReader<T> bucket_reader(this->bucket_filepath, DEFAULT_READ_CACHE_SIZE);

    // create a list of indices and randomly shuffle them
    std::vector<size_t> indices(bucket_reader.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937_64 gen(0);

    std::shuffle(indices.begin(), indices.end(), gen);

    // check whether the values in the dataset correspond to
    // those in the same positions in the bucket
    for (const auto& index: indices) {
        BOOST_CHECK(create_data<T>(index)==bucket_reader[index]);
    }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(random_tour_T, T, test_types, BucketFixture<T>)
{
    using namespace RACES::Archive;

    std::mt19937_64 gen(0);

    BucketReader<T> bucket_reader(this->bucket_filepath, DEFAULT_READ_CACHE_SIZE);

    // testing different tours with different generators
    std::set<T> last_values;
    for (size_t i=0; i<5; ++i) {
        std::mt19937_64 gen(i);

        auto tour = bucket_reader.random_tour(gen);

        last_values.insert(test_random_tour_on(tour, this->dataset));
    }
    BOOST_CHECK(last_values.size()>1);

    // testing different tours with the same generator
    last_values = std::set<T>();
    for (size_t i=0; i<5; ++i) {
        auto tour = bucket_reader.random_tour(gen);

        last_values.insert(test_random_tour_on(tour, this->dataset));
    }
    BOOST_CHECK(last_values.size()==1);

    // testing multiple time the same generator
    last_values = std::set<T>();
    auto tour = bucket_reader.random_tour(gen);
    for (size_t i=0; i<5; ++i) {
        last_values.insert(test_random_tour_on(tour, this->dataset));
    }

    BOOST_CHECK(last_values.size()==1);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(copy_bucket_T, T, test_types, BucketFixture<T>)
{
    using namespace RACES::Archive;

    BucketReader<T> bucket_reader(this->bucket_filepath);
    BOOST_CHECK_NO_THROW(BucketReader<T> bucket_reader2(bucket_reader));

    try {
        BucketReader<T> bucket_reader2(bucket_reader);

        BOOST_CHECK(bucket_reader.size() == bucket_reader2.size());

        auto it = bucket_reader2.begin();
        for (const T& value: bucket_reader) {
            BOOST_CHECK(value == *it);
            ++it;
        }
    } catch (std::exception& ex)
    {}
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(shuffle_bucket_with_split_T, T, test_types, BucketFixture<T>)
{
    shuffle_bucket(this->bucket_filepath, this->dataset, DEFAULT_WRITE_CACHE_SIZE);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(shuffle_bucket_without_split_T, T, test_types, BucketFixture<T>)
{
    shuffle_bucket(this->bucket_filepath, this->dataset, 2*sizeof(T));
}
