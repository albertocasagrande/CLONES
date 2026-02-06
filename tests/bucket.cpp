/**
 * @file bucket.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Bucket tests
 * @version 1.3
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

typedef boost::mpl::list<size_t, CLONES::Mutations::GenomicPosition> test_types;

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

template<typename TYPE>
struct BucketFixture
{
    CLONES::Archive::Bucket<TYPE> bucket;
    std::set<TYPE> dataset;

    BucketFixture():
        bucket{get_a_temporary_path(), DEFAULT_WRITE_CACHE_SIZE}
    {
        for (size_t i=0; i<DEFAULT_DATASET_SIZE; ++i) {
            TYPE data = create_data<TYPE>(i);
            BOOST_CHECK_NO_THROW(bucket.push_back(data));
            dataset.emplace(std::move(data));
        }

        bucket.flush();
    }

    ~BucketFixture()
    {
        bucket.flush();

        std::filesystem::remove(bucket.path());
    }
};

template<typename TYPE, typename RANDOM_GENERATOR>
TYPE
test_random_tour_on(const CLONES::Archive::BucketRandomTour<TYPE, RANDOM_GENERATOR>& tour,
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
void shuffle_bucket(CLONES::Archive::Bucket<TYPE>& bucket,
                    const std::set<TYPE>& dataset,
                    const size_t read_cache_size)
{
    using namespace CLONES::Archive;
    namespace fs = std::filesystem;

    std::mt19937_64 gen(0);

    {
        BOOST_CHECK_NO_THROW(bucket.shuffle(gen, read_cache_size,
                                            fs::temp_directory_path()));
    }

    {
        Bucket<TYPE> bucket2(bucket.path(), DEFAULT_READ_CACHE_SIZE);
        std::set<TYPE> local_dataset(dataset.begin(), dataset.end());

        BOOST_CHECK(bucket2.size()==dataset.size());

        for (const TYPE& data: bucket2) {
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
    using namespace CLONES::Archive;
    namespace fs = std::filesystem;

    fs::path bucket_filepath = get_a_temporary_path();

    BOOST_CHECK(!fs::exists(bucket_filepath));

    Bucket<T> bucket(bucket_filepath);

    BOOST_CHECK(fs::exists(bucket.path()));

    bucket.flush();

    fs::remove(bucket.path());

    BOOST_CHECK_THROW(Bucket<T>("/"), std::domain_error);

    BOOST_CHECK_THROW(Bucket<T>(get_a_temporary_path(), 0),
                      std::domain_error);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(load_bucket_T, T, test_types, BucketFixture<T>)
{
    using namespace CLONES::Archive;
    
    Bucket<T> load_bucket(this->bucket.path());
    BOOST_CHECK(load_bucket.size()==DEFAULT_DATASET_SIZE);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(sequential_bucket_T, T, test_types, BucketFixture<T>)
{
    using namespace CLONES::Archive;

    Bucket<T> reading_bucket(this->bucket.path(), DEFAULT_READ_CACHE_SIZE);

    BOOST_CHECK(reading_bucket.size()==DEFAULT_DATASET_SIZE);

    size_t i=0;
    for (const auto& value: reading_bucket) {
        BOOST_CHECK(create_data<T>(i)==value);
        ++i;
    }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(random_io_bucket_T, T, test_types, BucketFixture<T>)
{
    using namespace CLONES::Archive;

    // create a list of indices and randomly shuffle them
    std::vector<size_t> indices(this->bucket.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937_64 gen(0);

    std::shuffle(indices.begin(), indices.end(), gen);

    // check whether the values in the dataset correspond to
    // those in the same positions in the bucket
    for (const auto& index: indices) {
        BOOST_CHECK(create_data<T>(index)==this->bucket[index]);
    }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(random_tour_T, T, test_types, BucketFixture<T>)
{
    using namespace CLONES::Archive;

    std::mt19937_64 gen(0);

    // testing different tours with different generators
    std::set<T> last_values;
    for (size_t i=0; i<5; ++i) {
        std::mt19937_64 gen(i);

        auto tour = this->bucket.random_tour(gen, DEFAULT_READ_CACHE_SIZE);

        last_values.insert(test_random_tour_on(tour, this->dataset));
    }
    BOOST_CHECK(last_values.size()>1);

    // testing different tours with the same generator
    last_values = std::set<T>();
    for (size_t i=0; i<5; ++i) {
        auto tour = this->bucket.random_tour(gen, DEFAULT_READ_CACHE_SIZE);

        last_values.insert(test_random_tour_on(tour, this->dataset));
    }
    BOOST_CHECK(last_values.size()==1);

    // testing multiple time the same generator
    last_values = std::set<T>();
    auto tour = this->bucket.random_tour(gen, DEFAULT_READ_CACHE_SIZE);
    for (size_t i=0; i<5; ++i) {
        last_values.insert(test_random_tour_on(tour, this->dataset));
    }

    BOOST_CHECK(last_values.size()==1);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(copy_bucket_T, T, test_types, BucketFixture<T>)
{
    using namespace CLONES::Archive;

    BOOST_CHECK_NO_THROW(Bucket<T> bucket2(this->bucket));

    try {
        Bucket<T> bucket2(this->bucket);

        BOOST_CHECK(this->bucket.size() == bucket2.size());
        BOOST_CHECK(this->bucket.get_cache_size() == bucket2.get_cache_size());

        auto it = bucket2.begin();
        for (const T& value: this->bucket) {
            BOOST_CHECK(value == *it);
            ++it;
        }
    } catch (std::exception& ex)
    {}
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(shuffle_bucket_with_split_T, T, test_types, BucketFixture<T>)
{
    shuffle_bucket(this->bucket, this->dataset, DEFAULT_WRITE_CACHE_SIZE);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(shuffle_bucket_without_split_T, T, test_types, BucketFixture<T>)
{
    shuffle_bucket(this->bucket, this->dataset,
                   this->bucket.size()*sizeof(T));
}
