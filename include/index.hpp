/**
 * @file index.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines index
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

#ifndef __RACES_INDEX__
#define __RACES_INDEX__

#include <map>
#include <memory>
#include <type_traits>

#include "bucket.hpp"
#include "progress_bar.hpp"

namespace RACES
{

/**
 * @brief A structure to define the partitions of the values of a type
 *
 * A partition for a type is a set of non-overlapping sets
 * of values, named classes, whose union consists in all values
 * representable by the type.
 *
 * The default is partition of a type consists of singleton
 * classes.
 *
 * @tparam TYPE is the type for which this structure defines
 *      the partition
 */
template<typename TYPE>
struct partition
{
    /**
     * @brief Return a list of the values in a class
     *
     * @param[in] value is a representant of the class that is
     *      aimed
     * @return a list of the values in the class including
     *      `value`
     */
    static std::list<TYPE> get_class_of(const TYPE& value)
    {
        return std::list<TYPE>{value};
    }
};

namespace Archive
{

/**
 * @brief A base class for `IndexBuilder` and `IndexReader`
 *
 * An index is a data structure that associates a set of keys
 * to a set of values. Each key corresponds to a set of values
 * that are maintained in a bucket.
 * Any index is stored in a directory that contains the set
 * of key buckets and a map file containing the map key-bucket.
 * This class is a base class for both `IndexBuilder` and
 * `IndexReader`.
 *
 * @tparam KEY is a type of the index keys
 */
template<class KEY>
class IndexBase
{
    std::filesystem::path dir_path; //!< the path of the index directory

    size_t cache_size;      //!< the cache size in bytes

    std::string bucket_prefix;  //!< the bucket file prefixes

protected:
    /**
     * @brief Set the prefix of the index bucket filenames
     *
     * @param[in] bucket_prefix is the new prefix of bucket filename
     */
    inline void set_bucket_prefix(const std::string& bucket_prefix)
    {
        this->bucket_prefix = bucket_prefix;
    }

    /**
     * @brief Get the file format description string
     *
     * @return the file format description string
     */
    static constexpr std::string index_file_descr()
    {
        return "RACES index";
    }

    /**
     * @brief Get the file format version
     *
     * @return the file format version
     */
    static constexpr uint8_t index_file_version()
    {
        return 0;
    }
public:
    /**
     * @brief The empty constructor
     */
    IndexBase():
        cache_size{0}
    {}

    /**
     * @brief A constructor
     *
     * @param[in] index_path is the path of the index directory
     * @param[in] cache_size is the index cache size in bytes
     * @param[in] bucket_prefix is the prefix of the index bucket
     *      filenames
     */
    IndexBase(const std::filesystem::path index_path,
               size_t cache_size,
               const std::string bucket_prefix="bucket"):
        dir_path{index_path}, cache_size{cache_size},
        bucket_prefix{bucket_prefix}
    {
        if (cache_size==0) {
            throw std::domain_error("Index cache size must be greater than 0.");
        }
    }

    /**
     * @brief Get the path of a bucket file
     *
     * This method returns the path of the bucket file associated to
     * a key.
     *
     * @param[in] key is the key associated to the aimed bucket file path
     * @return the path of the bucket file associated to `key`
     */
    std::filesystem::path get_bucket_path(const KEY& key) const
    {
        std::ostringstream oss;

        oss << bucket_prefix << "_" << key << ".bin";

        return dir_path / oss.str();
    }

    /**
     * @brief Get the prefix of the bucket filenames
     *
     * @return the prefix of the bucket filenames
     */
    inline const std::string& get_bucket_prefix() const
    {
        return bucket_prefix;
    }

    /**
     * @brief Get the filename of the index map file
     *
     * Any index is saved in a directory as a file mapping
     * any key to a bucket, named *map file*, and the
     * corresponding buckets.
     * This method returns the filename of the map file.
     *
     * @return the filename of the index map file
     */
    inline static constexpr std::string map_filename()
    {
        return "map.bin";
    }

    /**
     * @brief Get the index map file path
     *
     * Any index is saved in a directory as a file mapping
     * any key to a bucket, named *map file*, and the
     * corresponding buckets.
     * This method returns the filename of the map file
     * path.
     *
     * @return the index map file path
     */
    inline std::filesystem::path get_map_path() const
    {
        return dir_path / map_filename();
    }

    /**
     * @brief Get the index directory path
     *
     * @return the index directory path
     */
    inline const std::filesystem::path& get_path() const
    {
        return dir_path;
    }

    /**
     * @brief Get the index cache size
     *
     * @return the index cache size in bytes
     */
    inline const size_t& get_cache_size() const
    {
        return cache_size;
    }

    /**
     * @brief Get the cache size per index bucket
     *
     * @param[in] num_of_keys is the number of index keys
     * @return the cache size per index bucket in bytes if the index
     *      has `num_of_keys` keys
     */
    inline size_t cache_size_per_bucket(const size_t num_of_keys) const
    {
        return cache_size/num_of_keys;
    }
};

/**
 * @brief The index writer base class
 *
 * An index is a data structure that associates a set of keys
 * to a set of values. Each key corresponds to a set of values
 * that are maintained in a bucket.
 * Any index is stored in a directory that contains the set
 * of key buckets and a map file containing the map key-bucket.
 * This class implements index writers. Its objects are used
 * to write data into indices.
 *
 * @tparam KEY the type of the index keys
 * @tparam VALUE the type of values
 * @tparam COMPARE the structure comparing KEY's values
 */
template<class KEY, class VALUE, class COMPARE = std::less<KEY>>
class IndexBuilder : private IndexBase<KEY>
{
public:
    /**
     * @brief The value type
     */
    using value_type = VALUE;

    /**
     * @brief The type of bucket writers
     *
     * A bucket is a container for values associated to the
     * same key.
     */
    using BucketWriterType = BucketWriter<VALUE>;

    /**
     * @brief The bucket map type
     */
    using BucketMapType = std::map<KEY, BucketWriterType>;

private:
    BucketMapType buckets;  //!< the key-bucket map

    /**
     * @brief Associated a new bucket to a key
     *
     * This method associated a new bucket to a key in the
     * key-bucket map. If the corresponding bucket file
     * already exists or new association cannot be added
     * to the index map, a `std::runtime_error` is thrown.
     *
     * @param[in] key is the key to which the new bucket is
     *      associated
     * @return an iterator of the key-bucket map referring
     *      to the just inserted pair `key`-`key`'s bucket
     */
    BucketMapType::iterator add_bucket_for(const KEY& key)
    {
        auto bucket_path = this->get_bucket_path(key);

        if (std::filesystem::exists(bucket_path)) {
            std::ostringstream oss;

            oss << "IndexBuilder::index: the bucket file "
                << "for \"" << key << "\" already exists.";
            throw std::runtime_error(oss.str());
        }

        const auto bucket_cache_size = this->cache_size_per_bucket(buckets.size()+1);

        for (auto& [key, bucket]: buckets) {
            bucket.set_cache_size(bucket_cache_size);
        }

        const auto result = buckets.emplace(key, BucketWriterType{bucket_path,
                                                                  bucket_cache_size});

        if (result.second) {
            return result.first;
        }

        std::ostringstream oss;

        oss << "The bucket for \"" << key << "\" cannot be added to the index.";
        throw std::runtime_error(oss.str());
    }

public:
    /**
     * @brief Construct a new index builder
     *
     * @param[in] index_path is the path to the directory storing the index
     * @param[in] cache_size is the write cache size in bytes
     * @param[in] bucket_prefix is the prefix of the bucket filenames
     */
    explicit IndexBuilder(const std::filesystem::path index_path,
                          const size_t cache_size = 10000000,
                          const std::string bucket_prefix="bucket"):
        IndexBase<KEY>{index_path, cache_size, bucket_prefix}
    {
        if (cache_size==0) {
            throw std::domain_error("IndexBuilder: The cache size must "
                                    "be greater than 0.");
        }

        if (std::filesystem::exists(index_path)) {
            std::ostringstream oss;

            oss << "IndexBuilder: \"" << to_string(index_path)
                << "\" already exists.";
            throw std::domain_error(oss.str());
        }

        std::filesystem::create_directory(index_path);
    }

    /**
     * @brief Insert a value in the index
     *
     * @param[in] key is the key of the value
     * @param[in] value is the value to be inserted
     */
    void insert(const KEY& key, const VALUE& value)
    {
        auto found = buckets.find(key);
        if (found == buckets.end()) {
            found = add_bucket_for(key);
        }

        found->second.push_back(value);
    }

    /**
     * @brief Insert a value in the index
     *
     * @param[in,out] key is the key of the value
     * @param[in,out] value is the value to be inserted
     */
    void insert(KEY&& key, VALUE&& value)
    {
        auto found = buckets.find(key);
        if (found == buckets.end()) {
            found = add_bucket_for(key);
        }

        found->second.push_back(std::move(value));
    }

    /**
     * @brief Insert a value in the index
     *
     * @param[in] key is the key of the value
     * @param[in,out] value is the value to be inserted
     */
    void insert(const KEY& key, VALUE&& value)
    {
        auto found = buckets.find(key);
        if (found == buckets.end()) {
            found = add_bucket_for(key);
        }

        found->second.push_back(std::move(value));
    }

    /**
     * @brief Shuffle the values in the index buckets
     *
     * @tparam ARGS are the types of the parameters
     * @param[in, out] args are the parameters
     */
    template <typename... ARGS>
    requires (!std::same_as<std::remove_cvref_t<
                typename std::tuple_element<sizeof...(ARGS) - 1, std::tuple<ARGS...>>::type
                >, UI::ProgressBar>)
    inline void shuffle(ARGS... args)
    {
        UI::ProgressBar progress_bar;

        shuffle(args..., progress_bar);
    }

    /**
     * @brief Shuffle the values in the index buckets
     *
     * @tparam RANDOM_GENERATOR is a random number generator type
     * @tparam ARGS are the types of the remaining parameters
     * @param[in] seed is the seed of a random number generator
     * @param[in, out] args are the remaining parameters
     */
    template<typename RANDOM_GENERATOR, typename... ARGS>
    inline void shuffle(const int seed, ARGS... args)
    {
        RANDOM_GENERATOR generator(seed);

        shuffle(generator, args...);
    }

    /**
     * @brief Shuffle the values in the index buckets
     *
     * @tparam RANDOM_GENERATOR is a random number generator type
     * @param[in, out] random_generator is a random number generator
     * @param[in, out] progress_bar is a progress bar
     */
    template<typename RANDOM_GENERATOR>
    inline void shuffle(RANDOM_GENERATOR& random_generator,
                        UI::ProgressBar& progress_bar)
    {
        shuffle(random_generator, std::filesystem::temp_directory_path(),
                progress_bar);
    }

    /**
     * @brief Shuffle the values in the index buckets
     *
     * @tparam RANDOM_GENERATOR is a random number generator type
     * @param[in, out] random_generator is a random number generator
     * @param[in] tmp_dir is the path of the temporary files
     * @param[in, out] progress_bar is a progress bar
     */
    template<typename RANDOM_GENERATOR>
    void shuffle(RANDOM_GENERATOR& random_generator,
                 const std::filesystem::path tmp_dir,
                 UI::ProgressBar& progress_bar)
    {
        progress_bar.set_message("Shuffling index");

        // flush the cache of all the buckets and free as much memory as possible
        for (auto& [key, bucket]: buckets) {
            progress_bar.update_elapsed_time();

            bucket.flush();
        }

        size_t progress{0};
        // use all the cache and sequentially shuffle each bucket
        for (auto& [key, bucket]: buckets) {
            bucket.shuffle(random_generator, this->get_cache_size(), tmp_dir,
                           progress_bar);

            progress_bar.set_progress(100*(++progress)
                                        /(buckets.size()+1));
        }

        progress_bar.set_progress(100, "Index shuffled");
    }

    /**
     * @brief Save the map file on disk
     *
     * This method saves the map file in the index directory
     */
    void save_map_on_disk()
    {
        for (auto& [key, bucket]: buckets) {
            bucket.flush();
        }

        Binary::Out archive(this->get_map_path());

        Binary::Out::write_header(archive, IndexBase<KEY>::index_file_descr(),
                                  IndexBase<KEY>::index_file_version());

        archive & this->get_bucket_prefix()
                & buckets.size();

        for (const auto& [key, bucket]: buckets) {
            archive & key;

            (void)bucket;
        }
    }

    /**
     * @brief Destroy the Index Builder
     *
     * This method saves the index map file in the index directory
     * and destroys the index builder object.
     */
    ~IndexBuilder()
    {
        save_map_on_disk();
    }
};

/**
 * @brief The index reader class
 *
 * An index is a data structure that associates a set of keys
 * to a set of values. Each key corresponds to a set of values
 * that are maintained in a bucket.
 * Any index is stored in a directory that contains the set
 * of key buckets and a map file containing the map key-bucket.
 * The objects of this class can read an index.
 *
 * @tparam KEY the type of the index keys
 * @tparam VALUE the type of values
 * @tparam RANDOM_GENERATOR the type of random number generator
 * @tparam COMPARE the structure comparing KEY's values
 */
template<class KEY, class VALUE, class RANDOM_GENERATOR,
         class COMPARE = std::less<KEY>>
class IndexReader : private IndexBase<KEY>
{
public:
    /**
     * @brief The key type
     */
    using key_type = KEY;

    /**
     * @brief The value type
     */
    using value_type = VALUE;

    /**
     * @brief The random number generator type
     */
    using random_generator_type = RANDOM_GENERATOR;

    /**
     * @brief The type of buckets
     *
     * A bucket is a container for values associated to the
     * same key.
     */
    using BucketType = Archive::BucketReader<VALUE>;

    /**
     * @brief The type of buckets
     *
     * A bucket is a container for values associated to the
     * same key.
     */
    using RandomTourType = Archive::BucketRandomTour<VALUE, RANDOM_GENERATOR>;

    /**
     * @brief The bucket map type
     */
    using BucketMapType = std::map<KEY, BucketType>;

    /**
     * @brief The bucket iterator map
     */
    using TourIteratorType = std::map<KEY, typename RandomTourType::const_iterator>;

private:
    BucketMapType buckets;      //!< the key-bucket map
    TourIteratorType bucket_iterators;     //!< the key-bucket iterator map

    /**
     * @brief Add a tour iterator to the key-bucket iterator map
     *
     * This method associates a key with a random tour iterator over the bucket
     * that corresponds to the key in the key-bucket iterator map. If the new
     * association cannot be added to the map, a `std::runtime_error` is thrown.
     *
     * @param[in] key is the key associated to the bucket whose random tour iterator
     *      must be added to the key-bucket iterator map
     * @return an iterator to the added association in the key-bucket iterator map
     */
    TourIteratorType::iterator add_tour_iterator_for(RANDOM_GENERATOR& generator,
                                                     const KEY& key)
    {
        const auto& bucket = operator[](key);

        const size_t bucket_cache_size = this->cache_size_per_bucket(buckets.size());
        const auto tour = bucket.random_tour(generator, bucket_cache_size);
        auto result = bucket_iterators.emplace(key, tour.begin());

        if (!result.second) {
            std::ostringstream oss;

            oss << "The iterator for " << key
                << "'s bucket cannot be initialised.";

            throw std::runtime_error(oss.str());
        }

        return result.first;
    }

public:
    /**
     * @brief The empty index reader constructor
     */
    IndexReader():
        IndexBase<KEY>{}, buckets{}, bucket_iterators{}
    {}

    /**
     * @brief A constructor
     *
     * @param[in] index_path is the path of the index directory
     * @param[in] cache_size is the read cache size in bytes
     */
    IndexReader(const std::filesystem::path index_path,
                const size_t cache_size = 10000000):
        IndexBase<KEY>{index_path, cache_size},
        buckets{}, bucket_iterators{}
    {
        if (!std::filesystem::exists(index_path)) {
            std::ostringstream oss;

            oss << "IndexReader: \"" << to_string(index_path)
                << "\" does not exist.";
            throw std::domain_error(oss.str());
        }

        if (!std::filesystem::is_directory(index_path)) {
            std::ostringstream oss;

            oss << "IndexReader: \"" << to_string(index_path)
                << "\" is not a directory.";
            throw std::domain_error(oss.str());
        }

        const auto map_path = this->get_map_path();
        if (!std::filesystem::exists(map_path)) {
            std::ostringstream oss;

            oss << "IndexReader: \"" << to_string(map_path)
                << "\" does not exist.";
            throw std::domain_error(oss.str());
        }

        Binary::In archive(map_path);

        Binary::In::read_header(archive, IndexBase<KEY>::index_file_descr(),
                                IndexBase<KEY>::index_file_version());

        {
            std::string buffer_prefix;

            archive & buffer_prefix;

            this->set_bucket_prefix(buffer_prefix);
        }

        size_t num_of_map_items;
        archive & num_of_map_items;

        const auto bucket_cache_size = this->cache_size_per_bucket(num_of_map_items);
        for (size_t i=0; i<num_of_map_items; ++i) {
            KEY key;

            archive & key;

            const auto bucket_path = this->get_bucket_path(key);

            buckets.emplace(key, BucketType{bucket_path,  bucket_cache_size});
        }
    }

    /**
     * @brief Get the index keys
     *
     * @return a vector containing the index keys
     */
    std::vector<KEY> get_keys() const
    {
        std::vector<KEY> keys;
        keys.reserve(buckets.size());

        for (const auto& [key, bucket] : buckets) {
            keys.push_back(key);
        }

        return keys;
    }

    /**
     * @brief The key-bucket operator
     *
     * This method returns the bucket of values associated to a key.
     *
     * @param[in] key is the searched bucket's key
     * @return a constant reference to the bucket of values
     *    associated to `key`.
     */
    inline const BucketType& operator[](const KEY& key) const
    {
        return buckets.at(key);
    }

    /**
     * @brief Count the extractable values for a key
     *
     * @param[in] key is a key
     * @return the number of extractable values associated to `key`
     */
    size_t extractable_for(const KEY& key) const
    {
        auto found_it = bucket_iterators.find(key);
        if (found_it != bucket_iterators.end()) {
            return found_it->second.remaining_values();
        }

        auto found = buckets.find(key);
        if (found != buckets.end()) {
            return found->second.size();
        }

        return 0;
    }

    /**
     * @brief Count the values associated to a key
     *
     * @param[in] key is the searched key
     * @return the number of values associated to `key`
     */
    size_t num_of_values(const KEY& key) const
    {
        auto found_it = buckets.find(key);
        if (found_it != buckets.end()) {
            return found_it->second.size();
        }

        return 0;
    }

    /**
     * @brief Return the number of keys in the index
     *
     * @return the number of keys in the index
     */
    inline size_t num_of_keys() const
    {
        return buckets.size();
    }

    /**
     * @brief Initialise bucket tours and caches some data
     *
     * This method boots up bucket tours randomizing the tour's
     * initial position and loading their first chunk.
     *
     * @param[in,out] generator is a random number generator
     * @param[in,out] progress_bar is a progress bar
     */
    void boot_up_tours(RANDOM_GENERATOR& generator,
                       UI::ProgressBar& progress_bar)
    {
        progress_bar.set_message("Booting context index");

        size_t booted_up{0};
        for (const auto& [key, bucket]: buckets) {
            auto found_it = bucket_iterators.find(key);

            if (found_it == bucket_iterators.end()) {
                found_it = add_tour_iterator_for(generator, key);
            }

            progress_bar.set_progress((100*(++booted_up))/buckets.size());
        }
    }

    /**
     * @brief Initialise bucket tours and caches some data
     *
     * This method boots up bucket tours randomizing the tour's
     * initial position and loading their first chunk.
     *
     * @tparam ARGS are the types of the remaining parameters
     * @param[in] seed is the seed of a random number generator
     * @param[in, out] args are the remaining parameters
     */
    template<typename... ARGS>
    inline void boot_up_tours(const int seed, ARGS... args)
    {
        RANDOM_GENERATOR generator(seed);

        boot_up_tours(generator, args...);
    }

    /**
     * @brief Initialise bucket tours and caches some data
     *
     * This method boots up bucket tours randomizing the tour's
     * initial position and loading their first chunk.
     *
     * @tparam ARGS are the types of the parameters
     * @param[in, out] args are the parameters
     */
    template <typename... ARGS>
    requires (!std::same_as<std::remove_cvref_t<
                typename std::tuple_element<sizeof...(ARGS) - 1, std::tuple<ARGS...>>::type
                >, UI::ProgressBar>)
    inline void boot_up_tours(ARGS... args)
    {
        UI::ProgressBar progress_bar;

        boot_up_tours(args..., progress_bar);
    }

    /**
     * @brief Extract a random value from the index
     *
     * This method extracts a random value associated to a key.
     * After the extraction the value is still present in the
     * index, but it is no more available for the extraction up
     * to the index reset.
     * If all values associated to the specified key have been
     * already extracted, a `std::runtime_error` is thrown.
     *
     * @param[in,out] generator is a random number generator
     * @param[in] key is a key
     * @return a random value associated to `key`
     */
    VALUE extract(RANDOM_GENERATOR& generator, const KEY& key)
    {
        auto found_it = bucket_iterators.find(key);

        if (found_it == bucket_iterators.end()) {
            found_it = add_tour_iterator_for(generator, key);
        }

        if (found_it->second.is_end()) {
            std::ostringstream oss;

            oss << "No more values available for " << key << ".";

            throw std::runtime_error(oss.str());
        }

        VALUE value{*(found_it->second)};
        ++(found_it->second);

        return value;
    }

    /**
     * @brief Choose a random value associated to a key
     *
     * This method chooses a random value associated to a key with
     * uniform distribution. The choosen value is *not* removed from
     * the index and can be choosen another time.
     *
     * @param[in,out] generator is the random number generator
     * @param[in] key is a key
     * @return a random value associated to `key`
     */
    template<typename VALUE2=VALUE,
                std::enable_if_t<std::is_same_v<VALUE2, VALUE>
                                 && uses_constant_space_on_disk<VALUE>::value, bool> = true>
    inline VALUE choose(RANDOM_GENERATOR& generator, const KEY& key) const
    {
        return buckets.at(key).choose(generator);
    }

    /**
     * @brief Count the values in associated to a key class
     *
     * @param[in] key is a representant of the class whose
     *      corresponding values are counted
     * @return the number of values associated to the class
     *      of `key`
     */
    template<typename KEY_PARTITION=partition<KEY>>
    size_t num_of_class_values(const KEY& key) const
    {
        size_t counter{0};
        for (const auto key : KEY_PARTITION::get_class_of(key)) {
            auto found_it = buckets.find(key);

            if (found_it != buckets.end()) {
                counter += found_it->second.size();
            }
        }

        return counter;
    }

    /**
     * @brief Count the extractable values from a key class
     *
     * @tparam KEY_PARTITION is a partition of key type
     * @param[in] key is the a representant of the searched key class
     * @return the number of extractable values associated to the
     *      class of `key`
     */
    template<typename KEY_PARTITION=partition<KEY>>
    size_t extractable_from_class(const KEY& key) const
    {
        size_t counter{0};
        for (const auto& class_key: KEY_PARTITION::get_class_of(key)) {
            counter += extractable_for(class_key);
        }

        return counter;
    }

    /**
     * @brief Extract a random value associated to a key class
     *
     * This method extracts a random value associated to a key class.
     * After the extraction the value is still present in the index,
     * but it is no more available for the extraction up to the
     * index reset.
     * If all values associated to the specified key class have been
     * already extracted, a `std::runtime_error` is thrown.
     *
     * @tparam KEY_PARTITION is a partition of key type
     * @param[in,out] generator is the random number generator
     * @param[in] key is a representative of the key class from which
     *      a value must be extracted
     * @return a pair whose first element is one of the keys in the
     *      class of `key`. The second element is one of the
     *      extractable values associated to the first element
     */
    template<typename KEY_PARTITION=partition<KEY>>
    std::pair<KEY, VALUE>
    extract_from_class(RANDOM_GENERATOR& generator, const KEY& key)
    {
        const auto key_class = KEY_PARTITION::get_class_of(key);

        size_t available_in_class{0};
        for (const auto& class_key: key_class) {
            available_in_class += extractable_for(class_key);
        }

        if (available_in_class==0) {
            throw std::runtime_error("No value available in the key class.");
        }

        std::uniform_int_distribution<size_t> dist(0, available_in_class-1);

        size_t pos = dist(generator);

        for (const auto& class_key: key_class) {
            const size_t available_for = extractable_for(class_key);

            if (pos < available_for) {
                return {class_key, extract(generator, class_key)};
            }

            pos -= available_for;
        }

        return {key, extract(generator, key)};
    }

    /**
     * @brief Choose a random value associated to a key class
     *
     * This method chooses a random value associated to a key class
     * with uniform distribution. The choosen value is *not* removed
     * from the index and can be choosen another time.
     *
     * @tparam KEY_PARTITION is a partition of key type
     * @param[in,out] generator is the random number generator
     * @param[in] key is a key
     * @return a pair whose first element is one of the keys in the
     *      class of `key`. The second element is one of the
     *      values associated to the first element
     */
    template<typename KEY_PARTITION=partition<KEY>, typename VALUE2=VALUE,
                std::enable_if_t<uses_constant_space_on_disk<VALUE>::value, bool> = true>
    std::pair<KEY, VALUE>
    choose_from_class(RANDOM_GENERATOR& generator, const KEY& key) const
    {
        const auto key_class = KEY_PARTITION::get_class_of(key);

        size_t available_in_class{0};
        for (const auto& class_key: key_class) {
            available_in_class += extractable_for(class_key);
        }

        if (available_in_class==0) {
            throw std::runtime_error("No value available in the key class.");
        }

        std::uniform_int_distribution<size_t> dist(0, available_in_class-1);

        size_t pos = dist(generator);

        for (const auto& class_key: key_class) {
            const size_t available_for = extractable_for(class_key);
            if (pos < available_for) {
                return {class_key, choose(generator, class_key)};
            }

            pos -= available_for;
        }

        return {key, choose(generator, key)};
    }

    /**
     * @brief Reset the index
     *
     * This method resets the index and labels all values as
     * extractable.
     */
    void reset()
    {
        bucket_iterators.clear();
    }
};

}   // Archive

}   // RACES

#endif // __RACES_INDEX__
