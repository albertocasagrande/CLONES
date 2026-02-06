/**
 * @file bucket.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines bucket
 * @version 1.4
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

#ifndef __CLONES_BUCKET__
#define __CLONES_BUCKET__

#include <vector>
#include <random>
#include <memory>
#include <filesystem>
#include <numeric> // iota

#include <algorithm> // shuffle

#include "archive.hpp"
#include "utils.hpp"

namespace CLONES
{

namespace Archive
{

template<typename VALUE, typename RANDOM_GENERATOR>
class BucketRandomTour;

/**
 * @brief A class for buckets
 *
 * A bucket is a container that stores a collection of values in
 * a file. These values are maintained in a specific order and can
 * be read sequentially from first to last. The order of values
 * within the bucket can also be shuffled if needed.
 *
 * @tparam VALUE is the type of the values contained in the bucket
 */
template<typename VALUE>
class Bucket
{
    std::filesystem::path filepath;   //!< the name of the file storing the bucket

    std::streampos size_pos;    //!< the position of the bucket size in the file
    std::streampos data_pos;    //!< the position of the first value in the file
    std::streampos file_size;   //!< the last position in the file

    size_t num_of_values;   //!< the number of values in the bucket (i.e., bucket size)

    size_t cacheable_values;    //!< the number of the cacheable values
    std::list<VALUE> write_cache; //!< the write cache

    /**
     * @brief Load a set of values from a file
     *
     * This static method loads a set of values from a file and fills
     * a buffer. The buffer must be large enough to store all values
     * in the file. If this is not the case, a `std::runtime_error`
     * is thrown.
     *
     * @param[in, out] buffer is the buffer to be filled
     * @param[in] filepath is the path of the file containing the values
     * @return an iterator referring to the position next to the last
     *      read value in the buffer.
     */
    static typename std::vector<VALUE>::iterator
    load_buffer(std::vector<VALUE>& buffer, const std::filesystem::path& filepath)
    {
        Archive::Binary::In archive(filepath);

        auto buffer_it = buffer.begin();
        while (!archive.eof() && buffer_it != buffer.end()) {
            archive & *buffer_it;

            ++buffer_it;
        }

        if (!archive.eof()) {
            throw std::runtime_error("load_buffer: The file is larger than the buffer");
        }

        return buffer_it;
    }

    /**
     * @brief Create a set of chunks
     *
     * A chunk is a file containing part of the values in the bucket.
     * This method creates a set of chunks and associates the archives
     * in `chucks` to them.
     *
     * @param[in, out] chunks is the vector of the created chucks
     * @param[in] prefix_name is the prefix name of the chunks
     * @param[in] tmp_dir is the directory of the chucks
     * @return the vector of the chunk paths
     */
    static std::vector<std::filesystem::path>
    create_chunks(std::vector<Archive::Binary::Out>& chunks,
                  const std::string& prefix_name,
                  const std::filesystem::path tmp_dir)
    {
        std::vector<std::filesystem::path> chunk_paths(chunks.size());

        size_t i{0}, name_num{0};
        for (auto& chunk_path: chunk_paths) {
            do {
                std::ostringstream oss;

                oss << prefix_name << ++name_num << ".tmp";

                chunk_path = tmp_dir / oss.str();
            } while (std::filesystem::exists(chunk_path));

            chunks[i].open(chunk_path);
            ++i;
        }

        return chunk_paths;
    }

    /**
     * @brief Distribute the bucket values into chunks
     *
     * This method distributes the bucket values into chunks having a
     * maximum specified size. A chunk is a file containing part of the
     * values in the bucket. The number of chunks depends on both bucket
     * and maximum chucks sizes. When a chunk reaches the maximum size,
     * it is full. Each value is randomly saved into one of the non-full
     * chunks according a uniform distribution over them.
     *
     * @tparam RANDOM_GENERATOR is the type of the random number generator
     * @param[in, out] random_generator is a random number generator
     * @param[in] max_chunk_size is the maximum size of chunks
     * @param[in] prefix_name is the prefix name of the chunks
     * @param[in] tmp_dir is the directory of the chucks
     * @return the vector of the chunk paths
     */
    template<typename RANDOM_GENERATOR>
    std::vector<std::filesystem::path>
    split_in_random_chunks(RANDOM_GENERATOR& random_generator,
                           const size_t& max_chunk_size,
                           const std::string& prefix_name = "tmp_chunk",
                           const std::filesystem::path& tmp_dir = "./")
    {
        const auto num_of_chunks = (size()-1)/max_chunk_size+1;
        auto last_chunk = num_of_chunks-1;

        std::uniform_int_distribution<size_t> chunk_dist(0, last_chunk);

        std::vector<Archive::Binary::Out> chunks(num_of_chunks);
        std::vector<size_t> chunk_sizes(num_of_chunks, 0);

        const auto chunk_paths = create_chunks(chunks, prefix_name, tmp_dir);
        std::vector<size_t> positions(num_of_chunks);
        std::iota(positions.begin(), positions.end(), 0);

        std::vector<VALUE> cache{max_chunk_size};
        std::streampos read_pos{0};

        size_t value_in_caches = load_buffer(cache, read_pos);
        auto cache_it = cache.begin();
        for (size_t i=0; i<num_of_values; ++i) {
            if (value_in_caches == 0) {
                value_in_caches = load_buffer(cache, read_pos);

                cache_it = cache.begin();
            }

            const auto pos = chunk_dist(random_generator);
            const auto index = positions[pos];

            ++(chunk_sizes[index]);
            chunks[index] & *cache_it;
            ++cache_it;
            --value_in_caches;

            if (chunk_sizes[index] == max_chunk_size) {
                if (pos != last_chunk) {
                    std::swap(positions[pos], positions[last_chunk]);
                }
                --last_chunk;

                chunk_dist = std::uniform_int_distribution<size_t>(0, last_chunk);
            }
        }

        return chunk_paths;
    }

    /**
     * @brief Load values into a buffer
     *
     * This method loads a set of values from a specified position of the bucket
     * file into a buffer and returns the number of values loaded into the buffer.
     * This value corresponds to the minimum among the buffer size and the number
     * of values in the buffer file from the specified file position.
     *
     * @param[in, out] buffer is the buffer that will store the read values
     * @param[in, out] read_pos is the position in the bucket file from which values
     *      are load
     * @return the number of values read from the bucket file
     */
    size_t load_buffer(std::vector<VALUE>& buffer, std::streampos& read_pos) const
    {
        Archive::Binary::In archive(filepath);

        const std::streampos final_pos{archive.size()};

        if (read_pos < data_pos) {
            read_pos = data_pos;
        } else {
            if (archive.size()==read_pos) {
                return 0;
            }
        }
        archive.seekg(read_pos);

        size_t read_values{0};
        for (auto& value: buffer) {
            if (final_pos==read_pos) {
                return read_values;
            }
            archive & value;

            read_pos = archive.tellg();
            ++read_values;
        }

        return read_values;
    }

    /**
     * @brief Load values into a buffer
     *
     * This method loads a set of values from a specified position of the bucket
     * file into a buffer and returns the number of values loaded into the buffer.
     * When the final position of the bucket file is reached the read position is
     * updated to read the first value in the file. The method stops to read
     * values from the bucket file when either the buffer has been filled or the
     * a specified final position has been reached. A Boolean flag allows the
     * method to proceed the first time the final position is reached.
     *
     * @param[in, out] buffer is the buffer that will store the read values
     * @param[in, out] read_pos is the position in the bucket file from which values
     *      are load
     * @param[in] final_pos is the reading final position
     * @param[in] init is a Boolean flag that must be set to `true` during the
     *      first buffer load
     * @return the number of values read from the bucket file
     */
    size_t load_buffer(std::vector<VALUE>& buffer, std::streampos& read_pos,
                       std::streampos final_pos, bool init=false) const
    {
        Archive::Binary::In archive(this->path());

        if (read_pos < this->get_data_pos()) {
            read_pos = this->get_data_pos();
        }
        if (final_pos < this->get_data_pos()) {
            final_pos = read_pos;
        }
        archive.seekg(read_pos);

        size_t read_values{0};
        for (auto& value: buffer) {
            if (archive.eof()) {
                read_pos = this->get_data_pos();
                archive.seekg(read_pos);
            }
            if (final_pos==read_pos) {
                if (!init) {
                    return read_values;
                }
                init = false;
            }
            archive & value;

            read_pos = archive.tellg();
            ++read_values;
        }

        return read_values;
    }

    /**
     * @brief The method that initializes Bucket objects
     *
     * This method is used to initialize Bucket objects during the construction
     * and the copy.
     */
    void init_bucket()
    {
        if (std::filesystem::exists(filepath)) {
            if (!std::filesystem::is_regular_file(filepath)) {
                std::ostringstream oss;

                oss << "\"" << to_string(filepath) << "\" is not a block file.";
                throw std::domain_error(oss.str());
            }

            Archive::Binary::In archive(filepath);

            Archive::Binary::In::read_header(archive, "CLONES Bucket", 0);

            size_pos = archive.tellg();
            archive & num_of_values;
            data_pos = archive.tellg();

            file_size = archive.size();
        } else {
            Archive::Binary::Out archive(filepath);

            Archive::Binary::Out::write_header(archive, "CLONES Bucket", 0);

            size_pos = archive.tellg();
            archive & num_of_values;
            data_pos = archive.tellg();

            archive.flush();
            file_size = data_pos;
        }
    }
protected:
    /**
     * @brief Get the position of values in bucket file
     *
     * @return the position of values in bucket file
     */
    inline std::streampos get_data_pos() const
    {
        return data_pos;
    }

    /**
     * @brief Compute the position of the i-th value in the bucket file
     *
     * @param[in] i is the index of the aimed value in the bucket file
     * @return the position of the `i`-th value in the bucket file or,
     *      whenever the bucket file contains less than i values, the
     *      end file position
     */
    template<typename VALUE2=VALUE, std::enable_if_t<std::is_same_v<VALUE2, VALUE>
                && uses_constant_space_on_disk<VALUE>::value, bool> = true>
    std::streampos get_value_pos(const size_t& i) const
    {
        const size_t values_in_disks = num_of_values-write_cache.size();

        if (i<values_in_disks) {
            const std::streamoff offset = i*(file_size-data_pos)/values_in_disks;
            return data_pos + offset;
        }

        return file_size;
    }

public:
    /**
     * @brief The type of stored values
     */
    using value_type = VALUE;

    /**
     * @brief A constant iterator for the values in the bucket
     */
    class const_iterator
    {
        Bucket<VALUE> const* bucket;    //!< a pointer to the bucket

        std::vector<VALUE> cache;   //!< the read cache
        std::streampos read_pos;    //!< the next position to be read in the bucket file
        size_t index;       //!< the position of the next value to be read in the cache
        size_t available_in_cache;  //!< the number of cached values

        /**
         * @brief A constructor
         *
         * @param[in] bucket is a pointer to the bucket over which iterate
         */
        explicit const_iterator(Bucket<VALUE> const *bucket):
            bucket{bucket}, cache{bucket->cacheable_values}, read_pos{0},
            index{0}, available_in_cache{0}
        {
            available_in_cache = bucket->load_buffer(cache, this->read_pos);
        }

    public:
        /**
         * @brief The empty constructor
         */
        const_iterator():
            bucket{nullptr}, cache{0}, read_pos{0},
            index{0}, available_in_cache{0}
        {}

        /**
         * @brief Move to the next value
         *
         * @return a reference to the updated object
         */
        const_iterator& operator++()
        {
            if (is_end() || bucket != nullptr) {
                ++index;

                if (index>=available_in_cache) {
                    index = 0;

                    available_in_cache = bucket->load_buffer(cache, read_pos);
                }
            }

            return *this;
        }

        /**
         * @brief Dereference the iterator
         *
         * @return a constant reference to the referenced value
         */
        inline const VALUE& operator*() const
        {
            if (is_end()) {
                throw std::runtime_error("No value is available.");
            }
            return cache[index];
        }

        /**
         * @brief Get the pointer to the object referenced by the iterator
         *
         * @return a constant pointer to the referenced value
         */
        inline const VALUE* operator->() const
        {
            return &(this->operator*());
        }

        /**
         * @brief Check whether the end of the iteration has been reached
         *
         * @return `true` if and only if the end of the iteration has
         *      been reached
         */
        bool is_end() const
        {
            return index == 0 && available_in_cache == 0;
        }

        /**
         * @brief Check whether two iterator refer to the same position
         *
         * @param[in] other is a constant iterator over a bucket
         * @return `true` if and only if the current iterator and `other`
         *           are referring to the same position in the bucket
         */
        bool operator==(const const_iterator& other) const
        {
            if (is_end() || other.is_end()) {
                return is_end() && other.is_end();
            }

            return bucket->path() == other.bucket->path()
                && read_pos == other.read_pos
                && index == other.index
                && available_in_cache == other.available_in_cache;
        }

        /**
         * @brief Check whether two iterator refer to different positions
         *
         * @param[in] other is a constant iterator over a bucket
         * @return `false` if and only if the current iterator and `other`
         *           are referring to the same position in the bucket
         */
        inline bool operator!=(const const_iterator& other) const
        {
            return !(*this == other);
        }

        friend class Bucket<VALUE>;
    };

    /**
     * @brief A bucket constructor
     *
     * This constructor creates a bucket object and associates it
     * with a file. If the specified file already exists, it loads
     * the bucket’s data --such as the number of stored values--
     * from the file. If the file does not exist, an empty bucket
     * is created and saved to the new file.
     *
     * @param[in] filepath is the path of the associated bucket file
     * @param[in] cache_size is the write cache size in bytes
     */
    Bucket(const std::filesystem::path filepath,
           const size_t cache_size=sizeof(VALUE)):
        filepath{filepath}, size_pos{0}, data_pos{0}, file_size{0},
        num_of_values{0}
    {
        set_max_cache_size(cache_size);

        init_bucket();
    }

    /**
     * @brief The copy constructor
     *
     * The copy constructor flushes the original object.
     *
     * @param[in, out] orig is the original version of the bucket object
     */
    Bucket(const Bucket& orig):
        filepath{orig.filepath}, size_pos{0}, data_pos{0}, file_size{0},
        num_of_values{0}, cacheable_values{orig.cacheable_values}, write_cache{}
    {
        const_cast<Bucket&>(orig).flush();

        init_bucket();
    }

    /**
     * @brief The copy operator
     *
     * The copy operator flushes the original object.
     *
     * @param[in, out] orig is the original version of the bucket object
     * @return A reference to the updated object
     */
    Bucket& operator=(const Bucket& orig)
    {
        const_cast<Bucket&>(orig).flush();

        filepath = orig.filepath;
        cacheable_values = orig.cacheable_values;

        init_bucket();

        return *this;
    }

    /**
     * @brief The move operator
     *
     * @param[in] orig is the original version of the bucket object
     * @return A reference to the updated object
     */
    Bucket& operator=(Bucket&& orig)
    {
        flush();

        std::swap(filepath, orig.filepath);
        std::swap(size_pos, orig.size_pos);
        std::swap(data_pos, orig.data_pos);
        std::swap(file_size, orig.file_size);
        std::swap(num_of_values, orig.num_of_values);
        std::swap(cacheable_values, orig.cacheable_values);
        std::swap(write_cache, orig.write_cache);

        return *this;
    }

    /**
     * @brief Insert a value in the bucket
     *
     * @param[in] value is the value to be inserted in the bucket
     */
    void push_back(VALUE&& value)
    {
        if (write_cache.size()>=cacheable_values) {
            flush();
        }

        write_cache.push_back(std::move(value));
        ++num_of_values;
    }

    /**
     * @brief Insert a value in the bucket
     *
     * @param[in] value is the value to be inserted in the bucket
     */
    void push_back(const VALUE& value)
    {
        if (write_cache.size()>=cacheable_values) {
            flush();
        }

        write_cache.push_back(value);
        ++num_of_values;
    }

    /**
     * @brief Shuffle the values in the bucket
     *
     * This method shuffles the bucket’s values randomly, ensuring a
     * uniform distribution. At any moment, the number of values held
     * in memory does not exceed the bucket’s write cache size. To
     * ensure this memory constraint, copies of the bucket values are
     * temporarily stored in a set of disk-based files.
     *
     * @tparam RANDOM_GENERATOR is a random number generator type
     * @param[in, out] random_generator is a random number generator
     * @param[in] tmp_dir is the path of the temporary files
     */
    template<typename RANDOM_GENERATOR>
    inline void shuffle(RANDOM_GENERATOR& random_generator,
                 const std::filesystem::path tmp_dir = std::filesystem::temp_directory_path())
    {
        shuffle(random_generator, get_cache_size(), tmp_dir);
    }

    /**
     * @brief Shuffle the values in the bucket
     *
     * This method shuffles the bucket’s values randomly, ensuring a
     * uniform distribution. At any moment, the number of values held
     * in memory does not exceed the specified cache size. To ensure
     * this memory constraint, copies of the bucket values are
     * temporarily stored in a set of disk-based files.
     *
     * @tparam RANDOM_GENERATOR is a random number generator type
     * @param[in, out] random_generator is a random number generator
     * @param[in] buffer_size is the buffer size in bytes
     * @param[in] tmp_dir is the path of the temporary files
     */
    template<typename RANDOM_GENERATOR>
    void shuffle(RANDOM_GENERATOR& random_generator,
                 size_t buffer_size,
                 const std::filesystem::path tmp_dir = std::filesystem::temp_directory_path())
    {
        flush();

        if (size() == 0) {
            return;
        }

        size_t buff_values = buffer_size/sizeof(VALUE);
        if (buff_values >= size()) {
            std::vector<VALUE> buffer(std::max(buff_values, size()));

            load_buffer(buffer, data_pos);

            const auto shuffled_path(get_a_temporary_path("CLONES_shuffled_tmp", tmp_dir));

            Binary::Out archive(shuffled_path);
            Binary::Out::write_header(archive, "CLONES Bucket", 0);

            archive & num_of_values;
            for (const auto& value: buffer) {
                archive & value;
            }

            std::filesystem::rename(shuffled_path, filepath);

            return;
        }

        // split buffer size between reader and writer
        buffer_size /= 2;
        buff_values = buffer_size/sizeof(VALUE);

        const auto chunk_paths = split_in_random_chunks(random_generator, buff_values,
                                                        "tmp_chunk", tmp_dir);

        std::vector<VALUE> buffer(buff_values);

        const auto shuffled_path(get_a_temporary_path("CLONES_shuffled_tmp", tmp_dir));
        Bucket shuffled_bucket(shuffled_path, buffer_size);

        for (const auto& chunk_path: chunk_paths) {
            const auto end_of_buffer = load_buffer(buffer, chunk_path);

            std::filesystem::remove(chunk_path);

            std::shuffle(buffer.begin(), end_of_buffer, random_generator);

            for (auto buffer_it=buffer.begin(); buffer_it != end_of_buffer; ++buffer_it) {
                shuffled_bucket.push_back(*buffer_it);
            }

            shuffled_bucket.flush();
        }

        std::filesystem::rename(shuffled_path, filepath);
    }

    /**
     * @brief Create an iterator referring to the bucket initial position
     *
     * This method also flushes the bucket.
     *
     * @return A constant iterator referring to the bucket initial position
     */
    inline Bucket<VALUE>::const_iterator begin() const
    {
        const_cast<Bucket<VALUE>*>(this)->flush();

        return Bucket<VALUE>::const_iterator(this);
    }

    /**
     * @brief Create an iterator referring to the bucket final position
     *
     * @return A constant iterator referring to the bucket final position
     */
    inline Bucket<VALUE>::const_iterator end() const
    {
        return Bucket<VALUE>::const_iterator();
    }

    /**
     * @brief Rename the bucket file
     *
     * This method renames the file associated to the bucket object.
     *
     * @param[in] new_filepath is the new path of the bucket file
     */
    void rename(const std::filesystem::path& new_filepath)
    {
        std::filesystem::rename(filepath, new_filepath);

        filepath = new_filepath;
    }

    /**
     * @brief Set the maximum bucket's read cache size
     *
     * @param[in] cache_size is the new maximum read cache size in bytes
     */
    void set_max_cache_size(const size_t& cache_size)
    {
        if (cache_size<sizeof(VALUE)) {
            std::ostringstream oss;

            oss << "Bucket: the minimum cache size is "
                << sizeof(VALUE) << ".";
            throw std::domain_error(oss.str());
        }

        this->cacheable_values = cache_size/sizeof(VALUE);
    }

    /**
     * @brief Get the bucket's read cache size
     *
     * This method returns the read cache size in bytes. The values
     * returned by this method may differ from those set by either
     * the constructor or `set_max_cache_size()`.
     *
     * @return the read cache size in bytes
     */
    inline size_t get_cache_size() const
    {
        return this->cacheable_values*sizeof(VALUE);
    }

    /**
     * @brief Get the number of values in the bucket
     *
     * @return the number of values in the bucket
     */
    inline const size_t& size() const
    {
        return num_of_values;
    }

    /**
     * @brief The random access operator
     *
     * This method returns the `i`-th value in the bucket. If the
     * bucket contains less than i+1 values, the method throws an
     * `std::out_of_range` exception.
     *
     * @param i is the index of value in the bucket order
     * @return the `i`-th value in the bucket
     */
    template<typename VALUE2=VALUE, std::enable_if_t<std::is_same_v<VALUE2, VALUE>
                && uses_constant_space_on_disk<VALUE>::value, bool> = true>
    VALUE operator[](const size_t& i) const
    {
        if (i>= size()) {
            throw std::out_of_range("The index is out of the bucket's boundaries.");
        }
        std::streamoff value_pos = get_value_pos(i);

        if (value_pos < file_size) {
            Binary::In archive(filepath);

            archive.seekg(value_pos);

            VALUE value;

            archive & value;

            return value;
        }

        auto it = write_cache.begin();
        std::advance(it, i-(num_of_values-write_cache.size()));

        return *it;
    }


    /**
     * @brief The random access operator
     *
     * This method returns the `i`-th value in the bucket. If the
     * bucket contains less than i+1 values, the method throws an
     * `std::out_of_range` exception.
     *
     * @param i is the index of value in the bucket order
     * @return the `i`-th value in the bucket
     */
    template<typename VALUE2=VALUE, std::enable_if_t<std::is_same_v<VALUE2, VALUE>
                && uses_constant_space_on_disk<VALUE>::value, bool> = true>
    VALUE operator[](size_t&& i) const
    {
        return operator[](i);
    }

    /**
     * @brief Choose a random value from the bucket
     *
     * This method chooses a random value from the bucket with uniform
     * distribution. The selected value is *not* removed from the bucket.
     * The selection is performed by using the associated random number
     * generator.
     *
     * @tparam RANDOM_GENERATOR is a type of random number generator
     * @return the choosen value
     */
    template<typename RANDOM_GENERATOR, typename VALUE2=VALUE,
                std::enable_if_t<std::is_same_v<VALUE2, VALUE>
                                 && uses_constant_space_on_disk<VALUE>::value, bool> = true>
    VALUE choose(RANDOM_GENERATOR& random_generator) const
    {
        if (this->size()==0) {
            throw std::runtime_error("No value in the bucket.");
        }

        std::uniform_int_distribution<size_t> dist(0, this->size()-1);

        return operator[](dist(random_generator));
    }

    /**
     * @brief Build a random tour for the flushed values in the bucket
     *
     * This method builds a random tour for the flushed values in the bucket.
     *
     * @tparam RANDOM_GENERATOR is the type of the random number generator
     * @param random_generator is the random number generator used for the tour
     * @param cache_size is the tour read cache size in bytes
     * @return a random tour for the buffer
     */
    template<typename RANDOM_GENERATOR>
    BucketRandomTour<VALUE, RANDOM_GENERATOR> random_tour(const RANDOM_GENERATOR& random_generator,
                                                          const size_t cache_size) const
    {
        return {*this, random_generator, cache_size};
    }

    /**
     * @brief Build a random tour for the flushed values in the bucket
     *
     * This method builds a random tour for the flushed values in the bucket.
     * The read cache size is set according to the bucket write cache size.
     *
     * @tparam RANDOM_GENERATOR is the type of the random number generator
     * @param random_generator is the random number generator used for the tour
     * @return a random tour for the buffer
     */
    template<typename RANDOM_GENERATOR>
    inline BucketRandomTour<VALUE, RANDOM_GENERATOR> random_tour(const RANDOM_GENERATOR& random_generator) const
    {
        return random_tour(random_generator, get_cache_size());
    }

    /**
     * @brief Get the bucket file path
     *
     * @return the bucket file path
     */
    inline std::filesystem::path path() const
    {
        return filepath;
    }

    /**
     * @brief Flush the write cache
     *
     * This method writes the values in the cache into the bucket
     * file on the disk.
     */
    void flush()
    {
        if (write_cache.empty()) {
            return;
        }

        {
            Archive::Binary::Out archive(filepath, std::ios::in);

            archive.seekp(size_pos);

            archive & num_of_values;
        }

        Archive::Binary::Out archive(filepath, std::ios::in | std::ios::app);

        archive.seekp(0, std::ios::end);

        for (const auto& value: write_cache) {
            archive & value;
        }

        file_size = archive.tellg();

        write_cache.clear();

        archive.flush();
    }

    /**
     * @brief Destroy the bucket object
     */
    ~Bucket()
    {
        flush();
    }

    friend class Bucket<VALUE>::const_iterator;

    template<typename VALUE2, typename RANDOM_GENERATOR>
    friend class BucketRandomTour;
};

/**
 * @brief A class to represent random tour over buckets
 *
 * This class implements random tours over flushed values in a bucket.
 *
 * Each tour visits every value in the bucket file exactly once. This
 * is achieved by using a constant amount of memory: at any time, a
 * fixed number of values are loaded into a buffer by reading the
 * bucket file in chunks. If all objects of the value type occupy equal
 * space on disk, the initial read position is chosen randomly with a
 * uniform distribution over the positions of the values in the bucket
 * file. If not, the first values loaded into the buffer come from the
 * beginning of the bucket file. The values in the buffer are then
 * iterated in random order, ensuring a uniform distribution over the
 * unvisited ones. Once all values in the buffer have been iterated,
 * the next chunk from the bucket file is loaded until the initial
 * reading position is reached.
 *
 * If the referred bucket changes, the behaviour of the objects of
 * `BucketRandomTour` is not defined.
 *
 * @tparam VALUE is the type of the values contained in the bucket
 * @tparam RANDOM_GENERATOR is the type of random number generator
 *      used to randomize the tour
 */
template<typename VALUE, typename RANDOM_GENERATOR>
class BucketRandomTour
{
    Bucket<VALUE> const& bucket;    //!< the bucket to be toured

    RANDOM_GENERATOR random_generator; //!< the random number generator

    size_t cacheable_values;    //!< the number of cached values

public:
    /**
     * @brief The randomized bucket random number generator
     */
    using random_generator_type = RANDOM_GENERATOR;

    /**
     * @brief A constant iterator for the values in the bucket
     *
     * The objects of this class use the bucket random number
     * generator to randomize the iteration over the bucket
     * values.
     */
    class const_iterator
    {
        Bucket<VALUE> const* bucket; //!< a pointer to the bucket

        RANDOM_GENERATOR random_generator;  //!< the random number generator

        std::vector<VALUE> cache;   //!< the read cache
        std::streampos initial_pos; //!< the position of the first read value from the file
        std::streampos read_pos;    //!< the position of the value to be read from the file
        size_t available_in_cache;  //!< the number of values available in the cache
        size_t iterated;            //!< the number of values already reached

        /**
         * @brief Randomly choose a value from the cache
         *
         * This method randomly chooses a value from the cache, moves it at the end of the
         * cache, and reduces the cache size by one.
         */
        void select_a_value_in_cache()
        {
            if (available_in_cache>0) {
                std::uniform_int_distribution<size_t> dist(0, available_in_cache-1);

                const size_t pos = dist(random_generator);

                std::swap(cache[pos], cache[available_in_cache-1]);

                ++iterated;
            }
        }

        /**
         * @brief A constructor
         *
         * @param[in, out] bucket is a pointer to the bucket whose values must be iterated
         * @param[in] initial_pos is the first position to be read from the bucket file
         * @param[in] cacheable_values is the maximum number of cached values
         */
        const_iterator(Bucket<VALUE> const* bucket, const std::streampos initial_pos,
                       const size_t cacheable_values):
            bucket{bucket}, cache{cacheable_values}, initial_pos{initial_pos},
            read_pos{initial_pos}, iterated{0}
        {
            available_in_cache = bucket->load_buffer(cache, this->read_pos,
                                                     this->initial_pos, true);

            select_a_value_in_cache();
        }

    public:
        /**
         * @brief The empty constructor
         */
        const_iterator():
            bucket{nullptr}, cache{0}, initial_pos{0}, read_pos{0}, available_in_cache{0},
            iterated{0}
        {}

        /**
         * @brief Move to the next value in the randomized order
         *
         * @return a reference to the updated object
         */
        const_iterator& operator++()
        {
            if (bucket != nullptr && !is_end()) {
                if (available_in_cache>0) {
                    --available_in_cache;
                }
                if (available_in_cache==0 && read_pos != initial_pos) {
                    available_in_cache = bucket->load_buffer(cache, read_pos, initial_pos);
                }

                if (is_end()) {
                    ++iterated;
                } else {
                    select_a_value_in_cache();
                }
            }

            return *this;
        }

        /**
         * @brief Dereference the iterator
         *
         * @return a constant reference to the referenced value
         */
        inline const VALUE& operator*() const
        {
            if (is_end()) {
                throw std::runtime_error("No more value available.");
            }
            return cache[available_in_cache-1];
        }

        /**
         * @brief Get the pointer to the object referenced by the iterator
         *
         * @return a constant pointer to the referenced value
         */
        inline const VALUE* operator->() const
        {
            return &(this->operator*());
        }

        /**
         * @brief Check whether the end of the iteration has been reached
         *
         * @return `true` if and only if the end of the iteration has
         *      been reached
         */
        bool is_end() const
        {
            return available_in_cache == 0 && read_pos == initial_pos;
        }

        /**
         * @brief Get the number of values in the tour already reached
         *
         * @return the number of values in the bucket that have already been
         *      reached by the tour
         */
        inline const size_t& reached_values() const
        {
            return iterated;
        }

        /**
         * @brief Get the number of values to complete the tour
         *
         * @return the number of values in the bucket that have not been
         *      reached by the tour yet
         */
        inline size_t remaining_values() const
        {
            return (bucket->size()+1)-reached_values();
        }

        /**
         * @brief Check whether two iterator refer to the same position
         *
         * @param[in] other is a constant iterator over a bucket
         * @return `true` if and only if the current iterator and `other`
         *           are referring to the same position in the bucket
         */
        bool operator==(const const_iterator& other) const
        {
            if (is_end() || other.is_end()) {
                return is_end() && other.is_end();
            }

            return bucket->path() == other.bucket->path()
                && read_pos == other.read_pos
                && initial_pos == other.initial_pos
                && available_in_cache == other.available_in_cache;
        }

        /**
         * @brief Check whether two iterator refer to different positions
         *
         * @param[in] other is a constant iterator over a bucket
         * @return `false` if and only if the current iterator and `other`
         *           are referring to the same position in the bucket
         */
        inline bool operator!=(const const_iterator& other) const
        {
            return !(*this == other);
        }

        friend class BucketRandomTour<VALUE, RANDOM_GENERATOR>;
    };

    /**
     * @brief A constructor for bucket random tours
     *
     * This constructor maintains a copy of a random number generator
     * passed as parameter.
     *
     * @param[in] bucket is bucket for which the random tour is
     *      requested
     * @param[in] random_generator is the random number generator that
     *      randomizes the tour
     * @param cache_size is the cache size in bytes
     */
    BucketRandomTour(const Bucket<VALUE>& bucket,
                     const RANDOM_GENERATOR& random_generator,
                     const size_t cache_size=sizeof(VALUE)):
        bucket{bucket}, random_generator{random_generator}
    {
        this->set_max_cache_size(cache_size);
    }

    /**
     * @brief Create an iterator referring to the random tour initial position
     *
     * The referred bucket is flushed.
     *
     * @return A constant iterator referring to the bucket initial position
     */
    inline BucketRandomTour<VALUE, RANDOM_GENERATOR>::const_iterator begin() const
    {
        const_cast<Bucket<VALUE>&>(bucket).flush();

        std::streampos begin_pos{bucket.get_data_pos()};
        if constexpr(uses_constant_space_on_disk<VALUE>::value) {
            if (bucket.size()>0) {
                std::uniform_int_distribution<size_t> dist(0, bucket.size()-1);

                RANDOM_GENERATOR random_gen_copy(random_generator);

                const size_t first_index = dist(random_gen_copy);

                begin_pos = bucket.get_value_pos(first_index);
            }
        }

        return const_iterator(&bucket, begin_pos, cacheable_values);
    }

    /**
     * @brief Create an iterator referring to the random tour final position
     *
     * @return A constant iterator referring to the random tour final position
     */
    BucketRandomTour<VALUE, RANDOM_GENERATOR>::const_iterator end() const
    {
        return BucketRandomTour<VALUE, RANDOM_GENERATOR>::const_iterator();
    }

    /**
     * @brief Set the maximum bucket's read cache size
     *
     * @param[in] cache_size is the new maximum read cache size in bytes
     */
    void set_max_cache_size(const size_t& cache_size)
    {
        if (cache_size<sizeof(VALUE)) {
            std::ostringstream oss;

            oss << "BucketRandomTour: the minimum cache size is "
                << sizeof(VALUE) << ".";
            throw std::domain_error(oss.str());
        }

        this->cacheable_values = cache_size/sizeof(VALUE);
    }

    /**
     * @brief Get the bucket's read cache size
     *
     * This method returns the read cache size in bytes. The values
     * returned by this method may differ from those set by either
     * the constructor or `set_max_cache_size()`.
     *
     * @return the read cache size in bytes
     */
    inline const size_t& get_cache_size() const
    {
        return this->cacheable_values*sizeof(VALUE);
    }

    /**
     * @brief Get the bucket
     *
     * @return a constant reference to the bucket
     */
    inline const Bucket<VALUE>& get_bucket() const
    {
        return bucket;
    }

    /**
     * @brief Get the random number generator
     *
     * @return a constant reference to random number generator
     */
    inline const RANDOM_GENERATOR& get_random_generator() const
    {
        return random_generator;
    }

    /**
     * @brief Get the random number generator
     *
     * @return a reference to random number generator
     */
    inline RANDOM_GENERATOR& get_random_generator()
    {
        return random_generator;
    }

    /**
     * @brief Set the random number generator
     *
     * @param random_generator is the new random number generator
     * @return a constant reference to the object's random number
     *      generator
     */
    const RANDOM_GENERATOR&
    set_random_generator(const RANDOM_GENERATOR& random_generator) const
    {
        this->random_generator = random_generator;

        return this->random_generator;
    }

    friend class Bucket<VALUE>;
    friend class BucketRandomTour<VALUE, RANDOM_GENERATOR>::const_iterator;
};

}   // Archive

}   // CLONES

#endif // __CLONES_BUCKET__