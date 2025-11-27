/**
 * @file bucket.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines bucket
 * @version 1.5
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

#ifndef __RACES_BUCKET__
#define __RACES_BUCKET__

#include <vector>
#include <random>
#include <memory>
#include <filesystem>
#include <numeric> // iota
#include <concepts>

#include <algorithm> // shuffle

#include "archive.hpp"
#include "utils.hpp"
#include "progress_bar.hpp"

namespace RACES
{

namespace Archive
{

template<typename VALUE, typename RANDOM_GENERATOR>
class BucketRandomTour;

/**
 * @brief A base class for buckets
 *
 * A bucket is a container that stores a collection of values in
 * a file. These values are maintained in a specific order and can
 * be read sequentially from first to last. The order of values
 * within the bucket can also be shuffled if needed.
 *
 * @tparam VALUE is the type of the values contained in the bucket
 */
template<typename VALUE>
class BucketBase
{
protected:
    std::filesystem::path filepath;   //!< the name of the file storing the bucket

    std::streampos size_pos;    //!< the position of the bucket size in the file
    std::streampos data_pos;    //!< the position of the first value in the file
    std::streampos final_pos;   //!< the last position in the file

    size_t num_of_values;   //!< the number of values in the bucket (i.e., bucket size)

    /**
     * @brief Read the bucket file header
     *
     * This method reads the bucket file header and sets the
     * position of the bucket size and that of the data in the
     * file.
     * @param archive is a reference to the input archive
     * @return a reference to the input archive
     */
    Archive::Binary::In& read_header(Archive::Binary::In& archive)
    {
        if (!archive.is_open()) {
            archive.open(filepath);
        }

        Archive::Binary::In::read_header(archive, "RACES Bucket", 0);

        size_pos = archive.tellg();
        archive & num_of_values;
        data_pos = archive.tellg();

        final_pos = archive.size();

        return archive;
    }

    /**
     * @brief Read the bucket file header
     *
     * This method reads the bucket file header and sets the
     * position of the bucket size and that of the data in the
     * file.
     */
    inline void read_header()
    {
        Archive::Binary::In archive(filepath);

        read_header(archive);
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
        if (i<num_of_values) {
            const std::streamoff offset = i*(final_pos-data_pos)/num_of_values;
            return data_pos + offset;
        }

        return final_pos;
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

public:
    BucketBase():
        filepath{}, size_pos{0}, data_pos{0}, final_pos{0},
        num_of_values{0}
    {}

    explicit BucketBase(const std::filesystem::path& filepath):
        BucketBase<VALUE>{}
    {
        this->filepath = filepath;
    }

    BucketBase(const BucketBase<VALUE>& orig_obj):
        filepath{orig_obj.filepath}, size_pos{orig_obj.size_pos},
        data_pos{orig_obj.data_pos}, final_pos{orig_obj.final_pos},
        num_of_values{orig_obj.num_of_values}
    {}

    BucketBase<VALUE>& operator=(BucketBase<VALUE>&& orig_obj)
    {
        swap(filepath, orig_obj.filepath);
        swap(size_pos, orig_obj.size_pos);
        swap(data_pos, orig_obj.data_pos);
        swap(final_pos, orig_obj.final_pos);
        swap(num_of_values, orig_obj.num_of_values);

        return *this;
    }

    /**
     * @brief Get the bucket size position in the file
     *
     * @return the bucket size position in the file
     */
    inline std::streampos get_size_pos() const
    {
        return size_pos;
    }

    /**
     * @brief Get the bucket file final position
     *
     * @return the bucket file final position
     */
    inline std::streampos get_final_pos() const
    {
        return final_pos;
    }

    /**
     * @brief Get the bucket data position in the file
     *
     * @return the bucket data position in the file
     */
    inline std::streampos get_data_pos() const
    {
        return data_pos;
    }

    /**
     * @brief Get the bucket file path
     *
     * @return the bucket file path
     */
    inline std::filesystem::path get_path() const
    {
        return filepath;
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
};

/**
 * @brief A bucket writer
 *
 * @tparam VALUE is the type of the values contained in the bucket
 */
template<typename VALUE>
class BucketWriter : public BucketBase<VALUE>
{
    std::vector<VALUE> cache; //!< the write cache

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
            throw std::runtime_error("BucketWriter::load_buffer(): The file is larger than the buffer");
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
     * @param[in, out] progress_bar is a progress bar
     * @return the vector of the chunk paths
     */
    template<typename RANDOM_GENERATOR>
    std::vector<std::filesystem::path>
    split_in_random_chunks(RANDOM_GENERATOR& random_generator,
                           const size_t& max_chunk_size,
                           const std::string& prefix_name,
                           const std::filesystem::path& tmp_dir,
                           UI::ProgressBar& progress_bar)
    {
        const auto num_of_chunks = (this->size()-1)/max_chunk_size+1;
        auto last_chunk = num_of_chunks-1;

        std::uniform_int_distribution<size_t> chunk_dist(0, last_chunk);

        std::vector<Archive::Binary::Out> chunks(num_of_chunks);
        std::vector<size_t> chunk_sizes(num_of_chunks, 0);

        const auto chunk_paths = create_chunks(chunks, prefix_name, tmp_dir);
        std::vector<size_t> positions(num_of_chunks);
        std::iota(positions.begin(), positions.end(), 0);

        std::vector<VALUE> cache{max_chunk_size};
        std::streampos read_pos{0};

        size_t value_in_caches = BucketBase<VALUE>::load_buffer(cache, read_pos);
        auto cache_it = cache.begin();
        for (size_t i=0; i<this->size(); ++i) {
            if (value_in_caches == 0) {
                value_in_caches = BucketBase<VALUE>::load_buffer(cache, read_pos);

                cache_it = cache.begin();

                progress_bar.update_elapsed_time();
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
     * @brief Write the bucket file header
     *
     * This method writes the bucket file header and sets the
     * position of the bucket size and that of the data in the
     * file.
     * @param archive is a reference to the output archive
     * @return a reference to the output archive
     */
    Archive::Binary::Out& write_header(Archive::Binary::Out& archive)
    {
        if (!archive.is_open()) {
            archive.open(this->filepath);
        }

        Archive::Binary::Out::write_header(archive, "RACES Bucket", 0);

        this->size_pos = archive.tellg();
        archive & this->num_of_values;
        this->data_pos = archive.tellg();

        archive.flush();

        this->final_pos = this->data_pos;

        return archive;
    }

    /**
     * @brief Write the buffer file header
     *
     * This method writes the buffer file header and sets the
     * position of the bucket size and that of the data in the
     * file.
     */
    inline void write_header()
    {
        Archive::Binary::Out archive(this->filepath);

        write_header(archive);
    }

    /**
     * @brief The method that initializes Bucket objects
     *
     * This method is used to initialize Bucket objects during the construction
     * and the copy.
     */
    void init_bucket()
    {
        if (std::filesystem::exists(this->filepath)) {
            if (!std::filesystem::is_regular_file(this->filepath)) {
                std::ostringstream oss;

                oss << "\"" << to_string(this->filepath)
                    << "\" is not a block file.";
                throw std::domain_error(oss.str());
            }

            this->read_header();
        } else {
            write_header();
        }
    }

    /**
     * @brief Shuffle the bucket values in memory
     *
     * This method shuffles the bucket’s values randomly, ensuring a
     * uniform distribution. The values are shuffled in memory and,
     * then, saved into the bucket.
     *
     * @tparam RANDOM_GENERATOR is a random number generator type
     * @param[in, out] random_generator is a random number generator
     * @param[in, out] progress_bar is a progress bar
     */
    template<typename RANDOM_GENERATOR>
    void shuffle_in_memory(RANDOM_GENERATOR& random_generator,
                           UI::ProgressBar& progress_bar)
    {
        flush();

        if (this->size() == 0) {
            return;
        }

        // temporary remove cache
        const size_t cacheable_values = cache.size();

        cache = std::vector<VALUE>(this->size());

        // load the values
        BucketBase<VALUE>::load_buffer(cache, this->data_pos);

        // shuffle the values
        std::shuffle(cache.begin(), cache.end(), random_generator);

        progress_bar.update_elapsed_time();

        std::filesystem::remove(this->filepath);

        // save the values
        Archive::Binary::Out archive{this->filepath};
        write_header(archive);
        for (const auto& value: cache) {
            archive & value;
        }
        this->final_pos = archive.tellg();

        progress_bar.update_elapsed_time();

        // restore cache
        cache = std::vector<VALUE>(cacheable_values);
    }


    /**
     * @brief Shuffle the bucket values on disk
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
     * @param[in, out] progress_bar is a progress bar
     */
    template<typename RANDOM_GENERATOR>
    void shuffle_on_disk(RANDOM_GENERATOR& random_generator,
                         size_t buffer_size,
                         const std::filesystem::path tmp_dir,
                         UI::ProgressBar& progress_bar)
    {
        flush();

        if (this->size() == 0) {
            return;
        }

        // temporary remove cache
        const size_t cacheable_values = cache.size();

        // split buffer size between reader and writer
        const size_t buff_values = (buffer_size / 2)/sizeof(VALUE);

        if (buff_values == 0) {
            throw std::domain_error("BucketWriter::shuffle_on_disk(): the minimum "
                                    "buffer size is " + std::to_string(2*sizeof(VALUE))
                                    + ". It has been set to "
                                    + std::to_string(buffer_size) + ".");
        }

        const auto chunk_paths = split_in_random_chunks(random_generator, buff_values,
                                                        "tmp_chunk", tmp_dir,
                                                        progress_bar);

        cache = std::vector<VALUE>(buff_values);

        std::filesystem::remove(this->filepath);

        Archive::Binary::Out archive(this->filepath);
        write_header(archive);

        for (const auto& chunk_path: chunk_paths) {
            const auto end_of_buffer = load_buffer(cache, chunk_path);

            std::filesystem::remove(chunk_path);

            std::shuffle(cache.begin(), end_of_buffer, random_generator);

            for (auto buffer_it=cache.begin(); buffer_it != end_of_buffer; ++buffer_it) {
                archive & *buffer_it;
            }
            this->final_pos = archive.tellg();

            progress_bar.update_elapsed_time();
        }

        // restore cache
        cache = std::vector<VALUE>(cacheable_values);
    }
public:

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
    BucketWriter(const std::filesystem::path filepath,
                 const size_t cache_size=1000*sizeof(VALUE)):
        BucketBase<VALUE>{filepath}
    {
        set_cache_size(cache_size);

        init_bucket();
    }

    /**
     * @brief The copy constructor
     *
     * The copy constructor flushes the original object.
     *
     * @param[in, out] orig is the original version of the bucket object
     */
    BucketWriter(const BucketWriter<VALUE>& orig):
        BucketBase<VALUE>{orig.filepath}
    {
        const_cast<BucketWriter&>(orig).flush();
        cache.reserve(orig.cache.capacity());

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
    BucketWriter<VALUE>& operator=(const BucketWriter<VALUE>& orig)
    {
        const_cast<BucketWriter&>(orig).flush();

        this->filepath = orig.filepath;
        cache.reserve(orig.cache.capacity());

        init_bucket();

        return *this;
    }

    /**
     * @brief The move operator
     *
     * @param[in] orig is the original version of the bucket object
     * @return A reference to the updated object
     */
    BucketWriter<VALUE>& operator=(BucketWriter<VALUE>&& orig)
    {
        flush();

        std::swap(cache, orig.cache);
        static_cast<BucketBase<VALUE> *>(this)->operator=(std::move(orig));

        return *this;
    }

    /**
     * @brief Insert a value in the bucket
     *
     * @param[in] value is the value to be inserted in the bucket
     */
    void push_back(VALUE&& value)
    {
        push_back(value);
    }

    /**
     * @brief Insert a value in the bucket
     *
     * @param[in] value is the value to be inserted in the bucket
     */
    void push_back(const VALUE& value)
    {
        if (cache.size()==cache.capacity()) {
            flush();
        }

        cache.push_back(value);

        ++this->num_of_values;
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
     * @tparam ARGS are the types of the remaining parameters
     * @param[in] seed is the seed of a random number generator
     * @param[in, out] args are the remaining parameters
     */
    template<typename RANDOM_GENERATOR, typename... ARGS>
    inline void shuffle(const int seed, ARGS... args)
    {
        RANDOM_GENERATOR random_generator(seed);

        shuffle(random_generator, args...);
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
     * @tparam ARGS are the types of the parameters
     * @param[in, out] args are the parameters
     */
    template <typename... ARGS>
    requires last_is_not<UI::ProgressBar, ARGS...>
    inline void shuffle(ARGS... args)
    {
        UI::ProgressBar progress_bar;

        shuffle(args..., progress_bar);
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
     * @param[in, out] progress_bar is a progress bar
     */
    template<typename RANDOM_GENERATOR>
    inline void shuffle(RANDOM_GENERATOR& random_generator,
                        const std::filesystem::path tmp_dir,
                        UI::ProgressBar& progress_bar)
    {
        shuffle(random_generator, get_cache_size(), tmp_dir, progress_bar);
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
     * @param[in, out] progress_bar is a progress bar
     */
    template<typename RANDOM_GENERATOR>
    inline void shuffle(RANDOM_GENERATOR& random_generator,
                        UI::ProgressBar& progress_bar)
    {
        shuffle(random_generator, get_cache_size(),
                std::filesystem::temp_directory_path(),
                progress_bar);
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
     * @param[in, out] progress_bar is a progress bar
     */
    template<typename RANDOM_GENERATOR>
    void shuffle(RANDOM_GENERATOR& random_generator,
                 const size_t buffer_size,
                 const std::filesystem::path tmp_dir,
                 UI::ProgressBar& progress_bar)
    {
        const size_t buff_values = buffer_size/sizeof(VALUE);
        if (buff_values >= this->size()) {
            shuffle_in_memory(random_generator, progress_bar);
        } else {
            shuffle_on_disk(random_generator, buffer_size,
                            tmp_dir, progress_bar);
        }
    }

    /**
     * @brief Set the maximum bucket's read cache size
     *
     * @param[in] cache_size is the new maximum read cache size in bytes
     */
    void set_cache_size(const size_t& cache_size)
    {
        if (cache_size<sizeof(VALUE)) {
            std::ostringstream oss;

            oss << "BucketWriter: the minimum cache size is "
                << sizeof(VALUE) << ".";
            throw std::domain_error(oss.str());
        }

        const size_t cacheable_values = cache_size/sizeof(VALUE);
        if (cache.size() > cacheable_values) {
            flush();
        }

        cache.reserve(cacheable_values);
    }

    /**
     * @brief Get the bucket's cache size
     *
     * This method returns the cache size in bytes. The values
     * returned by this method may differ from those set by either
     * the constructor or `set_cache_size()`. As a matter of
     * the fact, the returned value is evaluated as the number
     * of cacheable values multiplied by the value size.
     *
     * @return the cache size in bytes
     */
    inline size_t get_cache_size() const
    {
        return cache.capacity()*sizeof(VALUE);
    }

    /**
     * @brief Flush the write cache
     *
     * This method writes the values in the cache into the bucket
     * file on the disk.
     */
    void flush()
    {
        {
            Archive::Binary::Out archive(this->filepath, std::ios::in);

            archive.seekp(this->size_pos);

            archive & this->num_of_values;
        }

        if (cache.size()==0) {
            return;
        }

        Archive::Binary::Out archive(this->filepath, std::ios::in | std::ios::app);

        archive.seekp(0, std::ios::end);

        for (const auto& value : cache) {
            archive & value;
        }

        this->final_pos = archive.tellg();

        cache.resize(0);

        archive.flush();
    }

    /**
     * @brief Destroy the bucket object
     */
    ~BucketWriter()
    {
        flush();
    }
};

/**
 * @brief A bucket reader
 *
 * @tparam VALUE is the type of the values contained in the bucket
 */
template<typename VALUE>
class BucketReader : public BucketBase<VALUE>
{
    size_t cacheable_values;  //!< The number of values in the cache

    /**
     * @brief The method that initializes Bucket objects
     *
     * This method is used to initialize Bucket objects during the construction
     * and the copy.
     */
    void init_bucket()
    {
        if (!(std::filesystem::exists(this->filepath)
              &&std::filesystem::is_regular_file(this->filepath))) {
            std::ostringstream oss;

            oss << "BucketReader: \"" << to_string(this->filepath)
                << "\" is not a block file.";
            throw std::domain_error(oss.str());
        }

        this->read_header();
    }

    /**
     * @brief Load values into a buffer
     *
     * This method wraps `BucketBase<VALUE>::load_buffer()` to make it available to
     * `BucketReader::const_iterator`.
     *
     * @param[in, out] buffer is the buffer that will store the read values
     * @param[in, out] read_pos is the position in the bucket file from which values
     *      are load
     * @return the number of values read from the bucket file
     */
    inline size_t load_buffer(std::vector<VALUE>& buffer, std::streampos& read_pos) const
    {
        return BucketBase<VALUE>::load_buffer(buffer, read_pos);
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
        Archive::Binary::In archive(this->filepath);

        if (read_pos < this->data_pos) {
            read_pos = this->data_pos;
        }
        if (final_pos < this->data_pos) {
            final_pos = read_pos;
        }
        archive.seekg(read_pos);

        size_t read_values{0};
        for (auto& value: buffer) {
            if (archive.eof()) {
                read_pos = this->data_pos;
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
        BucketReader<VALUE> const* bucket_reader;    //!< a pointer to the bucket reader

        std::vector<VALUE> cache;   //!< the read cache
        std::streampos read_pos;    //!< the next position to be read in the bucket file
        size_t index;       //!< the position of the next value to be read in the cache
        size_t available_in_cache;  //!< the number of cached values

        /**
         * @brief A constructor
         *
         * @param[in] bucket_reader is a pointer to the bucket over which iterate
         */
        explicit const_iterator(BucketReader<VALUE> const *bucket_reader):
            bucket_reader{bucket_reader}, cache{bucket_reader->cacheable_values},
            read_pos{0}, index{0}, available_in_cache{0}
        {
            available_in_cache = bucket_reader->load_buffer(cache, this->read_pos);
        }

    public:
        /**
         * @brief The empty constructor
         */
        const_iterator():
            bucket_reader{nullptr}, cache{0}, read_pos{0},
            index{0}, available_in_cache{0}
        {}

        /**
         * @brief Move to the next value
         *
         * @return a reference to the updated object
         */
        const_iterator& operator++()
        {
            if (is_end() || bucket_reader != nullptr) {
                ++index;

                if (index>=available_in_cache) {
                    index = 0;

                    available_in_cache = bucket_reader->load_buffer(cache, read_pos);
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

            return bucket_reader->get_path() == other.bucket_reader->get_path()
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

        friend class BucketReader<VALUE>;
    };

    /**
     * @brief A bucket reader constructor
     *
     * This constructor creates a bucket reader object associated
     * with a file. If the specified file already exists, it loads
     * the bucket’s data --such as the number of stored values--
     * from the file. If the file does not exist, an
     * `std::domain_error` exception is thrown.
     *
     * @param[in] filepath is the path of the associated bucket file
     * @param[in] cache_size is the size of the read cache
     */
    explicit BucketReader(const std::filesystem::path filepath,
                          const size_t cache_size=1000 * sizeof(VALUE)):
        BucketBase<VALUE>{filepath},
        cacheable_values{cache_size/sizeof(VALUE)}
    {
        if (cacheable_values==0) {
            std::ostringstream oss;

            oss << "BucketReader: the minimum cache size is "
                << sizeof(VALUE) << ".";
            throw std::domain_error(oss.str());
        }

        init_bucket();
    }

    /**
     * @brief The copy constructor
     *
     * @param[in, out] orig is the original version of the bucket reader object
     */
    BucketReader(const BucketReader<VALUE>& orig):
        BucketBase<VALUE>{static_cast<const BucketBase<VALUE>&>(orig)},
        cacheable_values{orig.cacheable_values}
    {}

    /**
     * @brief The copy operator
     *
     * @param[in, out] orig is the original version of the bucket object
     * @return A reference to the updated object
     */
    BucketReader<VALUE>& operator=(const BucketReader<VALUE>& orig)
    {
        const_cast<BucketBase<VALUE> *>(this)->operator=(orig);
        cacheable_values = orig.cacheable_values;
    }

    /**
     * @brief Create an iterator referring to the bucket initial position
     *
     * This method also flushes the bucket.
     *
     * @return A constant iterator referring to the bucket initial position
     */
    inline BucketReader<VALUE>::const_iterator begin() const
    {
        return BucketReader<VALUE>::const_iterator(this);
    }

    /**
     * @brief Create an iterator referring to the bucket final position
     *
     * @return A constant iterator referring to the bucket final position
     */
    inline BucketReader<VALUE>::const_iterator end() const
    {
        return BucketReader<VALUE>::const_iterator();
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
        std::filesystem::rename(this->filepath, new_filepath);

        this->filepath = new_filepath;
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
        if (i>= this->size()) {
            throw std::out_of_range("The index is out of the bucket's boundaries.");
        }

        const std::streamoff value_pos = this->get_value_pos(i);

        Binary::In archive(this->filepath);

        archive.seekg(value_pos);

        VALUE value;

        archive & value;

        return value;
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
    BucketRandomTour<VALUE, RANDOM_GENERATOR>
    random_tour(const RANDOM_GENERATOR& random_generator,
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
    inline BucketRandomTour<VALUE, RANDOM_GENERATOR>
    random_tour(const RANDOM_GENERATOR& random_generator) const
    {
        return random_tour(random_generator, cacheable_values*sizeof(VALUE));
    }

    friend class BucketReader<VALUE>::const_iterator;

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
    BucketReader<VALUE> const& bucket_reader;    //!< the bucket reader to be toured

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
        BucketReader<VALUE> const* bucket_reader; //!< a pointer to the bucket reader

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
         * @param[in, out] bucket_reader is a pointer to the bucket reader whose values must be iterated
         * @param[in] initial_pos is the first position to be read from the bucket file
         * @param[in] cacheable_values is the maximum number of cached values
         */
        const_iterator(BucketReader<VALUE> const* bucket_reader, const std::streampos initial_pos,
                       const size_t cacheable_values):
            bucket_reader{bucket_reader}, cache{cacheable_values}, initial_pos{initial_pos},
            read_pos{initial_pos}, iterated{0}
        {
            available_in_cache = bucket_reader->load_buffer(cache, this->read_pos,
                                                            this->initial_pos, true);

            select_a_value_in_cache();
        }

    public:
        /**
         * @brief The empty constructor
         */
        const_iterator():
            bucket_reader{nullptr}, cache{0}, initial_pos{0}, read_pos{0}, available_in_cache{0},
            iterated{0}
        {}

        /**
         * @brief Move to the next value in the randomized order
         *
         * @return a reference to the updated object
         */
        const_iterator& operator++()
        {
            if (bucket_reader != nullptr && !is_end()) {
                if (available_in_cache>0) {
                    --available_in_cache;
                }
                if (available_in_cache==0 && read_pos != initial_pos) {
                    available_in_cache = bucket_reader->load_buffer(cache, read_pos, initial_pos);
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
            return (bucket_reader->size()+1)-reached_values();
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

            return bucket_reader->get_path() == other.bucket_reader->get_path()
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
     * @param[in] bucket_reader is bucket for which the random tour is
     *      requested
     * @param[in] random_generator is the random number generator that
     *      randomizes the tour
     * @param cache_size is the cache size in bytes
     */
    BucketRandomTour(const BucketReader<VALUE>& bucket_reader,
                     const RANDOM_GENERATOR& random_generator,
                     const size_t cache_size=sizeof(VALUE)):
        bucket_reader{bucket_reader}, random_generator{random_generator}
    {
        this->set_cache_size(cache_size);
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
        std::streampos begin_pos{bucket_reader.data_pos};
        if constexpr(uses_constant_space_on_disk<VALUE>::value) {
            if (bucket_reader.size()>0) {
                std::uniform_int_distribution<size_t> dist(0, bucket_reader.size()-1);

                RANDOM_GENERATOR random_gen_copy(random_generator);

                const size_t first_index = dist(random_gen_copy);

                begin_pos = bucket_reader.get_value_pos(first_index);
            }
        }

        return const_iterator(&bucket_reader, begin_pos, cacheable_values);
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
    void set_cache_size(const size_t& cache_size)
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
     * the constructor or `set_cache_size()`.
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
    inline const BucketReader<VALUE>& get_bucket_reader() const
    {
        return bucket_reader;
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

    friend class BucketReader<VALUE>;
    friend class BucketRandomTour<VALUE, RANDOM_GENERATOR>::const_iterator;
};

}   // Archive

}   // RACES

#endif // __RACES_BUCKET__