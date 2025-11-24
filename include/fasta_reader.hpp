/**
 * @file fasta_reader.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a FASTA file reader and support structures
 * @version 1.5
 * @date 2025-11-24
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

#ifndef __RACES_FASTA_READER__
#define __RACES_FASTA_READER__

#include <string>
#include <istream>
#include <algorithm>

#include <map>
#include <fstream>
#include <filesystem>
#include <sstream>

#include "utils.hpp"
#include "basic_IO.hpp"

#include "progress_bar.hpp"

namespace RACES
{

namespace IO
{

/**
 * @brief A name space for FASTA file format
 */
namespace FASTA
{

/**
 * @brief A class to represent FASTA sequence information
 */
struct SequenceInfo
{
    std::string name;   //!< The name of the FASTA sequence
    size_t length;      //!< The length of the FASTA sequence

    /**
     * @brief Get the identifier from the sequence header
     *
     * @param header is the header of the sequence
     * @return the identifier associated to the header name
     */
    static inline std::string get_id(const std::string& header)
    {
        return header.substr(1, header.find(" ")-1);
    }

    /**
     * @brief Check whether the header is valid
     *
     * @param header is the header of the sequence
     * @return `true` if and only if `header` is a valid
     *      sequence header
     */
    static inline bool is_valid(const std::string& header)
    {
        if (header.size() == 0) {
            return false;
        }

        return header[0] == '>';
    }
};

/**
 * @brief A class to represent FASTA sequence
 */
struct Sequence : public SequenceInfo
{
    std::string nucleotides;    //!< The string of sequence nucleotides
};

/**
 * @brief A class to read sequence information from a FASTA file
 *
 * @tparam READ_TYPE is the type of information read from the FASTA file
 */
template<typename READ_TYPE, std::enable_if_t<std::is_base_of_v<SequenceInfo, READ_TYPE>,bool> = true>
class Reader
{
protected:

    std::ifstream fasta_stream; //!< The FASTA file stream

    /**
     * @brief Read the next content of the stream up-to a new sequence
     *
     * @param progress_bar is a progress bar
     */
    void filter_remaining_sequence(RACES::UI::ProgressBar& progress_bar)
    {
        size_t counter{0};
        int c = fasta_stream.get();
        while (c != EOF && c != '>') {
            std::string line;
            std::getline(fasta_stream, line);

            if ((counter = (counter+1)%1000) == 0) {
                // let time elapse in progress bar
                progress_bar.set_progress(progress_bar.get_progress());
            }

            c = fasta_stream.get();
        }

        if (c == '>') {
            fasta_stream.unget();
        }
    }

    /**
     * @brief Read FASTA sequence nucleotide from a stream
     *
     * @param[out] nucleotides is a pointer to the string that will be filled by the
     *          sequence nucleotides
     * @param[in,out] progress_bar is a progress bar
     * @return The string of nucleotides up to the end of the stream or of the
     *          current sequence
     */
    inline size_t read_nucleotides(std::string* nucleotides, RACES::UI::ProgressBar& progress_bar)
    {
        return read_nucleotides(nucleotides, std::string().max_size(), progress_bar);
    }

    /**
     * @brief Read a nucleotide sequence from a FALSE stream
     *
     * @param[out] nucleotides is a pointer to the string that will be filled by the
     *          nucleotide sequence
     * @param[in] length is the length of the nucleotide sequence to be read
     * @param[in,out] progress_bar is a progress bar
     * @return The string of nucleotides up to the end of the stream or of the
     *          current sequence
     */
    size_t read_nucleotides(std::string* nucleotides, const size_t length,
                            RACES::UI::ProgressBar& progress_bar)
    {
        if (nucleotides != nullptr) {
            nucleotides->clear();
        }

        size_t counter{0}, read_length{0};
        int c;
        while ((c = fasta_stream.get()) != EOF) {
            std::string line;
            fasta_stream.unget();
            if (c == '>') {
                return read_length;
            }

            std::getline(fasta_stream, line);

            // let time elapse in progress bar
            if ((counter = (counter+1)%1000) == 0) {
                progress_bar.set_progress(progress_bar.get_progress());
            }

            line.erase(remove(line.begin(), line.end(), ' '), line.end());
            line.erase(remove(line.begin(), line.end(), '\t'), line.end());
            if (read_length + line.size() > length) {
                line.resize(length-read_length);
            }

            read_length += line.size();
            if (nucleotides != nullptr) {
                nucleotides->append(line);
            }

            if (read_length==length) {
                return length;
            }
        }

        return read_length;
    }

    /**
     * @brief Read FASTA sequence information from a stream
     *
     * @tparam FILTER is the type of sequence filter
     * @param[out] seq_info is the object that will be filled by the read information
     * @param[out] nucleotides is a pointer to the string that will be filled by the
     *          sequence nucleotides
     * @param[out] header is a pointer the string that will be filled by the sequence
     *          header
     * @param[in,out] progress_bar is a progress bar
     * @return `true` if and only if a sequence has been read from `fasta_stream`.
     *          If the method returns `true`, then `seq_info` and, whenever
     *          `nucleotides` or `header` are not `nullptr`, `*nucleotides` and
     *          `*hearer` are updated according to the read values
     */
    bool read(SequenceInfo& seq_info, std::string* nucleotides, std::string* header,
              RACES::UI::ProgressBar& progress_bar)
    {
        int c = fasta_stream.get();
        while (c != EOF) {
            fasta_stream.unget();
            if (c!='>') {
                std::string char_s(">");

                char_s[0] = c;
                throw std::runtime_error("expecting '>' got '"+char_s+"'");
            }
            std::string read_header;
            std::getline(fasta_stream, read_header);

            if (READ_TYPE::is_valid(read_header)) {
                auto seq_name = READ_TYPE::get_id(read_header);

                progress_bar.set_message("Reading " + seq_name);

                if (header != nullptr) {
                    *header = read_header;
                }

                seq_info.name = seq_name;
                seq_info.length = read_nucleotides(nucleotides, progress_bar);

                return true;
            } else {
                filter_remaining_sequence(progress_bar);

                c = fasta_stream.get();
            }
        }

        return false;
    }

public:

    /**
     * @brief The empty constructor
     */
    Reader()
    {}

    /**
     * @brief Build a new information reader
     *
     * @param fasta_filename is the name of the FASTA file to read
     */
    explicit Reader(const std::filesystem::path& fasta_filename)
    {
        open(fasta_filename);
    }

    /**
     * @brief Open a new FASTA file
     *
     * @param fasta_filename is the name of the FASTA file to open
     */
    void open(const std::filesystem::path& fasta_filename)
    {
        fasta_stream.open(fasta_filename, std::ios_base::in);

        if (!fasta_stream.good()) {
            std::ostringstream oss;

            oss << "\"" << to_string(fasta_filename) << "\" does not exist";
            throw std::runtime_error(oss.str());
        }
    }

    /**
     * @brief Close the current FASTA file
     */
    inline void close()
    {
        fasta_stream.close();
    }

    /**
     * @brief Read some data
     *
     * @param[out] read_object is a reference to an object that will be filled by the
     *          read data
     * @param[in,out] progress_bar is a progress bar
     * @return `true` if and only if a sequence has been read from `fasta_stream`.
     *           If the method returns `true`, then `read_object` is updated
     *           according to the read values
     */
    bool read(READ_TYPE& read_object, RACES::UI::ProgressBar& progress_bar)
    {
        return read(read_object, &(read_object.nucleotides), nullptr, progress_bar);
    }

    /**
     * @brief Read FASTA sequence information from a stream
     *
     * @tparam FILTER is the type of sequence filter
     * @param[out] read_object is a reference to an object that will be filled by the
     *          read data
     * @param progress_bar_stream is the output stream for the progress bar
     * @return `true` if and only if a sequence has been read from `fasta_stream`.
     *           If the method returns `true`, then `read_object` is updated
     *           according to the read values
     */
    bool read(READ_TYPE& read_object, std::ostream& progress_bar_stream)
    {
        RACES::UI::ProgressBar progress_bar(progress_bar_stream);

        return read(read_object, progress_bar);
    }

    /**
     * @brief Read FASTA sequence information from a stream
     *
     * @tparam FILTER is the type of sequence filter
     * @param[out] read_object is a reference to an object that will be filled by the
     *          read data
     * @return `true` if and only if a sequence has been read from `fasta_stream`.
     *           If the method returns `true`, then `read_object` is updated
     *           according to the read values
     */
    bool read(READ_TYPE& read_object)
    {
        RACES::UI::ProgressBar progress_bar(std::cout, true);

        return read(read_object, progress_bar);
    }

    /**
     * @brief Get the stream size
     *
     * @return The size of the stream size
     */
    inline size_t get_stream_size()
    {
        return IO::get_stream_size(fasta_stream);
    }

    /**
     * @brief Get the stream size
     *
     * @return The size of the stream size
     */
    inline std::streampos get_position()
    {
        return fasta_stream.tellg();
    }
};

template<>
inline bool Reader<SequenceInfo>::read(SequenceInfo& sequence, RACES::UI::ProgressBar& progress_bar)
{
    return read(sequence, nullptr, nullptr, progress_bar);
}

/**
 * @brief A sequence description in the index
 */
struct IndexEntry
{
    std::string name;   //!< Sequence name
    size_t length;      //!< Number of sequence bases
    size_t offset;      //!< File offset to the first base
    unsigned int linebases; //!< Number of bases in each line
    unsigned int linebytes; //!< Number of bytes per line

    /**
     * @brief The empty constructor
     */
    IndexEntry():
        name(""), length(0), offset(0), linebases(0), linebytes(0)
    {}

    /**
     * @brief Turn a base offset into the byte offset
     *
     * @param base_offset is the base offset in the sequence
     * @return the byte offset corresponding to `base_offset`
     */
    inline size_t get_byte_offset(const size_t base_offset) const
    {
        return ((base_offset / linebases) * linebytes
                    + base_offset % linebases);
    }

    /**
     * @brief Turn a byte offset into the base offset
     *
     * @param byte_offset is the byte offset in the sequence
     * @return the base offset corresponding to `byte_offset`
     */
    inline size_t get_base_offset(const size_t byte_offset) const
    {
        return ((byte_offset / linebytes) * linebases
                    + byte_offset % linebytes);
    }
};

/**
 * @brief The FASTA index class
 */
template<typename SEQUENCE_TYPE, std::enable_if_t<std::is_base_of_v<Sequence, SEQUENCE_TYPE>,bool> = true>
class Index
{
    std::map<std::string, IndexEntry> _map; //!< the key-entry map

public:
    /**
     * @brief Get the index file extension
     *
     * @return the index file extension
     */
    inline static const char* index_extension()
    {
        return ".fai";
    }

    /**
     * @brief Get the index filename
     *
     * @param fasta_filename is the FASTA filename
     * @return the index filename corresponding to `fasta_filename`
     */
    static std::filesystem::path
    get_index_filename(const std::filesystem::path& fasta_filename)
    {
        std::filesystem::path index_filename{fasta_filename};

        std::string ext = to_string(index_filename.extension())
                          + Index<SEQUENCE_TYPE>::index_extension();

        index_filename.replace_extension(ext);

        return index_filename;
    }

    /**
     * @brief The empty constructor
     *
     * This method builds an empty FASTA index
     */
    Index():
        _map()
    {}

    /**
     * @brief Check there the index accounts for a key
     *
     * @param key is the key whose entry is checked
     * @return `true` if and only if the index contains an
     *      entry whose key is `key`
     */
    inline bool has_key(const std::string& key) const
    {
        return _map.find(key) != _map.end();
    }

    /**
     * @brief Get the iterator to the begin of the index
     *
     * @return a constant iterator to the first element
     *      in the key-entry map
     */
    inline typename std::map<std::string, IndexEntry>::const_iterator
    begin() const
    {
        return _map.begin();
    }

    /**
     * @brief Get the iterator to the end of the index
     *
     * @return a constant iterator to the end of the
     *      key-entry map
     */
    inline typename std::map<std::string, IndexEntry>::const_iterator
    end() const
    {
        return _map.end();
    }

    /**
     * @brief Get an entry
     *
     * @param key is the key whose entry aimed
     * @return a constant reference to the index entry associated
     *      to `key`
     */
    const IndexEntry& operator[](const std::string& key) const
    {
        auto found = _map.find(key);

        if (found == _map.end()) {
            throw std::domain_error("The index has not entry having key \""
                                    + key +"\"");
        }
        return found->second;
    }

    /**
     * @brief Build a FASTA index
     *
     * @param fasta_filename is a FASTA file path
     * @param progress_bar is the progress bar
     * @return the FASTA index of `fasta_filename`
     */
    static Index build_index(const std::filesystem::path& fasta_filename,
                             RACES::UI::ProgressBar& progress_bar)
    {
        namespace fs = std::filesystem;

        Index index;

        size_t counter{0};

        std::ifstream fasta_stream(fasta_filename);
        const size_t file_size = IO::get_stream_size(fasta_stream);

        int c = fasta_stream.get();
        while (c != EOF) {
            fasta_stream.unget();

            std::string header;
            std::getline(fasta_stream, header);

            std::string seq_name;
            if (header.size()>0) {
                seq_name = Sequence::get_id(header);
                progress_bar.set_message("Found " + seq_name);
            }

            if (SEQUENCE_TYPE::is_valid(header)) {
                const auto key = SEQUENCE_TYPE::get_id(header);
                progress_bar.set_message("Reading " + key);

                IndexEntry entry;
                entry.name = seq_name;
                entry.offset = fasta_stream.tellg();

                c = fasta_stream.get();
                size_t last_progress{0}, line_number{0};
                bool byte_diff{false}, base_diff{false};
                while (c != EOF && c != '>') {
                    std::string line;
                    if (byte_diff) {
                        throw std::domain_error("\"" + key + "("+ entry.name
                                                + ")\"'s line "
                                                + std::to_string(line_number)
                                                + " differs from the previous "
                                                + "one by number of bytes.");
                    }

                    if (base_diff) {
                        throw std::domain_error("\"" + key + "("+ entry.name
                                                + ")\"'s line "
                                                + std::to_string(line_number)
                                                + " differs from the previous "
                                                + "one by length.");
                    }
                    fasta_stream.unget();

                    size_t first_line_pos = fasta_stream.tellg();
                    std::getline(fasta_stream, line);
                    size_t next_line_pos = fasta_stream.tellg();

                    if (entry.length==0) {
                        entry.linebytes = next_line_pos-first_line_pos;
                    } else {
                        byte_diff = (entry.linebytes != next_line_pos-first_line_pos);
                    }

                    if ((100*next_line_pos)/file_size > last_progress) {
                        last_progress = (100*next_line_pos)/file_size;
                        progress_bar.set_progress(last_progress);
                    } else {
                        // let time elapse in progress bar
                        if ((counter = (counter+1)%1000) == 0) {
                            progress_bar.set_progress(last_progress);
                        }
                    }

                    line.erase(remove(line.begin(), line.end(), ' '), line.end());
                    line.erase(remove(line.begin(), line.end(), '\t'), line.end());

                    if (entry.length==0) {
                        entry.linebases = line.size();
                    } else {
                        base_diff = (entry.linebases != line.size());
                    }

                    entry.length += line.size();
                    ++line_number;

                    c = fasta_stream.get();
                }

                index._map.emplace(key, entry);
            } else {
                c = fasta_stream.get();
                while (c != EOF && c != '>') {
                    std::string line;
                    std::getline(fasta_stream, line);

                    if ((counter = (counter+1)%1000) == 0) {
                        progress_bar.set_progress((100*fasta_stream.tellg())/file_size);
                    }

                    c = fasta_stream.get();
                }
            }
        }

        progress_bar.set_progress(100, "done");

        return index;
    }

    /**
     * @brief Build a FASTA index
     *
     * @param fasta_filename is a FASTA file path
     * @param progress_bar_stream is the progress bar stream
     * @return the FASTA index of `fasta_filename`
     */
    static Index build_index(const std::filesystem::path& fasta_filename,
                             std::ostream& progress_bar_stream)
    {
        RACES::UI::ProgressBar progress_bar(progress_bar_stream);

        return Index::build_index(fasta_filename, progress_bar);
    }

    /**
     * @brief Build a FASTA index
     *
     * @param fasta_filename is a FASTA file path
     * @return the FASTA index of `fasta_filename`
     */
    static Index build_index(const std::filesystem::path& fasta_filename)
    {
        RACES::UI::ProgressBar progress_bar(std::cout, true);

        return Index::build_index(fasta_filename, progress_bar);
    }

    /**
     * @brief Sort a container of sequence names according to their index offset
     *
     * @tparam CONTAINER is the type of the input container
     * @param sequence_names is a container of sequence names
     * @return a list of the elements in `sequence_names` sorted according to
     *      their positions in the index
     */
    template <typename CONTAINER, std::enable_if_t<std::is_same_v<typename CONTAINER::value_type, std::string>, bool> = true>
    std::list<std::string> sort_according_offset(CONTAINER&& sequence_names) const
    {
        return sort_according_offset(sequence_names);
    }

    /**
     * @brief Sort a container of sequence names according to their index offset
     *
     * @tparam CONTAINER is the type of the input container
     * @param sequence_names is a container of sequence names
     * @return a list of the elements in `sequence_names` sorted according to
     *      their positions in the index
     */
    template <typename CONTAINER, std::enable_if_t<std::is_same_v<typename CONTAINER::value_type, std::string>, bool> = true>
    std::list<std::string> sort_according_offset(const CONTAINER& sequence_names) const
    {
        std::map<size_t, std::string> offset_name_map;

        for (const auto& name : sequence_names) {
            auto found = _map.find(name);

            if (found == _map.end()) {
                throw std::domain_error("Unknown sequence \"" + name + "\"");
            }

            offset_name_map.emplace(found->second.offset, name);
        }

        std::list<std::string> sorted;

        for (const auto& [offset, name]: offset_name_map) {
            sorted.push_back(name);
        }

        return sorted;
    }

    /**
     * @brief Save a FASTA index
     *
     * @param output_stream is the output stream
     */
    void save(std::ostream& output_stream) const
    {
        bool first{true};
        for (const auto& map_it: _map) {
            const auto& entry = map_it.second;
            if (first) {
                first = false;
            } else {
                output_stream << std::endl;
            }
            output_stream << entry.name
                            << "\t" << entry.length
                            << "\t" << entry.offset
                            << "\t" << entry.linebases
                            << "\t" << entry.linebytes;
        }
    }

    /**
     * @brief Load a FASTA index
     *
     * @param input_stream is the output stream
     * @return the FASTA index
     */
    static Index<SEQUENCE_TYPE> load(std::istream& input_stream)
    {
        Index index;

        std::string line;
        while (std::getline(input_stream, line)) {
            std::stringstream ss(line);

            IndexEntry entry;

            ss >> entry.name >> entry.length >> entry.offset
               >> entry.linebases >> entry.linebytes;

            index._map.emplace(entry.name, entry);
        }

        return index;
    }
};

/**
 * @brief A class to read sequences from FASTA files
 *
 * This class exploits FASTA fai file to speed-up sequence reading.
 */
template<typename SEQUENCE_TYPE, std::enable_if_t<std::is_base_of_v<Sequence, SEQUENCE_TYPE>,bool> = true>
class IndexedReader : public Reader<SEQUENCE_TYPE>
{
    Index<SEQUENCE_TYPE> _index; //!< The index

    /**
     * @brief Read data
     *
     * @param[out] sequence is a reference to sequence to be filled by read
     *           data
     * @param[in] entry is the index entry of the aimed sequence
     * @param[in, out] progress_bar is the progress bar
     * @return `true` if and only if a sequence having name `sequence_name`
     *           can be read from the FASTA file. If the method returns `true`,
     *           then `read_object` is updated according to the read object
     */
    bool read(SEQUENCE_TYPE& sequence, const IndexEntry& entry,
              RACES::UI::ProgressBar& progress_bar)
    {
        sequence.length = entry.length;
        sequence.name = entry.name;

        // clear the error flag if EOF was reached
        this->fasta_stream.clear();

        this->fasta_stream.seekg(entry.offset);
        const size_t length = this->read_nucleotides(&(sequence.nucleotides),
                                                     progress_bar);

        if (length != entry.length) {
            throw std::runtime_error("The index length (" + std::to_string(entry.length)
                                     + ") and the actual length ("
                                     + std::to_string(length) + ") of \""
                                     + entry.name + "\" do not match");
        }

        return true;
    }
public:
    /**
     * @brief The empty constructor
     */
    IndexedReader():
        Reader<SEQUENCE_TYPE>()
    {}

    /**
     * @brief Build a new sequence reader
     *
     * @param progress_bar is the progress bar
     * @param fasta_filename is the filename of the FASTA to read
     */
    IndexedReader(const std::filesystem::path& fasta_filename,
                  RACES::UI::ProgressBar& progress_bar)
    {
        open(fasta_filename, progress_bar);
    }

    /**
     * @brief Build a new sequence reader
     *
     * @param fasta_filename is the filename of the FASTA to read
     */
    explicit IndexedReader(const std::filesystem::path& fasta_filename)
    {
        open(fasta_filename);
    }

    /**
     * @brief Get the index
     *
     * @return a constant reference to the index
     */
    inline const Index<SEQUENCE_TYPE>& get_index() const
    {
        return _index;
    }

    /**
     * @brief Open a new FASTA file
     *
     * @param fasta_filename is the name of the FASTA file to open
     * @param progress_bar is the progress bar
     */
    void open(const std::filesystem::path& fasta_filename,
              RACES::UI::ProgressBar& progress_bar)
    {
        const auto index_filename = _index.get_index_filename(fasta_filename);

        namespace fs = std::filesystem;
        if (fs::exists(index_filename)) {
            std::ifstream input_stream(index_filename);

            _index = Index<SEQUENCE_TYPE>::load(input_stream);
        } else {
            _index = Index<SEQUENCE_TYPE>::build_index(fasta_filename, progress_bar);

            std::ofstream output_stream(index_filename);

            _index.save(output_stream);
        }

        Reader<SEQUENCE_TYPE>::open(fasta_filename);
    }

    /**
     * @brief Open a new FASTA file
     *
     * @param fasta_filename is the name of the FASTA file to open
     */
    void open(const std::filesystem::path& fasta_filename)
    {
        RACES::UI::ProgressBar progress_bar(std::cout, true);

        open(fasta_filename, progress_bar);
    }

    /**
     * @brief Read data
     *
     * @param[out] sequence is a reference to sequence to be filled by read
     *           data
     * @param[in] sequence_name is the name of the sequence to be read
     * @param[in, out] progress_bar is the progress bar
     * @return `true` if and only if a sequence having name `sequence_name`
     *           can be read from the FASTA file. If the method returns `true`,
     *           then `read_object` is updated according to the read object
     */
    bool read(SEQUENCE_TYPE& sequence, const std::string& sequence_name,
              RACES::UI::ProgressBar& progress_bar)
    {
        if (!_index.has_key(sequence_name)) {
            return false;
        }

        return read(sequence, _index[sequence_name], progress_bar);
    }

    /**
     * @brief Read data
     *
     * @param[out] sequence is a reference to sequence to be filled by read
     *           data
     * @param[in] sequence_name is the name of the sequence to be read
     * @return `true` if and only if a sequence having name `sequence_name`
     *           can be read from the FASTA file. If the method returns `true`,
     *           then `read_object` is updated according to the read object
     */
    bool read(SEQUENCE_TYPE& sequence, const std::string& sequence_name)
    {
        RACES::UI::ProgressBar progress_bar(std::cout, true);

        const auto seq_read = read(sequence, sequence_name, progress_bar);

        progress_bar.set_message(sequence_name + " read")
                    .set_progress(100);

        return seq_read;
    }

    /**
     * @brief Read data
     *
     * @param[out] nucleotides is the string in which the read sequence
     *          fragment will be placed
     * @param[in] sequence_name is the name of the sequence to be read
     * @param[in] offset is the base offset from the beginning of the sequence
     * @param[in] length is the length of the fragment to be read
     * @param[in, out] progress_bar is the progress bar
     * @return `true` if and only if a sequence having name `sequence_name`
     *           and whose length is longer than `offset`+`length` can be
     *           read from the FASTA file. If the method returns `true`,
     *           then `nucleotides` is updated according to the read bases
     */
    bool read(std::string& nucleotides, const std::string& sequence_name,
              const size_t offset, const size_t length,
              RACES::UI::ProgressBar& progress_bar)
    {
        if (!_index.has_key(sequence_name)) {
            return false;
        }
        const auto& entry = _index[sequence_name];

        // clear the error flag if EOF was reached
        this->fasta_stream.clear();

        this->fasta_stream.seekg(entry.offset + entry.get_byte_offset(offset));

        this->read_nucleotides(&nucleotides, length, progress_bar);

        return true;
    }

    /**
     * @brief Read data
     *
     * @param[out] nucleotides is the string in which the read sequence
     *          fragment will be placed
     * @param[in] sequence_name is the name of the sequence to be read
     * @param[in] offset is the base offset from the beginning of the sequence
     * @param[in] length is the length of the fragment to be read
     * @return `true` if and only if a sequence having name `sequence_name`
     *           and whose length is longer than `offset`+`length` can be
     *           read from the FASTA file. If the method returns `true`,
     *           then `nucleotides` is updated according to the read bases
     */
    bool read(std::string& nucleotides, const std::string& sequence_name,
              const size_t offset, const size_t length)
    {
        RACES::UI::ProgressBar progress_bar(std::cout, true);

        return read(nucleotides, sequence_name, offset, length, progress_bar);
    }
};

}   // FASTA

}   // IO

}   // RACES

#endif // __RACES_FASTA_READER__