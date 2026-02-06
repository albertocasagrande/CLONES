/**
 * @file csv_reader.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to read CSV
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

#ifndef __CLONES_CSV_READER__
#define __CLONES_CSV_READER__

#include <map>
#include <vector>
#include <fstream>
#include <filesystem>

namespace CLONES
{

namespace IO
{

/**
 * @brief A class to read CSV
 */
class CSVReader
{
    const std::filesystem::path filename;  //!< The CSV filename

    const char col_sep;             //!< The column separator

    std::map<std::string, size_t> columns; //!< The column positions

    std::vector<std::string> header;    //!< The CSV header

    std::streampos first_row_pos;   //!< The position in the file of the first data row

public:
    /**
     * @brief A class to iterate over CSV rows
     */
    class const_iterator;

    /**
     * @brief A class representing a CSV row
     */
    class CSVRow
    {
        std::vector<std::string> fields;  //!< The row fields

        CSVRow(const std::string& row, const char& column_separator=',');
    public:
        /**
         * @brief The empty constructor
         */
        CSVRow();

        /**
         * @brief Get a field by index
         *
         * @param index is the index of the field
         * @return the value of the field
         */
        const std::string& get_field(const size_t& index) const;

        /**
         * @brief Get the vector of the row fields
         *
         * @return a constant reference to the vector of the fields
         */
        inline const std::vector<std::string>& get_fields() const
        {
            return fields;
        }

        inline size_t size() const
        {
            return fields.size();
        }

        friend class const_iterator;
    };

    class const_iterator
    {
        const CSVReader *reader;    //!< A pointer to the reader
        std::ifstream ifs;          //!< The CSV file stream
        CSVRow curr_row;            //!< The row of the current iterator
        bool eof;                   //!< A Boolean flag to denote when eof has been reached

        /**
         * @brief Read a CSV row
         *
         * @param row is the object in which the read row must be stored
         * @return `true` if and only if the operation succeeds
         */
        bool get_row(CSVRow& row);

        /**
         * @brief The constructor
         *
         * @param reader is a constant pointer to the reader
         * @param off is the offset from the `way` position
         * @param way is the file position from which the offset `off`
         *          must be accounted
         */
        const_iterator(const CSVReader* reader, const std::streamoff& off=0,
                       const std::ios_base::seekdir& way=std::ios_base::beg);

    public:

        /**
         * @brief Move to the next row
         *
         * @return a reference to the updated object
         */
        const_iterator& operator++();

        /**
         * @brief Dereference the iterator
         *
         * @return a constant reference to the referenced CSV row
         */
        inline const CSVRow& operator*() const
        {
            if (eof) {
                throw std::runtime_error("No rows available: EOF have been reached.");
            }
            return curr_row;
        }

        /**
         * @brief Get the pointer to the object referenced by the iterator
         *
         * @return a constant pointer to the referenced CSV row
         */
        inline const CSVRow* operator->() const
        {
            if (eof) {
                throw std::runtime_error("No rows available: EOF have been reached.");
            }
            return &curr_row;
        }

        /**
         * @brief Check whether two iterator refer to the same row
         *
         * @param other is a constant iterator over a CSV reader
         * @return `true` if and only if the current iterator and
         *           `other` are referring to the same CSV row
         */
        inline bool operator==(const const_iterator& other) const
        {
            return reader == other.reader && eof == other.eof
                    && const_cast<std::ifstream&>(ifs).tellg()
                        == const_cast<std::ifstream&>(other.ifs).tellg();
        }

        /**
         * @brief Check whether two iterator refer to different rows
         *
         * @param other is a constant iterator over a CSV reader
         * @return `false` if and only if the current iterator and
         *           `other` are referring to the same CSV row
         */
        inline bool operator!=(const const_iterator& other) const
        {
            return !(*this == other);
        }

        /**
         * @brief Get the reading position on the CSV file stream
         *
         * @return the reading position on the CSV file stream
         */
        inline std::streampos tellg()
        {
            return ifs.tellg();
        }

        friend class CSVReader;
    };

    /**
     * @brief A constructor
     *
     * @param filename is the path of the CSV file to be read
     * @param has_header is a Boolean flag to enable/disable header in CSV file
     * @param column_separator is the character used to separate columns
     */
    CSVReader(const std::filesystem::path& filename, const bool has_header=true,
              const char column_separator=',');

    /**
     * @brief A constructor
     *
     * This constructor assumes that the CSV file has a header
     *
     * @param filename is the path of the CSV file to be read
     * @param column_separator is the character used to separate columns
     */
    inline CSVReader(const std::filesystem::path& filename,
                     const char column_separator):
        CSVReader(filename, true, column_separator)
    {}

    /**
     * @brief Get an iterator referring to the first data row in the CSV
     *
     * @return an iterator referring to the first data row in the CSV
     */
    inline const_iterator begin() const
    {
        return const_iterator(this, first_row_pos);
    }

    /**
     * @brief Get an iterator referring to the end of the CSV
     *
     * @return an iterator referring to the end of the CSV
     */
    inline const_iterator end() const
    {
        return const_iterator(this, 0, std::ios_base::end);
    }

    /**
     * @brief Get the CSV header vector
     *
     * @return the CSV header vector
     */
    inline const std::vector<std::string>& get_header() const
    {
        return header;
    }

    /**
     * @brief Get the position of a column by name
     *
     * @param column_name is the name of the column whose position is request
     * @return the position of the column whose name is `column_name`
     */
    inline const size_t& get_column_position(const std::string& column_name) const
    {
        auto found = columns.find(column_name);

        if (found == columns.end()) {
            throw std::domain_error("Unknown column name \""
                                    + column_name + "\".");
        }

        return found->second;
    }

    /**
     * @brief Get the file size
     *
     * @return the file size
     */
    inline std::streampos size() const
    {
        return end().tellg();
    }

    friend class const_iterator;
};

}   // IO

}   // CLONES

#endif // __CLONES_CSV_READER__
