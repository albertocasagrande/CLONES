/**
 * @file driver_storage.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements class to load and store driver mutations
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

#include "driver_storage.hpp"

#include "csv_reader.hpp"

namespace CLONES
{

namespace Mutations
{

DriverStorage::DriverStorage()
{}

std::map<SID, std::string> DriverStorage::get_reverse_map() const
{
    std::map<SID, std::string> reverse_map;

    for (const auto& [code, mutation_entry] : mutation_map) {
        reverse_map[mutation_entry.mutation] = code;
    }

    return reverse_map;
}

std::list<GenomicPosition> DriverStorage::get_mutation_positions() const
{
    std::list<GenomicPosition> genomic_positions;

    for (const auto& [code, mutation_entry] : mutation_map) {
        genomic_positions.push_back(mutation_entry.mutation);
    }

    return genomic_positions;
}

DriverStorage DriverStorage::load(const std::filesystem::path& filename)
{
    DriverStorage storage;

    CLONES::IO::CSVReader csv_reader(filename, true, '\t');

    const size_t chr_col = csv_reader.get_column_position("chr");
    const size_t pos_col = csv_reader.get_column_position("from");
    const size_t alt_col = csv_reader.get_column_position("alt");
    const size_t ref_col = csv_reader.get_column_position("ref");
    const size_t type_col = csv_reader.get_column_position("mutation_type");
    const size_t code_col = csv_reader.get_column_position("driver_code");
    const size_t ttype_col = csv_reader.get_column_position("tumour_type");

    for (const auto& row : csv_reader) {
        const auto type = row.get_field(type_col);
        if (type=="SNV" || type=="indel") {
            const auto chr_str = row.get_field(chr_col);

            const auto chr_id = GenomicPosition::stochr(chr_str);
            const auto pos = static_cast<ChrPosition>(stoul(row.get_field(pos_col)));

            const SID mutation{chr_id, pos, row.get_field(ref_col),
                               row.get_field(alt_col), Mutations::Mutation::DRIVER};

            const std::string code{row.get_field(code_col)};
            const std::string ttype{row.get_field(ttype_col)};

            auto found = storage.mutation_map.find(code);
            if (found == storage.mutation_map.end()) {
                storage.mutation_map.insert({code, {std::move(mutation), {ttype}}});
            } else {
                if (found->second.mutation != mutation) {
                    std::ostringstream oss;

                    oss << "The code \"" << code << "\" is associated to both "
                        << mutation << " and " << found->second.mutation
                        << std::endl;

                    throw std::runtime_error(oss.str());
                }

                found->second.tumour_types.insert(ttype);
            }
        }
    }

    storage.source_path = filename;

    return storage;
}

}   // Mutations

}   // CLONES
