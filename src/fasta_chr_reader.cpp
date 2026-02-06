/**
 * @file fasta_chr_reader.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes to read chromosomes from FASTA streams
 * @version 1.1
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

#include <sstream>

#include "fasta_chr_reader.hpp"

namespace CLONES
{

namespace IO
{

namespace FASTA
{

template<>
bool Reader<ChromosomeData<SequenceInfo>>::read(ChromosomeData<SequenceInfo>& chr_info,
                                                CLONES::UI::ProgressBar& progress_bar)
{
    std::string header;

    while (read(chr_info, nullptr, &header, progress_bar)) {
        if (is_chromosome_header(header, chr_info.chr_id)) {
            return true;
        }
    }

    return false;
}

template<>
bool Reader<ChromosomeData<Sequence>>::read(ChromosomeData<Sequence>& chr,
                                            CLONES::UI::ProgressBar& progress_bar)
{
    std::string header;

    while (read(chr, &(chr.nucleotides), &header, progress_bar)) {
        if (is_chromosome_header(header, chr.chr_id)) {
            return true;
        }
    }

    return false;
}

template<>
void Index<ChromosomeData<Sequence>>::save(std::ostream& output_stream) const
{
    bool first{true};
    for (const auto& map_it: _map) {
        const auto& entry = map_it.second;
        if (first) {
            first = false;
        } else {
            output_stream << std::endl;
        }
        output_stream << map_it.first
                        << "\t" << entry.name
                        << "\t" << entry.length
                        << "\t" << entry.offset
                        << "\t" << entry.linebases
                        << "\t" << entry.linebytes;
    }
}

template<>
Index<ChromosomeData<Sequence>>
Index<ChromosomeData<Sequence>>::load(std::istream& input_stream)
{
    Index index;

    std::string line;
    while (std::getline(input_stream, line)) {
        std::stringstream ss(line);

        IndexEntry entry;
        std::string chr_name;

        ss >> chr_name >> entry.name >> entry.length >> entry.offset
           >> entry.linebases >> entry.linebytes;

        index._map.emplace(chr_name, entry);
    }

    return index;
}

template<>
bool IndexedReader<ChromosomeData<Sequence>>::read(ChromosomeData<Sequence>& chr,
                                                   const std::string& chr_name,
                                                   CLONES::UI::ProgressBar& progress_bar)
{
    if (!_index.has_key(chr_name)) {
        return false;
    }

    if (read(chr, _index[chr_name], progress_bar)) {
        chr.chr_id = Mutations::GenomicPosition::stochr(chr_name);

        return true;
    }

    return false;
}

}   // FASTA

}   // IO

}   // CLONES
