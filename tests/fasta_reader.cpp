/**
 * @file fasta_reader.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing CLONES::IO::FASTA::Reader class
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fasta_reader

#include <boost/test/unit_test.hpp>

#include <fstream>

#include "fasta_chr_reader.hpp"


BOOST_AUTO_TEST_CASE(reader_creation)
{
    using namespace CLONES::IO::FASTA;

    BOOST_CHECK_NO_THROW(Sequence sequence);
    BOOST_CHECK_NO_THROW(SequenceInfo seq_info);
}

struct NucleotideFixture
{
    std::string seq_name;
    size_t offset;
    size_t length;
    std::string nucleotides;
};

struct FASTAFixture
{
    std::list<std::pair<std::string, std::string>> sequences;
    std::list<NucleotideFixture> read_nucleotides;
    std::list<std::pair<std::string, std::string>> chromosomes;
    std::list<NucleotideFixture> chr_read_nucleotides;

    FASTAFixture():
        sequences({
        {"1", "AACCCTAACCCTAACCCTA"},
        {"7", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"},
        {"2", "ACCCTAACCCTAAGCCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"},
        {"Test3", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"},
        {"NC_2321.110", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"}
        }),
        read_nucleotides({
            {"7", 45, 8, "CTAACCCT"},
            {"7", 45, 14, "CTAACCCTAACCCT"},
            {"Test3", 45, 14, "CTAACCCTAACCCT"},
            {"Test3", 109, 20, "CCCTAA"},
            {"NC_2321.110", 109, 5, "CCCTA"},
            {"NC_2321.110", 109, 20, "CCCTAA"}
            }),
        chromosomes({
        {"7", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"},
        {"13", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"}
        }),
        chr_read_nucleotides({
            {"7", 45, 8, "CTAACCCT"},
            {"7", 45, 14, "CTAACCCTAACCCT"},
            {"13", 109, 5, "CCCTA"},
            {"13", 109, 20, "CCCTAA"}
            })
    {}
};

BOOST_FIXTURE_TEST_SUITE( FASTA_test, FASTAFixture )

BOOST_AUTO_TEST_CASE(read_sequence_infos)
{
    using namespace CLONES::IO::FASTA;

    Reader<SequenceInfo> fasta_reader(FASTA_FILE);
    SequenceInfo seq_info;
    for (const auto& [name, nucleotides]: sequences) {
        BOOST_CHECK(fasta_reader.read(seq_info));

        BOOST_CHECK_EQUAL(seq_info.name, name);
        BOOST_CHECK_EQUAL(seq_info.length, nucleotides.size());
    }

    BOOST_CHECK(!fasta_reader.read(seq_info));
}

BOOST_AUTO_TEST_CASE(read_sequences)
{
    using namespace CLONES::IO::FASTA;

    Reader<Sequence> fasta_reader(FASTA_FILE);

    Sequence sequence;
    for (const auto& [name, nucleotides]: sequences) {
        BOOST_CHECK(fasta_reader.read(sequence));

        BOOST_CHECK_EQUAL(sequence.name, name);
        BOOST_CHECK_EQUAL(sequence.length, nucleotides.size());
        BOOST_CHECK_EQUAL(sequence.nucleotides, nucleotides);
    }

    BOOST_CHECK(!fasta_reader.read(sequence));
}

BOOST_AUTO_TEST_CASE(index_reader)
{
    using namespace CLONES::IO::FASTA;

    std::filesystem::remove(Index<Sequence>::get_index_filename(FASTA_FILE));

    for (size_t i=0; i<2; ++i) {
        IndexedReader<Sequence> ireader(FASTA_FILE);

        Sequence sequence;
        for (auto it=sequences.rbegin(); it != sequences.rend(); ++it) {
            const auto& name = it->first;
            const auto& nucleotides = it->second;
            BOOST_CHECK(ireader.read(sequence, name));

            BOOST_CHECK_EQUAL(sequence.name, name);
            BOOST_CHECK_EQUAL(sequence.length, nucleotides.size());
            BOOST_CHECK_EQUAL(sequence.nucleotides, nucleotides);
        }
    }

    std::filesystem::remove(Index<Sequence>::get_index_filename(FASTA_FILE));
}

BOOST_AUTO_TEST_CASE(read_chromosomes)
{
    using namespace CLONES::IO::FASTA;

    Reader<ChromosomeData<Sequence>> fasta_reader(FASTA_FILE);

    ChromosomeData<Sequence> chromosome;
    for (const auto& [name, nucleotides]: chromosomes) {
        BOOST_CHECK(fasta_reader.read(chromosome));

        BOOST_CHECK_EQUAL(chromosome.name, name);
        BOOST_CHECK_EQUAL(chromosome.length, nucleotides.size());
        BOOST_CHECK_EQUAL(chromosome.nucleotides, nucleotides);
    }

    BOOST_CHECK(!fasta_reader.read(chromosome));
}

BOOST_AUTO_TEST_CASE(nucleotides_reader)
{
    using namespace CLONES::IO::FASTA;

    std::filesystem::remove(Index<Sequence>::get_index_filename(FASTA_FILE));

    for (size_t i=0; i<2; ++i) {
        IndexedReader<Sequence> ireader(FASTA_FILE);

        std::string nucleotides;
        for (auto it=read_nucleotides.rbegin(); it != read_nucleotides.rend(); ++it) {
            BOOST_CHECK(ireader.read(nucleotides, it->seq_name, it->offset, it->length));

            BOOST_CHECK_EQUAL(nucleotides, it->nucleotides);
        }
    }

    std::filesystem::remove(Index<Sequence>::get_index_filename(FASTA_FILE));
}

BOOST_AUTO_TEST_CASE(chr_index_reader)
{
    using namespace CLONES::IO::FASTA;

    std::filesystem::remove(Index<ChromosomeData<Sequence>>::get_index_filename(FASTA_FILE));

    for (size_t i=0; i<2; ++i) {
        IndexedReader<ChromosomeData<Sequence>> ireader(FASTA_FILE);

        ChromosomeData<Sequence> chr;
        for (auto it=chromosomes.rbegin(); it != chromosomes.rend(); ++it) {
            const auto& name = it->first;
            const auto& nucleotides = it->second;
            BOOST_CHECK(ireader.read(chr, name));

            BOOST_CHECK_EQUAL(CLONES::Mutations::GenomicPosition::chrtos(chr.chr_id), name);
            BOOST_CHECK_EQUAL(chr.length, nucleotides.size());
            BOOST_CHECK_EQUAL(chr.nucleotides, nucleotides);
        }
    }

    std::filesystem::remove(Index<ChromosomeData<Sequence>>::get_index_filename(FASTA_FILE));
}

BOOST_AUTO_TEST_CASE(build_index_reader_error)
{
    using namespace CLONES::IO::FASTA;

    BOOST_CHECK_NO_THROW(IndexedReader<Sequence> ireader);

    IndexedReader<Sequence> ireader;

    BOOST_CHECK_THROW(IndexedReader<Sequence>(FASTA_INDEX_ERR), std::domain_error);
    BOOST_CHECK_THROW(ireader.open(FASTA_INDEX_ERR), std::domain_error);
}

BOOST_AUTO_TEST_CASE(chr_nucleotides_reader)
{
    using namespace CLONES::IO::FASTA;

    std::filesystem::remove(Index<ChromosomeData<Sequence>>::get_index_filename(FASTA_FILE));

    for (size_t i=0; i<2; ++i) {
        IndexedReader<ChromosomeData<Sequence>> ireader(FASTA_FILE);

        std::string nucleotides;
        for (auto it=chr_read_nucleotides.rbegin(); it != chr_read_nucleotides.rend(); ++it) {
            BOOST_CHECK(ireader.read(nucleotides, it->seq_name, it->offset, it->length));

            BOOST_CHECK_EQUAL(nucleotides, it->nucleotides);
        }
    }

    std::filesystem::remove(Index<ChromosomeData<Sequence>>::get_index_filename(FASTA_FILE));
}

BOOST_AUTO_TEST_SUITE_END()
