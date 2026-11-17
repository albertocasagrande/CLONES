/**
 * @file genome_mutations.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing CLONES::Mutations::GenomeMutations class
 * @version 1.4
 * @date 2026-11-17
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
#define BOOST_TEST_MODULE genome_mutations

#include <boost/test/unit_test.hpp>

#include <fstream>

#include "genome_mutations.hpp"
#include "germline.hpp"
#include "driver_storage.hpp"
#include "context_index.hpp"
#include "fasta_chr_reader.hpp"
#include "utils.hpp"

#include <map>
#include <string>
#include <cctype>
#include <filesystem>

#include "error.hpp"

BOOST_AUTO_TEST_CASE(chromosome_mutations_creation)
{
    using namespace CLONES::Mutations;

    BOOST_CHECK_NO_THROW(ChromosomeMutations());
    BOOST_CHECK_NO_THROW(ChromosomeMutations(1, 100, 3));
    BOOST_CHECK_NO_THROW(ChromosomeMutations({{1, 1}, 100}, 3));
}

BOOST_AUTO_TEST_CASE(genome_mutations_creation)
{
    using namespace CLONES::Mutations;

    BOOST_CHECK_NO_THROW(GenomeMutations());
    BOOST_CHECK_NO_THROW(GenomeMutations({{1, 100, 3}, {2, 100, 2}}));

    auto index = ContextIndex<>::build_index(FASTA_FILE);

    auto chr_regions = index.get_chromosome_regions();

    BOOST_CHECK_NO_THROW(GenomeMutations(chr_regions, 2));

    BOOST_CHECK_THROW(GenomeMutations({{1, 100, 3}, {1, 100, 2}}),
                      CLONES::Error<std::domain_error>);
}

struct GenomeMutationsFixture
{
    CLONES::Mutations::GenomeMutations genome_mutations;

    GenomeMutationsFixture()
    {
        using namespace CLONES::Mutations;

        auto index = ContextIndex<>::build_index(FASTA_FILE);

        auto chr_regions = index.get_chromosome_regions();

        genome_mutations = GenomeMutations(chr_regions, 2);
    }
};

BOOST_FIXTURE_TEST_SUITE( genome_mutation_tests, GenomeMutationsFixture )

BOOST_AUTO_TEST_CASE(genome_sizes)
{
    BOOST_CHECK_EQUAL(2*genome_mutations.size(), genome_mutations.allelic_size());
}

BOOST_AUTO_TEST_CASE(genome_applied_SNVs)
{
    using namespace CLONES::Mutations;

    SID snv0(1, 32, 'A', 'G');  // chromosome 1, position 32

    SID snv1(snv0.chr_id, snv0.position-1, 'T', 'G'), // at snv0 5'
        snv2(snv0.chr_id, snv0.position+1, 'G', 'C'), // at snv0 3'
        snv3(snv0.chr_id+1, snv0.position, 'T', 'G'); // snv0 and snv3 differ in chromosome

    auto test_genome_mutations = genome_mutations;

    // apply snv0 in allele 0: ok
    BOOST_CHECK(test_genome_mutations.insert_in_object(snv0, 0));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        // snv in a wrong position
        SID snv_err(snv0.chr_id, chr.size()+1, 'A', 'G');

        BOOST_CHECK_THROW(test_genome_mutations.insert_in_object(snv_err, 0),
                          CLONES::Error<std::domain_error>);
    }

    // apply snv0 in allele 1: true (free context, non-free context in a different allele)
    BOOST_CHECK(test_genome_mutations.insert_in_object(snv0, 1));

    // apply snv0 in allele 1: false (non-free context)
    BOOST_CHECK(!test_genome_mutations.insert_in_object(snv0, 0));

    // apply snv1 in allele 1: false (non-free context)
    BOOST_CHECK(!test_genome_mutations.insert_in_object(snv1, 1));

    // apply snv2 in allele 1: false (non-free context)
    BOOST_CHECK(!test_genome_mutations.insert_in_object(snv2, 1));

    // apply snv3 in allele 7: ko (unknown allele)
    BOOST_CHECK_THROW(test_genome_mutations.insert_in_object(snv3, 7),
                      CLONES::Error<std::out_of_range>);

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv3.chr_id);

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_mutations().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_mutations().size(), 0);
    }

    // apply snv3 in allele 1: ok (same position, different chromosome)
    BOOST_CHECK(test_genome_mutations.insert_in_object(snv3, 1));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv3.chr_id);

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_mutations().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_mutations().size(), 1);
    }
}

BOOST_AUTO_TEST_CASE(genome_delete_SNVs)
{
    using namespace CLONES::Mutations;

    SID snv0(1, 32, 'A', 'G');

    SID snv1(snv0.chr_id, snv0.position-1, 'T', 'G');

    auto test_genome_mutations = genome_mutations;

    // remove snv0: false (no SNV in the position)
    BOOST_CHECK(!test_genome_mutations.remove_from_object_at(snv0));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        // snv in a wrong position
        SID snv_err(snv0.chr_id, chr.size()+1, 'A', 'G');

        BOOST_CHECK_THROW(test_genome_mutations.remove_from_object_at(snv_err),
                          CLONES::Error<std::domain_error>);
    }

    test_genome_mutations.insert_in_object(snv0, 0);

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_mutations().size(), 1);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_mutations().size(), 0);
    }

    // remove snv1: false (no SNV in the position)
    BOOST_CHECK(!test_genome_mutations.remove_from_object_at(snv1));

    // remove snv0: true (SNV removed)
    BOOST_CHECK(test_genome_mutations.remove_from_object_at(snv0));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_mutations().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_mutations().size(), 0);
    }
}

BOOST_AUTO_TEST_CASE(genome_amplify_region)
{
    using namespace CLONES::Mutations;

    SID snv0(1, 32, 'A', 'G');

    SID snv1(snv0.chr_id, 64, 'A', 'T'),
        snv2(snv0.chr_id, 110, 'T', 'G');

    GenomicRegion gr_32_64(snv0, 33), // from 32 to 64
                  gr_10_110({snv0.chr_id, 10}, 101),  // from 10 to 110
                  gr_60_80({snv0.chr_id, 60}, 21);    // from 60 to 80

    auto test_genome_mutations = genome_mutations;

    // apply snv0 (chromosome1, position 32)
    test_genome_mutations.insert_in_object(snv0, 0);

    // amplify allele 0 from 32 to 64: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_32_64, 0));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(gr_32_64.get_chromosome_id());

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 3);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_mutations().size(), 1);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_mutations().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(2).get_mutations().size(), 1);
        BOOST_CHECK(chr_allele.at(0).get_mutations()==chr_allele.at(2).get_mutations());
    }

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        // snv in a wrong position
        GenomicRegion gr_err(snv0, chr.size()+1);

        BOOST_CHECK_THROW(test_genome_mutations.amplify_region(gr_err, 0),
                          CLONES::Error<std::domain_error>);
    }

    // apply snv1 (chromosome1, position 64)
    test_genome_mutations.insert_in_object(snv1, 0);

    // amplify allele 0 from 10 to 110: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_10_110, 0));

    // apply snv2 (chromosome1, position 110)
    test_genome_mutations.insert_in_object(snv2, 0);

    // amplify allele 0 from 60 to 80: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_60_80, 0));

    {
        // now the chromosome contains 5 alleles:
        // - snv0 lays in alleles 0, 2, and 3
        // - snv1 lays in alleles 0, 3
        // - snv2 lays in alleles 0, 4

        const auto chr = test_genome_mutations.get_chromosomes().at(gr_32_64.get_chromosome_id());

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 5);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_mutations().size(), 3);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_mutations().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(2).get_mutations().size(), 1);
        BOOST_CHECK_EQUAL(chr_allele.at(3).get_mutations().size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(4).get_mutations().size(), 1);
        BOOST_CHECK(chr_allele.at(2).get_mutations()!=chr_allele.at(4).get_mutations());
    }
}

BOOST_AUTO_TEST_CASE(genome_remove_region)
{
    using namespace CLONES::Mutations;

    SID snv0(1, 32, 'A', 'G');

    SID snv1(snv0.chr_id, 64, 'A', 'T'),
        snv2(snv0.chr_id, 110, 'T', 'G');

    GenomicRegion gr_32_64(snv0, 33), // from 32 to 64
                  gr_10_110({snv0.chr_id, 10}, 101),  // from 10 to 110
                  gr_60_80({snv0.chr_id, 60}, 21);    // from 60 to 80

    auto test_genome_mutations = genome_mutations;

    // apply snv0 (chromosome1, position 32)
    test_genome_mutations.insert_in_object(snv0, 0);

    // amplify from 32 to 64: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_32_64, 0));

    // apply snv1 (chromosome1, position 64)
    BOOST_CHECK(test_genome_mutations.insert_in_object(snv1, 0));

    // amplify from 10 to 110: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_10_110, 0));

    // apply snv2 (chromosome1, position 110)
    test_genome_mutations.insert_in_object(snv2, 0);

    // amplify from 60 to 80: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_60_80, 0));

    // now the chromosome contains 5 alleles:
    // - snv0 lays in alleles 0, 2, and 3
    // - snv1 lays in alleles 0, 3
    // - snv2 lays in alleles 0, 4

    // remove allele 0 from 60 to 80: ok
    BOOST_CHECK(test_genome_mutations.remove_region(gr_60_80, 0));

    // now the chromosome contains 5 alleles:
    // - snv0 lays in alleles 0, 2, and 3
    // - snv1 lays in alleles 3
    // - snv2 lays in alleles 0, 4

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(gr_32_64.get_chromosome_id());

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 5);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_mutations().size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_mutations().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(2).get_mutations().size(), 1);
        BOOST_CHECK_EQUAL(chr_allele.at(3).get_mutations().size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(4).get_mutations().size(), 1);
        BOOST_CHECK(chr_allele.at(2).get_mutations()!=chr_allele.at(4).get_mutations());
    }

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        // snv in a wrong position
        GenomicRegion gr_err(snv0, chr.size()+1);

        BOOST_CHECK_THROW(test_genome_mutations.remove_region(gr_err, 0),
                          CLONES::Error<std::domain_error>);
    }

    // amplify allele 0 from 10 to 110: fail (allele 0 does not contain fragment 60-80 anymore)
    BOOST_CHECK(!test_genome_mutations.amplify_region(gr_10_110, 0));

    // remove allele 0 from 60 to 80: fail (allele 0 does not contain fragment 60-80 anymore)
    BOOST_CHECK(!test_genome_mutations.remove_region(gr_60_80, 0));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(germline_SNVs_match_the_reference)
{
    using namespace CLONES::Mutations;
    using namespace CLONES::IO::FASTA;

    // the test FASTA file contains both soft-masked (lower-case) and
    // unknown (N) bases
    std::map<ChromosomeId, std::string> sequences;
    {
        ChromosomeData<Sequence> chr_seq;
        Reader<ChromosomeData<Sequence>> reader(FASTA_FILE);
        while (reader.read(chr_seq)) {
            sequences[chr_seq.chr_id] = chr_seq.nucleotides;
        }
    }

    auto chr_regions = ContextIndex<>::build_index(FASTA_FILE).get_chromosome_regions();

    std::map<ChromosomeId, size_t> alleles_per_chromosome;
    for (const auto& chr_region : chr_regions) {
        alleles_per_chromosome[chr_region.get_chromosome_id()] = 2;
    }

    // a high rate makes many SNVs land on soft-masked bases
    GenomeMutations germline;
    BOOST_REQUIRE_NO_THROW(germline = GermlineMutations::generate(FASTA_FILE, chr_regions,
                                                                  DriverStorage(),
                                                                  alleles_per_chromosome,
                                                                  100, 0, std::cout, true));

    size_t num_of_SNVs{0};
    for (const auto& [chr_id, chr_mutations] : germline.get_chromosomes()) {
        const auto& sequence = sequences.at(chr_id);

        for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
            for (const auto& [position, sid] : allele.get_mutations()) {
                ++num_of_SNVs;

                // genomic positions are 1-based
                const char ref_base = static_cast<char>(std::toupper(sequence[position.position-1]));

                BOOST_CHECK_EQUAL(sid->ref.front(), ref_base);
                BOOST_CHECK(std::string("ACGT").find(sid->ref.front()) != std::string::npos);
                BOOST_CHECK(sid->alt.front() != sid->ref.front());
            }
        }
    }

    BOOST_CHECK_GT(num_of_SNVs, 0);
}

BOOST_AUTO_TEST_CASE(driver_storage_tumour_type_spellings)
{
    using namespace CLONES::Mutations;

    for (const std::string column : {"tumour_type", "tumor_type"}) {
        const auto path = get_a_temporary_path();
        {
            std::ofstream os(path);
            os << "chr\tfrom\tto\tref\talt\tmutation_type\t" << column
               << "\tdriver_gene\tdriver_code\n"
               << "1\t10\t10\tC\tT\tSNV\tTEST\tGENE\tGENE C10T\n";
        }

        BOOST_TEST_CONTEXT("column " << column) {
            BOOST_CHECK_NO_THROW(DriverStorage::load(path));
        }

        std::filesystem::remove(path);
    }
}
