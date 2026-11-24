/**
 * @file read_simulator.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing read simulator
 * @version 1.2
 * @date 2026-11-24
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
#define BOOST_TEST_MODULE read_simulator

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <random>
#include <algorithm>
#include <set>

#include <boost/test/unit_test.hpp>

#include "read_simulator.hpp"
#include "genome_fragment.hpp"
#include "tissue_simulation.hpp"
#include "ending_conditions.hpp"
#include "descendant_forest.hpp"
#include "mutation_engine.hpp"
#include "sequencer.hpp"
#include "utils.hpp"

#include "error.hpp"

CLONES::Mutations::SID build_SNV(const CLONES::Mutations::ChrPosition position)
{
    using namespace CLONES::Mutations;

    return SID(1, position, 'A', 'T', Mutation::PASSENGER);
}

void add_SNVs(CLONES::Mutations::SequencingSimulations::ChrSampleStatistics& stats,
              const std::vector<CLONES::Mutations::ChrPosition>& positions)
{
    using namespace CLONES::Mutations::SequencingSimulations;

    for (const auto& position : positions) {
        const auto snv = build_SNV(position);

        stats.get_data().insert({snv, SIDData(snv, 1)});
    }
}

BOOST_AUTO_TEST_CASE(chr_statistics_merge)
{
    using namespace CLONES::Mutations;
    using namespace CLONES::Mutations::SequencingSimulations;

    // the SIDs of the two statistics interleave: some SIDs of the second
    // statistics precede, follow, or lay between those of the first one
    ChrSampleStatistics stats_a, stats_b;
    add_SNVs(stats_a, {100, 300});
    add_SNVs(stats_b, {50, 100, 200, 400});

    stats_a += stats_b;

    const auto& data = stats_a.get_data();
    BOOST_REQUIRE_EQUAL(data.size(), 5);

    std::vector<ChrPosition> positions;
    for (const auto& [mutation, SID_data] : data) {
        positions.push_back(mutation.position);
    }
    BOOST_CHECK(positions == std::vector<ChrPosition>({50, 100, 200, 300, 400}));

    // the occurrences of the common SID are summed
    BOOST_CHECK_EQUAL(data.at(build_SNV(100)).num_of_occurrences, 2);
    BOOST_CHECK_EQUAL(data.at(build_SNV(50)).num_of_occurrences, 1);
}

BOOST_AUTO_TEST_CASE(genome_fragment_at_chromosome_end)
{
    using namespace CLONES;
    using namespace CLONES::Mutations;

    // positions:                1234567890
    const std::string reference("ACGTACGTAC");

    // a deletion of the positions 6-8 anchored on position 5
    CLONES::map<GenomicPosition, std::shared_ptr<SID>> germline, somatic;
    somatic.emplace(GenomicPosition(1, 5),
                    std::make_shared<SID>(1, 5, "ACGT", "A", Mutation::PASSENGER));

    // a 6-base fragment beginning in position 3 would require the
    // reference base in position 11, which does not exist
    const GenomeFragment fragment(reference, germline, somatic, GenomicPosition(1, 3), 6);

    // the fragment is truncated at the reference end
    BOOST_CHECK_EQUAL(fragment.get_sequence(), "GTAAC");

    const auto covered = fragment.get_covered_reference_region();
    BOOST_CHECK_EQUAL(covered.begin(), 3);
    BOOST_CHECK_EQUAL(covered.end(), 10);
}

/**
 * @brief A fixture providing a phylogenetic forest and a read simulator
 *
 * The forest is built on a reference genome with two chromosomes of 300 and
 * 400 bases. Its mutant has no passenger mutations and a single driver
 * deletion of the positions 150-249 of chromosome 2 on allele 0. Hence,
 * every sampled cell has two alleles per chromosome and allele 0 of
 * chromosome 2 consists of two fragments, i.e., [1,149] and [250,400].
 */
struct ReadSimulatorFixture
{
    TemporaryDirectory tmp_dir;     //!< The directory of the simulation data
    std::filesystem::path reference;    //!< The reference genome FASTA file
    CLONES::Mutations::PhylogeneticForest forest;   //!< The phylogenetic forest

    //! The deletion in chromosome 2
    static constexpr CLONES::Mutations::ChrPosition deletion_begin = 150,
                                                    deletion_length = 100;

    //! The read size
    static constexpr size_t read_size = 10;

    /**
     * @brief Write a random reference genome
     *
     * The indexed FASTA reader requires all the lines of a sequence, but
     * the last one, to have the same length.
     *
     * @param path is the path of the FASTA file
     * @param chr_sizes is the map associating chromosome names and sizes
     */
    static void write_reference(const std::filesystem::path& path,
                                const std::map<std::string, size_t>& chr_sizes)
    {
        const size_t line_size = 60;

        std::mt19937_64 generator(0);
        std::uniform_int_distribution<size_t> base_dist(0, 3);

        std::ofstream FASTA_stream(path);
        for (const auto& [name, size] : chr_sizes) {
            FASTA_stream << ">" << name << " dna:chromosome chromosome:GRCh38:"
                         << name << ":1:" << size << ":1 REF";
            for (size_t i=0; i<size; ++i) {
                if (i % line_size == 0) {
                    FASTA_stream << std::endl;
                }
                FASTA_stream << "ACGT"[base_dist(generator)];
            }
            FASTA_stream << std::endl;
        }
    }

    ReadSimulatorFixture():
        reference{tmp_dir.path() / "reference.fasta"}
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;
        using namespace CLONES::Mutants::Evolutions;
        using namespace CLONES::Mutations;

        write_reference(reference, {{"1", 300}, {"2", 400}});

        TissueSimulation simulation(tmp_dir.path() / "tissue", 0);

        simulation.set_tissue("Liver", {100, 100})
                  .add_mutant("A");

        auto& species = simulation.tissue().get_species("A");
        species.set_rate(CellEventType::DEATH, 0.01)
               .set_rate(CellEventType::DUPLICATION, 0.2);

        simulation.place_cell(species, {50, 50});
        simulation.death_activation_level = 100;

        TimeTest done(20);
        simulation.run(done);

        simulation.sample_tissue(SampleSpecification("S", RectangleSet({45, 45}, {55, 55})));

        const DescendantForest descendant_forest(simulation);

        auto context_index = ContextIndex<>::build_index(reference);
        auto rs_index = RSIndex::build_index(reference, 50, 1000000, 0);

        const GenomeMutations germline(context_index.get_chromosome_regions(), 2);

        MutationalProperties mutational_properties;
        mutational_properties.add_mutant("A", PassengerRates(0, 0, 0), {},
                                         {CNA(GenomicPosition(2, deletion_begin), deletion_length,
                                              CNA::Type::DELETION, 0)});

        // no SID is placed, but the engine requires the pre-neoplastic signatures
        const std::map<std::string, SBSSignature> SBS_signatures{{"SBS1", SBSSignature()}};
        const std::map<std::string, IDSignature> ID_signatures{{"ID1", IDSignature()}};

        MutationEngine<> engine(context_index, rs_index, SBS_signatures, ID_signatures,
                                mutational_properties, germline, DriverStorage());

        forest = engine.place_mutations(descendant_forest, 0, 0);
    }

    /**
     * @brief Get a read simulator that writes SAM files
     *
     * @param name is the name of the output subdirectory
     * @return a single-read simulator writing its SAM files in `name`
     */
    CLONES::Mutations::SequencingSimulations::ReadSimulator<>
    get_simulator(const std::string& name) const
    {
        using namespace CLONES::Mutations::SequencingSimulations;

        ReadSimulator<> simulator(tmp_dir.path() / name, reference, read_size);

        simulator.enable_SAM_writing(true);

        return simulator;
    }

    /**
     * @brief A read in a SAM file
     */
    struct SAMRead
    {
        CLONES::Mutations::ChromosomeId chr_id;     //!< The chromosome identifier
        CLONES::Mutations::ChrPosition begin;       //!< The first aligned position
        CLONES::Mutations::ChrPosition end;         //!< The last aligned position
    };

    /**
     * @brief Read all the reads in the SAM files of a directory
     *
     * @param name is the name of the output subdirectory
     * @return the reads in the SAM files of `name`
     */
    std::vector<SAMRead> get_reads(const std::string& name) const
    {
        using namespace CLONES::Mutations;

        std::vector<SAMRead> reads;
        for (const auto& entry : std::filesystem::directory_iterator(tmp_dir.path() / name)) {
            if (entry.path().extension() != ".sam") {
                continue;
            }

            std::ifstream SAM_stream(entry.path());
            std::string line;
            while (std::getline(SAM_stream, line)) {
                if (line.empty() || line[0] == '@') {
                    continue;
                }

                // QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ ...
                std::istringstream line_stream(line);
                std::string qname, flag, rname, mapq, cigar, rnext, pnext, tlen, seq;
                ChrPosition pos;
                line_stream >> qname >> flag >> rname >> pos >> mapq >> cigar
                            >> rnext >> pnext >> tlen >> seq;

                // the forest has no SIDs: the reads match the reference
                BOOST_REQUIRE_EQUAL(cigar, std::to_string(seq.size()) + "M");

                reads.push_back({GenomicPosition::stochr(rname), pos,
                                 static_cast<ChrPosition>(pos + seq.size() - 1)});
            }
        }

        return reads;
    }
};

/**
 * @brief Get an error-less sequencer
 *
 * @return a reference to an error-less sequencer
 */
CLONES::Sequencers::Illumina::ErrorLessSequencer<>& sequencer()
{
    static CLONES::Sequencers::Illumina::ErrorLessSequencer<> errorless_sequencer;

    return errorless_sequencer;
}

BOOST_FIXTURE_TEST_SUITE(target_regions, ReadSimulatorFixture)

BOOST_AUTO_TEST_CASE(default_target_regions)
{
    using namespace CLONES::Mutations::SequencingSimulations;

    std::ostringstream progress_stream;

    // by default, the whole chromosomes are targeted
    auto simulator = get_simulator("default");
    BOOST_REQUIRE_NO_THROW(simulator(sequencer(), forest, {1, 2}, 5, false, 1.0,
                                     true, true, "chr_", false, false,
                                     progress_stream, true));

    std::set<CLONES::Mutations::ChromosomeId> chr_ids;
    for (const auto& read : get_reads("default")) {
        chr_ids.insert(read.chr_id);
    }

    BOOST_CHECK(chr_ids.count(1) > 0);
    BOOST_CHECK(chr_ids.count(2) > 0);
}

BOOST_AUTO_TEST_CASE(reads_in_target_regions)
{
    using namespace CLONES::Mutations;

    // target regions partially beyond the chromosome end, or completely
    // beyond it, are trimmed or dropped
    const std::vector<GenomicRegion> target_regions{{1, 50, 120}, {1, 250, 400},
                                                    {2, 80, 200}, {2, 500, 600}};

    std::ostringstream progress_stream;

    auto simulator = get_simulator("targets");
    BOOST_REQUIRE_NO_THROW(simulator(sequencer(), forest, target_regions, {1, 2}, 50,
                                     false, 1.0, true, true, "chr_", false, false,
                                     progress_stream, true));

    const auto reads = get_reads("targets");

    BOOST_REQUIRE_GT(reads.size(), 0);

    for (const auto& read : reads) {
        const bool in_a_target = std::any_of(target_regions.begin(), target_regions.end(),
                                             [&read](const GenomicRegion& region) {
                                                 return region.get_chromosome_id() == read.chr_id
                                                        && region.begin() <= read.begin
                                                        && read.end <= region.end();
                                             });

        BOOST_TEST_CONTEXT("read on chr " << static_cast<int>(read.chr_id) << " ["
                                          << read.begin << "," << read.end << "]") {
            BOOST_CHECK(in_a_target);
        }
    }
}

BOOST_AUTO_TEST_CASE(target_regions_spanning_fragments)
{
    using namespace CLONES::Mutations;

    // the target region of chromosome 2 spans both the fragments of allele 0
    // and the deleted region between them
    const GenomicRegion chr2_target{2, 50, 350};

    std::ostringstream progress_stream;

    auto simulator = get_simulator("spanning");
    BOOST_REQUIRE_NO_THROW(simulator(sequencer(), forest, {chr2_target}, {2}, 200,
                                     false, 1.0, true, true, "chr_", false, false,
                                     progress_stream, true));

    std::map<ChrPosition, size_t> depth;
    for (const auto& read : get_reads("spanning")) {
        for (auto pos = read.begin; pos <= read.end; ++pos) {
            ++depth[pos];
        }
    }

    auto mean_depth = [&depth](const ChrPosition begin, const ChrPosition end) {
        double sum{0};
        for (auto pos = begin; pos <= end; ++pos) {
            sum += depth[pos];
        }
        return sum / (end - begin + 1);
    };

    // avoid the positions closer than a read to the fragment ends
    const auto deletion_end = deletion_begin + deletion_length - 1;
    const double flanking_depth = (mean_depth(chr2_target.begin() + read_size, deletion_begin - read_size)
                                   + mean_depth(deletion_end + read_size, chr2_target.end() - read_size)) / 2;
    const double deleted_depth = mean_depth(deletion_begin + read_size, deletion_end - read_size);

    // allele 1 is the only one covering the deleted region, while both the
    // alleles cover the flanking regions, so the depth ratio is about 1/2
    // (slightly less, because reads cannot begin in the last positions of
    // the shorter fragments of allele 0). Reads of allele 0 in the deleted
    // region would make the two depths equal.
    BOOST_REQUIRE_GT(flanking_depth, 0);
    BOOST_CHECK_GT(deleted_depth / flanking_depth, 0.4);
    BOOST_CHECK_LT(deleted_depth / flanking_depth, 0.6);
}

BOOST_AUTO_TEST_SUITE_END()
