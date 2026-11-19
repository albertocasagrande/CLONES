/**
 * @file mutations_sim.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the CLONES mutations simulator
 * @version 1.15
 * @date 2026-11-19
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <random>
#include <map>
#include <list>
#include <vector>
#include <charconv>
#include <cstdint>
#include <cctype>

#include <boost/program_options.hpp>

#include "common.hpp"

#include "tissue_simulation.hpp"
#include "descendant_forest.hpp"
#include "context_index.hpp"
#include "mutation_engine.hpp"
#include "germline.hpp"
#include "fasta_chr_reader.hpp"
#include "read_simulator.hpp"
#include "distributions.hpp"

#include "driver_storage.hpp"

#include "csv_reader.hpp"

#include "json_config.hpp"

#include "progress_bar.hpp"

class MutationsSimulator : public BasicExecutable
{
    std::filesystem::path simulation_filename;
    std::filesystem::path species_directory;
    std::filesystem::path snapshot_path;
    std::filesystem::path driver_mutations_filename;
    std::filesystem::path context_index_filename;
    std::filesystem::path rs_index_filename;
    std::filesystem::path ref_genome_filename;
    std::filesystem::path passenger_CNA_filename;
    std::filesystem::path SBS_filename;
    std::filesystem::path ID_filename;

    double coverage;
    double purity;
    std::string seq_output_directory;
    CLONES::Mutations::SequencingSimulations::ReadSimulator<>::Mode read_simulator_output_mode;
    bool paired_read;
    size_t read_size;
    double insert_size_mean;
    double insert_size_stddev;
    double sequencer_error_rate;
    std::string SIDs_csv_filename;
    std::string CNAs_csv_filename;
    std::string pre_neoplastic_SNV_signature_name;
    std::string pre_neoplastic_ID_signature_name;
    std::string germline_csv_filename;
    std::string germline_subject;
    bool epigenetic_FACS;
    bool write_SAM;

    int seed;
    size_t bytes_per_abs_position;
    bool quiet;

    static std::ifstream open_input_file(const std::filesystem::path& filename)
    {
        using namespace CLONES;

        std::ifstream is(filename);
        if (!is) {
            throw Error<std::runtime_error>("Cannot open \"" + to_string(filename)
                                            + "\" for reading.");
        }

        return is;
    }

    static std::ofstream open_output_file(const std::filesystem::path& filename)
    {
        using namespace CLONES;

        std::ofstream os(filename);
        if (!os) {
            throw Error<std::runtime_error>("Cannot open \"" + to_string(filename)
                                            + "\" for writing.");
        }

        return os;
    }

    /**
     * @brief Parse an integer CSV field
     *
     * The whole field, except leading and trailing blanks, must be an
     * integer representable by `INTEGER`.
     */
    template<typename INTEGER>
    static INTEGER parse_integer_field(const std::string& field, const std::string& field_name,
                                       const size_t& row_num)
    {
        using namespace CLONES;

        const auto first = field.find_first_not_of(" \t\r");
        const auto last = field.find_last_not_of(" \t\r");

        INTEGER value{0};
        if (first != std::string::npos) {
            const char* begin = field.data()+first;
            const char* end = field.data()+last+1;

            const auto [ptr, ec] = std::from_chars(begin, end, value);
            if (ec == std::errc() && ptr == end) {
                return value;
            }
        }

        throw Error<std::runtime_error>("Invalid " + field_name + " specification \""
                                        + field + "\" in row number "
                                        + std::to_string(row_num) + ".");
    }

    std::list<CLONES::Mutants::Evolutions::TissueSample>
    get_samples(CLONES::Mutants::Evolutions::TissueSimulation& simulation, const nlohmann::json& simulation_cfg) const
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;
        using namespace CLONES::Mutants::Evolutions;

        std::list<TissueSample> samples;

        if (simulation_cfg.contains("sample regions")) {
            const auto& sample_regions_json = simulation_cfg["sample regions"];
            if (!sample_regions_json.is_array()) {
                throw Error<std::runtime_error>("The \"sample regions\" field must be an array.");
            }

            for (const auto& sample_region_json: sample_regions_json) {
                auto sample_specification = ConfigReader::get_sample_specification(sample_region_json);
                samples.push_back(simulation.simulate_sampling(sample_specification));
            }

            return samples;
        }

        return simulation.get_tissue_samples();
    }

    template<typename ABSOLUTE_GENOTYPE_POSITION = uint32_t>
    CLONES::Mutations::ContextIndex<ABSOLUTE_GENOTYPE_POSITION> load_context_index(const std::string& filename) const
    {
        CLONES::Archive::Binary::In archive(filename);

        CLONES::Mutations::ContextIndex<ABSOLUTE_GENOTYPE_POSITION> context_index;

        if (quiet) {
            archive & context_index;
        } else {
            archive.load(context_index, "context index");
        }

        return context_index;
    }

    CLONES::Mutations::RSIndex load_rs_index(const std::string& filename) const
    {
        CLONES::Archive::Binary::In archive(filename);

        CLONES::Mutations::RSIndex rs_index;

        if (quiet) {
            archive & rs_index;
        } else {
            archive.load(rs_index, "repetition index");
        }

        return rs_index;
    }

    template<typename MUTATION_TYPE>
    static std::map<std::string, CLONES::Mutations::Signature<MUTATION_TYPE>>
    load_signatures(const std::filesystem::path& filename)
    {
        auto is = open_input_file(filename);

        return CLONES::Mutations::Signature<MUTATION_TYPE>::read_from_stream(is);
    }

    template<typename ABSOLUTE_GENOME_POSITION, typename RANDOM_GENERATOR>
    CLONES::Mutations::PhylogeneticForest
    place_mutations(CLONES::Mutations::MutationEngine<ABSOLUTE_GENOME_POSITION,RANDOM_GENERATOR>& engine,
                    const CLONES::Mutants::DescendantForest& forest,
                    const size_t& num_of_pre_neoplastic_SNVs,
                    const size_t& num_of_pre_neoplastic_IDs,
                    const std::string& pre_neoplastic_SNV_signature,
                    const std::string& pre_neoplastic_indel_signature) const
    {
        if (quiet) {
            return engine.place_mutations(forest, num_of_pre_neoplastic_SNVs,
                                          num_of_pre_neoplastic_IDs, seed,
                                          pre_neoplastic_SNV_signature,
                                          pre_neoplastic_indel_signature);
        }

        CLONES::UI::ProgressBar progress_bar(std::cout);

        progress_bar.set_message("Placing mutations");

        auto phylo_forest = engine.place_mutations(forest, num_of_pre_neoplastic_SNVs,
                                                   num_of_pre_neoplastic_IDs, progress_bar,
                                                   seed, pre_neoplastic_SNV_signature,
                                                   pre_neoplastic_indel_signature);

        progress_bar.set_message("Mutations placed");

        return phylo_forest;
    }

    nlohmann::json get_simulation_json() const
    {
        auto simulation_stream = open_input_file(simulation_filename);

        return nlohmann::json::parse(simulation_stream);
    }

    /**
     * @brief Complete and validate the references of the driver SIDs
     *
     * This method fills the missing "ref" fields of the driver SNVs in the
     * configuration by using the reference genome and checks that the
     * "ref" fields of all the driver SIDs match the reference genome.
     * It also rejects driver SNVs whose alternative base equals the
     * reference one and driver indels lacking the "ref" field.
     *
     * @param[in, out] simulation_cfg is the mutations simulation configuration
     */
    void complete_driver_references(nlohmann::json& simulation_cfg) const
    {
        using namespace CLONES;
        using namespace CLONES::Mutations;
        using namespace CLONES::IO::FASTA;

        if (!simulation_cfg.contains("mutant properties")
                || !simulation_cfg["mutant properties"].is_array()) {
            // `ConfigReader` reports the configuration errors
            return;
        }

        // collect the driver SIDs per chromosome
        std::map<ChromosomeId, std::list<nlohmann::json*>> chr_SIDs;
        for (auto& mutant_json : simulation_cfg["mutant properties"]) {
            if (!mutant_json.is_object() || !mutant_json.contains("driver mutations")
                    || !mutant_json["driver mutations"].is_array()) {
                continue;
            }

            for (auto& mutation_json : mutant_json["driver mutations"]) {
                if (!mutation_json.is_object() || !mutation_json.contains("type")
                        || !mutation_json.contains("chromosome")
                        || !mutation_json.contains("position")) {
                    continue;
                }

                const auto type = mutation_json["type"].template get<std::string>();
                if (type != "SNV" && type != "indel") {
                    continue;
                }

                if (type == "indel" && !mutation_json.contains("ref")) {
                    throw Error<std::runtime_error>("The driver indel " + mutation_json.dump()
                                                    + " lacks the \"ref\" field.");
                }

                const auto chr_name = mutation_json["chromosome"].template get<std::string>();
                chr_SIDs[GenomicPosition::stochr(chr_name)].push_back(&mutation_json);
            }
        }

        if (chr_SIDs.empty()) {
            return;
        }

        auto to_upper = [](std::string sequence) {
            for (auto& base : sequence) {
                base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
            }

            return sequence;
        };

        UI::ProgressBar progress_bar(std::cout, quiet);
        progress_bar.set_message("Checking driver references");

        ChromosomeData<Sequence> chr;
        Reader<ChromosomeData<Sequence>> reader(ref_genome_filename);
        while (!chr_SIDs.empty() && reader.read(chr, progress_bar)) {
            auto found = chr_SIDs.find(chr.chr_id);
            if (found == chr_SIDs.end()) {
                continue;
            }

            for (auto* mutation_json : found->second) {
                const auto position = (*mutation_json)["position"].template get<ChrPosition>();
                const auto description = mutation_json->dump();

                // genomic positions are 1-based
                if (position < 1 || position > chr.nucleotides.size()) {
                    throw Error<std::runtime_error>("The driver mutation " + description
                                                    + " lays outside the reference chromosome.");
                }

                if (!mutation_json->contains("ref")) {
                    // only driver SNVs may lack the "ref" field
                    const auto ref_base = to_upper(chr.nucleotides.substr(position-1, 1));
                    if (ref_base.find_first_not_of("ACGT") != std::string::npos) {
                        throw Error<std::runtime_error>("The reference base of the driver SNV "
                                                        + description + " is unknown (\""
                                                        + ref_base + "\").");
                    }

                    (*mutation_json)["ref"] = ref_base;
                } else {
                    const auto ref = to_upper((*mutation_json)["ref"].template get<std::string>());

                    // "-" and "" stand for an empty reference sequence (i.e., an insertion)
                    if (ref != "-" && ref != "") {
                        const auto genome_ref = to_upper(chr.nucleotides.substr(position-1,
                                                                                ref.size()));
                        if (ref != genome_ref) {
                            throw Error<std::runtime_error>("The reference of the driver mutation "
                                                            + description + " does not match the "
                                                            + "reference genome (\"" + genome_ref
                                                            + "\").");
                        }

                        // store the reference in upper case as the filled ones
                        (*mutation_json)["ref"] = ref;
                    }
                }

                const auto type = (*mutation_json)["type"].template get<std::string>();
                if (type == "SNV" && mutation_json->contains("alt")) {
                    const auto ref = to_upper((*mutation_json)["ref"].template get<std::string>());
                    const auto alt = to_upper((*mutation_json)["alt"].template get<std::string>());

                    if (ref == alt) {
                        throw Error<std::runtime_error>("The driver SNV " + mutation_json->dump()
                                                        + " does not change the reference base.");
                    }

                    (*mutation_json)["alt"] = alt;
                }
            }

            chr_SIDs.erase(found);
        }

        if (!chr_SIDs.empty()) {
            throw Error<std::runtime_error>("The reference genome lacks chromosome "
                                            + GenomicPosition::chrtos(chr_SIDs.begin()->first)
                                            + ", which contains some driver mutations.");
        }
    }

    void process_statistics(const CLONES::Mutations::PhylogeneticForest& forest) const
    {
        using namespace CLONES;
        using namespace CLONES::UI;
        using namespace CLONES::Mutations;

        if (SIDs_csv_filename.empty() && CNAs_csv_filename.empty()) {
            return;
        }

        ProgressBar bar(std::cout, quiet);

        MutationStatistics statistics;

        statistics.record(forest, bar);

        if (!SIDs_csv_filename.empty()) {
            auto os = open_output_file(SIDs_csv_filename);

            statistics.write_SIDs_table(os);

            if (!os) {
                throw Error<std::runtime_error>("Error writing \"" + SIDs_csv_filename + "\".");
            }
        }

        if (!CNAs_csv_filename.empty()) {
            auto os = open_output_file(CNAs_csv_filename);

            statistics.write_CNAs_table(os);

            if (!os) {
                throw Error<std::runtime_error>("Error writing \"" + CNAs_csv_filename + "\".");
            }
        }
    }

    template<typename TISSUE_SAMPLE,
             std::enable_if_t<std::is_base_of_v<CLONES::Mutants::Evolutions::TissueSample, TISSUE_SAMPLE>, bool> = true>
    static std::map<CLONES::Time, std::list<TISSUE_SAMPLE>>
    split_list_by_time(const std::list<TISSUE_SAMPLE>& sample_list)
    {
        std::map<CLONES::Time, std::list<TISSUE_SAMPLE>> time_map;

        for (const auto& sample : sample_list) {
            time_map[sample.get_time()].push_back(sample);
        }

        return time_map;
    }

    template<typename GENOME_WIDE_POSITION>
    static std::map<CLONES::Mutations::ChromosomeId, size_t>
    get_number_of_alleles(const CLONES::Mutations::ContextIndex<GENOME_WIDE_POSITION>&context_index,
                          const nlohmann::json& simulation_cfg)
    {
        using namespace CLONES::Mutations;

        auto chromosome_ids = context_index.get_chromosome_ids();
        return CLONES::ConfigReader::get_number_of_alleles(simulation_cfg, chromosome_ids);
    }

    static CLONES::Mutations::GenomicRegion get_CNA_region(const CLONES::IO::CSVReader::CSVRow& row, const size_t& row_num)
    {
        using namespace CLONES;
        using namespace CLONES::Mutations;

        if (row.size()<5) {
            throw Error<std::runtime_error>("The CNA CSV must contains at least 5 columns. "
                                            + std::to_string(row.size()) + " found in row number "
                                            + std::to_string(row_num) + ".");
        }
        ChromosomeId chr_id;

        try {
            chr_id = GenomicPosition::stochr(row.get_field(0));
        } catch (std::invalid_argument const&) {
            throw Error<std::runtime_error>("Unknown chromosome specification "
                                            + row.get_field(0) + " in row number "
                                            + std::to_string(row_num) + ".");
        }

        const auto begin_pos = parse_integer_field<ChrPosition>(row.get_field(1), "begin", row_num);
        const auto end_pos = parse_integer_field<ChrPosition>(row.get_field(2), "end", row_num);

        if (begin_pos>end_pos) {
            throw Error<std::runtime_error>("The CNA begin (i.e., "
                                            + std::to_string(begin_pos)
                                            + ") lays after its end (i.e., "
                                            + std::to_string(end_pos)
                                            + ") in row number "
                                            + std::to_string(row_num) + ".");
        }

        GenomicPosition pos(chr_id, begin_pos);

        return {pos, end_pos+1-begin_pos};
    }

    void saving_statistics_data_and_images(const CLONES::Mutations::SequencingSimulations::SampleSetStatistics& statistics,
                                           const std::string& base_name="chr_") const
    {
        statistics.save_VAF_CSVs(base_name, std::cout, quiet);
#if WITH_MATPLOT
        statistics.save_coverage_images(base_name, std::cout, quiet);
        statistics.save_SID_histograms(base_name, std::cout, quiet);
#endif // WITH_MATPLOT
    }

    template<typename GENOME_WIDE_POSITION>
    static void add_exposures(CLONES::Mutations::MutationEngine<GENOME_WIDE_POSITION>& engine,
                              const nlohmann::json& exposures_json, const std::string& mutation_name)
    {
        using namespace CLONES;

        ConfigReader::expecting(mutation_name, exposures_json,
                                "The elements of the \"exposures\" field");
        const auto& type_exposures_json = exposures_json[mutation_name];

        auto default_exposure = ConfigReader::get_default_exposure(mutation_name,
                                                                   type_exposures_json);
        engine.add(default_exposure);

        auto timed_exposures = ConfigReader::get_timed_exposures(mutation_name,
                                                                 type_exposures_json);

        for (const auto& [time, exposure] : timed_exposures) {
            engine.add(time, exposure);
        }
    }

    class EpigeneticLabelling
    {
        CLONES::Mutants::DescendantForest const* forest;
    public:
        EpigeneticLabelling(const CLONES::Mutants::DescendantForest& forest):
            forest{&forest}
        {}

        std::string operator()(const CLONES::Mutants::CellId& cell_id) const
        {
            const auto node = forest->get_node(cell_id);

            return node.get_epistate_name();
        }
    };

    /**
     * @brief Build the read simulator
     *
     * The read simulator creates, or wipes when overwriting, the
     * sequencing output directory. Hence, it must be built only when
     * the reads are about to be simulated.
     */
    CLONES::Mutations::SequencingSimulations::ReadSimulator<> build_read_simulator() const
    {
        using namespace CLONES::Mutations::SequencingSimulations;

        if (paired_read) {
            std::normal_distribution<double> insert_size_dist(insert_size_mean,
                                                              insert_size_stddev);
            return ReadSimulator<>(seq_output_directory, ref_genome_filename, read_size,
                                   insert_size_dist, read_simulator_output_mode, true,
                                   "r", seed);
        }

        return ReadSimulator<>(seq_output_directory, ref_genome_filename, read_size,
                               read_simulator_output_mode, true, "r", seed);
    }

    template<typename GENOME_WIDE_POSITION>
    void run_abs_position() const
    {
        using namespace CLONES;
        using namespace CLONES::Mutants;
        using namespace CLONES::Mutations;
        using namespace CLONES::Mutations::SequencingSimulations;

        nlohmann::json simulation_cfg = get_simulation_json();

        if (!simulation_cfg.contains("exposures")) {
            throw Error<std::runtime_error>("The passengers simulation configuration "
                                            "must contain an \"exposures\" field");
        }

        complete_driver_references(simulation_cfg);

        if (!quiet) {
            if (simulation_cfg.contains("sample regions")) {
                std::cout << "Using \"sample regions\" configuration field as sample regions." << std::endl;
            } else {
                std::cout << "Using pre-collected samples." << std::endl;
            }
        }

        Mutants::DescendantForest descendant_forest;
        MutationalProperties mutational_properties;

        {
            auto species_simulation = load_species_simulation(snapshot_path, quiet);

            auto samples = get_samples(species_simulation, simulation_cfg);

            descendant_forest = Mutants::DescendantForest(species_simulation, samples);

            mutational_properties = ConfigReader::get_mutational_properties(simulation_cfg);
        }

        auto SBS_signatures = load_signatures<SBSType>(SBS_filename);
        auto ID_signatures = load_signatures<IDType>(ID_filename);

        auto context_index = load_context_index<GENOME_WIDE_POSITION>(context_index_filename);
        auto rs_index = load_rs_index(rs_index_filename);

        auto num_of_alleles = get_number_of_alleles(context_index, simulation_cfg);

        const auto passenger_CNAs = MutationsSimulator::load_passenger_CNAs(passenger_CNA_filename);

        const auto driver_storage = DriverStorage::load(driver_mutations_filename);

        auto germline_mutations_per_kbase = ConfigReader::get_number_of_germline_mutations_per_kbase(simulation_cfg);

        auto chromosome_regions = context_index.get_chromosome_regions();

        GenomeMutations germline;

        if (germline_csv_filename.size()==0) {
            germline = GermlineMutations::generate(ref_genome_filename, chromosome_regions, driver_storage,
                                                   num_of_alleles, germline_mutations_per_kbase,
                                                   seed, std::cout, quiet);
        } else {
            germline = GermlineMutations::load(germline_csv_filename, num_of_alleles, germline_subject);
        }

        MutationEngine<GENOME_WIDE_POSITION> engine(context_index, rs_index,
                                                    SBS_signatures, ID_signatures,
                                                    mutational_properties, germline, driver_storage,
                                                    passenger_CNAs);

        const auto& exposures_json = simulation_cfg["exposures"];

        add_exposures(engine, exposures_json, "SBS");
        add_exposures(engine, exposures_json, "indel");

        auto num_of_pnp_SNV = ConfigReader::get_number_of_pre_neoplastic_mutations(simulation_cfg, "SNV");
        auto num_of_pnp_indel = ConfigReader::get_number_of_pre_neoplastic_mutations(simulation_cfg, "indel");

        auto forest = place_mutations(engine, descendant_forest, num_of_pnp_SNV,
                                      num_of_pnp_indel, pre_neoplastic_SNV_signature_name,
                                      pre_neoplastic_ID_signature_name);

        if (epigenetic_FACS) {
            EpigeneticLabelling labelling{forest};

            forest.partition_samples(labelling);
        }

        if (coverage>0) {
            auto read_simulator = build_read_simulator();

            read_simulator.enable_SAM_writing(write_SAM);

            CLONES::Sequencers::Illumina::BasicSequencer sequencer(sequencer_error_rate);

            auto statistics = read_simulator(sequencer, forest, coverage, purity,
                                             true, true, "chr_", false, false,
                                             std::cout, quiet);

            // VAF tables, coverage images, and SID histograms require all the
            // samples to report the very same SIDs
            statistics.canonize();

            saving_statistics_data_and_images(statistics);
        }

        process_statistics(forest);
    }

    void test_file_existence(const boost::program_options::variables_map& vm, const std::string& name,
                             const std::filesystem::path& file, const std::string& producer="") const
    {
        if (!vm.count(name)) {
            std::string msg = "The " + name + " file is mandatory.";
            if (producer.size()>0) {
                msg = msg + " You can produce it by using `" + producer + "`.";
            }
            print_help_and_exit(msg, 1);
        }

        if (!std::filesystem::exists(file)) {
            print_help_and_exit("\"" + to_string(file) + "\"  does not exist.", 1);
        }

        if (!std::filesystem::is_regular_file(file)) {
            print_help_and_exit("\"" + to_string(file) + "\"  is not a regular file.", 1);
        }
    }

    void validate(const boost::program_options::variables_map& vm) const
    {
        namespace fs = std::filesystem;

        test_file_existence(vm, "mutation", simulation_filename);

        if (!vm.count("species simulation")) {
            print_help_and_exit("The species simulation directory is mandatory. "
                                "You can produce it by using `mutants_sim`.", 1);
        }

        if (!fs::exists(species_directory)) {
            print_help_and_exit("\"" + to_string(species_directory) + "\"  does not exist.", 1);
        }

        if (!fs::is_directory(species_directory)) {
            print_help_and_exit("\"" + to_string(species_directory) + "\"  is not a directory.", 1);
        }

        test_file_existence(vm, "driver mutations", driver_mutations_filename);
        test_file_existence(vm, "context index", context_index_filename, "build_context_index");
        test_file_existence(vm, "repetition index", rs_index_filename, "build_repetition_index");
        test_file_existence(vm, "SBS signature", SBS_filename);
        test_file_existence(vm, "ID signature", ID_filename);
        test_file_existence(vm, "reference genome", ref_genome_filename);
        test_file_existence(vm, "passenger CNAs", passenger_CNA_filename);

        if (vm.count("germline-file") != vm.count("germline-subject")) {
            if (vm.count("germline-file")>0) {
                print_help_and_exit("A germline file has been specified. Please, "
                                    "specify a subject among those in the IGSR "
                                    "VCF files.", 1);
            }

            print_help_and_exit("A germline subject has been specified. Please, "
                                "specify a germline file.", 1);
        }

        if (!SIDs_csv_filename.empty() && fs::exists(SIDs_csv_filename)) {
            print_help_and_exit("\"" + SIDs_csv_filename + "\" already exists.", 1);
        }

        if (!CNAs_csv_filename.empty() && fs::exists(CNAs_csv_filename)) {
            print_help_and_exit("\"" + CNAs_csv_filename + "\" already exists.", 1);
        }

        if (coverage<0) {
            print_help_and_exit("The coverage must be non-negative.", 1);
        }

        if (purity<0 || purity>1) {
            print_help_and_exit("The purity must be in the interval [0,1].", 1);
        }

        if (sequencer_error_rate<0 || sequencer_error_rate>1) {
            print_help_and_exit("The sequencer error rate must be in the interval [0,1].", 1);
        }

        if (read_size==0) {
            print_help_and_exit("The read size must be positive.", 1);
        }

        if (vm.count("insert-size-stddev")>0 && !paired_read) {
            print_help_and_exit("\"--insert-size-stddev\" requires \"--insert-size-mean\".", 1);
        }

        if (paired_read) {
            // check the insert size distribution feasibility before the
            // (possibly long) mutation placement
            try {
                using dist_type = CLONES::left_truncated_normal_distribution<double>;
    
                dist_type::fit_untruncated(insert_size_mean, insert_size_stddev,
                                           static_cast<double>(read_size));
            } catch (const std::exception& e) {
                print_help_and_exit("Invalid insert size distribution for read size "
                                    + std::to_string(read_size) + ": " + e.what(), 1);
            }
        }

        if (coverage>0) {
            if (vm.count("output-directory")==0) {
                print_help_and_exit("\"--output-directory\" is mandatory when "
                                    "coverage differs from 0.", 1);
            }

            if (fs::exists(seq_output_directory)) {
                if (vm.count("overwrite")==0) {
                    print_help_and_exit("The output directory \"" + seq_output_directory
                                        + "\" already exists.", 1);
                }

                if (!fs::is_directory(seq_output_directory)) {
                    print_help_and_exit("\"" + seq_output_directory + "\" is not a directory.", 1);
                }
            }
        }
    }

public:
    static std::vector<CLONES::Mutations::CNA> load_passenger_CNAs(const std::filesystem::path& CNAs_csv)
    {
        using namespace CLONES;
        using namespace CLONES::Mutations;

        std::vector<CNA> CNAs;

        CLONES::IO::CSVReader csv_reader(CNAs_csv, true, '\t');

        size_t row_num{2};
        for (const auto& row : csv_reader) {
            const auto region = get_CNA_region(row, row_num);

            const auto major = row.get_field(3);
            if (major=="NA" || parse_integer_field<int>(major, "major", row_num)>1) {
                CNAs.emplace_back(region.get_initial_position(), region.size(), CNA::Type::AMPLIFICATION);
            }

            const auto minor = row.get_field(4);
            if (minor=="NA" || parse_integer_field<int>(minor, "minor", row_num)<1) {
                CNAs.emplace_back(region.get_initial_position(), region.size(), CNA::Type::DELETION);
            }

            ++row_num;
        }

        return CNAs;
    }

    MutationsSimulator(int argc, char* argv[]):
        BasicExecutable(argv[0], {{"mutations", "Mutation simulation options"},
                                  {"sequencing", "Sequencing simulation options"},
                                  {"mutants", "Mutants evolution related options"},
                                  {"generic", "Generic options"}}),
        paired_read(false),
        insert_size_mean(0), insert_size_stddev(0),
        sequencer_error_rate(0), SIDs_csv_filename(""),
        CNAs_csv_filename(""), epigenetic_FACS(false)
    {
        namespace po = boost::program_options;

        using namespace CLONES::Mutations::SequencingSimulations;

        visible_options.at("mutations").add_options()
            ("germline-file,G", po::value<std::string>(&germline_csv_filename),
             "CSV file reporting the IGSR VCF file of each chromosome")
            ("germline-subject,J", po::value<std::string>(&germline_subject),
             "name of the germline subject" )
            ("pnp-SNV,S",
             po::value<std::string>(&pre_neoplastic_SNV_signature_name)->default_value("SBS1"),
             "name of the pre-neoplastic SNV signature" )
            ("pnp-indel,I",
             po::value<std::string>(&pre_neoplastic_ID_signature_name)->default_value("ID1"),
             "name of the pre-neoplastic indel signature" )
            ("mutations-CSV,M", po::value<std::string>(&SIDs_csv_filename),
             "mutations CSV output file")
            ("CNAs-CSV,C", po::value<std::string>(&CNAs_csv_filename),
             "CNAs CSV output file")
        ;

        visible_options.at("sequencing").add_options()
            ("coverage,c", po::value<double>(&coverage)->default_value(0.0),
             "coverage for the simulated reads (>0 to simulate sequencing)")
            ("purity,p", po::value<double>(&purity)->default_value(1.0),
             "purity of the sample (a value in [0,1])")
            ("output-directory,d", po::value<std::string>(&seq_output_directory),
             "sequencing output directory")
            ("overwrite,o", "overwrite the output directory")
            ("read-size,r", po::value<size_t>(&read_size)->default_value(150),
             "simulated read size")
            ("insert-size-mean,i", po::value<double>(&insert_size_mean),
             "insert size mean (enable paired-read)")
            ("insert-size-stddev,t", po::value<double>(&insert_size_stddev),
             "insert size standard deviation")
            ("sequencer-error-rate,e", po::value<double>(&sequencer_error_rate),
             "sequencer error rate per base")
            ("write-SAM,w", "write the simulated read SAM files")
        ;

        visible_options.at("mutants").add_options()
            ("epigenetic-FACS,F", "distinguish between epigenetic states")
        ;

        visible_options.at("generic").add_options()
            ("seed,s", po::value<int>(&seed)->default_value(0),
             "random generator seed")
            ("quiet,q", "avoid output messages")
            ("help,h", "get the help")
        ;

        hidden_options.add_options()
            ("mutation", po::value<std::filesystem::path>(&simulation_filename),
             "the name of the file describing the mutations simulation")
            ("species simulation", po::value<std::filesystem::path>(&species_directory),
             "the species simulation directory")
            ("context index", po::value<std::filesystem::path>(&context_index_filename),
             "the genome context index")
            ("repetition index", po::value<std::filesystem::path>(&rs_index_filename),
             "the genome repetition index")
            ("SBS signature", po::value<std::filesystem::path>(&SBS_filename),
             "the SBS signature")
            ("ID signature", po::value<std::filesystem::path>(&ID_filename),
             "the ID signature")
            ("reference genome", po::value<std::filesystem::path>(&ref_genome_filename),
             "the filename of the reference genome FASTA file")
            ("driver mutations",  po::value<std::filesystem::path>(&driver_mutations_filename),
             "the name of file containing the driver mutations")
            ("passenger CNAs", po::value<std::filesystem::path>(&passenger_CNA_filename),
             "the filename of the passenger CNAs")
        ;

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("mutation", 1);
        positional_options.add("species simulation", 1);
        positional_options.add("context index", 1);
        positional_options.add("repetition index", 1);
        positional_options.add("SBS signature", 1);
        positional_options.add("ID signature", 1);
        positional_options.add("reference genome", 1);
        positional_options.add("driver mutations", 1);
        positional_options.add("passenger CNAs", 1);

        po::variables_map vm;

        try {
            po::command_line_parser parser{argc, argv};
            po::store(parser.options(program_options).positional(positional_options).run(), vm);
            po::notify(vm);
        } catch (const po::error& except) {
            print_help_and_exit(except.what(), 1);
        }

        if (vm.count("help")) {
            print_help_and_exit(0);
        }

        quiet = vm.count("quiet")>0;

        read_simulator_output_mode = ((vm.count("overwrite")>0)?
                                      ReadSimulator<>::Mode::OVERWRITE:
                                      ReadSimulator<>::Mode::CREATE);

        paired_read = vm.count("insert-size-mean")>0;

        epigenetic_FACS = vm.count("epigenetic-FACS")>0;
        write_SAM = vm.count("write-SAM")>0;

        validate(vm);

        try {
            using namespace CLONES::Mutations;

            CLONES::Archive::Binary::In archive(context_index_filename);

            bytes_per_abs_position = ContextIndex<>::read_bytes_per_absolute_position(archive);

            snapshot_path = get_last_snapshot_path(species_directory, "species simulation");
        } catch (const std::exception& except) {
            print_help_and_exit(except.what(), 1);
        }
    }

    void run() const
    {
        try {
            switch (bytes_per_abs_position) {
                case 2:
                    run_abs_position<uint16_t>();
                    break;
                case 4:
                    run_abs_position<uint32_t>();
                    break;
                case 8:
                    run_abs_position<uint64_t>();
                    break;
                default:
                    print_help_and_exit("Unsupported bytes per absolute position ("
                                        + std::to_string(bytes_per_abs_position)
                                        + "). Supported values are 2, 4, or 8.", 1);
            }
        } catch (const std::exception& except) {
            print_help_and_exit(except.what(), 1);
        }
    }
};

int main(int argc, char* argv[])
{
    MutationsSimulator passenger_sim(argc, argv);

    passenger_sim.run();

    return 0;
}
