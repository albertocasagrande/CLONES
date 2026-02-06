/**
 * @file germline.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements the functions to generate and load germline mutations
 * @version 1.6
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

#include "germline.hpp"

#include <fstream>
#include <random>
#include <iostream>
#include <regex>
#include <cctype>

#include "csv_reader.hpp"
#include "fasta_chr_reader.hpp"

#include "utils.hpp"

namespace CLONES
{

namespace Mutations
{

GermlineMutations::GermlineMutations(const size_t& genome_size, const double& mutations_per_ref_kbase,
                                     const int& seed):
    bases{'A','C','G','T'}, base_dist(1,3), rand_gen(seed), genome_size{genome_size},
    expected_mutations{static_cast<size_t>(genome_size*mutations_per_ref_kbase/1000)},
    not_processed_size{genome_size}
{
    mutations_not_placed = expected_mutations;

    uint8_t i{0};
    for (const char& base : bases) {
        bases_pos[base] = i;

        ++i;
    }
}

SID GermlineMutations::get_candidate_germline_SNV(const Chromosome& chromosome,
                                                  const ChrPosition& position)
{
    const std::string& sequence = chromosome.nucleotides;

    const uint8_t orig_base_pos = bases_pos.at(sequence[position]);
    const char alt_base = bases[(orig_base_pos+base_dist(rand_gen))%4];

    return SID(chromosome.chr_id, position, sequence[position],
               alt_base, Mutation::GERMINAL);
}


size_t GermlineMutations::get_mutations_in(const ChromosomeMutations& chr_mutations)
{
    const double in_chr_prob = static_cast<double>(chr_mutations.size())/not_processed_size;

    std::binomial_distribution<size_t> mut_dist(mutations_not_placed, in_chr_prob);
    return mut_dist(rand_gen);
}


/**
 * @brief Check whether a genomic position is at least one base away other positions
 *
 * @param positions_to_avoid is a set of genomic positions
 * @param genomic_position is a genomic position
 * @return `true` if and only if `genomic_position` is at least one base away
 *          from any position in `positions_to_avoid`
 */
bool is_context_free_from(const GenomicPosition& genomic_position,
                          const std::set<GenomicPosition>& positions_to_avoid)
{
    auto found = positions_to_avoid.lower_bound(genomic_position);

    if (found != positions_to_avoid.begin()) {
        --found;

        if (genomic_position.chr_id != found->chr_id) {
            ++found;
        }
    }

    while (found != positions_to_avoid.end() && genomic_position.chr_id == found->chr_id
            && found->position <= genomic_position.position+1) {

        if (genomic_position.position <= found->position+1
            && genomic_position.position+1 >= found->position) {
                return false;
        }

        ++found;
    }

    return true;
}

inline
std::set<GenomicPosition>
get_mutation_position_set(const DriverStorage& driver_storage)
{
    auto SID_pos_list = driver_storage.get_mutation_positions();
    return {SID_pos_list.begin(), SID_pos_list.end()};
}

GenomeMutations
GermlineMutations::generate(const std::filesystem::path& reference_fasta_filename,
                            const std::list<GenomicRegion>& chromosome_regions,
                            const DriverStorage& driver_storage,
                            const std::map<ChromosomeId, size_t>& alleles_per_chromosome,
                            const double& mutations_per_ref_kbase,
                            const int& seed, std::ostream& progress_bar_stream,
                            const bool quiet)
{
    CLONES::UI::ProgressBar progress_bar(progress_bar_stream, quiet);

    progress_bar.set_message("Generating germline");

    GenomeMutations germline(chromosome_regions, alleles_per_chromosome);

    GermlineMutations generator(germline.size(), mutations_per_ref_kbase,
                                seed);

    using namespace CLONES::IO::FASTA;

    auto driver_positions = get_mutation_position_set(driver_storage);

    ChromosomeData<Sequence> chr_seq;
    Reader<ChromosomeData<Sequence>> chr_reader(reference_fasta_filename);
    while (chr_reader.read(chr_seq, progress_bar)) {
        auto& chr_mutations = germline.get_chromosome(chr_seq.chr_id);
        auto mutations_in_chr = generator.get_mutations_in(chr_mutations);

        if (chr_mutations.size()>2) {
            std::uniform_int_distribution<ChrPosition> pos_dist(1, chr_mutations.size()-2);

            size_t num_of_alleles = chr_mutations.get_alleles().size();
            while (mutations_in_chr > 0) {
                ChrPosition pos = pos_dist(generator.rand_gen);
                if (chr_seq.nucleotides[pos-1]!='N' && chr_seq.nucleotides[pos]!='N'
                        && chr_seq.nucleotides[pos+1]!='N') {

                    auto snv = generator.get_candidate_germline_SNV(chr_seq, pos);

                    if (is_context_free_from(snv, driver_positions)) {

                        // This assume no deletion
                        AlleleId allele_id = generator.rand_gen()%num_of_alleles;

                        if (chr_mutations.apply(snv, allele_id)) {
                            --mutations_in_chr;
                            --generator.mutations_not_placed;
                        }
                    }
                }

                const auto percentage = generator.get_process_percentage();
                if (progress_bar.get_progress() < percentage) {
                    progress_bar.set_progress(percentage);
                }
            }

            generator.not_processed_size -= chr_mutations.size();
        }
    }

    progress_bar.set_progress(100, "Germline generated");

    return germline;
}

std::map<ChromosomeId, std::filesystem::path>
load_germline_chr_map(const std::filesystem::path& germline_data_file)
{
    std::map<ChromosomeId, std::filesystem::path> germline_chr_map;

    CLONES::IO::CSVReader germline_data(germline_data_file, true, '\t');

    auto germline_path = germline_data_file.parent_path();

    for (const auto& row : germline_data) {
        ChromosomeId chr_id = GenomicPosition::stochr(row.get_field(0));

        germline_chr_map.insert({chr_id, germline_path/row.get_field(1)});
    }

    return germline_chr_map;
}

bool read_up_to_header(std::string& header, std::ifstream& VCF_stream)
{
    const std::regex header_regex = build_regex("^#CHROM.*");
    std::smatch matches;

    std::string line;
    while (std::getline(VCF_stream, line)) {
        if (std::regex_match(line, matches, header_regex)) {
            header = line;
            return true;
        }
        if (line.size()>0 && line[0] != '#') {
            return false;
        }
    }

    return false;
}

bool is_male(const std::map<ChromosomeId, std::filesystem::path>& germline_chr_map,
             const std::string& subject)
{
    const auto chr_Y = GenomicPosition::stochr("Y");

    auto chr_it = germline_chr_map.find(chr_Y);
    if (chr_it == germline_chr_map.end()) {
        return false;
    }

    std::ifstream VCF_X_stream(chr_it->second);
    if (!VCF_X_stream.good()) {
        throw std::runtime_error("The file \"" + to_string(chr_it->second)
                                 + "\" is not readable.");
    }

    std::string header;
    if (read_up_to_header(header, VCF_X_stream)) {
        const std::regex subject_regex = build_regex(".*"+subject+".*");
        std::smatch matches;

        return std::regex_match(header, matches, subject_regex);
    }

    return false;
}

std::list<GenomicRegion>
get_chromosome_regions_from_VCF(const std::map<ChromosomeId, std::filesystem::path>& germline_chr_map)
{
    if (germline_chr_map.size() == 0) {
        return {};
    }

    std::ifstream VCF_stream(germline_chr_map.begin()->second);

    if (!VCF_stream.good()) {
        throw std::runtime_error("The file \"" + to_string(germline_chr_map.begin()->second)
                                 + "\" is not readable.");
    }

    const std::regex chr_regex = build_regex("^##contig=<([0-9a-zA-Z=]+,)*ID=([0-9]+|X|Y).*");
    const std::regex chr_length = build_regex("^.*(<|,)length=([0-9]+)(,|>).*");
    std::smatch matches;

    std::list<GenomicRegion> chr_regions;

    std::string line;
    while (std::getline(VCF_stream, line) && line.size()>0 && line[0] == '#') {
        if (std::regex_match(line, matches, chr_regex)) {
            auto chr_id = GenomicPosition::stochr(matches[2].str());

            if (germline_chr_map.count(chr_id)>0
                    && std::regex_match(line, matches, chr_length)) {

                GenomicRegion::Length length = std::stoul(matches[2].str());
                chr_regions.push_back({chr_id, length});
            }
        }
    }

    return chr_regions;
}

GenomeMutations
init_genome_mutations_from_VCF(const std::map<ChromosomeId, std::filesystem::path>& germline_chr_map,
                               const std::map<ChromosomeId, size_t>& alleles_per_chromosome)
{

    if (germline_chr_map.size()==0) {
        return {};
    }

    auto chr_regions = get_chromosome_regions_from_VCF(germline_chr_map);

    return {chr_regions, alleles_per_chromosome};
}

std::vector<size_t>
find_occurrences(const std::string& line, const std::string& pattern, size_t pos=0)
{
    std::vector<size_t> occurrences;

    while ((pos=line.find(pattern, pos+1))!=std::string::npos) {
        occurrences.push_back(pos);
    }

    return occurrences;
}

std::vector<std::string>
tokenize_by(const std::string& line, const char& sep='\t')
{
    std::vector<std::string> tokens;

    std::istringstream iss(line);
    std::string token;
    while(std::getline(iss, token, sep)) {
        tokens.push_back(token);
    }

    return tokens;
}

size_t
get_subject_column(const std::string& header, const std::string& subject)
{
    auto tokens = tokenize_by(header, '\t');
    size_t pos{0};
    for (const auto& token : tokens) {
        if (token == subject) {
            return pos;
        }

        ++pos;
    }

    throw std::runtime_error(subject + " not present.");
}

MutationType::Type get_mutation_type(const std::string& line, size_t pos=0)
{
    pos = line.find("VT=", pos);

    if (pos != std::string::npos) {
        pos += 3;

        if (line.size() < pos+3) {
            return MutationType::Type::UNKNOWN;
        }
        switch (line[pos]) {
        case 'S':
            if (line[pos+1]!='N' || line[pos+2]!='P') {
                return MutationType::Type::UNKNOWN;
            }

            pos += 3;
            if (line.size() == pos || line[pos] == '\t'
                    || line[pos] == ';') {
                return MutationType::Type::SBS;
            }

            break;
        case 'I':
            if (line.find("INDEL", pos) == pos) {
                return MutationType::Type::INDEL;
            }
            break;
        case 'M':
            if (line.find("MNP", pos) == pos) {
                return MutationType::Type::MNP;
            }
            break;
        default:
            return MutationType::Type::UNKNOWN;
        }
    }

    return MutationType::Type::UNKNOWN;
}

template<typename INTEGER_TYPE>
bool read_int(INTEGER_TYPE& value, const std::string& line, size_t& pos, size_t end)
{
    if (end > line.size()) {
        end = line.size();
    }

    const size_t begin{pos};

    value = 0;
    while (pos < end && isdigit(line[pos])) {
        value = value*10+(line[pos]-'0');
        ++pos;
    }

    return pos != begin;
}

std::string get_elem(const std::string& line,
                     const std::vector<size_t>& column_separators,
                     const size_t& column_index)
{
    if (column_index==0) {
        return line.substr(0, column_separators[0]);
    }

    const size_t elem_pos = column_separators[column_index-1]+1;
    const size_t elem_size = column_separators[column_index]-elem_pos;
    return line.substr(elem_pos, elem_size);
}

SID get_SID_from_line(const std::string& line,
                      const std::vector<size_t>& column_separators)
{
    auto chr_str = line.substr(0, column_separators[0]);

    ChromosomeId chr_id = GenomicPosition::stochr(chr_str);

    ChrPosition chr_pos;

    size_t pos = column_separators[0]+1;
    read_int(chr_pos, line, pos, column_separators[1]);

    const std::string ref = get_elem(line, column_separators, 3);
    std::string alt = get_elem(line, column_separators, 4);

    pos = alt.find(',');
    if (pos != std::string::npos) {
        alt.resize(pos);
    }

    return {chr_id, chr_pos, ref, alt, Mutation::GERMINAL};
}

std::vector<AlleleId>
get_alleles_ids(const std::string& line, size_t pos, size_t end)
{
    std::vector<AlleleId> alleles;
    while (pos < end) {
        AlleleId allele;

        if (read_int(allele, line, pos, end)) {
            alleles.push_back(allele);
            ++pos;
        } else {
            return alleles;
        }
    }

    return alleles;
}

void add_SID(GenomeMutations& mutations, const std::string& line,
             const std::vector<size_t>& column_separators,
             const size_t& num_of_alleles,
             const std::vector<AlleleId>& allele_ids)
{
    const auto mutation = get_SID_from_line(line, column_separators);

    if (allele_ids[0] == 1) {
        mutations.apply(mutation, 0);
    }

    if (num_of_alleles>1 && allele_ids.size()>1 && allele_ids[1]==1) {
        mutations.apply(mutation, 1);
    }
}

void add_mutation(GenomeMutations& mutations, const std::string& line,
                  const std::vector<size_t>& column_separators,
                  const size_t& num_of_alleles, const size_t& subject_column)
{
    size_t last_char_position;

    if (subject_column<column_separators.size()) {
        last_char_position = column_separators[subject_column];
    } else {
        last_char_position = line.size();
    }

    auto allele_ids = get_alleles_ids(line, column_separators[subject_column-1]+1,
                                      last_char_position);

    if (allele_ids.size()==0) {
        return;
    }

    bool in_first_allele = (allele_ids[0] == 1);
    bool in_second_allele = (num_of_alleles>1 && allele_ids.size()>1
                                && allele_ids[1]==1);

    if (in_first_allele || in_second_allele) {
        switch (get_mutation_type(line, column_separators[6])) {
            case MutationType::Type::SBS:
            case MutationType::Type::INDEL:
                add_SID(mutations, line, column_separators,
                        num_of_alleles, allele_ids);
                break;
            default:
                break;
        }
    }
}

inline size_t get_file_size(const std::filesystem::path& filename)
{
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);

    return ifs.tellg();
}

void add_VCF_mutations(GenomeMutations& mutations, const std::filesystem::path& VCF_file,
                       const size_t& num_of_alleles, const std::string& subject,
                       const size_t& total_VCF_size, size_t& processed_VCF_size,
                       CLONES::UI::ProgressBar& progress_bar)
{
    if (num_of_alleles==0) {
        return;
    }

    std::ifstream VCF_stream(VCF_file);

    if (!VCF_stream.good()) {
        throw std::runtime_error("The file \"" + to_string(VCF_file)
                                 + "\" is not readable.");
    }

    std::string header;
    if (!read_up_to_header(header, VCF_stream)) {
        throw std::runtime_error("The VCF file \"" + to_string(VCF_file)
                                 + "\" misses the standard header line.");
    }

    size_t subject_column;

    try {
        subject_column = get_subject_column(header, subject);
    } catch(std::runtime_error& ex) {
        throw std::runtime_error("\"" + to_string(VCF_file)
                                 + "\" does not contains "
                                 + subject + "'s mutations.");
    }

    std::string line;
    std::size_t line_num{0};
    std::streampos last_pos{0};
    while (std::getline(VCF_stream, line)) {
        ++line_num;

        if (line.size()>0 && line[0] != '#') {
            auto column_separators = find_occurrences(line, "\t");

            if (column_separators.size()+1 < subject_column) {
                throw std::runtime_error("Wrong number of columns in file \""
                                         + to_string(VCF_file)
                                         + "\" line number "
                                         + std::to_string(line_num) + ".");
            }

            add_mutation(mutations, line, column_separators,
                         num_of_alleles, subject_column);
        }

        //const size_t percentage = (100*(processed_VCF_size+VCF_stream.tellg()))/total_VCF_size;
        std::streampos file_pos = VCF_stream.tellg();
        if (file_pos>=0 && file_pos - last_pos > 1000000) {
            last_pos = file_pos;
            const size_t percentage = (100*(static_cast<size_t>(last_pos)+processed_VCF_size))/total_VCF_size;
            progress_bar.set_progress(std::min(static_cast<size_t>(100), percentage));
        }
    }

    processed_VCF_size += get_file_size(VCF_file);

    const size_t percentage = (100*processed_VCF_size)/total_VCF_size;
    progress_bar.set_progress(std::min(static_cast<size_t>(100), percentage));
}

size_t get_total_VCF_size(const std::map<ChromosomeId, std::filesystem::path>& germline_chr_map,
                          const std::map<ChromosomeId, size_t>& alleles_per_chromosome)
{
    size_t total_size{0};

    for (const auto& [chr_id, VCF_path]: germline_chr_map) {
        auto found = alleles_per_chromosome.find(chr_id);
        if (found != alleles_per_chromosome.end()) {
            total_size += get_file_size(VCF_path);
        }
    }

    return total_size;
}

GenomeMutations
GermlineMutations::load(const std::filesystem::path& germline_data_file,
                        const std::map<ChromosomeId, size_t>& alleles_per_chromosome,
                        const std::string& subject, std::ostream& progress_bar_stream,
                        const bool quiet)
{
    CLONES::UI::ProgressBar progress_bar(progress_bar_stream, quiet);

    progress_bar.set_message("Loading germline");

    auto germline_chr_map = load_germline_chr_map(germline_data_file);
    auto mutations = init_genome_mutations_from_VCF(germline_chr_map,
                                                    alleles_per_chromosome);

    const size_t total_VCF_size = get_total_VCF_size(germline_chr_map,
                                                     alleles_per_chromosome);

    size_t processed_VCF{0};
    for (const auto& [chr_id, VCF_path]: germline_chr_map) {
        auto found = alleles_per_chromosome.find(chr_id);
        if (found != alleles_per_chromosome.end()) {
            add_VCF_mutations(mutations, VCF_path, found->second, subject,
                              total_VCF_size, processed_VCF, progress_bar);
        }
    }
    progress_bar.set_progress(100, "Germline loaded");

    return mutations;
}

}   // Mutations

}   // CLONES
