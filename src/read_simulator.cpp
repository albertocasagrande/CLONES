/**
 * @file read_simulator.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes to simulate sequencing
 * @version 1.4
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

#include <string>
#include <regex>

#include <algorithm>

#include "read_simulator.hpp"

#include "utils.hpp"

namespace CLONES
{

namespace Mutations
{

namespace SequencingSimulations
{

ChrCoverage::ChrCoverage()
{}

ChrCoverage::ChrCoverage(const ChromosomeId& chromosome_id, const GenomicRegion::Length& size):
    chr_id(chromosome_id), coverage(size)
{
    if (size == 0) {
        throw std::domain_error("A chromosome cannot have length 0.");
    }
}


void ChrCoverage::increase_coverage(const ChrPosition& begin_pos, const size_t& read_size)
{
    if (begin_pos+read_size-2 > coverage.size()) {
        using namespace CLONES::Mutants;

        throw std::runtime_error("The chromosome " + GenomicPosition::chrtos(chr_id) + " has length "
                                 + std::to_string(coverage.size())+". Coverage cannot be increased "
                                 + "up to position " + std::to_string(begin_pos) + "+"
                                 + std::to_string(read_size) + " as requested.");
    }

    for (auto base=coverage.begin()+begin_pos-1;
            base != coverage.begin()+begin_pos+read_size-1; ++base) {
        ++(*base);
    }
}


void check_in(const GenomicPosition& pos, const ChromosomeId& chr_id)
{
    if (pos.chr_id != chr_id) {
        using namespace CLONES::Mutants;

        std::ostringstream oss;

        oss << pos << " does not belong to the chromosome "
            << GenomicPosition::chrtos(chr_id) << ".";

        throw std::domain_error(oss.str());
    }
}

const BaseCoverage& ChrCoverage::get_coverage(const GenomicPosition& position) const
{
    check_in(position, chr_id);

    if (position.position == 0 && coverage.size() < position.position) {
        throw std::out_of_range("Position must lay in the interval [1,"
                                + std::to_string(coverage.size())
                                + "]: requested coverage of position "
                                + std::to_string(position.position)
                                + ".");
    }

    return coverage[position.position-1];
}

void check_same_chr(const ChrCoverage& a, const ChrCoverage& b)
{
    if (a.get_chr_id() != b.get_chr_id()) {
        using namespace CLONES::Mutants;

        std::ostringstream oss;

        oss << "The two sequencing chromosome statistics refer "
            << "to different chromosomes: " << GenomicPosition::chrtos(a.get_chr_id())
            << " and " << GenomicPosition::chrtos(b.get_chr_id()) << ".";

        throw std::domain_error(oss.str());
    }

    if (a.size() != b.size()) {
        using namespace CLONES::Mutants;

        std::ostringstream oss;

        oss << "The two coverages have different lengths: " << a.size()
            << " and " << b.size() << ".";

        throw std::domain_error(oss.str());
    }
}

void update_coverage(std::vector<BaseCoverage>& a, const std::vector<BaseCoverage>& b)
{
    auto b_it=b.begin();
    for (auto a_it=a.begin(); a_it != a.end(); ++a_it, ++b_it) {
        *a_it += *b_it;
    }
}

ChrCoverage& ChrCoverage::operator+=(const ChrCoverage& chr_coverage)
{
    if (coverage.size()==0) {
        if (chr_coverage.coverage.size()!=0) {
            *this = chr_coverage;
        }

        return *this;
    }

    if (chr_coverage.coverage.size()==0) {
        return *this;
    }

    check_same_chr(*this, chr_coverage);

    update_coverage(coverage, chr_coverage.coverage);

    return *this;
}

ChrCoverage& ChrCoverage::operator+=(ChrCoverage&& chr_coverage)
{
    if (coverage.size()==0) {
        // `*this` has been initialized yet
        if (chr_coverage.coverage.size()!=0) {
            *this = std::move(chr_coverage);
        }

        return *this;
    }

    // `chr_coverage` has been initialized yet
    if (chr_coverage.coverage.size()==0) {
        return *this;
    }

    check_same_chr(*this, chr_coverage);

    update_coverage(coverage, chr_coverage.coverage);

    return *this;
}

SIDData::SIDData():
    num_of_occurrences(0)
{}

SIDData::SIDData(const SID& mutation, const BaseCoverage num_of_occurrences):
    num_of_occurrences(num_of_occurrences)
{
    causes.insert(mutation.cause);
    nature_set.insert(mutation.nature);
}

void SIDData::update(const SID& mutation)
{
    ++num_of_occurrences;

    causes.insert(mutation.cause);
    nature_set.insert(mutation.nature);
}

template<typename T>
void update_set(std::set<T>& S, const std::set<T>& P)
{
    std::set<T> new_S;

    for (const auto& value : P) {
        S.insert(value);
    }
}

void SIDData::update(const SIDData& mutation_data)
{
    num_of_occurrences += mutation_data.num_of_occurrences;

    update_set(causes, mutation_data.causes);
    update_set(nature_set, mutation_data.nature_set);
}

ChrSampleStatistics::ChrSampleStatistics():
    ChrCoverage()
{}

ChrSampleStatistics::ChrSampleStatistics(const ChromosomeId& chromosome_id,
                                         const GenomicRegion::Length& size,
                                         const std::list<SampleGenomeMutations>& mutations_list):
    ChrCoverage(chromosome_id, size)
{
    for (const SampleGenomeMutations& mutations : mutations_list) {
        account_for(mutations.germline_mutations->get_chromosome(chromosome_id));
        for (const auto& cell_genome_ptr : mutations.mutations) {
            account_for(cell_genome_ptr->get_chromosome(chromosome_id));
        }
    }
}

void ChrSampleStatistics::account_for(const ChromosomeMutations& chromosome)
{
    for (const auto& [allele_id, allele]: chromosome.get_alleles()) {
        for (const auto& [frag_pos, fragment]: allele.get_fragments()) {
            for (const auto& [SID_pos, SID_ptr] : fragment.get_mutations()) {
                auto SID_data_it = SID_data.find(*SID_ptr);

                if (SID_data_it == SID_data.end()) {
                    SID_data.insert({*SID_ptr, {*SID_ptr, 0}});
                }
            }
        }
    }
}

void check_in(const SID& mutation, const ChromosomeId& chr_id)
{
    if (mutation.chr_id != chr_id) {
        using namespace CLONES::Mutants;

        std::ostringstream oss;

        oss << mutation << " does not lay in the chromosome "
            << GenomicPosition::chrtos(chr_id) << ".";

        throw std::domain_error(oss.str());
    }
}

void ChrSampleStatistics::update(const SID& mutation)
{
    check_in(mutation, get_chr_id());

    auto SID_data_it = SID_data.find(mutation);

    if (SID_data_it == SID_data.end()) {
        SID_data.insert({mutation, {mutation, 1}});
    } else {
        SID_data_it->second.update(mutation);
    }
}

void ChrSampleStatistics::add(const Read& read)
{
    increase_coverage(read.get_covered_reference_region());

    for (const auto& mutation: read.get_mutations()) {
        update(mutation);
    }
}

SIDData ChrSampleStatistics::get_data(const SID& mutation) const
{
    check_in(mutation, get_chr_id());

    auto SID_data_it = SID_data.find(mutation);

    if (SID_data_it == SID_data.end()) {
        return {mutation, 0};
    }

    return SID_data_it->second;
}

std::map<SID, SIDData>::iterator
find_not_before(const SID& mutation, const std::map<SID, SIDData>& data,
                std::map<SID, SIDData>::iterator it)
{
    std::less<SID> before;
    while (it != data.end() && before(it->first, mutation)) {
        ++it;
    }

    return it;
}

void update_data(std::map<SID, SIDData>& a, const std::map<SID, SIDData>& b)
{
    auto a_it = a.begin();
    auto b_it = b.begin();

    std::less<SID> before;
    while (b_it != b.end()) {

        // find occurrence of b_it->first in a
        a_it = find_not_before(b_it->first, a, a_it);

        if (a_it != a.end()) {
            // found an SID not laying before b_it->first

            if (!before(b_it->first, a_it->first)) {
                // if the found SID does not lay after b_it->first,
                // then it is in the same SID

                a_it->second.update(b_it->second);
            }
            ++b_it;
        } else {
            for (; b_it != b.end(); ++b_it) {
                a.insert(*b_it);
            }
        }
    }
}

ChrSampleStatistics& ChrSampleStatistics::operator+=(const ChrSampleStatistics& chr_stats)
{
    static_cast<ChrCoverage&>(*this) += chr_stats;

    update_data(SID_data, chr_stats.SID_data);

    return *this;
}

ChrSampleStatistics& ChrSampleStatistics::operator+=(ChrSampleStatistics&& chr_stats)
{
    static_cast<ChrCoverage*>(this)->operator+=(std::move(chr_stats));

    update_data(SID_data, chr_stats.SID_data);

    return *this;
}

std::filesystem::path SampleStatistics::get_coverage_filename(const ChromosomeId& chr_id) const
{
    return get_data_directory()/("chr_"+GenomicPosition::chrtos(chr_id)+"_coverage.dat");
}

SampleStatistics::SampleStatistics(const std::filesystem::path& data_directory,
                                   const std::string& sample_name,
                                   const bool save_coverage):
    data_dir(data_directory/sample_name), sample_name(sample_name),
    save_coverage(save_coverage)
{}

void create_dir(const std::filesystem::path& directory)
{
    using namespace std::filesystem;

    if (!exists(directory)) {
        create_directory(directory);

        return;
    }

    if (!is_directory(directory)) {
        throw std::runtime_error("\"" + to_string(directory)
                                 + "\" exists, but it is not a directory.");
    }
}

void SampleStatistics::add_chr_statistics(const ChrSampleStatistics& chr_stats)
{
    const auto& chr_id = chr_stats.get_chr_id();

    if (chr_ids.count(chr_id)>0) {
        using namespace CLONES::Mutants;

        throw std::runtime_error("SampleStatistics " + sample_name + " already "
                                 + "includes statistics about chromosome "
                                 + GenomicPosition::chrtos(chr_id) + ".");
    }

    for (const auto&[mutation, data] : chr_stats.get_data()) {
        SID_data.insert({mutation, data});
        locus_coverage[mutation] = chr_stats.get_coverage(mutation);
    }

    chr_ids.insert(chr_id);

    if (save_coverage) {
        create_dir(get_data_directory());

        const auto chr_filename = get_coverage_filename(chr_id);
        Archive::Binary::Out out_archive(chr_filename);
        static_cast<const ChrCoverage&>(chr_stats).save(out_archive);
    }
}

ChrCoverage SampleStatistics::get_chr_coverage(const ChromosomeId& chr_id) const
{
    if (!save_coverage) {
        throw std::runtime_error("Coverage data is not available.");
    }

    if (chr_ids.count(chr_id)==0) {
        using namespace CLONES::Mutants;

        throw std::runtime_error("SampleStatistics " + sample_name + " does not "
                                 + "include statistics about chromosome "
                                 + GenomicPosition::chrtos(chr_id) + ".");
    }

    const auto chr_filename = get_coverage_filename(chr_id);
    Archive::Binary::In in_archive(chr_filename);
    return ChrCoverage::load(in_archive);
}

SIDData SampleStatistics::get_data(const SID& mutation) const
{
    auto found = SID_data.find(mutation);
    if (found == SID_data.end()) {
        return {mutation, 0};
    }

    return found->second;
}

const BaseCoverage& SampleStatistics::get_coverage(const GenomicPosition& pos) const
{
    auto found = locus_coverage.find(pos);
    if (found == locus_coverage.end()) {
        std::ostringstream oss;

        oss << "\"" << sample_name
            << "\" does not contain any SID in position "
            << pos << ".";
        throw std::runtime_error(oss.str());
    }

    return found->second;
}

void print_SID_data(std::ofstream& os, const size_t& coverage,
                    const SIDData& mutation_data, const char separator='\t')
{
    if (mutation_data.num_of_occurrences != 0) {
        os << separator << mutation_data.num_of_occurrences
           << separator << coverage
           << separator << (static_cast<double>(mutation_data.num_of_occurrences)/coverage);
    } else {
        os << separator << "0" << separator << coverage << separator << "0";
    }
}

std::ofstream& SampleSetStatistics::stream_VAF_csv_header(std::ofstream& os, const char& separator) const
{
    os << "chr" << separator << "from" << separator << "to" << separator
        << "ref" << separator << "alt" << separator
        << "type" << separator << "causes" << separator << "classes";

    if (stats_map.size()>1) {
        for (const auto& [sample_name, sample_stats]: stats_map) {
            os << separator << sample_name << " occurrences"
               << separator << sample_name << " coverage"
               << separator << sample_name << " VAF";
        }
        os << separator << "total occurrences"
           << separator << "total coverage"
           << separator << "total VAF" << std::endl;
    } else {
        os << separator << "occurrences"
           << separator << "coverage"
           << separator << "VAF" << std::endl;
    }

    return os;
}

void check_in(const std::map<std::string, SampleStatistics>& stats_map,
              const std::string& sample_name)
{
    if (stats_map.count(sample_name)>0) {
        throw std::runtime_error("The statistics of the sample \"" + sample_name
                                 + "\" are already contained by the object.");
    }
}

SampleSetStatistics::SampleSetStatistics(const std::filesystem::path& data_directory):
    data_dir(data_directory)
{
    create_dir(data_dir);
}

std::map<std::string, std::map<GenomicPosition, BaseCoverage>>
get_sample_coverages_for_SIDs(const std::map<std::string, SampleStatistics>& stats_map,
                              const std::map<SID, BaseCoverage>& SID_occurrences)
{
    using SIDCoverage = std::map<GenomicPosition, BaseCoverage>;
    std::map<std::string, SIDCoverage> sample_locus_coverage;

    for (auto& [sample_name, sample_stats]: stats_map) {
        ChrCoverage chr_coverage;
        auto it = sample_locus_coverage.insert({sample_name,
                                                sample_stats.get_coverage()}).first;

        for (auto& [mutation, occurrences] : SID_occurrences) {
            auto& sample_coverage = it->second;
            if (sample_coverage.count(mutation)==0) {
                if (mutation.chr_id!=chr_coverage.get_chr_id()
                        || !chr_coverage.is_initialized()) {
                    chr_coverage = sample_stats.get_chr_coverage(mutation.chr_id);
                }
                sample_coverage.insert({GenomicPosition(mutation),
                                        chr_coverage.get_coverage(mutation)});
            }
        }
    }

    return sample_locus_coverage;
}

template<typename K, typename V>
bool have_different_keys(const std::map<K, V>& map_a, const std::map<K, V>& map_b)
{
    if (map_a.size() != map_b.size()) {
        return true;
    }

    auto it_a = map_a.begin();
    auto it_b = map_b.begin();

    while (it_a != map_a.end()) {
        if (it_a->first != it_b->first) {
            return true;
        }

        ++it_a; ++it_b;
    }

    return false;
}

bool SampleSetStatistics::is_canonical() const
{
    if (stats_map.size()<2) {
        return true;
    }

    auto it = stats_map.begin();
    const auto& first_SID_data = it->second.get_data();

    for (++it; it != stats_map.end(); ++it) {
        if (have_different_keys(first_SID_data, it->second.get_data())) {
            return false;
        }
    }

    return true;
}

const SampleStatistics& SampleSetStatistics::add_chr_statistics(const std::string& sample_name,
                                                                const ChrSampleStatistics& chr_statistics,
                                                                const bool coverage_save)
{
    auto found = stats_map.find(sample_name);

    if (found == stats_map.end()) {
        found = stats_map.insert({sample_name, SampleStatistics(data_dir, sample_name, coverage_save)}).first;
    }

    found->second.add_chr_statistics(chr_statistics);

    return found->second;
}

const SampleStatistics& SampleSetStatistics::add_statistics(const SampleStatistics& sample_statistics)
{
    const auto& name = sample_statistics.get_sample_name();
    auto found = stats_map.find(name);

    if (found != stats_map.end()) {
        throw std::domain_error("The current object already contains a statistics for sample.");
    }

    return stats_map.insert({name, sample_statistics}).first->second;
}

const SampleStatistics& SampleSetStatistics::operator[](const std::string& sample_name) const
{
    auto found = find(sample_name);
    if (found == end()) {
        throw std::runtime_error("This object does not contain data about of "
                                 "the sample \"" + sample_name + "\".");
    }

    return found->second;
}

std::list<std::string> SampleSetStatistics::get_sample_names() const
{
    std::list<std::string> sample_names;

    for (const auto& [sample_name, sample_statistics] : stats_map) {
        sample_names.push_back(sample_name);
    }

    return sample_names;
}

void SampleSetStatistics::save_VAF_CSVs(const std::string& base_name,
                                        std::ostream& progress_bar_stream,
                                        const bool& quiet) const
{
    CLONES::UI::ProgressBar progress_bar(progress_bar_stream, quiet);

    size_t chr_processes{0};
    for (const auto& chr_id : repr_chr_ids) {
        std::string chr_str = GenomicPosition::chrtos(chr_id);

        progress_bar.set_progress(100*chr_processes/repr_chr_ids.size(),
                                  "Saving chr. " + chr_str + " VAFs");

        save_VAF_CSV(base_name + chr_str + ".csv", chr_id);

        ++chr_processes;
    }
    progress_bar.set_progress(100, "VAFs saved");
}

std::map<SID, SIDData>
get_total_SID_data(const std::map<std::string, SampleStatistics>& stats_map)
{
    std::map<SID, SIDData> total_SID_data;

    for (const auto& [sample_name, sample_stats]: stats_map) {
        update_data(total_SID_data, sample_stats.get_data());
    }

    return total_SID_data;
}

std::map<GenomicPosition, BaseCoverage>
get_total_SID_coverage(const std::map<std::string, SampleStatistics>& stats_map,
                       const std::map<SID, SIDData>& total_SID_data)
{
    std::map<GenomicPosition, BaseCoverage> total_locus_coverage;

    for (const auto& [mutation, total_mutation_data]: total_SID_data) {
        auto it = total_locus_coverage.insert({mutation,0}).first;
        auto& coverage = it->second;
        for (const auto& [sample_name, sample_stats]: stats_map) {
            coverage += sample_stats.get_coverage(mutation);
        }
    }

    return total_locus_coverage;
}

void SampleSetStatistics::canonize()
{
    auto total_SID_data = get_total_SID_data(stats_map);

    for (auto& [sample_name, sample_stats]: stats_map) {
        ChrCoverage chr_coverage;

        for (auto& [mutation, total_data] : total_SID_data) {
            auto& SID_data = sample_stats.SID_data;

            // if the SID is not mentioned in the sample statistics
            if (SID_data.count(mutation)==0) {
                // add the SID among those mentioned by the sample
                // statistics
                SID_data.insert({mutation, {mutation, 0}});

                // if the coverage of the chromosome containing the SID
                // has not been loaded yet
                if (mutation.chr_id!=chr_coverage.get_chr_id()
                        || !chr_coverage.is_initialized()) {

                    // load it
                    chr_coverage = sample_stats.get_chr_coverage(mutation.chr_id);
                }

                // add the coverage for the added SID
                auto& locus_coverage = sample_stats.locus_coverage;
                locus_coverage.insert({GenomicPosition(mutation),
                                       chr_coverage.get_coverage(mutation)});
            }
        }
    }
}

std::ostream& print_join(std::ostream& os, const std::set<std::string>& S, const char& sep=';')
{
    if (S.size()>0) {
        auto S_it = S.begin();
        os << *S_it;
        while (++S_it != S.end()) {
            os << sep << *S_it;
        }
    }

    return os;
}

std::set<std::string> get_descriptions(const std::set<Mutation::Nature>& types)
{
    std::set<std::string> string_types;

    for (const auto& type: types) {
        string_types.insert(SID::get_nature_description(type));
    }

    return string_types;
}

void SampleSetStatistics::save_VAF_CSV(const std::filesystem::path& filename,
                                       const ChromosomeId& chr_id) const
{
    if (!is_canonical()) {
        throw std::runtime_error("Only canonical statistics can save CAV CSV. Canonize it.");
    }

    const char separator='\t';

    std::ofstream os(data_dir/filename);

    stream_VAF_csv_header(os, separator);

    if (stats_map.size()==0) {
        return;
    }

    auto total_SID_data = get_total_SID_data(stats_map);

    for (const auto& [mutation, total_mutation_data]: total_SID_data) {
        os << GenomicPosition::chrtos(chr_id) << separator << mutation.position
           << separator << mutation.position
           << separator << mutation.ref  << separator << mutation.alt
           << separator << (mutation.is_SBS()?"SNV":"indel") << separator;

        print_join(os, total_mutation_data.causes, ';');

        os << separator;

        auto descriptions = get_descriptions(total_mutation_data.nature_set);

        print_join(os, descriptions, ';');

        size_t total_mutation_coverage{0};
        for (const auto& [sample_name, sample_stats]: stats_map) {
            const auto& mutation_coverage = sample_stats.get_coverage(mutation);
            auto mutation_data = sample_stats.get_data(mutation);

            print_SID_data(os, mutation_coverage, mutation_data, separator);

            total_mutation_coverage += mutation_coverage;
        }

        if (stats_map.size()>1) {
            print_SID_data(os, total_mutation_coverage, total_mutation_data, separator);
        }

        os << std::endl;
    }
}

#if WITH_MATPLOT

std::vector<double> down_sample_coverage(const std::vector<BaseCoverage>& coverage, const size_t& new_size=300)
{
    std::vector<double> down_sampled;

    down_sampled.reserve(new_size);

    auto bin_size = coverage.size()/new_size;
    size_t base_in_bins{0};
    double sum{0};
    for (const auto& base_coverage: coverage) {
        sum += base_coverage;
        ++base_in_bins;

        if (base_in_bins == bin_size) {
            down_sampled.push_back(sum/bin_size);

            base_in_bins = 0;
            sum = 0;
        }
    }

    return down_sampled;
}

void SampleSetStatistics::save_coverage_images(const std::string& base_name,
                                               std::ostream& progress_bar_stream,
                                               const bool& quiet) const
{
    CLONES::UI::ProgressBar progress_bar(progress_bar_stream, quiet);

    size_t chr_processes{0};
    for (const auto& chr_id : repr_chr_ids) {
        std::string chr_str = GenomicPosition::chrtos(chr_id);

        progress_bar.set_progress(100*chr_processes/repr_chr_ids.size(),
                                  "Saving chr. " + chr_str + " coverage");

        save_coverage_image( base_name + chr_str + "_coverage.jpg", chr_id);

        ++chr_processes;
    }
    progress_bar.set_progress(100, "Coverage images saved");
}

void SampleSetStatistics::save_coverage_image(const std::filesystem::path& filename,
                                     const ChromosomeId& chromosome_id) const
{
    if (!is_canonical()) {
        throw std::runtime_error("Only canonical statistics can save coverage image. Canonize it.");
    }

    using namespace matplot;

    size_t num_of_plots = stats_map.size();
    const size_t new_size = 300;

    auto f = figure(true);
    f->size(1500, 1500);

    if (num_of_plots==1) {
        const auto& chr_stats = stats_map.begin()->second;

        const auto coverage = chr_stats.get_chr_coverage(chromosome_id).get_coverage_vector();

        std::vector<std::vector<double>> d_coverage = {down_sample_coverage(coverage, new_size)};

        area(d_coverage);

        long chr_size = coverage.size();
        xticks({1, new_size/2, new_size});
        xticklabels({"0",std::to_string(chr_size/2),std::to_string(chr_size)});

        xlabel("Pos. in sequence");
        ylabel("Coverage");
    } else {
        std::vector<std::vector<double>> total(1);

        tiledlayout(++num_of_plots, 1);

        long chr_size;
        for (const auto& [sample_name, sample_stats]: stats_map) {
            nexttile();

            const auto coverage = sample_stats.get_chr_coverage(chromosome_id).get_coverage_vector();

            chr_size = coverage.size();

            std::vector<std::vector<double>> d_coverage = {down_sample_coverage(coverage, new_size)};

            area(d_coverage);

            const auto latex_name = std::regex_replace(sample_name, std::regex("_"), "\\\\_");
            title("Sample \""+latex_name+"\"");

            xticks({1, new_size/2, new_size});
            xticklabels({"0",std::to_string(chr_size/2),std::to_string(chr_size)});

            ylabel("Coverage");
            if (total[0].size() != d_coverage[0].size()) {
                total[0] = d_coverage[0];
            } else {
                auto total_it = total[0].begin();
                auto d_coverage_it = d_coverage[0].begin();

                for (; total_it != total[0].end(); ++total_it,++d_coverage_it) {
                    *total_it += *d_coverage_it;
                }
            }
        }

        nexttile();

        area(total);
        title("Overall");

        xticks({1, new_size/2, new_size});
        xticklabels({"0",std::to_string(chr_size/2),std::to_string(chr_size)});

        xlabel("Pos. in sequence");
        ylabel("Coverage");
    }

    sgtitle("Chromosome " + GenomicPosition::chrtos(chromosome_id) + "'s coverage");

    f->save(data_dir/filename);
}

std::vector<double> get_VAF_data(const std::map<SID, SIDData>& SID_data,
                                 const std::map<GenomicPosition, BaseCoverage>& locus_coverage,
                                 const ChromosomeId& chr_id,
                                 const double& threshold=0.0)
{
    std::vector<double> VAF;

    for (const auto& [mutation, data]: SID_data) {
        if (mutation.chr_id == chr_id && data.num_of_occurrences > 0) {
            auto value = static_cast<double>(data.num_of_occurrences)/locus_coverage.at(mutation);

            if (value > threshold) {
                VAF.push_back(value);
            }
        }
    }

    return VAF;
}

inline
std::vector<double>
get_VAF_data(const SampleStatistics& sample_statistics, const ChromosomeId& chr_id,
             const double& threshold=0.0)
{
    return get_VAF_data(sample_statistics.get_data(), sample_statistics.get_coverage(),
                        chr_id, threshold);
}

void SampleSetStatistics::save_SID_histograms(const std::string& base_name,
                                              std::ostream& progress_bar_stream,
                                              const bool& quiet) const
{
    CLONES::UI::ProgressBar progress_bar(progress_bar_stream, quiet);

    size_t chr_processes{0};
    for (const auto& chr_id : repr_chr_ids) {
        std::string chr_str = GenomicPosition::chrtos(chr_id);

        progress_bar.set_progress(100*chr_processes/repr_chr_ids.size(),
                                  "Saving chr. " + chr_str + " histogram");

        save_SID_histogram( base_name + chr_str + "_hist.jpg", chr_id);

        ++chr_processes;
    }
    progress_bar.set_progress(100, "Histogram images saved");
}

void SampleSetStatistics::save_SID_histogram(const std::filesystem::path& filename,
                                             const ChromosomeId& chromosome_id) const
{
    if (!is_canonical()) {
        throw std::runtime_error("Only canonical statistics can save SID histogram. Canonize it.");
    }

    using namespace matplot;

    size_t num_of_bins = 50;
    size_t num_of_plots = stats_map.size();

    auto f = figure(true);
    f->size(1500, 1500);

    if (num_of_plots==1) {
        const auto& sample_stats = (stats_map.begin()->second);

        std::vector<double> VAF = get_VAF_data(sample_stats, chromosome_id, 0.15);

        if (VAF.size()>0) {
            hist(VAF,num_of_bins);
        }

        xlabel("VAF");
        ylabel("Num. of SIDs");
    } else {
        auto total_SID_data = get_total_SID_data(stats_map);

        auto total_locus_coverage = get_total_SID_coverage(stats_map,
                                                           total_SID_data);

        tiledlayout(++num_of_plots, 1);
        for (const auto& [sample_name, sample_stats]: stats_map) {
            nexttile();

            std::vector<double> VAF = get_VAF_data(sample_stats, chromosome_id, 0.15);

            if (VAF.size()>0) {
                hist(VAF,num_of_bins);
            }

            const auto latex_name = std::regex_replace(sample_name, std::regex("_"), "\\\\_");
            title("Sample \""+latex_name+"\"");

            ylabel("Num. of SIDs");
        }

        nexttile();

        std::vector<double> VAF = get_VAF_data(total_SID_data, total_locus_coverage,
                                               chromosome_id, 0.15);

        if (VAF.size()>0) {
            hist(VAF,num_of_bins);
        }
        title("Overall");

        xlabel("VAF");
        ylabel("Num. of SIDs");
    }

    sgtitle("Chromosome " + GenomicPosition::chrtos(chromosome_id)+"'s SID VAF");

    f->save(data_dir/filename);
}

#endif // WITH_MATPLOT

} // SequencingSimulator

} // Mutations

} // CLONES
