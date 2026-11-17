/**
 * @file binary_logger.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a binary simulation logger
 * @version 1.9
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

#include <sstream>
#include <ctime>
#include <limits>
#include <cmath> // log10 and ceil
#include <regex>

#include "binary_logger.hpp"
#include "tissue.hpp"

#include "tissue_simulation.hpp"
#include "utils.hpp"

#include "error.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

std::string snapshot_prefix = "snapshot";

BinaryLogger::BinaryLogger():
    BinaryLogger("clones_" + BasicLogger::get_time_string())
{
}

BinaryLogger::BinaryLogger(const std::filesystem::path& simulation_dir, const size_t& cells_per_file):
    BasicLogger(simulation_dir), cell_archive(),
    cells_per_file(cells_per_file), cell_in_current_file(0), current_file_number(0)
{
}

BinaryLogger::BinaryLogger(const BinaryLogger& orig):
    BasicLogger(orig),  cell_archive(),
    cells_per_file(orig.cells_per_file), cell_in_current_file(orig.cell_in_current_file),
    current_file_number(orig.current_file_number)
{
}

std::filesystem::path BinaryLogger::get_cell_archive_path(const std::filesystem::path& directory, const uint16_t& file_number)
{
    using namespace std;

    std::ostringstream oss;

    size_t digits = static_cast<size_t>(ceil(log10(numeric_limits<uint16_t>::max())));

    oss << "cells_" << std::setfill('0') << std::setw(digits) << file_number;

    return directory / (oss.str()+".dat") ;
}

std::filesystem::path BinaryLogger::get_snapshot_path() const
{
    using namespace std::chrono;

    const auto curr_time = system_clock::now().time_since_epoch();

    const auto epoch_ms = duration_cast<microseconds>(curr_time).count();

    std::ostringstream oss;

    oss << snapshot_prefix << "_" << std::hex << epoch_ms << ".dat";

    return directory / oss.str();
}

std::filesystem::path BinaryLogger::find_last_snapshot_in(const std::filesystem::path& directory)
{
    namespace fs = std::filesystem;

    if (!fs::exists(directory)) {
        throw Error<std::runtime_error>("The path \"" + to_string(directory)
                                        + "\" does not exist.");
    }

    if (!fs::is_directory(directory)) {
        throw Error<std::runtime_error>("\"" + to_string(directory)
                                        + "\" is not a directory.");
    }

    std::regex re = build_regex(to_string(directory / (snapshot_prefix+"_[a-zA-Z0-9]+.dat")));

    bool found{false};
    std::string last;
    for (const auto & entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            const auto e_string = to_string(entry.path());
            if (std::regex_match(e_string, re)) {
                if (!found || e_string.compare(last)>0) {
                    last = e_string;
                    found = true;
                }
            }
        }
    }

    if (!found) {
        throw Error<std::runtime_error>("No CLONES simulation snapshot in \""
                                        + to_string(directory) + "\".");
    }

    return last;
}

BinaryLogger& BinaryLogger::operator=(const BinaryLogger& orig)
{
    close();

    directory = orig.directory;
    cells_per_file = orig.cells_per_file;
    cell_in_current_file = orig.cell_in_current_file;
    current_file_number = orig.current_file_number;

    return *this;
}

inline std::streampos compute_bytes_per_cell()
{
    Cell cell;

    return Archive::Binary::ByteCounter::bytes_required_by(cell);
}

BinaryLogger::CellReader::CellReader(const std::filesystem::path& directory):
    directory(directory), cells_per_archive(0), bytes_per_cell(compute_bytes_per_cell()),
    last_file_number(0), number_of_cells(0), cell_archive_id(0), cell_archive()
{
    namespace fs = std::filesystem;

    if (!fs::is_directory(directory)) {
        throw Error<std::runtime_error>("The simulation directory \"" + to_string(directory)
                                        + "\" does not exist.");
    }

    const auto info_path = BinaryLogger::get_cell_file_info_path(directory);

    // if no cell was saved
    if (!fs::exists(info_path)) {
        if (fs::exists(get_cell_archive_path(directory, 1))) {
            throw Error<std::runtime_error>("Missing cell file information \""
                                            + to_string(info_path) + "\".");
        }

        return;
    }

    std::ifstream cell_file(info_path, std::fstream::in);
    if (!(cell_file >> cells_per_archive >> last_file_number) || cells_per_archive == 0) {
        throw Error<std::runtime_error>("Wrong cell file information in \""
                                        + to_string(info_path) + "\".");
    }

    // cell files are numbered from 1: 0 means that no cell file was created
    if (last_file_number == 0) {
        return;
    }

    {   // compute the total number of cells: all files except the last are full
        const auto filename = get_cell_archive_path(directory, last_file_number);

        if (!fs::is_regular_file(filename)) {
            throw Error<std::runtime_error>("Missing cell file \"" + to_string(filename) + "\".");
        }

        const auto lastfile_bytes = fs::file_size(filename);

        number_of_cells = static_cast<uint64_t>(lastfile_bytes/bytes_per_cell)
                            + cells_per_archive * (last_file_number-1);
    }
}

Cell BinaryLogger::CellReader::operator[](const CellId& cell_id)
{
    if (cell_id < Cell::first_tumour_cell_id()
            || cell_id - Cell::first_tumour_cell_id() >= number_of_cells) {
        throw Error<std::out_of_range>("The cell " + std::to_string(cell_id)
                                       + " has not been created.");
    }

    const auto file_number = get_cell_file_number(cell_id);

    if (cell_archive == nullptr || cell_archive_id != file_number) {
        cell_archive = std::make_unique<Archive::Binary::In>(get_cell_archive_path(directory,
                                                                                   file_number));
        cell_archive_id = file_number;
    }

    cell_archive->seekg(get_cell_pos(cell_id));

    auto cell = Cell::load(*cell_archive);

    if (cell.get_id()!=cell_id) {
        std::ostringstream oss;

        oss << "Wrong cell file format. Aiming cell id "
            << cell_id << " read " << cell << ".";
        throw Error<std::runtime_error>(oss.str());
    }

    return cell;
}

void BinaryLogger::save_cell_file_info() const
{
    // save the number of cells per file and the next file number to be used
    std::ofstream cell_file(get_cell_file_info_path(), std::fstream::out);
    cell_file << cells_per_file << " " << current_file_number << std::endl;
}

void BinaryLogger::rotate_cell_file()
{
    if (cell_archive.is_open()) {
        cell_archive.close();
    }

    ++current_file_number;
    save_cell_file_info();

    auto filename = get_cell_archive_path(directory, current_file_number);

    cell_in_current_file=0;
    cell_archive.open(filename, std::fstream::binary);
}

void BinaryLogger::record_cell(const CellInTissue& cell)
{
    if (!cell_archive.is_open()) {
        // cell files are numbered from 1: a logger that has not created any
        // cell file yet (`current_file_number==0`) must start a new one even
        // if the directory already exists
        if (std::filesystem::exists(directory) && current_file_number > 0) {
            auto filename = get_cell_archive_path(directory, current_file_number);

            cell_archive.open(filename, std::fstream::binary | std::fstream::app);

            if (!std::filesystem::exists(get_cell_file_info_path())) {
                save_cell_file_info();
            }
        } else {
            std::filesystem::create_directory(directory);

            rotate_cell_file();
        }
    }

    cell_archive & static_cast<const Cell&>(cell);

    if (++cell_in_current_file>=cells_per_file) {
        rotate_cell_file();
    }
}

void BinaryLogger::record(const CellEventType& type, const CellInTissue& cell, const Time& time)
{
    (void)time;

    if (type==CellEventType::DUPLICATION ||
            type==CellEventType::DUP_AND_EPI_SWITCH ||
            type==CellEventType::MUTATION) {

        record_cell(cell);
    }
}

void BinaryLogger::record_initial_cell(const CellInTissue& cell)
{
    if (cell.get_id() != cell.get_parent_id()) {
        throw Error<std::domain_error>("The provided cell is not an initial cell.");
    }

    record_cell(cell);
}

std::filesystem::path BinaryLogger::snapshot(TissueSimulation& simulation)
{
    flush_archives();

    const auto snapshot_path = get_snapshot_path();

    {
        Archive::Binary::Out archive(snapshot_path);

        simulation.snapshot_info.emplace_back(simulation, snapshot_path);

        archive & simulation;

        archive.flush();
    }

    return snapshot_path;
}

void BinaryLogger::reset(const std::filesystem::path& output_directory)
{
    close();

    directory = output_directory;
    cell_in_current_file = 0;
    current_file_number = 0;
}

BinaryLogger::~BinaryLogger()
{
    close();
}

}   // Evolutions

}   // Mutants

}   // CLONES
