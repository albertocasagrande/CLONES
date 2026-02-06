/**
 * @file common.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines auxiliary classes and functions for executables
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

#include <map>
#include <string>
#include <filesystem>

#include <boost/program_options.hpp>

#include "simulation.hpp"
#include "json_config.hpp"

/**
 * @brief A basic class for CLONES executables
 */
class BasicExecutable
{
protected:
    const std::string program_name;         //!< The program name
    std::map<std::string, boost::program_options::options_description> visible_options; //!< the program visible option sections

    boost::program_options::options_description hidden_options;  //!< the program hidden options
    boost::program_options::positional_options_description positional_options;  //!< the program positional options
public:
    /**
     * @brief Print the program help message and exit
     *
     * @param message is the message to print before the help message
     * @param err_code is the exit code
     */
    void print_help_and_exit(const std::string& message, const int& err_code) const;

    /**
     * @brief Print the program help message and exit
     *
     * @param err_code is the exit code
     */
    void print_help_and_exit(const int& err_code) const;

    /**
     * @brief Construct a new BasicExecutable object
     *
     * @param program_name is the program name
     * @param option_section_names is a list of pairs section name-description
     */
    BasicExecutable(const std::string& program_name,
                    const std::vector<std::pair<std::string,std::string>>& option_sections);

    /**
     * @brief Find the last species simulation snapshot path
     *
     * @param simulation_dir is the simulation path provided as an executable parameter
     * @param parameter_name is the executable parameter name
     * @return the last simulation snapshot in the directory `simulation_dir`
     */
    std::filesystem::path get_last_snapshot_path(const std::string& simulation_dir,
                                                 const std::string& parameter_name);

    /**
     * @brief Load a species simulation snapshot
     *
     * @param snapshot_path is the path to the snapshot to be loaded
     * @param quiet is a Boolean flag to switch on/off (false/true) descriptive messages
     * @return the loaded simulation snapshot
     */
    static CLONES::Mutants::Evolutions::Simulation
    load_species_simulation(const std::filesystem::path snapshot_path, const bool quiet=false);
};
