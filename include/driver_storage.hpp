/**
 * @file driver_storage.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines class to load and store driver mutations
 * @version 1.3
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

#ifndef __CLONES_DRIVER_STORAGE__
#define __CLONES_DRIVER_STORAGE__

#include <map>
#include <set>
#include <list>
#include <string>
#include <filesystem>

#include "sid.hpp"

namespace CLONES
{

namespace Mutations
{
class DriverStorage
{
public:
    /**
     * @brief The mutation entry type in the code-mutation map
     * 
     * This class represent the co-domain in the code-mutation map.
     * It maintains the mutation and the set of all the tumour type
     * for which it appears as a driver mutation.
     */
    struct MutationEntry
    {
        SID mutation;                       //!< A mutation
        std::set<std::string> tumour_types; //!< The set of tumours for which the mutation is a driver
    };

    /**
     * @brief The empty constructor
     */
    DriverStorage();

    /**
     * @brief Get the driver SID mutations
     *
     * @return a constant reference to a map associating a driver mutation
     *      code to the corresponding driver mutation entry
     */
    inline const std::map<std::string, MutationEntry>& get_code2mutation_map() const
    {
        return mutation_map;
    }

    /**
     * @brief Get the reverse map
     *
     * @return a map from the known driver mutations to their codes
     */
    std::map<SID, std::string> get_reverse_map() const;

    /**
     * @brief Get the SNV positions
     *
     * @return A list containing the genomic position of the SID mutations
     */
    std::list<GenomicPosition> get_mutation_positions() const;

    /**
     * @brief Get the source path of the drivers
     *
     * @return a constant reference to the driver source path
     */
    const std::filesystem::path& get_source_path() const
    {
        return source_path;
    }

    /**
     * @brief Load the driver SID mutations
     *
     * @param filename is the driver mutation filename
     * @return a map associating a driver code to the corresponding
     *      driver mutation
     */
    static DriverStorage load(const std::filesystem::path& filename);

private:
    std::map<std::string, MutationEntry> mutation_map;    //!< The code-mutation map

    std::filesystem::path source_path;  //!< The driver filename
};

}   // Mutations

}   // CLONES

#endif // __CLONES_DRIVER_MUTATION_STORAGE__