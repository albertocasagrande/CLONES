/**
 * @file sample_specification.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines sample specification
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

#ifndef __CLONES_SAMPLE_SPECIFICATION__
#define __CLONES_SAMPLE_SPECIFICATION__

#include "position_set.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

/**
 * @brief Representation for a sample specification
 */
class SampleSpecification
{
public:
    /**
     * @brief A constructor
     *
     * @param bounding_box is the sample bounding box
     */
    explicit SampleSpecification(const RectangleSet& bounding_box);

    /**
     * @brief A constructor
     *
     * @param bounding_box is the sample bounding box
     * @param num_of_cells is the maximum number of cells to sample
     */
    SampleSpecification(const RectangleSet& bounding_box, const size_t num_of_cells);

    /**
     * @brief A constructor
     *
     * @param name is the specification name
     * @param bounding_box is the sample bounding box
     */
    SampleSpecification(const std::string& name, const RectangleSet& bounding_box);

    /**
     * @brief A constructor
     *
     * @param name is the specification name
     * @param bounding_box is the sample bounding box
     * @param num_of_cells is the maximum number of cells to sample
     */
    SampleSpecification(const std::string& name, const RectangleSet& bounding_box,
                        const size_t num_of_cells);

    /**
     * @brief Set the specification name
     *
     * @param name is the new specification name
     * @return a constant reference to the updated specification name
     */
    inline const std::string& set_name(const std::string& name)
    {
        return (this->name = name);
    }

    /**
     * @brief Get the sample specification name
     *
     * @return a constant reference to the sample specification name
     */
    inline const std::string& get_name() const
    {
        return name;
    }

    /**
     * @brief Get the sample specification bounding box
     *
     * @return a constant reference to the sample specification bounding box
     */
    inline const RectangleSet& get_bounding_box() const
    {
        return bbox;
    }

    /**
     * @brief Get the number of cells to be sampled
     *
     * @return a constant reference to the he number of cells
     *      to be sampled
     */
    inline const size_t& get_num_of_cells() const
    {
        return num_of_cells;
    }

    /**
     * @brief Save a sample specification in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & name
                & bbox
                & num_of_cells;
    }

    /**
     * @brief Load a sample specification from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded sample specification
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SampleSpecification load(ARCHIVE& archive)
    {
        std::string name;
        RectangleSet bbox;
        size_t num_of_cells;

        archive & name
                & bbox
                & num_of_cells;

        return {name, bbox, num_of_cells};
    }

protected:
    std::string name;       //!< The name of the sample
    RectangleSet bbox;      //!< The set of tissue region to be sampled
    size_t num_of_cells;    //!< The number of cells to sample

    /**
     * @brief Get the default sample specification name for a bounding box
     *
     * @param bounding_box is the sample bounding box
     * @param num_of_cells is the maximum number of cells to sample
     * @return the default name for the sample specification
     */
    static std::string get_default_name(const RectangleSet& bounding_box,
                                        const size_t num_of_cells);
};

}   // Evolutions

}   // Mutants

}   // CLONES


/**
 * @brief Test the equivalence between two sample specifications
 *
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two sample specifications
 *      deal with the same name, the same bounding box, and the same number
 *      of cells to be sampled
 */
inline
bool operator==(const CLONES::Mutants::Evolutions::SampleSpecification& lhs,
                const CLONES::Mutants::Evolutions::SampleSpecification& rhs)
{
    return (lhs.get_name() == rhs.get_name() && lhs.get_bounding_box() == rhs.get_bounding_box()
                && lhs.get_num_of_cells() == rhs.get_num_of_cells());
}

#endif // __CLONES_RATE_UPDATE__
