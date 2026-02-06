/**
 * @file tissue_sample.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines tissue samples
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

#ifndef __CLONES_TISSUE_SAMPLE__
#define __CLONES_TISSUE_SAMPLE__

#include <iostream>

#include <limits>
#include <list>
#include <memory>

#include "position_set.hpp"
#include "cell.hpp"
#include "time.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

using TissueSampleId = uint32_t;

#define NORMAL_SAMPLE_ID std::numeric_limits<Mutants::Evolutions::TissueSampleId>::max()

/**
 * @brief A class to represent tissue samples
 */
class TissueSample
{
    static TissueSampleId counter;   //!< A counter for the non-null samples

    TissueSampleId  id;         //!< The sample identifier
    CLONES::Time time;           //!< The sampling time
    RectangleSet bounding_box;  //!< The sample bounding box
    size_t tumour_cells_in_bbox; //!< The number of tumour cells in the bounding box
    std::shared_ptr<std::list<CellId>> cell_ids; //!< The list of cell identifier

    std::string name;           //!< The name of the sample
public:

    /**
     * @brief An empty constructor
     */
    TissueSample();

    /**
     * @brief Construct a new cell sample
     *
     * @param time is the sampling time
     * @param bounding_box is the sample bounding box
     * @param tumour_cells_in_bbox is the number of tumour cells in
     *      the bounding box
     */
    TissueSample(const CLONES::Time& time,
                 const CLONES::Mutants::RectangleSet& bounding_box,
                 const size_t& tumour_cells_in_bbox);

    /**
     * @brief Construct a new cell sample
     *
     * @param name is the sample name
     * @param time is the sampling time
     * @param bounding_box is the sample bounding box
     * @param tumour_cells_in_bbox is the number of tumour cells in
     *      the bounding box
     */
    TissueSample(const std::string& name, const CLONES::Time& time,
                 const CLONES::Mutants::RectangleSet& bounding_box,
                 const size_t& tumour_cells_in_bbox);

    /**
     * @brief Construct a new cell sample
     *
     * @param time is the sampling time
     * @param bounding_box is the sample bounding box
     * @param tumour_cells_in_bbox is the number of tumour cells in
     *      the bounding box
     * @param cell_ids is a list of the identifiers of the sample cells
     */
    TissueSample(const CLONES::Time& time,
                 const CLONES::Mutants::RectangleSet& bounding_box,
                 const size_t& tumour_cells_in_bbox,
                 const std::list<CLONES::Mutants::CellId>& cell_ids);

    /**
     * @brief Construct a new cell sample
     *
     * @param time is the sampling time
     * @param bounding_box is the sample bounding box
     * @param tumour_cells_in_bbox is the number of tumour cells in
     *      the bounding box
     * @param cell_ids is a list of the identifiers of the sample cells
     */
    TissueSample(const CLONES::Time& time,
                 const CLONES::Mutants::RectangleSet& bounding_box,
                 const size_t& tumour_cells_in_bbox,
                 std::list<CLONES::Mutants::CellId>&& cell_ids);

    /**
     * @brief Construct a new cell sample
     *
     * @param name is the sample name
     * @param time is the sampling time
     * @param bounding_box is the sample bounding box
     * @param tumour_cells_in_bbox is the number of tumour cells in
     *      the bounding box
     * @param cell_ids is a list of the identifiers of the sample cells
     */
    TissueSample(const std::string& name, const CLONES::Time& time,
                 const CLONES::Mutants::RectangleSet& bounding_box,
                 const size_t& tumour_cells_in_bbox,
                 const std::list<CLONES::Mutants::CellId>& cell_ids);

    /**
     * @brief Construct a new cell sample
     *
     * @param name is the sample name
     * @param time is the sampling time
     * @param bounding_box is the sample bounding box
     * @param tumour_cells_in_bbox is the number of tumour cells in
     *      the bounding box
     * @param cell_ids is a list of the identifiers of the sample cells
     */
    TissueSample(const std::string& name, const CLONES::Time& time,
                 const CLONES::Mutants::RectangleSet& bounding_box,
                 const size_t& tumour_cells_in_bbox,
                 std::list<CLONES::Mutants::CellId>&& cell_ids);

    /**
     * @brief Add a cell id among those in the sample
     *
     * @param cell_id is the cell id to be added in the sample
     */
    void add_cell_id(const CLONES::Mutants::CellId& cell_id);

    /**
     * @brief Set the sample name
     *
     * @param name is the new sample name
     * @return a constant reference to the updated sample name
     */
    inline const std::string& set_name(const std::string& name)
    {
        return (this->name = name);
    }

    /**
     * @brief Get the sample name
     *
     * @return a constant reference to the sample name
     */
    inline const std::string& get_name() const
    {
        return name;
    }

    /**
     * @brief Get the sample id
     *
     * @return a constant reference to the sample id
     */
    inline const TissueSampleId& get_id() const
    {
        return id;
    }

    /**
     * @brief Get the sampling time
     *
     * @return a constant reference to the sampling time
     */
    inline const CLONES::Time& get_time() const
    {
        return time;
    }

    /**
     * @brief Get the sample bounding box
     *
     * @return a constant reference to the sample bounding box
     */
    inline const CLONES::Mutants::RectangleSet& get_bounding_box() const
    {
        return bounding_box;
    }

    /**
     * @brief Get the number of tumour cells in the bounding box
     *
     * @return a constant reference to the number of tumour cells
     *      in the bounding box
     */
    inline const size_t& get_tumour_cells_in_bbox() const
    {
        return tumour_cells_in_bbox;
    }

    /**
     * @brief Get the list of the sampled cell identifiers
     *
     * @return a constant reference to the list of the sampled cell identifiers
     */
    inline const std::list<CLONES::Mutants::CellId>& get_cell_ids() const
    {
        return *cell_ids;
    }

    /**
     * @brief Save a tissue sample in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & time
                & bounding_box
                & tumour_cells_in_bbox
                & id
                & *cell_ids
                & name;
    }

    /**
     * @brief Load a tissue sample from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the tissue sample
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static TissueSample load(ARCHIVE& archive)
    {
        Time time;

        archive & time;

        RectangleSet bounding_box = RectangleSet::load(archive);

        size_t tumour_cells_in_bbox;
        archive & tumour_cells_in_bbox;

        TissueSample tissue_sample(time, bounding_box, tumour_cells_in_bbox);

        archive & tissue_sample.id;

        if (tissue_sample.id+1>TissueSample::counter) {
            TissueSample::counter = tissue_sample.id+1;
        }

        archive & *(tissue_sample.cell_ids)
                & tissue_sample.name;

        return tissue_sample;
    }
};

}   // Evolutions

}   // Mutants

}   // CLONES

/**
 * @brief Write the tissue sample in an output stream
 *
 * @param os is the output stream
 * @param tissue_sample is the cell sample to be streamed
 * @return a reference to the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::TissueSample& tissue_sample);

#endif // __CLONES_TISSUE_SAMPLE__
