/**
 * @file tissue_sample.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements tissue samples
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

#include "tissue_sample.hpp"

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{

TissueSampleId TissueSample::counter = 0;

TissueSample::TissueSample():
    id(0), time(0), bounding_box({{0,0},{0,0}}), name("")
{}

TissueSample::TissueSample(const Time& time, const RectangleSet& bounding_box,
                           const size_t& tumour_cells_in_bbox):
    TissueSample(time, bounding_box, tumour_cells_in_bbox, {})
{}

TissueSample::TissueSample(const std::string& name, const CLONES::Time& time,
                           const CLONES::Mutants::RectangleSet& bounding_box,
                           const size_t& tumour_cells_in_bbox):
    TissueSample(name, time, bounding_box, tumour_cells_in_bbox, {})
{}

TissueSample::TissueSample(const Time& time, const RectangleSet& bounding_box,
                           const size_t& tumour_cells_in_bbox,
                           const std::list<CLONES::Mutants::CellId>& cell_ids):
    TissueSample("", time, bounding_box, tumour_cells_in_bbox, cell_ids)
{
    name = "S_"+std::to_string(id);
}

TissueSample::TissueSample(const std::string& name, const Time& time,
                           const RectangleSet& bounding_box,
                           const size_t& tumour_cells_in_bbox,
                           const std::list<CLONES::Mutants::CellId>& cell_ids):
    id(counter++), time(time), bounding_box(bounding_box),
    tumour_cells_in_bbox(tumour_cells_in_bbox), cell_ids(cell_ids),
    name(name)
{}

void TissueSample::add_cell_id(const CLONES::Mutants::CellId& cell_id)
{
    if (bounding_box.size() == cell_ids.size()) {
        throw std::domain_error("The sample already contains all the cell ids");
    }

    cell_ids.push_back(cell_id);
}

}   // Evolutions

}   // Mutants

}   // CLONES

std::ostream& operator<<(std::ostream& os, const CLONES::Mutants::Evolutions::TissueSample& tissue_sample)
{
    const CLONES::Mutants::RectangleSet& bounding_box =  tissue_sample.get_bounding_box();

    os << "# " << tissue_sample.get_time()
       << " " << tissue_sample.get_cell_ids().size()
       << " " << bounding_box.size()
       << " "<< bounding_box.lower_corner
       << " "<< bounding_box.upper_corner << std::endl;

    for (const auto& cell_id: tissue_sample.get_cell_ids()) {
        os << cell_id << std::endl;
    }

    return os;
}
