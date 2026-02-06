/**
 * @file genomic_region.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements genomic region
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

#include <iostream>

#include "genomic_region.hpp"

#include "sid.hpp"

namespace CLONES
{

namespace Mutations
{

using Length = GenomicRegion::Length;

GenomicRegion::GenomicRegion():
    initial_pos(0), length(0)
{}

GenomicRegion::GenomicRegion(const ChromosomeId chromosome_id, const Length length):
    initial_pos(chromosome_id, 1), length(length)
{
    if (length == 0) {
        throw std::domain_error("the genomic region length must greater than 0");
    }
}

GenomicRegion::GenomicRegion(const GenomicPosition initial_pos, const Length length):
    initial_pos(initial_pos), length(length)
{
    if (length == 0) {
        throw std::domain_error("the genomic region length must greater than 0");
    }
}

bool GenomicRegion::overlaps(const GenomicRegion& genomic_region) const
{
    return (this->contains(genomic_region.get_initial_position())
            || genomic_region.contains(get_initial_position())
            || this->contains(genomic_region.get_final_position())
            || genomic_region.contains(get_final_position()));
}

GenomicRegion GenomicRegion::split(const GenomicPosition& split_point)
{
    if (!contains(split_point)) {
        throw std::domain_error("the fragment does not contain the split point");
    }

    if (split_point==initial_pos) {
        throw std::domain_error("the fragment initial position and the split point are the same");
    }

    GenomicRegion new_region(split_point, length-(split_point.position-initial_pos.position));

    length = split_point.position - initial_pos.position;

    return new_region;
}

GenomicRegion& GenomicRegion::join(GenomicRegion& contiguous_region)
{
    if (follows(contiguous_region)) {
        std::swap(contiguous_region, *this);
    }

    if (!precedes(contiguous_region)) {
        throw std::domain_error("the two genomic regions are not contiguous");
    }

    length += contiguous_region.length;

    // reset contiguous_region length
    contiguous_region.length = 0;

    return *this;
}

bool GenomicRegion::contains(const SID& mutation) const
{
    const auto mutation_region = mutation.get_region();

    return contains(mutation_region);
}

bool GenomicRegion::strictly_contains(const SID& mutation) const
{
    const auto mutation_region = mutation.get_region();

    return strictly_contains(mutation_region);
}

}   // Mutations

}   // CLONES

namespace std
{

std::ostream& operator<<(std::ostream& out, const CLONES::Mutations::GenomicRegion& genomic_region)
{
    auto begin_pos = genomic_region.get_initial_position();

    out << "GenomicRegion(chr: " << CLONES::Mutations::GenomicPosition::chrtos(begin_pos.chr_id)
        << ", begin: " << begin_pos.position
        << ", length: " << genomic_region.size() << ")";

    return out;
}

}
