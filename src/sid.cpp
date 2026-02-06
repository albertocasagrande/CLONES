/**
 * @file sid.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements SNV, Insertion, and Deletion mutations
 * @version 1.2
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

#include "sid.hpp"

namespace CLONES
{

namespace Mutations
{

SID::SID():
    Mutation()
{}

inline void decode_empty_sequence(std::string& sequence)
{
    if (sequence == "-") {
        sequence = "";
    }
}

inline void decode_empty_sequences(std::string& ref, std::string& alt)
{
    decode_empty_sequence(ref);
    decode_empty_sequence(alt);
}

SID::SID(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const char ref_base, const char alt_base, const Mutation::Nature& nature):
    SID(chromosome_id, chromosomic_position, ref_base, alt_base, "", nature)
{}

SID::SID(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const char ref_base, const char alt_base, const std::string& cause,
         const Mutation::Nature& nature):
    Mutation(chromosome_id, chromosomic_position, nature, cause),
    ref(1, ref_base), alt(1, alt_base)
{
    decode_empty_sequences(ref, alt);
}

SID::SID(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const std::string& ref, const std::string& alt,
         const Mutation::Nature& nature):
    Mutation(chromosome_id, chromosomic_position, nature), ref(ref), alt(alt)
{
    decode_empty_sequences(this->ref, this->alt);
}

SID::SID(const GenomicPosition& genomic_position,
         const std::string& ref, const std::string& alt,
         const Mutation::Nature& nature):
    Mutation(genomic_position, nature), ref(ref), alt(alt)
{
    decode_empty_sequences(this->ref, this->alt);
}

SID::SID(GenomicPosition&& genomic_position,
         const std::string& ref, const std::string& alt,
         const Mutation::Nature& nature):
    Mutation(std::move(genomic_position), nature), ref(ref), alt(alt)
{
    decode_empty_sequences(this->ref, this->alt);
}

SID::SID(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const std::string& ref, const std::string& alt,
         const std::string& cause, const Mutation::Nature& nature):
    Mutation(chromosome_id, chromosomic_position, nature, cause), ref(ref), alt(alt)
{
    decode_empty_sequences(this->ref, this->alt);
}

SID::SID(const GenomicPosition& genomic_position,
         const std::string& ref, const std::string& alt,
         const std::string& cause, const Mutation::Nature& nature):
    Mutation(genomic_position, nature, cause), ref(ref), alt(alt)
{
    decode_empty_sequences(this->ref, this->alt);
}

SID::SID(GenomicPosition&& genomic_position,
         const std::string& ref, const std::string& alt,
         const std::string& cause, const Mutation::Nature& nature):
    Mutation(std::move(genomic_position), nature, cause), ref(ref), alt(alt)
{
    decode_empty_sequences(this->ref, this->alt);
}

}   // Mutations

}   // CLONES

namespace std
{

bool less<CLONES::Mutations::SID>::operator()(const CLONES::Mutations::SID &lhs,
                                             const CLONES::Mutations::SID &rhs) const
{
    {
        less<CLONES::Mutations::GenomicPosition> gp_op;

        if (gp_op(lhs, rhs)) {
            return true;
        }

        if (gp_op(rhs, lhs)) {
            return false;
        }
    }

    const auto ref_cmp = lhs.ref.compare(rhs.ref);

    if (ref_cmp < 0) {
        return true;
    }

    if (ref_cmp > 0) {
        return false;
    }

    return lhs.alt.compare(rhs.alt) < 0;
}

bool operator==(const CLONES::Mutations::SID &lhs,
                const CLONES::Mutations::SID &rhs)
{
    return lhs.chr_id==rhs.chr_id && lhs.position==rhs.position
            && lhs.nature==rhs.nature && lhs.cause==rhs.cause
            && lhs.ref==rhs.ref && lhs.alt==rhs.alt;
}

std::ostream& operator<<(std::ostream& out, const CLONES::Mutations::SID& sid)
{
      out << static_cast<CLONES::Mutations::GenomicPosition>(sid)
        << "[" << sid.ref << ">" <<  sid.alt << "]";

    return out;
}

}   // std
