/**
 * @file id_signature.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines indel signature
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

#ifndef __CLONES_ID_SIGNATURE__
#define __CLONES_ID_SIGNATURE__

#include <string>
#include <functional> // std::less
#include <iostream>
#include <sstream>

#include "mutation.hpp"
#include "signature.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief A class to represent indel type
 *
 * An indel type consists in a fragment type (either HOMOPOLYMER, HETEROPOLYMER,
 * or MICROHOMOLOGY), a first level index that is either a base (for HOMOPOLYMER),
 * or the size of the unit (for HETEROPOLYMER and MICROHOMOLOGY) and a second
 * level index that is either the number of repetitions (for HOMOPOLYMER and
 * HETEROPOLYMER) or the homology size (for MICROHOMOLOGY), and a Boolean flag
 * establishing whether the mutation is an insertion or a deletion.
 */
class IDType : public MutationType
{
public:
    using RepetitionType = uint16_t;

    /**
     * @brief The type of the fragment involved in the mutation
     */
    enum class FragmentType
    {
        HOMOPOLYMER,    //!< A repeated sequence whose nucleotides are the same
        HETEROPOLYMER,  //!< A repeated sequence whose nucleotides may differ
        MICROHOMOLOGY  //!< A fragment followed by a sequence that matches its prefix
    };

    FragmentType ftype;         //!< The type of mutated fragment

    uint8_t fl_index;           //!< The first level index
    RepetitionType sl_index;    //!< The second level index

    bool insertion;     //!< A Boolean flag establishing whether the mutation is an insertion

    /**
     * @brief The empty constructor
     */
    IDType();

    /**
     * @brief A constructor
     *
     * @param fragment_type is the fragment type
     * @param first_level_index is the first level index
     * @param second_level_index is the second level index
     * @param insertion is a Boolean flag establishing whether the mutation is an insertion
     */
    IDType(const FragmentType& fragment_type, const uint8_t& first_level_index,
           const RepetitionType& second_level_index, const bool& insertion);

    /**
     * @brief A constructor
     *
     * An indel type is conventionally represented by a string in the
     * form `{number}:{"Del"|"Ins"}:{'C':'T':'R':'M'}:{number}`. The
     * first number is the unit size for both HOMOPOLYMER and
     * HETEROPOLYMER, and the fragment size for MICROHOMOLOGY. When
     * the string represents a HETEROPOLYMER or a MICROHOMOLOGY it
     * is the first level index.
     * The strings "Del" and "Ins" represent deletions and insertions,
     * respectively. The character between the second and the third
     * ':' denotes the fragment type: 'R' stands for HETEROPOLYMER,
     * 'M' for MICROHOMOLOGY, and both 'C' and 'T' for HOMOPOLYMER.
     * Dealing with HOMOPOLYMER, the character is the first level
     * index.
     * The last number is the size of the MICROHOMOLOGY (for
     * microhomologies), the number of repetitions in the reference
     * sequence (for polymer insertions), or the number of repetitions in
     * the altered sequence (for polymer deletions).
     *
     * @param type is the textual representation of an ID type
     */
    explicit IDType(const std::string& type);

    /**
     * @brief Get the type of the mutation type
     *
     * @return the type of the mutation type (i.e.,
     *      `MutationType::Type::INDEL`)
     */
    static inline constexpr Type type()
    {
        return Type::INDEL;
    }

    /**
     * @brief Get the mutation type name
     *
     * @return the mutation type name
     */
    static inline const std::string name()
    {
        return "indel";
    }

    /**
     * @brief Test whether two ID types are equivalent
     *
     * @param type is the ID type to compare
     * @return `true` if and only if the two ID types are equivalent
     */
    inline bool operator==(const IDType& type) const
    {
        return (ftype == type.ftype) && (fl_index == type.fl_index)
                && (sl_index == type.sl_index)
                && (insertion == type.insertion);
    }

    /**
     * @brief Test whether two ID types differ
     *
     * @param type is the ID type to compare
     * @return `true` if and only if the two ID types differ
     */
    inline bool operator!=(const IDType& type) const
    {
        return !(*this == type);
    }
private:
    /**
     * @brief Read a size from a string
     *
     * @tparam TYPE is the aimed type
     * @param num_str is the string representation of the string
     * @param type is the string ID type representation from which `num_str` comes from
     * @return a `TYPE` interpretation of the string `num_str`
     */
    template<typename TYPE>
    static TYPE read_size(const std::string& num_str, const std::string& type="")
    {
        try {
            int num = stoi(num_str);

            if (num >= 0 && num <= std::numeric_limits<TYPE>::max()) {
                return static_cast<TYPE>(num);
            }
        } catch (std::exception& e) {
        }

        std::ostringstream oss;
        if (type != "") {
            oss << "\"" << type << "\" does not represent an ID type: ";
        }
        oss << "\"" << num_str << "\" should be a number in the interval [0,"
            << std::numeric_limits<TYPE>::max() << "].";
        throw std::domain_error(oss.str());
    }
};

}   // Mutations

}   // CLONES


namespace std
{

template<>
struct less<CLONES::Mutations::IDType>
{
    bool operator()(const CLONES::Mutations::IDType &lhs,
                    const CLONES::Mutations::IDType &rhs) const;
};

/**
 * @brief Stream the SBS type in a stream
 *
 * @param out is the output stream
 * @param type is the ID type to stream
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const CLONES::Mutations::IDType& type);

/**
 * @brief Stream the ID type from a stream
 *
 * @param in is the input stream
 * @param type is the object where the streamed ID type will be placed
 * @return a reference to the input stream
 */
std::istream& operator>>(std::istream& in, CLONES::Mutations::IDType& type);

}  // std

namespace CLONES
{

namespace Mutations
{

using IDSignature = Signature<IDType>;

}   // Mutations

}   // CLONES


#endif // __CLONES_ID_SIGNATURE__