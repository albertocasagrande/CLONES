/**
 * @file id_signature.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines indel signature
 * @version 1.1
 * @date 2025-11-14
 *
 * @copyright Copyright (c) 2023-2025
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

#ifndef __RACES_ID_SIGNATURE__
#define __RACES_ID_SIGNATURE__

#include <string>
#include <functional> // std::less
#include <iostream>
#include <sstream>

#include "mutation.hpp"
#include "signature.hpp"

#include "id_context.hpp"

namespace RACES
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
class IDType : public MutationType, public IDContext
{
    bool insertion;     //!< A Boolean flag establishing whether the mutation is an insertion

public:

    /**
     * @brief The empty constructor
     */
    IDType();

    /**
     * @brief A constructor
     *
     * @param fragment_type is the fragment type
     * @param first_level_code is the first level code
     * @param first_level_code is the second level code
     * @param insertion is a Boolean flag establishing whether the mutation is an insertion
     */
    IDType(const FragmentType& fragment_type, const FirstLevelType& first_level_code,
           const SecondLevelType& second_level_code, const bool& insertion);

    /**
     * @brief A constructor
     *
     * An indel context is conventionally represented by a string in the form
     * `{number}:{"Del"|"Ins"}:{'A','C','G','T','R','M'}:{number}`. The first
     * number is the unit size for both homopolymer and heteropolymer, and the
     * fragment size for microhomology. When the string represents a
     * heteropolymer or a microhomology it is the first level index.
     * The strings "Del" and "Ins" represent deletions and insertions,
     * respectively. The character between the second and the third ':' denotes
     * the fragment type: 'R' stands for heteropolymer, 'M' for microhomology,
     * and 'A', 'C', 'G', and 'T' for homopolymer. Dealing with homopolymer,
     * the character is the first level index.
     * The last number is the size of the microhomology (for microhomologies)
     * and the number of repetitions in the reference sequence (for polymer
     * insertions).
     *
     * @param type is the textual representation of an ID type
     */
    explicit IDType(const std::string& type);

    inline bool is_insertion() const
    {
        return insertion;
    }

    inline bool is_deletion() const
    {
        return !insertion;
    }

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
        return static_cast<const IDContext&>(*this)
                    == static_cast<const IDContext&>(type)
                && insertion == type.insertion;
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

}   // RACES


namespace std
{

template<>
struct less<RACES::Mutations::IDType>
{
    bool operator()(const RACES::Mutations::IDType &lhs,
                    const RACES::Mutations::IDType &rhs) const;
};

/**
 * @brief Stream the SBS type in a stream
 *
 * @param out is the output stream
 * @param type is the ID type to stream
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const RACES::Mutations::IDType& type);

/**
 * @brief Stream the ID type from a stream
 *
 * @param in is the input stream
 * @param type is the object where the streamed ID type will be placed
 * @return a reference to the input stream
 */
std::istream& operator>>(std::istream& in, RACES::Mutations::IDType& type);

}  // std

namespace RACES
{

namespace Mutations
{

using IDSignature = Signature<IDType>;

}   // Mutations

}   // RACES


#endif // __RACES_ID_SIGNATURE__