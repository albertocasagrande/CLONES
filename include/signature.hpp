/**
 * @file signature.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines mutational signatures
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

#ifndef __CLONES_SIGNATURE__
#define __CLONES_SIGNATURE__

#include <string>
#include <map>
#include <set>
#include <istream>
#include <sstream>

namespace CLONES
{

namespace Mutations
{


template<typename MUTATION_TYPE>
class Signature;

/**
 * @brief A class to represent the result of a signature expression
 *
 * This class is meant to represent temporary object that are evaluated during
 * the computation of expressions of the kind:
 *
 * \f$\alpha_1 * \beta_1 + alpha_2 * beta_2 + \ldots\f$
 *
 * where the \f$\alpha_i\f$'s are real values in the interval \f$[0,1]\f$ and
 * the \f$\beta_i\f$'s are `Signature` objects.
 *
 * Even if the final result of above expression is a signature, the partial
 * results may be different from a probability distribution and that is why
 * this class is needed.
 *
 * @tparam MUTATION_TYPE is the mutation type
 */
template<typename MUTATION_TYPE>
class SignatureExprResult
{
    std::map<MUTATION_TYPE, double> value_map; //!< the mutation type-value map

    /**
     * @brief The constructor
     *
     * This constructor is private and it is meant to be exclusively called by
     * `Signature`'s methods.
     *
     * @param value_map is a mutation type-value map
     */
    SignatureExprResult(const std::map<MUTATION_TYPE, double>& value_map):
            value_map(value_map)
    {}
public:
    /**
     * @brief The empty constructor
     */
    SignatureExprResult():
            value_map()
    {}

    /**
     * @brief Cast to `Signature<MUTATION_TYPE>`
     *
     * This method tries to cast a signature expression to a
     * signature. When the expression does not represent a
     * probability distribution a `std::domain_error` is thrown.
     *
     * @return the corresponding `Signature<MUTATION_TYPE>` object
     */
    inline operator Signature<MUTATION_TYPE>()
    {
        return Signature<MUTATION_TYPE>(value_map);
    }

    /**
     * @brief Inplace multiply by an arithmetic value
     *
     * @tparam T is the type of the multiplicand
     * @param value is the multiplicand
     * @return a reference to the updated object
     */
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    SignatureExprResult<MUTATION_TYPE>& operator*(const T& value)
    {
        if (value>1 || value<0) {
            std::ostringstream oss;

            oss << "the multiplying value must be in the real interval [0,1]. "
                << std::to_string(value) << " passed.";

            throw std::domain_error(oss.str());
        }

        for (auto& [type, type_value]: value_map) {
            type_value *= value;
        }

        return *this;
    }

    /**
     * @brief Inplace add a signature expression value
     *
     * @param expression_value is a signature expression value
     * @return a reference to the updated object
     */
    SignatureExprResult<MUTATION_TYPE>&
    operator+(SignatureExprResult<MUTATION_TYPE>&& expression_value)
    {
        for (auto& [type, prob]: value_map) {
            auto it = expression_value.value_map.find(type);
            if (it != expression_value.value_map.end()) {
                prob += expression_value.value_map[type];
            }
        }

        for (const auto& [type, prob]: expression_value.value_map) {
            if (value_map.count(type)==0) {
                value_map[type] = prob;
            }
        }

        return *this;
    }

    /**
     * @brief Inplace add a signature
     *
     * @param signature is a signature
     * @return a reference to the updated object
     */
    SignatureExprResult<MUTATION_TYPE>& operator+(const Signature<MUTATION_TYPE>& signature)
    {
        for (auto& [type, prob]: value_map) {
            prob += signature(type);
        }

        for (const auto& [type, prob]: signature) {
            if (value_map.count(type)==0) {
                value_map[type] = prob;
            }
        }

        return *this;
    }

    template<typename MUTATION_TYPE2>
    friend class Signature;
};

/**
 * @brief A class to represent a signature
 *
 * A signature is a probability distribution on
 * the set of mutation types.
 *
 * @tparam MUTATION_TYPE is the mutation type
 */
template<typename MUTATION_TYPE>
class Signature
{
    std::map<MUTATION_TYPE, double> dist_map; //!< the signature probability distribution map

    /**
     * @brief Check whether two value are about the same
     *
     * @param x is a value
     * @param y is a value
     * @return `true` if and only if `x` and `y` differ for an establish value at most
     */
    static inline bool is_about(const double& x, const double& y)
    {
        return std::fabs(x - y) <= 1e-6;
    }

    /**
     * @brief Read a row of a file representing signatures
     *
     * @param in is the input stream
     * @param delimiter is the character delimiting the columns
     * @return the vector of the field in the read row
     */
    static std::vector<std::string> read_row(std::istream& in, const char& delimiter)
    {
        std::string line;
        getline(in, line);

        std::istringstream iss(line);

        std::vector<std::string> row;
        std::string cell;

        while (getline(iss, cell, delimiter)) {
            if (cell.back()=='\r') {
                cell.pop_back();
            }
            row.push_back(cell);
        }

        return row;
    }

    /**
     * @brief Read a map name-(mutation type-probability map) from a stream
     *
     * This method read a map name-(mutation type-probability map) from a
     * stream.
     *
     * @param in is the input stream
     * @param delimiter is the character delimiting the columns
     * @return a map name-(mutation type-probability map)
     */
    static std::map<std::string, std::map<MUTATION_TYPE, double>>
    read_map_from_stream(std::istream& in, const char& delimiter)
    {
        std::map<std::string, std::map<MUTATION_TYPE, double>> result;
        std::vector<std::string> name_vector = read_row(in, delimiter);

        unsigned int row_number = 2;
        std::vector<std::string> row = read_row(in, delimiter);
        while (row.size() != 0) {
            if (row.size() != name_vector.size()) {
                std::ostringstream oss;

                oss << "The header and the row number " << row_number << " differ in size.";
                throw std::runtime_error(oss.str());
            }

            MUTATION_TYPE type(row.front());

            auto row_it = row.begin();
            auto name_it = name_vector.begin();
            for (++row_it, ++name_it;row_it != row.end(); ++row_it, ++name_it) {
                result[*name_it][type] = std::stod(*row_it);
            }

            ++row_number;
            row = read_row(in,delimiter);
        }

        return result;
    }
public:
    using const_iterator = typename std::map<MUTATION_TYPE, double>::const_iterator;

    /**
     * @brief The empty constructor
     */
    Signature():
        dist_map()
    {
        dist_map[MUTATION_TYPE()] = 1;
    }

    /**
     * @brief A constructor
     *
     * @param distribution is a mutation type-value map representing a distribution
     */
    explicit Signature(const std::map<MUTATION_TYPE, double>& distribution):
        dist_map(distribution)
    {
        double partial = 0;
        for (const auto& [type, prob]: dist_map) {
            if (prob<0 || prob>1) {
                std::ostringstream oss;

                oss << type << ": " << prob;
                throw std::domain_error("The parameter is not a probability distribution: "
                                        + oss.str());
            }
            partial += prob;
        }
        if (!is_about(partial,1)) {
            std::ostringstream oss;
            oss << "The parameter is not a probability distribution: 1 minus the sum of"
                << " the probabilities is " << (1-partial);

            throw std::domain_error(oss.str());
        }
    }

    /**
     * @brief Get the initial constant iterator
     *
     * @return the initial constant iterator
     */
    inline const_iterator begin() const
    {
        return dist_map.begin();
    }

    /**
     * @brief Get the final constant iterator
     *
     * @return the final constant iterator
     */
    inline const_iterator end() const
    {
        return dist_map.end();
    }

    /**
     * @brief Get the probability associated to a mutation type
     *
     * @param type is the mutation type whose probability is aimed
     * @return the probability of `type`
     */
    double operator()(const MUTATION_TYPE& type) const
    {
        auto it = dist_map.find(type);
        if (it != dist_map.end()) {
            return it->second;
        }

        return 0;
    }

    /**
     * @brief Multiply by an arithmetic value
     *
     * @tparam T is the type of the multiplicand
     * @param value is the multiplicand
     * @return the resulting signature expression value
     */
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    inline SignatureExprResult<MUTATION_TYPE> operator*(const T& value) const
    {
        return SignatureExprResult<MUTATION_TYPE>(dist_map) * value;
    }

    /**
     * @brief Add a signature
     *
     * @param signature is a signature
     * @return the resulting signature expression value
     */
    inline SignatureExprResult<MUTATION_TYPE> operator+(const Signature& signature) const
    {
        return SignatureExprResult<MUTATION_TYPE>(dist_map) + signature;
    }

    /**
     * @brief Read signature from a input stream
     *
     * @param in is the input stream
     * @return a map that associates the name of the signatures in the file
     *         and the corresponding signature.
     */
    static std::map<std::string, Signature<MUTATION_TYPE>>
    read_from_stream(std::istream& in)
    {
        std::map<std::string, Signature<MUTATION_TYPE>> result;

        std::string curLocale = setlocale(LC_ALL, nullptr);
        setlocale(LC_ALL,"C");
        auto name_signature_map = read_map_from_stream(in,'\t');
        setlocale(LC_ALL, curLocale.c_str());
        for (const auto& [name, signature]: name_signature_map) {
            try {
                result.emplace(std::string(name), Signature<MUTATION_TYPE>(signature));
            } catch (std::domain_error& ex) {
                std::ostringstream oss;

                oss << "Column \"" << name << "\" is not a signature: "
                    << ex.what();
                throw std::runtime_error(oss.str());
            }
        }

        return result;
    }

    /**
     * @brief Read signature from a input stream
     *
     * @param in is the input stream
     * @param signature_names is the set of the requested signature
     * @return a map that associates the name of the signatures in the file that
     *         match `signature_names` and the corresponding signature.
     */
    static std::map<std::string, Signature<MUTATION_TYPE>>
    read_from_stream(std::istream& in, const std::set<std::string>& signature_names)
    {
        auto result = read_from_stream(in);

        auto it = result.begin();

        while (it != result.end()) {
            if (signature_names.count(it->first)==0) {
                result.erase(it++);
            } else {
                ++it;
            }
        }

        return result;
    }
};

/**
 * @brief Multiply an arithmetic value and a signature
 *
 * @tparam MUTATION_TYPE is the mutation type
 * @tparam T is the type of the arithmetic value
 * @param value is the arithmetic value
 * @param signature
 * @return a `SignatureExprResult<MUTATION_TYPE>` object representing the
 *         the multiplication result
 */
template<typename MUTATION_TYPE, typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
inline SignatureExprResult<MUTATION_TYPE> operator*(const T& value, const Signature<MUTATION_TYPE>& signature)
{
    return signature * value;
}

}   // Mutations

}   // CLONES

#endif // __CLONES_SIGNATURE__