/**
 * @file logics.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a logic about the simulation
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

#ifndef __CLONES_LOGIC__
#define __CLONES_LOGIC__

#include <string>
#include <memory>
#include <ostream>

#include "cell_event.hpp"

struct Context;

namespace CLONES
{

namespace Mutants
{

namespace Evolutions
{
class Tissue;
class Simulation;
}

namespace Logics
{

/**
 * @brief The variable that represent a simulation quantity
 */
class Variable
{
public:
    /**
     * @brief Variable type
     */
    enum class Type {
        CARDINALITY,    //!< The variable represents a species cardinality
        EVENT,          //!< The variable represents a species event
        TIME            //!< The variable represents the simulation time
    };

    /**
     * @brief The empty constructor
     */
    Variable();

    /**
     * @brief Evaluate a variable in a context
     *
     * @tparam CONTEXT is the type of the context in which the variable
     *   must be evaluated
     * @param context is the context in which the variable must be
     *   evaluated
     * @return The variable evalutation, i.e., the cardinality of the
     *   species in the context
     */
    template<typename CONTEXT>
    inline double evaluate(const CONTEXT& context) const
    {
        return context.evaluate(*this);
    }

    /**
     * @brief Get the identified of the species
     *
     * @return the identified of the species
     */
    inline const SpeciesId& get_species_id() const
    {
        return species_id;
    }

    /**
     * @brief Get the name of the variable species
     *
     * @return the name of the variable species
     */
    inline const std::string& get_name() const
    {
        return species_name;
    }

    /**
     * @brief Get the variable type
     *
     * @return the variable type
     */
    inline const Type& get_type() const
    {
        return type;
    }

    friend class CLONES::Mutants::Evolutions::Tissue;
    friend class CLONES::Mutants::Evolutions::Simulation;
    friend std::ostream& operator<<(std::ostream& os, const Variable& variable);
private:

    Type type;                  //!< The variable type
    CellEventType event_type;   //!< The event type
    SpeciesId species_id;       //!< The identifier of the variable species
    std::string species_name;           //!< The name of the variable species

    /**
     * @brief Build a variable representing a species cardinality
     *
     * @param species_id is the identifier of the species
     * @param species_name is the name of the species
     */
    Variable(const SpeciesId& species_id, const std::string& species_name);

    /**
     * @brief Build a variable representing a species event number
     *
     * @param species_id is the identifier of the species
     * @param species_name is the name of the species
     */
    Variable(const CellEventType& event_type, const SpeciesId& species_id,
             const std::string& species_name);
};

/**
 * @brief Write in a stream a species variable representation
 *
 * @param os is the output stream
 * @param variable is the species variable to write
 * @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Variable& variable);

/**
 * @brief A class to represent expressions
 */
class Expression
{
public:
    /**
     * @brief The types of expressions
     */
    enum class Type
    {
        VARIABLE,   // variable
        VALUE,      // value
        SUM,        // sum of expressions
        SUB,        // subtraction of expressions
        MULT        // multiplication of expressions
    };

    /**
     * @brief A constructor for expressions consisting in a variable
     *
     * @param variable is the expression variable
     */
    Expression(const Variable& variable);

    /**
     * @brief A constructor for expressions consisting in a value
     *
     * @param value is the expression value
     */
    Expression(const double& value);

    /**
     * @brief Evaluate the expression in a context
     *
     * @tparam CONTEXT is the context type
     * @param context is the context in which the expression must be
     *   evaluated
     * @return the value of the expression in `context`
     */
    template<typename CONTEXT>
    double evaluate(const CONTEXT& context) const
    {
        switch(type) {
            case Type::VALUE:
                return value;
            case Type::VARIABLE:
                return variable.evaluate(context);
            case Type::SUM:
                return lhs->evaluate(context)+rhs->evaluate(context);
            case Type::SUB:
                return lhs->evaluate(context)-rhs->evaluate(context);
            case Type::MULT:
                return (lhs->evaluate(context))*(rhs->evaluate(context));
            default:
                throw std::domain_error("Unsupported expression type code "
                                        + std::to_string(static_cast<unsigned int>(type)));
        }
    }

    /**
     * @brief Get the expression type
     *
     * @return the expression type
     */
    inline const Type& get_type() const
    {
        return type;
    }

    friend Expression operator+(const Expression& lhs, const Expression& rhs);
    friend Expression operator-(const Expression& lhs, const Expression& rhs);
    friend Expression operator*(const Expression& lhs, const Expression& rhs);

    friend Expression operator+(Expression&& lhs, Expression&& rhs);
    friend Expression operator-(Expression&& lhs, Expression&& rhs);
    friend Expression operator*(Expression&& lhs, Expression&& rhs);

    friend std::ostream& operator<<(std::ostream& os, const Expression& expression);
private:
    Type type;          //!< the type of expressions
    Variable variable;  //!< the variable when the expression is a variable
    double value;       //!< the value when the expression is a value
    std::shared_ptr<Expression> lhs;    //!< the left hand-side expression
    std::shared_ptr<Expression> rhs;    //!< the right hand-side expression

    /**
     * @brief Build an expression involving two sub-expressions
     *
     * @param type is the type of the new expression
     * @param lhs is the left hand-side of the new expression
     * @param rhs is the left hand-side of the new expression
     */
    Expression(const Type& type, const Expression& lhs, const Expression& rhs);

    /**
     * @brief Build an expression involving two sub-expressions
     *
     * @param type is the type of the new expression
     * @param lhs is the left hand-side of the new expression
     * @param rhs is the left hand-side of the new expression
     */
    Expression(const Type& type, Expression&& lhs, Expression&& rhs);
};

/**
 * @brief Add two expressions
 *
 * @param lhs is the left hand-side of the sum
 * @param rhs is the left hand-side of the sum
 * @return the expression that is the sum of the two parameters
 */
inline Expression operator+(const Expression& lhs, const Expression& rhs)
{
    return Expression(Expression::Type::SUM, lhs, rhs);
}

/**
 * @brief Add two expressions
 *
 * @param lhs is the left hand-side of the sum
 * @param rhs is the left hand-side of the sum
 * @return the expression that is the sum of the two parameters
 */
inline Expression operator+(Expression&& lhs, Expression&& rhs)
{
    return Expression(Expression::Type::SUM, std::move(lhs), std::move(rhs));
}

/**
 * @brief Subtract two expressions
 *
 * @param lhs is the left hand-side of the difference
 * @param rhs is the left hand-side of the difference
 * @return the expression that is the difference of the two parameters
 */
inline Expression operator-(const Expression& lhs, const Expression& rhs)
{
    return Expression(Expression::Type::SUB, lhs, rhs);
}

/**
 * @brief Subtract two expressions
 *
 * @param lhs is the left hand-side of the difference
 * @param rhs is the left hand-side of the difference
 * @return the expression that is the difference of the two parameters
 */
inline Expression operator-(Expression&& lhs, Expression&& rhs)
{
    return Expression(Expression::Type::SUB, std::move(lhs), std::move(rhs));
}

/**
 * @brief Multiply two expressions
 *
 * @param lhs is the left hand-side of the product
 * @param rhs is the left hand-side of the product
 * @return the expression that is the product of the two parameters
 */
inline Expression operator*(const Expression& lhs, const Expression& rhs)
{
    return Expression(Expression::Type::MULT, lhs, rhs);
}

/**
 * @brief Multiply two expressions
 *
 * @param lhs is the left hand-side of the product
 * @param rhs is the left hand-side of the product
 * @return the expression that is the product of the two parameters
 */
inline Expression operator*(Expression&& lhs, Expression&& rhs)
{
    return Expression(Expression::Type::MULT, std::move(lhs), std::move(rhs));
}

/**
 * @brief Write in a stream an expression representation
 *
 * @param os is the output stream
 * @param expression is the expression to write
 * @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Expression& expression);

/**
 * @brief A class to represent expression relations
 */
class Relation
{
public:
    /**
     * @brief The types of relations
     */
    enum class Type
    {
        GT, // >
        GE, // >=
        EQ, // ==
        NE  // !=
    };

    /**
     * @brief Evaluate the relation in a context
     *
     * @tparam CONTEXT is the context type
     * @param context is the context in which the relation must be
     *   evaluated
     * @return the value of the relation in `context`
     */
    template<typename CONTEXT>
    bool evaluate(const CONTEXT& context) const
    {
        auto l_value = lhs.evaluate(context);
        auto r_value = rhs.evaluate(context);

        switch(type) {
            case Type::GT:
                return l_value > r_value;
            case Type::GE:
                return l_value >= r_value;
            case Type::EQ:
                return l_value == r_value;
            case Type::NE:
                return l_value != r_value;
            default:
                throw std::runtime_error("Unknown relation type code "
                                         + std::to_string(static_cast<unsigned int>(type)));
        }
    }

    /**
     * @brief Measure the distance of the relation from being satisfied
     *
     * @tparam CONTEXT is the context type
     * @param context is the context in which the relation must be
     *   evaluated
     * @return a value that measure the distance of the relation from
     *   being satisfied in the context
     */
    template<typename CONTEXT>
    double sat_distance(const CONTEXT& context) const
    {
        auto l_dist = lhs.evaluate(context);
        auto r_dist = rhs.evaluate(context);

        switch(type) {
            case Type::GT:
                return (l_dist > r_dist? 0:r_dist-l_dist);
            case Type::GE:
                return (l_dist >= r_dist? 0:r_dist-l_dist+1);
            case Type::EQ:
                return (l_dist > r_dist? l_dist-r_dist:r_dist-l_dist);
            case Type::NE:
                return (l_dist != r_dist? 0: 1);
            default:
                throw std::runtime_error("Unknown relation type code "
                                         + std::to_string(static_cast<unsigned int>(type)));
        }
    }

    /**
     * @brief Get the relation type
     *
     * @return the relation type
     */
    inline const Type& get_type() const
    {
        return type;
    }

    friend Relation operator>(const Expression& lhs, const Expression& rhs);
    friend Relation operator>=(const Expression& lhs, const Expression& rhs);
    friend Relation operator==(const Expression& lhs, const Expression& rhs);
    friend Relation operator!=(const Expression& lhs, const Expression& rhs);
    friend Relation operator<=(const Expression& lhs, const Expression& rhs);
    friend Relation operator<(const Expression& lhs, const Expression& rhs);

    friend Relation operator>(Expression&& lhs, Expression&& rhs);
    friend Relation operator>=(Expression&& lhs, Expression&& rhs);
    friend Relation operator==(Expression&& lhs, Expression&& rhs);
    friend Relation operator!=(Expression&& lhs, Expression&& rhs);
    friend Relation operator<=(Expression&& lhs, Expression&& rhs);
    friend Relation operator<(Expression&& lhs, Expression&& rhs);

    friend std::ostream& operator<<(std::ostream& os, const Relation& relation);
private:
    Type type;          //!< The relation type
    Expression lhs;     //!< The relation left hand-side operator
    Expression rhs;     //!< The relation right hand-side operator

    /**
     * @brief A constructor
     *
     * @param type is the type of the relation
     * @param lhs is the left hand-side of the relation
     * @param rhs is the left hand-side of the relation
     */
    Relation(const Type& type, Expression&& lhs, Expression&& rhs);

    /**
     * @brief A constructor
     *
     * @param type is the type of the relation
     * @param lhs is the left hand-side of the relation
     * @param rhs is the left hand-side of the relation
     */
    Relation(const Type& type, const Expression& lhs, const Expression& rhs);
};

/**
 * @brief Write in a stream a relation representation
 *
 * @param os is the output stream
 * @param relation is the relation to write
 * @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Relation& relation);

/**
 * @brief Build the relation ">"
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `lhs > rhs`
 */
inline Relation operator>(const Expression& lhs, const Expression& rhs)
{
    return Relation(Relation::Type::GT, lhs, rhs);
}

/**
 * @brief Build the relation ">"
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `lhs > rhs`
 */
inline Relation operator>(Expression&& lhs, Expression&& rhs)
{
    return Relation(Relation::Type::GT, std::move(lhs), std::move(rhs));
}

/**
 * @brief Build the relation ">="
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `lhs >= rhs`
 */
inline Relation operator>=(const Expression& lhs, const Expression& rhs)
{
    return Relation(Relation::Type::GE, lhs, rhs);
}

/**
 * @brief Build the relation ">="
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `lhs >= rhs`
 */
inline Relation operator>=(Expression&& lhs, Expression&& rhs)
{
    return Relation(Relation::Type::GE, std::move(lhs), std::move(rhs));
}

/**
 * @brief Build the relation "=="
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `lhs == rhs`
 */
inline Relation operator==(const Expression& lhs, const Expression& rhs)
{
    return Relation(Relation::Type::EQ, lhs, rhs);
}

/**
 * @brief Build the relation "=="
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `lhs == rhs`
 */
inline Relation operator==(Expression&& lhs, Expression&& rhs)
{
    return Relation(Relation::Type::EQ, std::move(lhs), std::move(rhs));
}

/**
 * @brief Build the relation "!="
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `lhs != rhs`
 */
inline Relation operator!=(const Expression& lhs, const Expression& rhs)
{
    return Relation(Relation::Type::NE, lhs, rhs);
}

/**
 * @brief Build the relation "!="
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `lhs != rhs`
 */
inline Relation operator!=(Expression&& lhs, Expression&& rhs)
{
    return Relation(Relation::Type::NE, std::move(lhs), std::move(rhs));
}

/**
 * @brief Build the relation "<="
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `rhs >= lhs`
 */
inline Relation operator<=(const Expression& lhs, const Expression& rhs)
{
    return Relation(Relation::Type::GE, rhs, lhs);
}

/**
 * @brief Build the relation "<="
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `rhs >= lhs`
 */
inline Relation operator<=(Expression&& lhs, Expression&& rhs)
{
    return Relation(Relation::Type::GE, std::move(rhs), std::move(lhs));
}

/**
 * @brief Build the relation "<"
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `rhs < lhs`
 */
inline Relation operator<(const Expression& lhs, const Expression& rhs)
{
    return Relation(Relation::Type::GT, rhs, lhs);
}

/**
 * @brief Build the relation "<"
 *
 * @param lhs is the left hand-side of the relation
 * @param rhs is the left hand-side of the relation
 * @return the relation `rhs < lhs`
 */
inline Relation operator<(Expression&& lhs, Expression&& rhs)
{
    return Relation(Relation::Type::GT, std::move(rhs), std::move(lhs));
}

/**
 * @brief A class to represent formulas
 */
struct Formula
{
    /**
     * @brief The types of formulas
     */
    enum class Type
    {
        RELATION,   // relation
        NEG,        // negation of a formula
        AND,        // conjunction of two formulas
        OR          // disjunction of two formulas
    };

    Type type;  //!< The type of the formulas
    std::shared_ptr<Relation> rel;  //!< The relation is this formula is a relation
    std::shared_ptr<Formula> lhs;   //!< The left hand-side subformula
    std::shared_ptr<Formula> rhs;   //!< The right hand-side subformula

    /**
     * @brief The empty constructor
     */
    Formula();

    /**
     * @brief Build a formula that consists in a relation
     *
     * @param relation is the relation
     */
    Formula(const Relation& relation);

    /**
     * @brief Build a formula that consists in a relation
     *
     * @param relation is the relation
     */
    Formula(Relation&& relation);

    /**
     * @brief Evaluate the formula in a context
     *
     * @tparam CONTEXT is the context type
     * @param context is the context in which the formula must be
     *   evaluated
     * @return the value of the formula in `context`
     */
    template<typename CONTEXT>
    bool evaluate(const CONTEXT& context) const
    {
        switch(type) {
            case Type::RELATION:
                return rel->evaluate(context);
            case Type::NEG:
                return !(lhs->evaluate(context));
            case Type::AND:
                return (lhs->evaluate(context) && rhs->evaluate(context));
            case Type::OR:
                return (lhs->evaluate(context) || rhs->evaluate(context));
            default:
                throw std::runtime_error("Unknown formula type code "
                                         + std::to_string(static_cast<unsigned int>(type)));
        }
    }

    /**
     * @brief Measure the distance of the formula from being satisfied
     *
     * @tparam CONTEXT is the context type
     * @param context is the context in which the formula must be
     *   evaluated
     * @return a value that measure the distance of the formula from
     *   being satisfied in the context
     */
    template<typename CONTEXT>
    double sat_distance(const CONTEXT& context) const
    {
        switch(type) {
            case Type::RELATION:
                return rel->sat_distance(context);
            case Type::NEG:
                return (lhs->sat_distance(context)?1:0);
            case Type::AND:
                return lhs->sat_distance(context) + rhs->sat_distance(context);
            case Type::OR:
                return std::min(lhs->sat_distance(context), rhs->sat_distance(context));
            default:
                throw std::runtime_error("Unknown formula type code "
                                         + std::to_string(static_cast<unsigned int>(type)));
        }
    }

    /**
     * @brief Get the formula type
     *
     * @return the formula type
     */
    inline const Type& get_type() const
    {
        return type;
    }

    friend Formula operator!(const Formula& subformula);
    friend Formula operator&&(const Formula& lhs, const Formula& rhs);
    friend Formula operator||(const Formula& lhs, const Formula& rhs);
    friend Formula operator!(Formula&& subformula);
    friend Formula operator&&(Formula&& lhs, Formula&& rhs);
    friend Formula operator||(Formula&& lhs, Formula&& rhs);

    friend std::ostream& operator<<(std::ostream& os, const Formula& formula);
private:

    /**
     * @brief A constructor for a formula having two sub-fomulas
     *
     * @param type is the type of the new formula
     * @param lhs is the left hand-side of the new formula
     * @param rhs is the left hand-side of the new formula
     */
    Formula(const Type& type, const Formula& lhs, const Formula& rhs);

    /**
     * @brief A constructor for a formula having two sub-fomulas
     *
     * @param type is the type of the new formula
     * @param lhs is the left hand-side of the new formula
     * @param rhs is the left hand-side of the new formula
     */
    Formula(const Type& type, Formula&& lhs, Formula&& rhs);

    /**
     * @brief A constructor for a formula having one sub-fomula
     *
     * @param type is the type of the new formula
     * @param subformula is the subformula
     */
    Formula(const Type& type, const Formula& subformula);

    /**
     * @brief A constructor for a formula having one sub-fomula
     *
     * @param type is the type of the new formula
     * @param subformula is the subformula
     */
    Formula(const Type& type, Formula&& subformula);
};

/**
 * @brief Write in a stream a formula representation
 *
 * @param os is the output stream
 * @param formula is the formula to write
 * @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Formula& formula);

/**
 * @brief Negate a formula
 *
 * @param subformula is the formula to be negated
 * @return the formula that negate `subformula`
 */
inline Formula operator!(const Formula& subformula)
{
    return Formula(Formula::Type::NEG, subformula);
}

/**
 * @brief Negate a formula
 *
 * @param subformula is the formula to be negated
 * @return the formula that negate `subformula`
 */
inline Formula operator!(Formula&& subformula)
{
    return Formula(Formula::Type::NEG, std::move(subformula));
}

/**
 * @brief Get the conjunction of two formulas
 *
 * @param lhs is the left hand-side of the conjunction
 * @param rhs is the right hand-side of the conjunction
 * @return the conjunction of `lhs` and `rhs`
 */
inline Formula operator&&(const Formula& lhs, const Formula& rhs)
{
    return Formula(Formula::Type::AND, lhs, rhs);
}

/**
 * @brief Get the conjunction of two formulas
 *
 * @param lhs is the left hand-side of the conjunction
 * @param rhs is the right hand-side of the conjunction
 * @return the conjunction of `lhs` and `rhs`
 */
inline Formula operator&&(Formula&& lhs, Formula&& rhs)
{
    return Formula(Formula::Type::AND, std::move(lhs), std::move(rhs));
}

/**
 * @brief Get the disjunction of two formulas
 *
 * @param lhs is the left hand-side of the disjunction
 * @param rhs is the right hand-side of the disjunction
 * @return the disjunction of `lhs` and `rhs`
 */
inline Formula operator||(const Formula& lhs, const Formula& rhs)
{
    return Formula(Formula::Type::OR, lhs, rhs);
}

/**
 * @brief Get the disjunction of two formulas
 *
 * @param lhs is the left hand-side of the disjunction
 * @param rhs is the right hand-side of the disjunction
 * @return the disjunction of `lhs` and `rhs`
 */
inline Formula operator||(Formula&& lhs, Formula&& rhs)
{
    return Formula(Formula::Type::OR, std::move(lhs), std::move(rhs));
}

}   // Logics

}   // Mutants

}   // CLONES

#endif // __CLONES_LOGIC__