/**
 * @file label_tour.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines the label tour class
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

#ifndef __CLONES_LABEL_TOUR__
#define __CLONES_LABEL_TOUR__

#include <vector>

#include <concepts>
#include <type_traits>

#include "cell.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief The concept of node
 */
template <typename NODE>
concept isNode = requires(NODE node) {
    { node.parent() } -> std::same_as<NODE>;
    { node.children() } -> std::same_as<std::vector<NODE>>;
    { node.is_leaf() } -> std::same_as<bool>;
    { node.is_root() } -> std::same_as<bool>;
    { node.get_id() } -> std::convertible_to<Mutants::CellId>;
};

/**
 * @brief The concept of forest
 */
template <typename FOREST>
concept isForest = requires(const FOREST forest) {
    typename FOREST::const_node;

    { forest.get_roots() } -> std::same_as<std::vector<typename FOREST::const_node>>;
} && isNode<typename FOREST::const_node>;

/**
 * @brief The concept of labelling functor
 */
template <typename F, typename NODE_TYPE>
concept isLabellingFunctor = requires(F f, const typename F::label_type& parent_label, const NODE_TYPE& node)
{
    typename F::label_type;

    { f(parent_label, node) } -> std::convertible_to<typename F::label_type>;
};

/**
 * @brief Tours of the forest node labels
 *
 * This class implements tours of the labels of the nodes of a forest. It
 * implements a constant iterator that allows the on-line generation labels
 * of each node in the forest.
 *
 * @tparam FOREST is the type of the forest whose nodes must be labelled
 * @tparam LABELLING_FUNCTOR is the type of the labelling functor
 */
template<typename FOREST, typename LABELLING_FUNCTOR>
    requires isForest<FOREST>
        && isLabellingFunctor<LABELLING_FUNCTOR, typename FOREST::const_node>
class LabelTour
{
    FOREST const* forest;   //!< A pointer to the forest
    bool only_leaves;       //!< A Boolean flag to enable/disable internal node visit

    LABELLING_FUNCTOR l_functor;
    typename LABELLING_FUNCTOR::label_type init_label;

public:
    using forest_type = FOREST;
    using const_node = typename FOREST::const_node;
    using labelling_functor_type = LABELLING_FUNCTOR;
    using label_type = typename LABELLING_FUNCTOR::label_type;

    /**
     * @brief A constant iterator for the label tour
     *
     * This class implements a constant iterator that visits all nodes or
     * leaves in the tour forest and generates the corresponding labels.
     * The asymptotic complexity of the successor operator is linear in the
     * number of forest nodes in the worst case. However, this asymptotic
     * complexity also holds for the full tour.
     * In the worst case, the memory required by each iterator is
     * logarithmic in the number of forest nodes.
     */
    class const_iterator
    {
    public:
        /**
         * @brief The dereferenced value type
         */
        using value_type = std::pair<Mutants::CellId, label_type>;
    private:
        forest_type const* forest;   //!< A pointer to the forest

        labelling_functor_type l_functor;

        bool only_leaves;   //!< A Boolean flag to enable/disable internal node visit
        bool tour_end;      //!< A Boolean flag to mark the end of the tour

        std::stack<value_type> iterator_stack; //!< The recursion stack
        value_type node_label; //!< The current node label

        /**
         * @brief A constructor
         *
         * @param[in] forest is a constant pointer to the forest
         * @param[in] labelling_functor is the labelling functor
         * @param[in] init_label is the initialization label
         * @param[in] only_leaves is a Boolean flag to enable/disable
         *      internal node visit
         */
        const_iterator(const forest_type* forest,
                       const labelling_functor_type& labelling_functor,
                       const label_type& init_label,
                       const bool& only_leaves):
            forest{forest}, l_functor{labelling_functor},
            only_leaves{only_leaves}, tour_end{false}
        {
            if (forest != nullptr) {
                auto forest_roots = forest->get_roots();
                for (auto root_it = forest_roots.rbegin();
                        root_it != forest_roots.rend(); ++root_it) {
                    label_type label = labelling_functor(init_label, *root_it);

                    iterator_stack.emplace(root_it->get_id(), std::move(label));
                }

                std::swap(node_label, iterator_stack.top());

                iterator_stack.pop();

                if (only_leaves) {
                    this->operator++();
                }
            } else {
                tour_end = true;
            }
        }

        const_iterator(const forest_type* forest, const bool& only_leaves):
            forest{forest}, only_leaves{only_leaves}, tour_end{true}
        {}
    public:
        /**
         * @brief The empty constructor
         */
        const_iterator():
            const_iterator{nullptr, false}
        {}

        /**
         * @brief The successor operator
         *
         * This method moves the iterator to the next node of the tour
         * and, at the same time, computes the label of the reached node.
         * The asymptotic complexity of this method is linear in the
         * number of the forest nodes. However, the full tour has the
         * same complexity.
         *
         * @return a reference to the updated object
         */
        const_iterator& operator++()
        {
            if (forest == nullptr) {
                return *this;
            }

            const_node node{forest, node_label.first};

            if (node.is_leaf()) {
                if (iterator_stack.empty()) {
                    tour_end = true;

                    return *this;
                }

                // take a new cell genome mutations object from the stack
                std::swap(node_label, iterator_stack.top());
                iterator_stack.pop();

                if (!only_leaves) {
                    return *this;
                }

                // update the node
                node = const_node{forest, node_label.first};
            }

            bool next_node_found = node.is_leaf();
            while (!next_node_found) {
                auto children = node.children();

                // place all children's cell genome mutations, but the first one, into the stack
                for (auto child_it = children.rbegin();
                        child_it != children.rend()-1; ++child_it) {

                    label_type child_label = l_functor(node_label.second, *child_it);

                    iterator_stack.emplace(child_it->get_id(), std::move(child_label));
                }

                // apply the first children mutations to the current cell genome mutations
                std::swap(node, children.front());
                node_label.second = l_functor(node_label.second, node);

                next_node_found = node.is_leaf() || !only_leaves;
            }

            // update the node cell id
            node_label.first = node.get_id();

            return *this;
        }

        /**
         * @brief The dereference operator
         *
         * This method returns a constant reference to the pair
         * identifier/label of the current node in the tour. It takes
         * constant time.
         *
         * @return a constant reference to the pair identifier/label
         *      of the current node in the tour
         */
        inline const value_type& operator*() const
        {
            return node_label;
        }

        /**
         * @brief Check whether the tour end has been reached
         *
         * @return `true` if and only if the tour end has been
         *      reached
         */
        inline const bool& is_end() const
        {
            return tour_end;
        }

        /**
         * @brief Equality operator
         *
         * This method checks whether two `const_iterator` are
         * the same.
         *
         * @param[in] rhs is the right-hand side of the equality
         * @return `true` if and only if this object and `rhs`
         *      iterates over the same tour and refer to the
         *      same node
         */
        bool operator==(const const_iterator& rhs) const
        {
            if (this->forest != rhs.forest) {
                return false;
            }

            if (this->iterator_stack.size() != rhs.iterator_stack.size()) {
                return false;
            }

            if (this->iterator_stack.size()==0) {
                if (this->tour_end != rhs.tour_end) {
                    return false;
                }

                if (this->tour_end) {
                    return true;
                }
            }

            return this->iterator_stack.top() == rhs.iterator_stack.top();
        }

        /**
         * @brief inequality operator
         *
         * This method checks whether two `const_iterator` differ.
         *
         * @param[in] rhs is the right-hand side of the inequality
         * @return `true` if and only if this object and `rhs`
         *      iterates over different tours or refer to different
         *      nodes
         */
        inline bool operator!=(const const_iterator& rhs) const
        {
            return !(*this == rhs);
        }

        friend class LabelTour<FOREST, LABELLING_FUNCTOR>;
    };

    /**
     * @brief The empty constructor
     */
    LabelTour():
        forest{nullptr}, only_leaves{false}
    {}

    /**
     * @brief A constructor
     *
     * @param[in] forest is a constant reference to the forest
     * @param[in] labelling_functor is the labelling functor
     * @param[in] init_label is the initialization label
     * @param[in] only_leaves is a Boolean flag to enable/disable internal node visit
     */
    LabelTour(const forest_type& forest,
              const labelling_functor_type& labelling_functor,
              const label_type& init_label,
              const bool only_leaves):
        forest{&forest}, only_leaves{only_leaves},
        l_functor{labelling_functor}, init_label{init_label}
    {}

    /**
     * @brief Get the constant iterator to the tour begin
     *
     * @return the constant iterator to the tour begin
     */
    const_iterator begin() const
    {
        if (forest == nullptr) {
            return const_iterator{forest, only_leaves};
        }

        if (forest->num_of_nodes()==0) {
            return const_iterator{forest, only_leaves};
        }

        return const_iterator{forest, l_functor,
                              init_label, only_leaves};
    }

    /**
     * @brief Get the constant iterator to the tour end
     *
     * @return the constant iterator to the tour end
     */
    inline const_iterator end() const
    {
        return const_iterator{forest, only_leaves};
    }

    /**
     * @brief Get the associated forest
     *
     * @return a constant reference to the tour forest
     */
    inline const forest_type& get_forest() const
    {
        return *forest;
    }

    friend forest_type;
};

}

}

#endif // __CLONES_LABEL_TOUR__