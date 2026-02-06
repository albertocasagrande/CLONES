/**
 * @file irb_tree.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines indexed red-black trees
 * @version 1.2
 * @date 2026-02-06
 *
 * @copyright Copyright (c) 2024
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

#ifndef __CLONES_IRB_TREE__
#define __CLONES_IRB_TREE__

#include <memory>
#include <utility>
#include <functional>

namespace CLONES
{

/**
 * @brief A Red-Black Tree implementation that supports indexing
 *
 * This implementation of the Red-Black trees (e.g., see Chapter 13 [1]) supports
 * indexing by adding the field `subtree_size` to each node in the tree. This
 * field accounts the number of nodes in the subtree rooted in the specific node.
 * Finding the i-th key in the tree can be achieved by using the `subtree_size`
 * fields in time \f$O(\log \texttt{size()})\f$.
 * 
 * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, Clifford Stein:
 *     Introduction to Algorithms, 3rd Edition. MIT Press 2009,
 *     ISBN 978-0-262-03384-8, pp. I-XIX, 1-1292
 *
 * @tparam KEY is the type of the RB tree keys
 * @tparam COMPARE is the order among keys
 * @tparam ALLOCATOR is the key allocator
 */
template<class KEY, class COMPARE=std::less<KEY>, class ALLOCATOR=std::allocator<KEY>>
class IRBTree
{
    /**
     * @brief The tree node
     */
    struct IRBNode
    {
        /**
         * @brief The node color type
         */
        enum Color {
            BLACK,
            RED
        };

        /**
         * @brief The type of the children side
         */
        enum Side {
            LEFT = 0,
            RIGHT = 1
        };

        KEY* key;       //!< The node's key
        Color color;    //!< The node's color
        size_t subtree_size;    //!< The size of the subtree rooted in the node

        IRBNode* parent;        //!< The node's parent
        IRBNode* child[2];   //!< The node's children

        /**
         * @brief Construct a new indexed reb-black tree
         *
         * @param key is the node key
         */
        explicit IRBNode(const KEY& key):
            color{RED}, subtree_size{1},
            parent{nullptr}, child{nullptr, nullptr}
        {
            ALLOCATOR allocator;

            this->key = allocator.allocate(1);

            // build key
            using traits = std::allocator_traits<ALLOCATOR>;

            traits::construct(allocator, this->key, key);
        }

        /**
         * @brief Update the size of the subtree rooted on the node
         *
         * This method updates the size of the subtree rooted on a
         * node and all its ancestors. The complexity of this method
         * is \f$O(h)\f$ where \f$h\f$ is the height of the tree.
         *
         * @param delta is the difference between the new size
         *      and the previous size of the subtree rooted on
         *      the current node
         */
        void delta_size(const int delta)
        {
            IRBNode* node = this;
            while (node != nullptr) {
                node->subtree_size = static_cast<size_t>(node->subtree_size + delta);

                node = node->parent;
            }
        }

        /**
         * @brief Get a string representation of a side
         *
         * @param side is a side
         * @return a string representation of `side`
         */
        inline static std::string side_str(const Side& side)
        {
            if (side == LEFT) return "left";

            return "right";
        }

        /**
         * @brief Set a node as one of the children of a node
         *
         * The complexity of the method is \f$O(1)\f$.
         *
         * @param new_child is a pointer to the new child
         * @param side is the side of the new child
         * @return a pointer to the new child
         */
        IRBNode* set_child(IRBNode* new_child, const Side side)
        {
            if (child[side] != nullptr) {
                throw std::runtime_error("The node already have a "
                                          + side_str(side) + " child.");
            }
            if (new_child == nullptr) {
                throw std::runtime_error("The new child must differ from nullptr.");
            }
            child[side] = new_child;
            new_child->parent = this;

            delta_size(new_child->subtree_size);

            return new_child;
        }

        /**
         * @brief Test whether this node is a root
         *
         * @return `true` if and only if this node does not
         *      have a parent
         */
        inline bool is_root() const
        {
            return parent == nullptr;
        }

        /**
         * @brief Test whether this node is a child on a specific side
         *
         * @param side is a side
         * @return `true` if and only if this is node is the `side`
         *      child of a node
         */
        bool is_child(const Side& side) const
        {
            if (parent == nullptr) {
                return false;
            }

            return parent->child[side] == this;
        }

        /**
         * @brief Get the node side
         *
         * @param node is the node whose side is aimed
         * @return `IRBNode::LEFT` if the node is a left child.
         *      `IRBNode::RIGHT` if the node is a right child.
         */
        inline static Side get_side(const IRBNode* node)
        {
            if (node->is_child(LEFT)) return LEFT;

            return RIGHT;
        }

        /**
         * @brief Get the node side
         *
         * @return `IRBNode::LEFT` if the node is a left child.
         *      `IRBNode::RIGHT` if the node is a right child.
         */
        inline Side get_side() const
        {
            return get_side(this);
        }

        /**
         * @brief Get the opposite side
         *
         * @param side is a side
         * @return `IRBNode::LEFT`, if `side` is `IRBNode::RIGHT`, or
         *      `IRBNode::RIGHT`, if `side` is `IRBNode::LEFT`
         */
        static inline Side get_opposite_side(const Side side)
        {
            if (side == LEFT) return RIGHT;

            return LEFT;
        }

        /**
         * @brief Get the node opposite side
         *
         * @param node is the node whose opposite side is aimed
         * @return `IRBNode::LEFT` if the node is a right child.
         *      `IRBNode::RIGHT` if the node is a left child.
         */
        static inline Side get_opposite_side(const IRBNode* node)
        {
            return get_opposite_side(get_side(node));
        }

        /**
         * @brief Get the opposite side
         *
         * @return `IRBNode::LEFT` if the node is a right child.
         *      `IRBNode::RIGHT` if the node is a left child.
         */
        inline Side get_opposite_side() const
        {
            return get_opposite_side(this);
        }

        /**
         * @brief Test whether a node is black
         *
         * @param node is a pointer to a node
         * @return `true` if and only if `node` is black or
         *       `node` is `nullptr`
         */
        static inline bool is_black(const IRBNode* node)
        {
            if (node == nullptr) {
                return true;
            }

            return node->color == BLACK;
        }

        /**
         * @brief Get the left/right-most node in the subtree
         *
         * The complexity of this method is \f$O(h)\f$ where \f$h\f$
         * is the height of the tree.
         *
         * @param side is a side
         * @return a constant pointer to a `side`-most node in the
         *      subtree rooted on the current node
         */
        const IRBNode* get_all_on(const Side& side) const
        {
            IRBNode const* node = this;
            while (node->child[side] != nullptr) {
                node = node->child[side];
            }

            return node;
        }

        /**
         * @brief Get the maximum-key node on the subtree
         *
         * The complexity of this method is \f$O(h)\f$ where \f$h\f$
         * is the height of the tree.
         *
         * @return a constant pointer to the node containing
         *      the maximum key among those in the subtree
         *      rooted on the current node
         */
        inline const IRBNode* get_subtree_max() const
        {
            return get_all_on(RIGHT);
        }

        /**
         * @brief Get the minimum-key node on the subtree
         *
         * The complexity of this method is \f$O(h)\f$ where \f$h\f$
         * is the height of the tree.
         *
         * @return a constant pointer to the node containing
         *      the minimum key among those in the subtree
         *      rooted on the current node
         */
        inline const IRBNode* get_subtree_min() const
        {
            return get_all_on(LEFT);
        }

        /**
         * @brief Get the left/right-most node in the subtree
         *
         * The complexity of this method is \f$O(h)\f$ where \f$h\f$
         * is the height of the tree.
         *
         * @param side is a side
         * @return a pointer to a `side`-most node in the
         *      subtree rooted on the current node
         */
        IRBNode* get_all_on(const Side& side)
        {
            IRBNode* node = this;
            while (node->child[side] != nullptr) {
                node = node->child[side];
            }

            return node;
        }

        /**
         * @brief Get the maximum-key node on the subtree
         *
         * The complexity of this method is \f$O(h)\f$ where \f$h\f$
         * is the height of the tree.
         *
         * @return a pointer to the node containing the maximum
         *      key among those in the subtree rooted on the
         *      current node
         */
        inline IRBNode* get_subtree_max()
        {
            return get_all_on(RIGHT);
        }

        /**
         * @brief Get the minumum-key node on the subtree
         *
         * The complexity of this method is \f$O(h)\f$ where \f$h\f$
         * is the height of the tree.
         *
         * @return a pointer to the node containing the minumum
         *      key among those in the subtree rooted on the
         *      current node
         */
        inline IRBNode* get_subtree_min()
        {
            return get_all_on(LEFT);
        }

        /**
         * @brief Get the next node in the visit
         *
         * This method generalizes the algorithm `Tree-Successor` Chapter 12,
         * pp 292, [1]. Its complexity is \f$O(h)\f$ where \f$h\f$ is the height
         * of the tree.
         *
         * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
         *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
         *     (3rd. ed.). The MIT Press.
         *
         * @param side is a side
         * @return a constant pointer to the next node in the in-order or
         *      reversed in-order visits when `side` is `IRBNode::LEFT` or
         *      `IRBNode::RIGHT`, respectively
         */
        const IRBNode* torwards(const Side& side) const
        {
            if (child[side] != nullptr) {
                const Side opposite_side = get_opposite_side(side);
                return child[side]->get_all_on(opposite_side);
            }

            IRBNode const* node = this;
            while (node->is_child(side)) {
                node = node->parent;
            }

            return node->parent;
        }

        /**
         * @brief Get the successor
         *
         * This method implements the algorithm `Tree-Successor` Chapter 12,
         * pp 292, [1]. Its complexity is \f$O(h)\f$ where \f$h\f$ is the height
         * of the tree.
         *
         * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
         *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
         *     (3rd. ed.). The MIT Press.
         *
         * @return a constant pointer to the next node in the in-order visit
         */
        inline const IRBNode* successor() const
        {
            return torwards(RIGHT);
        }

        /**
         * @brief Get the predecessor
         *
         * This method generalizes the algorithm `Tree-Successor` Chapter 12,
         * pp 292, [1]. Its complexity is \f$O(h)\f$ where \f$h\f$ is the height
         * of the tree.
         *
         * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
         *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
         *     (3rd. ed.). The MIT Press.
         *
         * @return a constant pointer to the predecessor node in the in-order
         *      visit
         */
        inline const IRBNode* predecessor() const
        {
            return torwards(LEFT);
        }

        /**
         * @brief Get the next node in the visit
         *
         * This method generalizes the algorithm `Tree-Successor` Chapter 12,
         * pp 292, [1]. Its complexity is \f$O(h)\f$ where \f$h\f$ is the height
         * of the tree.
         *
         * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
         *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
         *     (3rd. ed.). The MIT Press.
         *
         * @param side is a side
         * @return a constant pointer to the next node in the in-order or
         *      reversed in-order visits when `side` is `IRBNode::LEFT` or
         *      `IRBNode::RIGHT`, respectively
         */
        IRBNode* torwards(const Side& side)
        {
            if (child[side] != nullptr) {
                const Side sibling = get_opposite_side(side);
                return child[side]->get_all_on(sibling);
            }

            IRBNode* node = this;
            while (node->is_child(side)) {
                node = node->parent;
            }

            return node->parent;
        }

        /**
         * @brief Get the successor
         *
         * This method implements the algorithm `Tree-Successor` Chapter 12,
         * pp 292, [1]. Its complexity is \f$O(h)\f$ where \f$h\f$ is the height
         * of the tree.
         *
         * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
         *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
         *     (3rd. ed.). The MIT Press.
         *
         * @return a pointer to the next node in the in-order visit
         */
        inline IRBNode* successor()
        {
            return torwards(RIGHT);
        }

        /**
         * @brief Get the predecessor
         *
         * This method generalizes the algorithm `Tree-Successor` Chapter 12,
         * pp 292, [1]. Its complexity is \f$O(h)\f$ where \f$h\f$ is the height
         * of the tree.
         *
         * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
         *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
         *     (3rd. ed.). The MIT Press.
         *
         * @return a pointer to the predecessor node in the in-order visit
         */
        inline IRBNode* predecessor()
        {
            return torwards(LEFT);
        }

        /**
         * @brief Clone a subtree rooted in the this node
         *
         * The complexity of this method is \f$\Theta(n)\f$ where
         * \f$n\f$ is the size of the subtree to be cloned.
         *
         * @param parent is the parent of the new subtree root
         * @return the root of the clone of the subtree rooted
         *      in this node
         */
        IRBNode* clone(IRBNode* parent=nullptr) const
        {
            auto cloned = new IRBNode(*key);
            cloned->parent = parent;
            cloned->color = color;
            cloned->subtree_size = subtree_size;

            for (size_t i=0; i<2; ++i) {
                if (child[i] != nullptr) {
                    cloned->child[i] = child[i]->clone(cloned);
                }
            }

            return cloned;
        }

        /**
         * @brief Destroy the subtree rooted on this node
         *
         * The complexity of this method is \f$\Theta(n)\f$ where
         * \f$n\f$ is the number of nodes in the subtree.
         */
        void destroy_subtree()
        {
            parent = nullptr;

            for (size_t i=0; i<2; ++i) {
                if (child[i] != nullptr) {
                    child[i]->destroy_subtree();
                }
            }

            delete this;
        }

        /**
         * @brief The destroyer
         */
        ~IRBNode()
        {
            if (key != nullptr) {
                ALLOCATOR allocator;

                allocator.deallocate(key, 1);
            }
        }
    };

    IRBNode* root; //!< The tree's root

    /**
     * @brief Find a node by key
     *
     * This method implements the algorithm `Iterative-Tree-Search`
     * Chapter 12, pp 291, [1]. Its complexity is \f$O(h)\f$ where \f$h\f$ is
     * the height of the tree.
     *
     * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
     *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
     *     (3rd. ed.). The MIT Press.
     *
     * @param key is a key
     * @return a constant pointer to the tree node whose key is `key` if
     *      such a node exists. `nullptr` if no node in the tree has key
     *      `key`
     */
    const IRBNode* find_by_key(const KEY& key) const
    {
        IRBNode* node = root;

        COMPARE cmp;
        while (node != nullptr) {
            if (cmp(key, *(node->key))) {
                node = node->child[IRBNode::LEFT];
            } else {
                if (!cmp(*(node->key), key)) {
                    return node;
                }

                node = node->child[IRBNode::RIGHT];
            }
        }

        return node;
    }

    /**
     * @brief Find the \f$i\f$-th node in the tree
     *
     * This method find the \f$i\f$-th node in the in-order visit by using the
     * node's field `subtree_size`.
     *
     * Let \f$n.\textrm{T}\f$ be the subtree rooted in the node \f$n\f$. Moreover, let
     * \f$n.\textrm{left}\f$ and \f$n.\textrm{right}\f$ be the left and right children
     * of \f$n\f$, respectively. Finally, let \f$|n|\f$ be the number of nodes in the
     * subtree rooted in \f$n\f$.
     * If \f$i<|n|\f$, then \f$i\f$-th node of \f$n.T\f$ is:
     * - the \f$i\f$-th node of \f$n.\textrm{left}.\textrm{T}\f$ if and only if that
     *   \f$i < |n.\textrm{left}|\f$;
     * - \f$n\f$ if and only if \f$i = |n.\textrm{left}|\f$;
     * - the \f$(i-|n.\textrm{left}|)\f$-th node of \f$n.\textrm{right}.\textrm{T}\f$
     *   if and only if \f$i > |n.\textrm{left}|\f$.
     *
     * We can identify the \f$i\f$-th node of the tree by successively iteratively
     * the subtree containing it and, possibly, updating the index \f$i\f$. Since
     * each subtree selection takes time \f$O(1)\f$ by using the `subtree_size`,
     * the overall time costs of this method is \f$O(h)\f$ where \f$h\f$ is the
     * height of the tree.
     *
     * @param index is an index
     * @return a constant pointer to the node whose key is the `index`-th in the
     *      sorted set of the tree keys. If `index` is greater than or equal to
     *      the size of the tree, an `std::out_of_range` exception is thrown
     */
    const IRBNode* find_by_index(size_t index) const
    {
        if (index < size()) {
            IRBNode* node = root;
            while (node != nullptr) {
                if (node->child[IRBNode::LEFT] != nullptr) {
                    if (node->child[IRBNode::LEFT]->subtree_size > index) {
                        node = node->child[IRBNode::LEFT];
                    } else {
                        index -= node->child[IRBNode::LEFT]->subtree_size;
                        if (index==0) {
                            return node;
                        }
                        --index;
                        node = node->child[IRBNode::RIGHT];
                    }
                } else {
                    if (index==0) {
                        return node;
                    }
                    --index;
                    node = node->child[IRBNode::RIGHT];
                }
            }
        }

        throw std::out_of_range("The container has " + std::to_string(size())
                                + " keys; requested the element in position "
                                + std::to_string(index) + ".");
    }

    /**
     * @brief Find a node by key
     *
     * This method implements the algorithm `Iterative-Tree-Search`
     * Chapter 12, pp 291, [1]. Its complexity is \f$O(h)\f$ where \f$h\f$ is
     * the height of the tree.
     *
     * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
     *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
     *     (3rd. ed.). The MIT Press.
     *
     * @param key is a key
     * @return a constant pointer to the tree node whose key is `key` if
     *      such a node exists. `nullptr` if no node in the tree has key
     *      `key`
     */
    IRBNode* find_by_key(const KEY& key)
    {
        IRBNode* node = root;

        COMPARE cmp;
        while (node != nullptr) {
            if (cmp(key, *(node->key))) {
                node = node->child[IRBNode::LEFT];
            } else {
                if (!cmp(*(node->key), key)) {
                    return node;
                }

                node = node->child[IRBNode::RIGHT];
            }
        }

        return node;
    }

    /**
     * @brief Find the \f$i\f$-th node in the tree
     *
     * This method find the \f$i\f$-th node in the in-order visit by using the
     * node's field `subtree_size`.
     *
     * Let \f$n.\textrm{T}\f$ be the subtree rooted in the node \f$n\f$. Moreover, let
     * \f$n.\textrm{left}\f$ and \f$n.\textrm{right}\f$ be the left and right children
     * of \f$n\f$, respectively. Finally, let \f$|n|\f$ be the number of nodes in the
     * subtree rooted in \f$n\f$.
     * If \f$i<|n|\f$, then \f$i\f$-th node of \f$n.T\f$ is:
     * - the \f$i\f$-th node of \f$n.\textrm{left}.\textrm{T}\f$ if and only if that
     *   \f$i < |n.\textrm{left}|\f$;
     * - \f$m\f$ if and only if \f$i = |n.\textrm{left}|\f$;
     * - the \f$(i-|n.\textrm{left}|-1)\f$-th node of \f$n.\textrm{right}.\textrm{T}\f$
     *   if and only if \f$i > |n.\textrm{left}|\f$.
     *
     * We can identify the \f$i\f$-th node of the tree by successively iteratively
     * the subtree containing it and, possibly, updating the index \f$i\f$. Since
     * each subtree selection takes time \f$O(1)\f$ by using the `subtree_size`,
     * the overall time costs of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param index is an index
     * @return a constant pointer to the node whose key is the `index`-th in the
     *      sorted set of the tree keys. If `index` is greater than or equal to
     *      the size of the tree, an `std::out_of_range` exception is thrown
     */
    IRBNode* find_by_index(size_t index)
    {
        if (index < size()) {
            IRBNode* node = root;
            while (node != nullptr) {
                if (node->child[IRBNode::LEFT] != nullptr) {
                    if (node->child[IRBNode::LEFT]->subtree_size > index) {
                        node = node->child[IRBNode::LEFT];
                    } else {
                        index -= node->child[IRBNode::LEFT]->subtree_size;
                        if (index==0) {
                            return node;
                        }
                        --index;
                        node = node->child[IRBNode::RIGHT];
                    }
                } else {
                    if (index==0) {
                        return node;
                    }
                    --index;
                    node = node->child[IRBNode::RIGHT];
                }
            }
        }

        throw std::out_of_range("The container has " + std::to_string(size())
                                + " keys; requested the element in position "
                                + std::to_string(index) + ".");
    }

    /**
     * @brief Rotate a subtree on a side
     *
     * This method extends the algorithm `Left-Rotate` Chapter 13, pp 313,
     * [1]. Behind the parametrization of the rotation side, this method
     * also takes care of updating the `subtree_size` fields of the involved
     * nodes.
     *
     * The complexity of this method is \f$O(1)\f$.
     *
     * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
     *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
     *     (3rd. ed.). The MIT Press.
     *
     * @param x is the root of the subtree
     * @param side is the side on which the subtree is rotated
     */
    void rotate_on(IRBNode* x, const typename IRBNode::Side side)
    {
        const auto opposite_side = IRBNode::get_opposite_side(side);
        IRBNode* y = x->child[opposite_side];

        if (y == nullptr) {
            throw std::domain_error("Cannot perform rotation on the "
                                    + IRBNode::side_str(side) +".");
        }
        x->child[opposite_side] = y->child[side];

        if (y->child[side] != nullptr) {
            y->child[side]->parent = x;
        }

        y->parent = x->parent;
        if (x->parent == nullptr) {
            root = y;
        } else {
            replace_by(x, y);
        }

        y->child[side] = x;
        x->parent = y;

        // update subtree sizes
        y->subtree_size = x->subtree_size;
        if (y->child[opposite_side] != nullptr) {
            x->subtree_size -= (1 + y->child[opposite_side]->subtree_size);
        } else {
            --x->subtree_size;
        }
    }

    //!< @private
    inline IRBNode* insert_fixup_case1(IRBNode*& parent, IRBNode*& uncle)
    {
        // case 1
        parent->color = IRBNode::BLACK;
        uncle->color = IRBNode::BLACK;
        parent->parent->color = IRBNode::RED;

        return parent->parent;
    }

    //!< @private
    inline IRBNode* insert_fixup_case3(IRBNode*& parent, const typename IRBNode::Side& uncle_side)
    {
        parent->color = IRBNode::BLACK;
        parent->parent->color = IRBNode::RED;

        rotate_on(parent->parent, uncle_side);

        return parent;
    }

    /**
     * @brief Fix the coloring of the tree after an insertion
     *
     * This method implements the algorithm `RB-Insert-Fixup` Chapter 13,
     * pp 316, [1]. The complexity of this method is
     * \f$O(\log \texttt{size()})\f$.
     *
     * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
     *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
     *     (3rd. ed.). The MIT Press.
     *
     * @param node is the inserted node
     */
    void insert_fixup(IRBNode* node)
    {
        IRBNode* parent = node->parent;

        while (parent != nullptr) {
            if (parent->color == IRBNode::BLACK) {
                return;
            }

            auto parent_side = parent->get_side();
            auto uncle_side = IRBNode::get_opposite_side(parent_side);

            IRBNode* uncle = parent->parent->child[uncle_side];

            auto uncle_color = IRBNode::BLACK;
            if (uncle != nullptr) {
                uncle_color = uncle->color;
            }

            if (uncle_color == IRBNode::RED) {
                // case 1
                node = insert_fixup_case1(parent, uncle);
            } else {
                if (node->is_child(uncle_side)) {
                    // case 2
                    node = node->parent;
                    rotate_on(node, parent_side);
                }

                // case 3
                insert_fixup_case3(node->parent, uncle_side);
            }

            parent = node->parent;
        }

        root->color = IRBNode::BLACK;
    }

    /**
     * @brief Replace a node by another node
     *
     * This method takes time \f$O(1)\f$.
     *
     * @param node is the node to be replaced
     * @param alternative is the node that will replace
     *      `node`
     */
    void replace_by(IRBNode* node, IRBNode* alternative)
    {
        IRBNode* parent = node->parent;

        if (parent == nullptr) {
            root = alternative;
        } else {
            if (node == parent->child[IRBNode::LEFT]) {
                parent->child[IRBNode::LEFT] = alternative;
            } else {
                parent->child[IRBNode::RIGHT] = alternative;
            }
        }
    }

    /**
     * @brief Delete a leaf of the tree
     *
     * This method deletes a leaf from the tree and updates the
     * fields `subtree_size` of all its ancestors.
     *
     * The complexity of this method is \f$\Theta(\log \texttt{size()})$.
     *
     * @param leaf is the node to be removed
     */
    void remove_leaf(IRBNode* leaf)
    {
        if (leaf->parent == nullptr) {
            root = nullptr;
        } else {
            leaf->parent->delta_size(-1);
            replace_by(leaf, nullptr);
        }

        delete leaf;
    }

    /**
     * @brief Delete a node and replace it by one of its children
     *
     * This method deletes a node having only one child replacing it
     * by its child. Moreover, the method updates the fields
     * `subtree_size` of all its ancestors.
     *
     * The time complexity of this method \f$O(\log \texttt{size()})\f$.
     *
     * @param node is the node to be delete
     * @param side is the side of the child
     * @return a pointer to the `node`'s child that replaced `node`
     */
    IRBNode* transplant(IRBNode* node, const typename IRBNode::Side side)
    {
        IRBNode* child = node->child[side];
        if (child != nullptr) {
            child->parent = node->parent;
            replace_by(node, child);
            child->parent->delta_size(-1);
        }

        if (node->parent == nullptr) {
            root = child;
        }

        delete node;

        return child;
    }

    /**
     * @brief Delete a node, replace it by one of its children, and fix the tree coloring
     *
     * This method deletes a node having only one child replacing it by its child.
     * Moreover, the method fixes the coloring of the tree nodes. Finally, it updates
     * the fields `subtree_size` of all its ancestors.
     *
     * The time complexity of this method \f$O(\log \texttt{size()})\f$.
     *
     * @param node is the node to be delete
     * @param side is the side of the child
     */
    void transplant_and_fix(IRBNode* node, const typename IRBNode::Side side)
    {
        bool to_fix = IRBNode::is_black(node);

        node = transplant(node, side);

        if (node != nullptr && to_fix) {
            remove_fixup(node);
        }
    }

    //!< @private
    inline void remove_fixup_case1(IRBNode*& sibling,
                                   const typename IRBNode::Side& node_side,
                                   const typename IRBNode::Side& sibling_side)
    {
        sibling->color = IRBNode::BLACK;
        auto parent = sibling->parent;
        parent->color = IRBNode::RED;
        rotate_on(parent, node_side);
        sibling = parent->child[sibling_side];
    }

    //!< @private
    inline void remove_fixup_case3(IRBNode*& node, IRBNode*& sibling,
                                   const typename IRBNode::Side& node_side,
                                   const typename IRBNode::Side& sibling_side)
    {
        sibling->child[node_side]->color = IRBNode::BLACK;
        sibling->color = IRBNode::RED;
        rotate_on(sibling, sibling_side);
        sibling = node->parent->child[sibling_side];
    }

    //!< @private
    inline void remove_fixup_case4(IRBNode*& node, IRBNode*& sibling,
                                   const typename IRBNode::Side& node_side,
                                   const typename IRBNode::Side& sibling_side)
    {
        sibling->color = node->parent->color;
        node->parent->color = IRBNode::BLACK;
        sibling->child[sibling_side]->color = IRBNode::BLACK;
        rotate_on(node->parent, node_side);
    }

    /**
     * @brief Fix the tree coloring after a node removal
     *
     * This method extends the algorithm `RB-Delete-Fixup` Chapter 13,
     * pp 326, [1] by supporting the `subtree_size` field.
     * Please, notice that this method is invocated before the actual
     * removal when the node to be removed is a leaf.
     *
     * The time complexity of this method \f$O(\log \texttt{size()})\f$.
     *
     * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
     *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
     *     (3rd. ed.). The MIT Press.
     *
     * @param node is the pointer of the node to be removed
     */
    void remove_fixup(IRBNode* node)
    {
        if (node != nullptr) {
            while (!node->is_root() && node->color == IRBNode::BLACK) {
                const auto node_side = node->get_side();
                const auto sibling_side = IRBNode::get_opposite_side(node_side);

                IRBNode* sibling = node->parent->child[sibling_side];
                if (sibling->color == IRBNode::RED) {
                    // case 1
                    remove_fixup_case1(sibling, node_side, sibling_side);
                }

                if (IRBNode::is_black(sibling->child[node_side])
                        && IRBNode::is_black(sibling->child[sibling_side])) {

                    // case 2
                    sibling->color = IRBNode::RED;
                    node = node->parent;
                } else {
                    if (IRBNode::is_black(sibling->child[sibling_side])) {
                        // case 3
                        remove_fixup_case3(node, sibling, node_side, sibling_side);
                    }

                    // case 4
                    remove_fixup_case4(node, sibling, node_side, sibling_side);
                    node = root;
                }
            }

            node->color = IRBNode::BLACK;
        }
    }

    /**
     * @brief Remove a node from the tree and delete it
     *
     * This method extends the algorithm `RB-Insert` Chapter 13, pp 315,
     * [1] by supporting the `subtree_size` field.
     *
     * The time complexity of this method \f$O(\log \texttt{size()})\f$.
     *
     * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
     *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
     *     (3rd. ed.). The MIT Press.
     *
     * @param node is the node to be removed from the tree and deleted
     */
    void remove(IRBNode* node)
    {
        if (node->child[IRBNode::LEFT] == nullptr) {
            // the node does not have the left child

            if (node->child[IRBNode::RIGHT] == nullptr) {
                // the node is a leaf
                if (node->color == IRBNode::BLACK) {
                    remove_fixup(node);
                }

                remove_leaf(node);
            } else {
                // the node has only the right child
                transplant_and_fix(node, IRBNode::RIGHT);
            }
        } else {
            // the node has the left child
            if (node->child[IRBNode::RIGHT] == nullptr) {
                // the node has only the left child
                transplant_and_fix(node, IRBNode::LEFT);
            } else {
                // the node has both the children
                auto succ_ptr = node->successor();

                std::swap(node->key, succ_ptr->key);

                remove(succ_ptr);
            }
        }
    }

public:

    /**
     * @brief A constant iterator type for the tree nodes
     */
    class const_iterator
    {
        IRBNode const* root;    //!< The tree root
        IRBNode const* node;    //!< The referenced node

        /**
         * @brief Construct a new constant iterator
         *
         * @param root is the tree root
         * @param node is the referenced node
         */
        explicit const_iterator(IRBNode const *root, IRBNode const* node):
            root{root}, node{node}
        {}
    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   KEY;
        using pointer           =   const KEY*;
        using reference         =   const KEY&;
        using iterator_category =   std::bidirectional_iterator_tag;

        /**
         * @brief An empty constructor
         */
        const_iterator():
            root{nullptr}, node{nullptr}
        {}

        /**
         * @brief Reference operator
         *
         * @return a reference to the species pointer by the iterator
         */
        inline reference operator*() const
        {
            return *(node->key);
        }

        /**
         * @brief Pointer operator
         *
         * @return a pointer to the species pointer by the iterator
         */
        inline pointer operator->() const
        {
            return node->key;
        }

        /**
         * @brief The prefix increment
         *
         * @return a reference to the updated object
         */
        const_iterator& operator++()
        {
            if (node != nullptr) {
                node = node->successor();
            } else {
                if (root != nullptr) {
                    node = root->get_subtree_min();
                }
            }

            return *this;
        }

        /**
         * @brief The postfix increment
         *
         * @return a copy of the original object
         */
        inline const_iterator operator++(int)
        {
            const_iterator copy(*this);

            this->operator++();

            return copy;
        }

        /**
         * @brief The prefix decrement
         *
         * @return a reference to the updated object
         */
        const_iterator& operator--()
        {
            if (node != nullptr) {
                node = node->predecessor();
            } else {
                if (root != nullptr) {
                    node = root->get_subtree_max();
                }
            }

            return *this;
        }

        /**
         * @brief The postfix decrement
         *
         * @return a copy of the original object
         */
        inline const_iterator operator--(int)
        {
            const_iterator copy(*this);

            this->operator--();

            return copy;
        }

        /**
         * @brief Test whether two iterators are the same
         *
         * @param a is the first iterator to compare
         * @param b is the second iterator to compare
         * @return `true` if and only if the two iterators
         *      refer to the same object
         */
        friend inline bool operator==(const const_iterator& a, const const_iterator& b)
        {
            return (a.node == b.node);
        }

        /**
         * @brief Test whether two iterators are the same
         *
         * @param a is the first iterator to compare
         * @param b is the second iterator to compare
         * @return `true` if and only if the two iterators
         *      refer to the same object
         */
        friend inline bool operator!=(const const_iterator& a, const const_iterator& b)
        {
            return (a.node != b.node);
        }

        template<typename KEY2, typename VALUE2, typename COMPARE2> friend class IRBTree;
    };


    /**
     * @brief An iterator type for the tree nodes
     */
    class iterator
    {
        IRBNode* root;    //!< The tree root
        IRBNode* node;    //!< The referenced node

        /**
         * @brief Construct a new constant iterator
         *
         * @param root is the tree root
         * @param node is the referenced node
         */
        explicit iterator(IRBNode *root, IRBNode* node):
            root{root}, node{node}
        {}
    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   KEY;
        using pointer           =   KEY*;
        using reference         =   KEY&;
        using iterator_category =   std::bidirectional_iterator_tag;

        /**
         * @brief An empty constructor
         */
        iterator():
            root{nullptr}, node{nullptr}
        {}

        /**
         * @brief Reference operator
         *
         * @return a reference to the species pointer by the iterator
         */
        inline reference operator*()
        {
            return *(node->key);
        }

        /**
         * @brief Pointer operator
         *
         * @return a pointer to the species pointer by the iterator
         */
        inline pointer operator->()
        {
            return node->key;
        }

        /**
         * @brief The prefix increment
         *
         * @return a reference to the updated object
         */
        iterator& operator++()
        {
            if (node != nullptr) {
                node = node->successor();
            } else {
                if (root != nullptr) {
                    node = root->get_subtree_min();
                }
            }

            return *this;
        }

        /**
         * @brief The postfix increment
         *
         * @return a copy of the original object
         */
        inline iterator operator++(int)
        {
            iterator copy(*this);

            this->operator++();

            return copy;
        }

        /**
         * @brief The prefix decrement
         *
         * @return a reference to the updated object
         */
        iterator& operator--()
        {
            if (node != nullptr) {
                node = node->predecessor();
            } else {
                if (root != nullptr) {
                    node = root->get_subtree_max();
                }
            }

            return *this;
        }

        /**
         * @brief The postfix decrement
         *
         * @return a copy of the original object
         */
        inline iterator operator--(int)
        {
            iterator copy(*this);

            this->operator--();

            return copy;
        }

        /**
         * @brief Test whether two iterators are the same
         *
         * @param a is the first iterator to compare
         * @param b is the second iterator to compare
         * @return `true` if and only if the two iterators
         *      refer to the same object
         */
        friend inline bool operator==(const iterator& a, const iterator& b)
        {
            return (a.node == b.node);
        }

        /**
         * @brief Test whether two iterators are the same
         *
         * @param a is the first iterator to compare
         * @param b is the second iterator to compare
         * @return `true` if and only if the two iterators
         *      refer to the same object
         */
        friend inline bool operator!=(const iterator& a, const iterator& b)
        {
            return (a.node != b.node);
        }

        template<typename KEY2, typename VALUE2, typename COMPARE2> friend class IRBTree;
    };

protected:

    /**
     * @brief Insert a new key in the tree
     *
     * This method implements the algorithm `RB-Insert` Chapter 13, pp 315,
     * [1]. The complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
     *     Clifford Stein. 2009. Introduction to Algorithms, Third Edition
     *     (3rd. ed.). The MIT Press.
     *
     * @param key is the to-be-inserted key
     * @return a pointer to the inserted node
     */
    std::pair<iterator, bool> insert_node(KEY&& key)
    {
        if (root == nullptr) {
            root = new IRBNode(std::move(key));
            root->color = IRBNode::BLACK;

            return std::make_pair(iterator(root, root), true);
        }

        IRBNode* node = root;
        IRBNode* parent = root;

        COMPARE cmp;
        while (node != nullptr) {
            parent = node;
            if (cmp(key, *(node->key))) {
                node = node->child[IRBNode::LEFT];
            } else {
                if (!cmp(*(node->key), key)) {
                    return std::make_pair(iterator(root, node), false);
                }

                node = node->child[IRBNode::RIGHT];
            }
        }

        auto new_child = new IRBNode(std::move(key));
        if (cmp(key, *(parent->key))) {
            parent->set_child(new_child, IRBNode::LEFT);
        } else {
            parent->set_child(new_child, IRBNode::RIGHT);
        }

        insert_fixup(new_child);

        return std::make_pair(iterator(root, new_child), true);
    }

public:

    /**
     * @brief The empty constructor
     */
    IRBTree():
        root{nullptr}
    {}

    /**
     * @brief The move constructor
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @param orig is the original tree
     */
    IRBTree(IRBTree&& orig):
        root{orig.root}
    {
        orig.root = nullptr;
    }

    /**
     * @brief The copy constructor
     *
     * The time complexity of this method is \f$O(\texttt{size()})\f$.
     *
     * @param orig is the original tree
     */
    IRBTree(const IRBTree& orig):
        IRBTree()
    {
        if (orig.root != nullptr) {
            root = orig.root->clone();
        }
    }

    /**
     * @brief The move operator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @param orig is the original tree
     * @return a reference to the updated object
     */
    inline IRBTree& operator=(IRBTree&& orig)
    {
        std::swap(orig.root, root);

        return *this;
    }

    /**
     * @brief The copy operator
     *
     * The time complexity of this method is \f$O(\texttt{size()})\f$.
     *
     * @param orig is the original tree
     * @return a reference to the updated object
     */
    IRBTree& operator=(const IRBTree& orig)
    {
        if (root != nullptr) {
            root->destroy_subtree();
        }

        if (orig.root != nullptr) {
            root = orig.root->clone();
        }

        return *this;
    }

    /**
     * @brief The number of nodes in the tree
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return the number of nodes in the tree
     */
    size_t size() const
    {
        if (root==nullptr) {
            return 0;
        }

        return root->subtree_size;
    }

    /**
     * @brief Count how many nodes contains a key
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is the key whose node number is aimed
     * @return the number of nodes whose key is `key`
     */
    size_t count(const KEY& key) const
    {
        if (find_by_key(key)==nullptr) {
            return 0;
        }

        return 1;
    }

    /**
     * @brief Insert a key
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is the key to be inserted
     * @return a pair whose first value is an iterator to a node
     *      whose key is `key` and whose second value is a Boolean
     *      value: `true` if the insertion succeeds; `false`
     *      otherwise.
     */
    inline std::pair<iterator, bool> insert(KEY&& key)
    {
        return insert_node(std::move(key));
    }

    /**
     * @brief Insert a key
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is the key to be inserted
     * @return a pair whose first value is an iterator to a node
     *      whose key is `key` and whose second value is a Boolean
     *      value: `true` if the insertion succeeds; `false`
     *      otherwise.
     */
    inline std::pair<iterator, bool> insert(KEY key)
    {
        return insert_node(std::move(key));
    }

    /**
     * @brief Find a key in the tree
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is a key value
     * @return a constant iterator refering to the
     *      node whose key is `key` if it exists or
     *      `end()` otherwise.
     */
    inline const_iterator find(const KEY& key) const
    {
        return const_iterator(root, find_by_key(key));
    }

    /**
     * @brief Find a key in the tree
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is a key value
     * @return an iterator refering to the node whose
     *      key is `key` if it exists or `end()`
     *      otherwise.
     */
    inline iterator find(const KEY& key)
    {
        return iterator(root, find_by_key(key));
    }

    /**
     * @brief Find the i-th node in the in-order visit
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param index is the index
     * @return a reference to the `index`-th node in the in-order visit
     *      of the tree
     */
    inline KEY& operator[](const size_t index)
    {
        return *(find_by_index(index)->key);
    }

    /**
     * @brief Find the i-th node in the in-order visit
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param index is the index
     * @return a constant reference to the `index`-th node in the in-order
     *      visit of the tree
     */
    inline const KEY& operator[](const size_t index) const
    {
        return *(find_by_index(index)->key);
    }

    /**
     * @brief Find the i-th node in the in-order visit
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param index is the index
     * @return an iterator refering to the `index`-th node in the in-order visit
     *      of the tree
     */
    inline iterator get(const size_t index)
    {
        return iterator(root, find_by_index(index));
    }

    /**
     * @brief Find the i-th node in the in-order visit
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param index is the index
     * @return a constant iterator refering to the `index`-th node in
     *      the in-order visit of the tree
     */
    inline const_iterator get(const size_t index) const
    {
        return const_iterator(root, find_by_index(index));
    }

    /**
     * @brief Delete a key from the tree
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is the key to be removed
     * @return the number of nodes removed: 1 if a node has been removed,
     *      0 if no node in the tree has `key` as the key
     */
    size_t erase(const KEY& key)
    {
        auto node = find_by_key(key);
        if (node != nullptr) {
            remove(node);

            return 1;
        }

        return 0;
    }

    /**
     * @brief Delete a node from the tree
     *
     * @param it is an iterator of the tree
     * @return the iterator following `it`
     */
    iterator erase(const iterator& it)
    {
        if (it.root != root) {
            throw std::domain_error("The iterator does not refer to the tree");
        }

        iterator next_it = it;
        ++next_it;
        if (it.node != nullptr) {
            remove(it.node);
        }

        return next_it;
    }

    /**
     * @brief Get the first in-order constant iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return an in-order constant iterator refering to the minimal value node
     */
    const_iterator begin() const noexcept
    {
        if (root == nullptr) {
            return const_iterator();
        }

        return const_iterator(root, root->get_subtree_min());
    }

    /**
     * @brief Get the first in-order iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return an in-order iterator refering to the minimal value node
     */
    iterator begin() noexcept
    {
        if (root == nullptr) {
            return iterator();
        }

        return iterator(root, root->get_subtree_min());
    }

    /**
     * @brief Get the first in-order constant iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return an in-order constant iterator refering to the minimal value node
     */
    inline const_iterator cbegin() const noexcept
    {
        return begin();
    }

    /**
     * @brief Get the last in-order iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return the last in-order iterator
     */
    inline iterator end() noexcept
    {
        return iterator(root, nullptr);
    }

    /**
     * @brief Get the last in-order constant iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return the last in-order constant iterator
     */
    inline const_iterator end() const noexcept
    {
        return const_iterator(root, nullptr);
    }

    /**
     * @brief Get the last in-order constant iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return the last in-order constant iterator
     */
    inline const_iterator cend() const noexcept
    {
        return end();
    }

    /**
     * @brief Remove all the node from the tree
     *
     * The time complexity of this method is \f$\Theta(\texttt{size()})\f$.
     */
    void clear()
    {
        if (root != nullptr) {
            root->destroy_subtree();

            root = nullptr;
        }
    }

    /**
     * @brief Destroy the IRBTree object
     *
     * The time complexity of this method is \f$\Theta(\texttt{size()})\f$.
     */
    ~IRBTree()
    {
        clear();
    }
};

}   // CLONES

#endif // __CLONES_IRB_TREE__