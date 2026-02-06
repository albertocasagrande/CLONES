/**
 * @file imap.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines indexed map
 * @version 1.2
 * @date 2026-02-06
 *
 * @copyright Copyright (c) 2024-2025
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

#ifndef __CLONES_IMAP__
#define __CLONES_IMAP__

#include "irb_tree.hpp"

namespace CLONES
{

/**
 * @brief An indexed map type
 *
 * @tparam KEY is the type of the RB tree keys
 * @tparam COMPARE is the order among keys
 * @tparam ALLOCATOR is the key allocator
 */
template<class KEY, class VALUE, class COMPARE=std::less<KEY>, class ALLOCATOR=std::allocator<std::pair<const KEY,VALUE>>>
class imap
{
    using tree_key = std::pair<const KEY, VALUE>;

    struct tree_key_less
    {
        constexpr bool operator()(const tree_key& a, const tree_key& b) const
        {
            COMPARE cmp;

            return cmp(a.first, b.first);
        }
    };

    using tree_type = IRBTree<tree_key, tree_key_less, ALLOCATOR>;

    /**
     * @brief Build a tree key by using a map key
     *
     * @param key
     * @return tree_key
     */
    static tree_key build_tree_key(const KEY& key)
    {
        std::pair<KEY, VALUE> p;
        p.first = key;

        return static_cast<tree_key>(p);
    }

    tree_type tree; //!< The tree containing the map association
public:
    using iterator = typename tree_type::iterator;
    using const_iterator = typename tree_type::const_iterator;

    /**
     * @brief The empty constructor
     */
    imap():
        tree()
    {}

    /**
     * @brief The copy constructor
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @param orig is the original map
     */
    imap(const imap& orig):
        tree(orig.tree)
    {}

    /**
     * @brief The move constructor
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @param orig is the original map
     */
    imap(imap&& orig):
        tree(std::move(orig.tree))
    {}

    /**
     * @brief The copy operator
     *
     * The time complexity of this method is \f$O(\texttt{size()})\f$.
     *
     * @param orig is the original map
     * @return a reference to the updated object
     */
    inline imap& operator=(const imap& orig)
    {
        tree = orig.tree;

        return *this;
    }

    /**
     * @brief The move operator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @param orig is the original map
     * @return a reference to the updated object
     */
    inline imap& operator=(imap&& orig)
    {
        tree = std::move(orig.tree);

        return *this;
    }

    /**
     * @brief The size of the map
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return The number of keys in the map
     */
    inline size_t size() const
    {
        return tree.size();
    }

    /**
     * @brief Count how many key are stored in the map
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is the key whose node number is aimed
     * @return the number of nodes whose key is `key`
     */
    inline size_t count(const KEY& key) const
    {
        const auto tkey = build_tree_key(key);

        return tree.count(tkey);
    }

    /**
     * @brief Insert a pair in the map
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param p is the pair to be inserted
     * @return a pair whose first value is an iterator to a node whose
     *      first value is `p->first` and whose second value is a
     *      Boolean value: `true` if the insertion succeeds; `false`
     *      otherwise.
     */
    inline std::pair<iterator, bool> insert(const std::pair<KEY,VALUE>& p)
    {
        auto const_p = static_cast<tree_key>(p);
        return tree.insert(const_p);
    }

    /**
     * @brief Find a the pair whose first value is the searched key
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is a key
     * @return a constant iterator referring to the node whose first
     *      value is `key` if it exists or `end()` otherwise.
     */
    inline const_iterator find(const KEY& key) const
    {
        const auto tkey = build_tree_key(key);

        return tree.find(tkey);
    }

    /**
     * @brief Find a the pair whose first value is the searched key
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is a key
     * @return an iterator referring to the pair whose first value is
     *      `key` if it exists or `end()` otherwise.
     */
    inline iterator find(const KEY& key)
    {
        const auto tkey = build_tree_key(key);

        return tree.find(tkey);
    }

    /**
     * @brief Get the value associated to a key
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is a key
     * @return a reference to the value associated to `key`
     */
    VALUE& operator[](const KEY& key)
    {
        const auto tkey = build_tree_key(key);

        auto found = tree.find(tkey);
        if (found == nullptr) {
           found = tree.insert(tkey)->first;
        }

        return found->second;
    }

    /**
     * @brief Get an iterator to the n-th pair in the map
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param index is an index
     * @return an iterator to the n-th pair in key-in-order
     *      visit of the map pairs
     */
    inline iterator get(const size_t index)
    {
        return tree.get(index);
    }

    /**
     * @brief Get an iterator to the n-th pair in the map
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param index is an index
     * @return a constant iterator to the n-th pair in key-in-order
     *      visit of the map pairs
     */
    inline const_iterator get(const size_t index) const
    {
        return tree.get(index);
    }

    /**
     * @brief Delete a pair from the map
     *
     * The time complexity of this method is \f$O(\log \texttt{size()})\f$.
     *
     * @param key is the first value in the pair to be removed
     * @return the number of key removed: 1 if key has been removed,
     *      0 if no pair in the map has `key` as the first value
     */
    inline size_t erase(const KEY& key)
    {
        const auto tkey = build_tree_key(key);

        return tree.erase(tkey);
    }

    /**
     * @brief Delete a pair from the map
     *
     * @param it is an iterator of the map
     * @return the iterator following `it`
     */
    inline iterator erase(const iterator& it)
    {
        return tree.erase(it);
    }

    /**
     * @brief Get the first in-order iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return a constant in-order iterator referring to the
     *      first pair in the map
     */
    inline const_iterator begin() const noexcept
    {
        return tree.begin();
    }

    /**
     * @brief Get the first in-order iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return a constant in-order iterator referring to the
     *      first pair in the map
     */
    inline const_iterator cbegin() const noexcept
    {
        return tree.cbegin();
    }

    /**
     * @brief Get the first in-order iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return an in-order iterator referring to the first
     *      pair in the map
     */
    inline iterator begin() noexcept
    {
        return tree.begin();
    }

    /**
     * @brief Get the last in-order iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return the last in-order iterator
     */
    inline const_iterator end() const noexcept
    {
        return tree.end();
    }

    /**
     * @brief Get the last in-order iterator
     *
     * The time complexity of this method is \f$O(1)\f$.
     *
     * @return the last in-order iterator
     */
    inline const_iterator cend() const noexcept
    {
        return tree.cend();
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
        return tree.end();
    }

    /**
     * @brief Clean the map
     *
     * The time complexity of this method is \f$O(\texttt{size()})\f$.
     */
    inline void clear()
    {
        tree.clear();
    }

    /**
     * @brief Destroy the imap object
     *
     * The time complexity of this method is \f$\Theta(\texttt{size()})\f$.
     */
    ~imap()
    {
        clear();
    }
};

}

#endif // __CLONES_IMAP__
