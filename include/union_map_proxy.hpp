/**
 * @file union_map_proxy.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines union map representation
 * @version 1.1
 * @date 2026-07-14
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

#ifndef __CLONES_UNION_MAP_PROXY__
#define __CLONES_UNION_MAP_PROXY__

#include <map>


namespace CLONES
{

/**
 * @brief A template class to represent the union of two maps
 *
 * The objects of this class class mimic the union of two maps.
 *
 * @tparam KEY the key type
 * @tparam VALUE the value type
 * @tparam COMPARE the compare class
 * @tparam ALLOCATOR the allocator class
 */
template<class KEY, class VALUE, class COMPARE = std::less<KEY>,
         class ALLOCATOR = std::allocator<std::pair<const KEY, VALUE>>>
class union_map_proxy
{
    using map_type = std::map<KEY, VALUE, COMPARE, ALLOCATOR>;

    map_type const& map1;   //!< The first map in the union
    map_type const& map2;   //!< The second map in the union

public:
    /**
     * @brief A constant iterator for the union map proxy
     */
    class const_iterator
    {
    public:
        using difference_type   =   std::ptrdiff_t;
        using map_type = std::map<KEY, VALUE, COMPARE, ALLOCATOR>;
        using value_type        =   map_type::value_type;
        using allocator_type	=   ALLOCATOR;
        using const_map_iterator	=  map_type::const_iterator;
        using const_pointer     =   std::allocator_traits<ALLOCATOR>::const_pointer;
        using reference         =   const value_type&;
        using iterator_category =   std::bidirectional_iterator_tag;

        /**
         * @brief The empty constructor
         */
        const_iterator():
            union_obj{nullptr}
        {}

        /**
         * @brief The reference operator
         *
         * @return The reference to the pointed object
         */
        reference operator*() const
        {
            return *(this->operator->());
        }

        /**
         * @brief Get the pointer operator
         *
         * @return the pointer operator
         */
        inline const_pointer operator->() const
        {
            return get_current_it().operator->();
        }

        /**
         * @brief Prefix increment
         *
         * @return a reference to the updated
         *      constant iterator
         */
        const_iterator& operator++()
        {
            ++(get_current_it());

            return *this;
        }

        /**
         * @brief Postfix increment
         *
         * @return a copy of the original iterator
         */
        const_iterator operator++(int)
        {
            const_iterator orig{union_obj, map1_it, map2_it};

            this->operator++();

            return orig;
        }

        /**
         * @brief Prefix decrement
         *
         * @return a reference to the updated
         *      constant iterator
         */
        const_iterator& operator--()
        {
            if (map1_it == union_obj->map1.begin()) {
                --map2_it;

                return *this;
            }

            if (map2_it == union_obj->map2.begin()) {
                --map1_it;

                return *this;
            }

            auto map1_it_copy{map1_it};
            auto map2_it_copy{map2_it};

            --map1_it_copy;
            --map2_it_copy;

            COMPARE cmp;
            if (cmp(map2_it_copy->first, map1_it_copy->first)) {
                --map1_it;
            } else {
                --map2_it;
            }

            return *this;
        }

        /**
         * @brief Postfix decrement
         *
         * @return a copy of the original iterator
         */
        const_iterator operator--(int)
        {
            const_iterator orig{union_obj, map1_it, map2_it};

            this->operator--();

            return orig;
        }

        /**
         * @brief Test whether two constant iterators are equal
         *
         * @param obj a constant iterator
         * @return `true` if and only if the two constant iterators
         *   point to the same object
         */
        bool operator==(const const_iterator& obj) const
        {
            return (union_obj == obj.union_obj
                    && map1_it == obj.map1_it
                    && map2_it == obj.map2_it);
        }

        /**
         * @brief Test whether two constant iterators differ
         *
         * @param obj a constant iterator
         * @return `true` if and only if the two constant iterators
         *   point to different objects
         */
        inline bool operator!=(const const_iterator& obj) const
        {
            return !(this->operator==(obj));
        }

    private:
        union_map_proxy<KEY, VALUE, COMPARE, ALLOCATOR> const* union_obj;

        map_type::const_iterator map1_it;    //!< The first map iterator
        map_type::const_iterator map2_it;    //!< The second map iterator

        /**
         * @brief A constructor
         *
         * @param union_obj is a pointer to the iterated union map
         * @param map1_it is a constant iterator to the first map in the union
         * @param map2_it is a constant iterator to the second map in the union
         */
        const_iterator(const union_map_proxy<KEY, VALUE, COMPARE, ALLOCATOR>* union_obj,
                       std::map<KEY, VALUE, COMPARE, ALLOCATOR>::const_iterator map1_it,
                       std::map<KEY, VALUE, COMPARE, ALLOCATOR>::const_iterator map2_it):
            union_obj{union_obj}, map1_it{map1_it}, map2_it{map2_it}
        {}

        /**
         * @brief Get the current map iterator
         *
         * @return A constant reference to the current map iterator
         */
        const const_map_iterator& get_current_it() const
        {
            if (map1_it == union_obj->map1.end()) {
                return map2_it;
            }

            if (map2_it == union_obj->map2.end()) {
                return map1_it;
            }

            COMPARE cmp;
            if (cmp(map2_it->first, map1_it->first)) {
                return map2_it;
            }
            return map1_it;
        }

        /**
         * @brief Get the current map iterator
         *
         * @return A non-constant reference to the current map iterator
         */
        inline const_map_iterator& get_current_it()
        {
            return const_cast<const_map_iterator &>(
                static_cast<const const_iterator*>(this)->get_current_it()
            );
        }

        friend class union_map_proxy<KEY, VALUE>;
    };

    /**
     * @brief A constructor
     *
     * @param map1 the first map in the union
     * @param map2 the second map in the union
     */
    union_map_proxy(const std::map<KEY, VALUE>& map1, const std::map<KEY, VALUE>& map2):
        map1{map1}, map2{map2}
    {}

    /**
     * @brief Get a constant iterator pointing to the first element in the union
     *
     * @return the constant iterator pointing to the first element in the union
     */
    const_iterator begin() const
    {
        return const_iterator{this, map1.begin(), map2.begin()};
    }

    /**
     * @brief Get a constant iterator pointing to the first element in the union
     *
     * @return the constant iterator pointing to the first element in the union
     */
    const_iterator end() const
    {
        return const_iterator{this, map1.end(), map2.end()};
    }

    /**
     * @brief Get an iterator to the first element greater than the given key
     *
     * @param key the key of the aimed value
     * @return a constant iterator to the first element greater than the given key
     */
    inline const_iterator upper_bound(KEY&& key) const
    {
        return upper_bound(key);
    }

    /**
     * @brief Get an iterator to the first element greater than the given key
     *
     * @param key the key of the aimed value
     * @return a constant iterator to the first element greater than the given key
     */
    const_iterator upper_bound(const KEY& key) const
    {
        return {this, map1.upper_bound(key),
                map2.upper_bound(key)};
    }

    /**
     * @brief Get an iterator to the first element not less than the given key
     *
     * @param key the key of the aimed value
     * @return a constant iterator to the first element not less than the given key
     */
    inline const_iterator lower_bound(KEY&& key) const
    {
        return lower_bound(key);
    }

    /**
     * @brief Get an iterator to the first element not less than the given key
     *
     * @param key the key of the aimed value
     * @return a constant iterator to the first element not less than the given key
     */
    const_iterator lower_bound(const KEY& key) const
    {
        return {this, map1.lower_bound(key),
                map2.lower_bound(key)};
    }

    /**
     * @brief Finds element with specific key
     *
     * @param key the key of the aimed value
     * @return a constant iterator to the aimed value, if one of the two maps
     *   contains `key`, or to the end of the union, otherwise.
     */
    const_iterator find(KEY&& key) const
    {
        return find(key);
    }

    /**
     * @brief Finds element with specific key
     *
     * @param key the key of the aimed value
     * @return a constant iterator to the aimed value, if one of the two maps
     *   contains `key`, or to the end of the union, otherwise.
     */
    const_iterator find(const KEY& key) const
    {
        auto map1_it = map1.lower_bound(key);
        auto map2_it = map2.lower_bound(key);

        if (map1_it != map1.end()) {
            if (map1_it->first == key) {
                return {this, map1_it, map2_it};
            }
        }

        if (map2_it != map2.end()) {
            if (map2_it->first == key) {
                return {this, map1_it, map2_it};
            }
        }

        return end();
    }
};

} // CLONES

#endif // __CLONES_UNION_MAP_PROXY__