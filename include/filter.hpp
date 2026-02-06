/**
 * @file filter.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines filters
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

#ifndef __CLONES_FILTER__
#define __CLONES_FILTER__

#include <set>
#include <type_traits>

namespace CLONES
{

/**
 * @brief A base filter class
 *
 * This class is meant to represent set filters. It provides the
 * `filtered` method which returns `true` when an object it is
 * filtered by the class object. The `filtered` method of this class
 * always returns `false`.
 *
 * @tparam OBJECT is the type of the objects in the set to be filtered
 */
template<typename OBJECT>
struct BaseFilter
{
    /**
     * @brief The empty constructor
     */
    BaseFilter()
    {}

    /**
     * @brief Establish whether an object must be filtered by the set
     *
     * @param object is the object to be tested
     * @return `false`
     */
    inline bool filtered(const OBJECT& object) const
    {
        (void)object;

        return false;
    }
};

/**
 * @brief A filter class for object not belonging to a container
 *
 * This class instances filter objects that do not belong to a
 * given container.
 *
 * @tparam OBJECT is the type of the objects in the set to be filtered
 */
template<typename OBJECT>
class FilterNotIn : public BaseFilter<OBJECT>
{
    std::set<OBJECT> not_filtered;  //!< Set of the object that should not be filtered
public:

    /**
     * @brief A constructor
     *
     * @tparam CONTAINER is the type of the container storing the object that should not be filtered
     * @param container is the container storing the object that should not be filtered
     */
    template<typename CONTAINER, std::enable_if_t<std::is_same_v<OBJECT, typename CONTAINER::value_type>, bool> = true>
    FilterNotIn(const CONTAINER& container):
        BaseFilter<OBJECT>()
    {
        for (const auto& value : container) {
            not_filtered.insert(value);
        }
    }

    /**
     * @brief Establish whether an object must be filtered by the set
     *
     * @param object is the object to be tested
     * @return `true` if and only if `object` must be filtered
     */
    inline bool filtered(const OBJECT& object) const
    {
        return not_filtered.count(object)==0;
    }
};

}   // CLONES

#endif // __CLONES_FILTER__