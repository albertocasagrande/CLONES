/**
 * @file mutation_list.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to represent mutation lists
 * @version 1.3
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

#ifndef __CLONES_MUTATION_LIST__
#define __CLONES_MUTATION_LIST__

#include <list>
#include <string>

#include "simulation.hpp"
#include "sid.hpp"
#include "cna.hpp"
#include "mutation_spec.hpp"

namespace CLONES
{

namespace Mutations
{

/**
 * @brief The genomic characterization of a mutant
 */
struct MutationList
{
    /**
     * @brief Define a choice among mutation type
     */
    enum MutationType {
        SID_TURN,   //!< A SID mutation must be applied
        CNA_TURN,   //!< A CNA mutation must be applied
        WGD_TURN    //!< A whole genome doubling must be applied
    };

    std::list<MutationSpec<SID>> SIDs;  //!< The mutant SIDs
    std::list<CNA> CNAs;                //!< The mutant CNAs

    std::list<MutationType> application_order;  //!< The mutation application order

    /**
     * @brief A constant iterator for the mutation list class
     */
    struct const_iterator {

    protected:
        std::list<MutationSpec<SID>>::const_iterator SID_it;    //!< The SID list iterator
        std::list<CNA>::const_iterator CNA_it;                  //!< The CNA list iterator

        std::list<MutationType>::const_iterator order_it;       //!< The order list iterator

        /**
         * @brief Get the iterator referring to the mutation list begin
         *
         * @param mutations is a mutation list
         * @return the constant iterator referring to the first element in `mutations`
         */
        static const_iterator get_begin(const MutationList& mutations);

        /**
         * @brief Get the iterator referring to the mutation list end
         *
         * @param mutations is a mutation list
         * @return the constant iterator referring to the element following
         *   the last element in `mutations`
         */
        static const_iterator get_end(const MutationList& mutations);
    public:
        /**
         * @brief The empty constructor
         */
        const_iterator();

        /**
         * @brief The pre-increment operator
         *
         * @return the updated constant iterator
         */
        const_iterator& operator++();

        /**
         * @brief The pre-decrement operator
         *
         * @return the updated constant iterator
         */
        const_iterator& operator--();

        /**
         * @brief Get the type of the referred element in the list
         *
         * @return the type of the referred element in the list
         */
        inline const MutationType& get_type() const
        {
            return *order_it;
        }

        /**
         * @brief Get the last CNA not following the referred element in the list
         *
         * @return the last CNA not following the referred element in the list
         */
        inline const CNA& get_last_CNA() const
        {
            return *CNA_it;
        }

        /**
         * @brief Get the last SID not following the referred element in the list
         *
         * @return the last SID not following the referred element in the list
         */
        inline const MutationSpec<SID>& get_last_SID() const
        {
            return *SID_it;
        }

        /**
         * @brief Test whether to iterator refer to the same element
         *
         * @param iterator is an iterator
         * @return `true` if and only if the object and `iterator` refer
         *   to the same object
         */
        inline bool operator==(const const_iterator& iterator) const
        {
            return order_it==iterator.order_it;
        }

        /**
         * @brief Test whether to iterator refer to different elements
         *
         * @param iterator is an iterator
         * @return `true` if and only if the object and `iterator` refer
         *   to different objects
         */
        inline bool operator!=(const const_iterator& iterator) const
        {
            return order_it!=iterator.order_it;
        }

        friend struct MutationList;
    };

    /**
     * @brief The empty constructor
     */
    MutationList();

    /**
     * @brief A constructor
     *
     * This method chooses an arbitrary order for the driver mutations:
     * first all the SIDs, then all the CNAs, and, finally, the whole genome
     * doubling (WGD).
     *
     * @param SIDs is the vector of species specific SIDs
     * @param CNAs is the vector of species specific CNAs
     * @param wg_doubling is a Boolean flag to enable whole genome doubling
     */
    MutationList(const std::list<MutationSpec<SID>>& SIDs,
                 const std::list<CNA>& CNAs,
                 const bool& wg_doubling=false);

    /**
     * @brief A constructor
     *
     * This method takes as a parameter a list that defines the order
     * in which SIDs, CNAs, and WGD must be applied. The relative
     * orders of SIDs and CNAs are defined by the respective lists, while
     * the application order list defines the order of the mutation kind.
     *
     * @param SIDs is the vector of species specific SIDs
     * @param CNAs is the vector of species specific CNAs
     * @param application_order is the list of application order
     */
    MutationList(const std::list<MutationSpec<SID>>& SIDs,
                 const std::list<CNA>& CNAs,
                 const std::list<MutationType>& application_order);

    /**
     * @brief Insert a CNA in the list
     *
     * @param cna is the CNA to be inserted
     * @return a reference to the updated object
     */
    MutationList& insert(const CNA& cna);

    /**
     * @brief Insert a SID specification in the list
     *
     * @param sid is the SID specification to be inserted
     * @return a reference to the updated object
     */
    MutationList& insert(const MutationSpec<SID>& sid);

    /**
     * @brief Insert a whole genome doubling in the list
     * @return a reference to the updated object
     */
    MutationList& insert_WGD();

    /**
     * @brief Get the default order for driver mutations
     *
     * @param SIDs is a list of SID mutations
     * @param CNAs is a list CNA mutations
     * @param wg_doubling is a Boolean flag to enable whole genome doubling
     * @return the application order list for mutations that requests to
     *   apply all the SIDs, the CNAs, and, when `wg_doubling` is set to
     *   `true`, the whole genome doubling in this order.
     */
    static std::list<MutationType>
    get_default_order(const std::list<MutationSpec<SID>>& SIDs,
                      const std::list<CNA>& CNAs, const bool& wg_doubling);

    /**
     * @brief Get the iterator referring to the mutation list begin
     *
     * @return the constant iterator to the first element in the list
     */
    inline const_iterator begin() const
    {
        return MutationList::const_iterator::get_begin(*this);
    }

    /**
     * @brief Get the iterator referring to the mutation list end
     *
     * @return the constant iterator to the element following the last
     *   in the list
     */
    inline const_iterator end() const
    {
        return MutationList::const_iterator::get_end(*this);
    }

    /**
     * @brief Save a mutation list in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & SIDs
                & CNAs
                & application_order;
    }

    /**
     * @brief Load novel mutations from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load novel mutations
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static MutationList load(ARCHIVE& archive)
    {
        MutationList mutations;

        archive & mutations.SIDs
                & mutations.CNAs
                & mutations.application_order;

        return mutations;
    }

    /**
     * @brief Get the number of mutations
     * 
     * @return the number of mutations
     */
    inline size_t size() const
    {
        return application_order.size();
    }
};

struct DriverMutations : public MutationList
{
    std::string name;                   //!< The mutant name

    /**
     * @brief The empty constructor
     */
    DriverMutations();

    /**
     * @brief A constructor
     *
     * This method chooses an arbitrary order for the driver mutations:
     * first all the SIDs, then all the CNAs, and, finally, the whole genome
     * doubling (WGD).
     *
     * @param mutant_name is the name of the mutant
     * @param SIDs is the vector of species specific SIDs
     * @param CNAs is the vector of species specific CNAs
     * @param wg_doubling is a Boolean flag to enable whole genome doubling
     */
    DriverMutations(const std::string& mutant_name,
                    const std::list<MutationSpec<SID>>& SIDs,
                    const std::list<CNA>& CNAs,
                    const bool& wg_doubling=false);

    /**
     * @brief A constructor
     *
     * This method takes as a parameter a list that defines the order
     * in which SIDs, CNAs, and WGD must be applied. The relative
     * orders of SIDs and CNAs are defined by the respective lists, while
     * the application order list defines the order of the mutation kind.
     *
     * @param mutant_name is the name of the mutant
     * @param SIDs is the vector of species specific SIDs
     * @param CNAs is the vector of species specific CNAs
     * @param application_order is the list of application order
     */
    DriverMutations(const std::string& mutant_name,
                    const std::list<MutationSpec<SID>>& SIDs,
                    const std::list<CNA>& CNAs,
                    const std::list<MutationType>& application_order);

    /**
     * @brief Save driver mutations
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        static_cast<const MutationList*>(this)->save(archive);
    
        archive & name;
    }

    /**
     * @brief Load driver mutations
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load driver mutations object
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static DriverMutations load(ARCHIVE& archive)
    {
        DriverMutations d_mutations; 
        
        static_cast<MutationList&>(d_mutations) = MutationList::load(archive);
        archive & d_mutations.name;

        return d_mutations;
    }
};

}   // Mutations

}   // CLONES

#endif // __CLONES_MUTATION_LIST__
