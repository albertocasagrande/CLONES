/**
 * @file mutation_list.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to represent mutation lists
 * @version 1.2
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

#include <sstream>

#include "mutation_list.hpp"

namespace CLONES
{

namespace Mutations
{

MutationList::const_iterator
MutationList::const_iterator::get_begin(const MutationList& mutations)
{
    const_iterator begin_it;

    begin_it.SID_it = mutations.SIDs.begin();
    begin_it.CNA_it = mutations.CNAs.begin();
    begin_it.order_it = mutations.application_order.begin();

    return begin_it;
}

MutationList::const_iterator
MutationList::const_iterator::get_end(const MutationList& mutations)
{
    const_iterator end_it;

    end_it.SID_it = mutations.SIDs.end();
    end_it.CNA_it = mutations.CNAs.end();
    end_it.order_it = mutations.application_order.end();

    return end_it;
}

MutationList::const_iterator::const_iterator()
{}

MutationList::const_iterator&
MutationList::const_iterator::operator++()
{
    switch(*order_it) {
        case SID_TURN:
            ++SID_it;
            break;
        case CNA_TURN:
            ++CNA_it;
            break;
        default:
            break;
    }

    ++order_it;

    return *this;
}

MutationList::const_iterator&
MutationList::const_iterator::operator--()
{
    --order_it;

    switch(*order_it) {
        case SID_TURN:
            --SID_it;
            break;
        case CNA_TURN:
            --CNA_it;
            break;
        default:
            break;
    }

    return *this;
}

MutationList::MutationList()
{}

//! @private
template<typename MUTATION>
std::set<MUTATION> build_mutation_set(const std::string& mutant_name,
                                      const std::list<MUTATION>& mutations)
{
    std::set<MUTATION> mutation_set;

    for (const auto& mutation : mutations) {
        if (mutation_set.count(mutation)>0) {
            std::ostringstream oss;

            oss << mutation << " added twice among "
                << mutant_name << "'s driver mutations.";

            throw std::domain_error(oss.str());
        }

        mutation_set.insert(mutation);
    }

    return mutation_set;
}

std::list<DriverMutations::MutationType>
MutationList::get_default_order(const std::list<MutationSpec<SID>>& SIDs,
                                const std::list<CNA>& CNAs,
                                const bool& wg_doubling)
{
    std::list<DriverMutations::MutationType> application_order;

    for (size_t i=0; i<SIDs.size(); ++i) {
        application_order.push_back(DriverMutations::MutationType::SID_TURN);
    }

    for (size_t i=0; i<CNAs.size(); ++i) {
        application_order.push_back(DriverMutations::MutationType::CNA_TURN);
    }

    if (wg_doubling) {
        application_order.push_back(DriverMutations::MutationType::WGD_TURN);
    }

    return application_order;
}


MutationList::MutationList(const std::list<MutationSpec<SID>>& SIDs,
                           const std::list<CNA>& CNAs,
                           const bool& wg_doubling):
    SIDs{SIDs}, CNAs{CNAs},
    application_order{DriverMutations::get_default_order(SIDs,CNAs,wg_doubling)}
{}

MutationList::MutationList(const std::list<MutationSpec<SID>>& SIDs,
                           const std::list<CNA>& CNAs,
                           const std::list<MutationType>& application_order):
    SIDs{SIDs}, CNAs{CNAs}, application_order{application_order}
{
    size_t SID_count{0}, CNA_count{0};
    for (auto order_it = application_order.begin();
         order_it != application_order.end(); ++order_it) {

        switch(*order_it) {
            case SID_TURN:
                ++SID_count;
                break;
            case CNA_TURN:
                ++CNA_count;
                break;
            case WGD_TURN:
                break;
            default:
                throw std::domain_error("Unsupported driver mutation type.");
        }
    }

    if (SID_count != SIDs.size()) {
        throw std::domain_error("The number of SNVs/indels differs from "
                                "that of the same kind of mutations "
                                "in the application order list.");
    }

    if (CNA_count != CNAs.size()) {
        throw std::domain_error("The number of CNAs differs from "
                                "that of the same kind of mutations "
                                "in the application order list.");
    }
}

MutationList& MutationList::insert(const CNA& cna)
{
    CNAs.push_back(cna);
    application_order.push_back(CNA_TURN);

    return *this;
}

MutationList& MutationList::insert(const MutationSpec<SID>& sid)
{
    SIDs.push_back(sid);
    application_order.push_back(SID_TURN);

    return *this;
}

MutationList& MutationList::insert_WGD()
{
    application_order.push_back(WGD_TURN);

    return *this;
}

DriverMutations::DriverMutations()
{}

DriverMutations::DriverMutations(const std::string& mutant_name,
                                 const std::list<MutationSpec<SID>>& SIDs,
                                 const std::list<CNA>& CNAs,
                                 const bool& wg_doubling):
    MutationList(SIDs, CNAs, wg_doubling), name{mutant_name}
{}

DriverMutations::DriverMutations(const std::string& mutant_name,
                                 const std::list<MutationSpec<SID>>& SIDs,
                                 const std::list<CNA>& CNAs,
                                 const std::list<MutationType>& application_order):
    MutationList(SIDs, CNAs, application_order), name{mutant_name}
{}

}   // Mutations

}   // CLONES
