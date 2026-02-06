/**
 * @file phyloXML.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements phyloXML stream
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

#include <ostream>

#include "phyloXML.hpp"

#include "descendant_forest.hpp"
#include "palette.hpp"

namespace CLONES
{

namespace Mutants
{
namespace IO
{

phyloXMLStream::phyloXMLStream(std::ostream& os, const std::string& indentation_symbols):
    os(os), closed(false), indent_level(1), indent_symbols(indentation_symbols)
{
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl
       << "<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
       << "xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" "
       << "xmlns=\"http://www.phyloxml.org\">" << std::endl;
}

void phyloXMLStream::change_indentation_level(const size_t level)
{
    indent_level = level;

    std::ostringstream oss;
    for (size_t i=0; i<indent_level; ++i) {
        oss << indent_symbols;
    }

    indent_string = oss.str();
}

phyloXMLStream& phyloXMLStream::operator<<(const DescendantsForest& forest)
{
    if (closed) {
        throw std::runtime_error("The stream has been already closed");
    }

    const auto roots = forest.get_roots();
    os << indent_string << "<phylogeny rooted=\"true\">" << std::endl;

    change_indentation_level(indent_level+1);

    for (const auto& root: roots) {
        *this << root;
    }

    change_indentation_level(indent_level-1);

    os << indent_string << "</phylogeny>" << std::endl;

    return *this;
}

phyloXMLStream& phyloXMLStream::operator<<(const UI::Color& color)
{
    os << indent_string << "<color>" << std::endl;
    change_indentation_level(indent_level+1);

       os << indent_string << "<red>" << static_cast<unsigned int>(color.red) << "</red>" << std::endl
          << indent_string << "<green>" << static_cast<unsigned int>(color.green) << "</green>" << std::endl
           << indent_string << "<blue>" << static_cast<unsigned int>(color.blue) << "</blue>" << std::endl;

    change_indentation_level(indent_level-1);
    os << indent_string << "</color>" << std::endl;

    return *this;
}


phyloXMLStream& phyloXMLStream::operator<<(const SpeciesId& species_id)
{
    os << indent_string << "<taxonomy>" << std::endl;

    change_indentation_level(indent_level+1);
    os << indent_string << "<id>" << species_id << "</id>" << std::endl;

    change_indentation_level(indent_level-1);
    os  << indent_string << "</taxonomy>" << std::endl;

    return *this;
}

phyloXMLStream& phyloXMLStream::operator<<(const DescendantsForest::const_node& node)
{
    if (closed) {
        throw std::runtime_error("The stream has been already closed");
    }

    const Cell& cell = node;

    if (node.is_root()) {
        // this is the root node
        os << indent_string << "<clade>" << std::endl;
    } else {
        // this is not the root
        const Cell parent = node.parent();

        os << indent_string << "<clade branch_length=\""
           << (cell.get_birth_time()-parent.get_birth_time())  << "\">" << std::endl;
    }

    change_indentation_level(indent_level+1);

    *this << cell.get_species_id() << UI::palette[cell.get_species_id()];

    os << indent_string << "<name>" << cell.get_id() << "</name>" << std::endl;

    const auto children = node.children();

    for (const auto& child : children) {
        *this << child;
    }

    change_indentation_level(indent_level-1);

    os << indent_string << "</clade>" << std::endl;

    return *this;
}

void phyloXMLStream::close()
{
    if (!closed) {
        os << "</phyloxml>" << std::endl;

        closed = true;
    }
}

phyloXMLStream::~phyloXMLStream()
{
    close();
}

}  // IO

}  // Mutants

}  // Race
