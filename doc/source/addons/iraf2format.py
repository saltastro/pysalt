################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
# Redistribution and use in source and binary forms, with or without       #
# modification, are permitted provided that the following conditions       #
# are met:                                                                 #
#                                                                          #
#     * Redistributions of source code must retain the above copyright     #
#       notice, this list of conditions and the following disclaimer.      #
#     * Redistributions in binary form must reproduce the above copyright  #
#       notice, this list of conditions and the following disclaimer       #
#       in the documentation and/or other materials provided with the      #
#       distribution.                                                      #
#     * Neither the name of the South African Astronomical Observatory     #
#       (SAAO) nor the names of its contributors may be used to endorse    #
#       or promote products derived from this software without specific    #
#       prior written permission.                                          #
#                                                                          #
# THIS SOFTWARE IS PROVIDED BY THE SAAO ''AS IS'' AND ANY EXPRESS OR       #
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED           #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE   #
# DISCLAIMED. IN NO EVENT SHALL THE SAAO BE LIABLE FOR ANY                 #
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL       #
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS  #
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    #
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE          #
# POSSIBILITY OF SUCH DAMAGE.                                              #
############################################################################

#! /usr/bin/env python

"""Tool to convert IRAF style help files to ReStructuredText or XML format."""

from __future__ import with_statement

import re
import xml.dom.minidom

class IrafHelpParser:
    """Parser used to convert IRAF style help files to ReStructuredText or XML format."""

    def add_element(self,name,attribute=None):
        """Add element to DOM tree."""

        element=self.dom.createElement(name)

        if attribute is not None:
            element.setAttribute('name',str(attribute))

        self.current_node=self.current_node.appendChild(element)

    def add_reference(self,key=None):
        self.add_element('reference',key)

    def add_section(self,name=None):
        """Add section to DOM tree."""

        if self.in_section:
            self.move_up()

        self.add_element('section',name)

        self.in_section=True

    def move_up(self):
        """Move up one level in the DOM tree."""

        p=self.current_node.tagName

        self.current_node=self.current_node.parentNode

    def add_text(self,s):
        """Add a text node to the DOM tree."""

        self.current_node.appendChild(self.dom.createTextNode(s))

    def __init__(self,file):
        """Read file and parse markup into DOM."""

        # Read file into list and strip whitespace
        with open(file) as f:
            lines=[l.strip() for l in f.readlines()]

        # Generate a empty DOM and set node to root element
        implement=xml.dom.minidom.getDOMImplementation()

        self.dom=implement.createDocument(None, "document", None)

        self.current_node=self.dom.documentElement

        # Helper variables
        self.in_document=False
        self.in_section=False

        while True:
            try:
                l=lines.pop(0)

                if re.match('\.help',l):
                    self.in_document=True

                    name=re.split('\.help',l)[1].split()[0]
                    self.add_element('help',name)

                elif re.match('\.ih',l):
                    name=lines.pop(0)
                    self.add_section(name.capitalize())

                    # Check for references
                    if re.match('SEE ALSO',name.strip()):
                        references=lines.pop(0)
                        for r in references.split():
                            self.add_reference(r.strip(',').strip())

                elif re.match('\.nf',l):
                    self.add_element('source')

                elif re.match('\.fi',l):
                    self.move_up()

                elif re.match('\.ls',l):
                    name=re.split('\.ls',l)[1]
                    self.add_element('parameter',name)

                elif re.match('\.le',l):
                    self.move_up()

                elif re.match('\.ju',l):
                    pass

                elif re.match('\.re',l):
                    pass

                elif re.match('\.sh',l):
                    pass

                elif re.match('\.br',l):
                    pass

                elif re.match('\.ce',l):
                    pass

                elif re.match('\.sp',l):
                    pass

                elif re.match('\.in',l):
                    pass

                elif re.match('\.endhelp',l):
                    break

                elif self.in_document:
                    self.add_text(l)

            except IndexError:
                break

        # Clean some redundant whitespace nodes
        self.clean_tree()

    def clean_tree(self):
        """Removes redundant blank lines."""
        source_list=self.dom.getElementsByTagName('source')

        # List of nodes to be deleted
        delete=[]

        def isempty(node):
            if node is None:
                return False
            elif node.nodeType==node.TEXT_NODE and node.data.strip()=='':
                return True
            else:
                return False

        for node in source_list:
            # Travel up the tree
            sibling=node.previousSibling
            while isempty(sibling):
                delete.append(sibling)
                sibling=sibling.previousSibling

            # Travel down the tree
            sibling=node.nextSibling
            while isempty(sibling):
                delete.append(sibling)
                sibling=sibling.nextSibling

        # Delete all marked nodes
        for node in delete:
            node.parentNode.removeChild(node)

    def element_to_rst(self,node,indent=0):
        """Recursively move trough all child nodes of *node* and output ReStructuredText."""

        rst=""

        # Check if we should append a newline becouse previous block was a
        # source example
        if node.previousSibling is not None and node.previousSibling.nodeName=='source':
            rst+='\n'

        if node.nodeType == node.TEXT_NODE:
            if node.parentNode.nodeName=='source':
                rst+=' '*indent+node.data+'\n'
            else:
                rst+=' '*indent+node.data.replace('"','``')+'\n'
        elif node.nodeName == 'help':
            name=node.getAttribute('name')
            rst+='.. _'+name.strip()+':'+'\n\n'
            rst+='*'*len(name)+'\n'
            rst+=name+'\n'
            rst+='*'*len(name)+'\n\n'
        elif node.nodeName == 'section':
            name=node.getAttribute('name')
            rst+='\n'
            rst+=name+'\n'
            rst+='='*len(name)+'\n\n'
        elif node.nodeName == 'parameter':
            name=node.getAttribute('name')
            rst+='\n'+'*'+name.strip()+'*'+'\n'
            indent+=4
        elif node.nodeName == 'source':
            indent+=4
        elif node.nodeName == 'reference':
            name=node.getAttribute('name')
            rst+=' :ref:`'+name.strip()+'`'

        # Check if we should append :: because next non empty block is
        # a source code example.
        if node.nextSibling is not None and node.nextSibling.nodeName=='source':
            rst=rst.rstrip().rstrip('.').rstrip(':')+'::'+'\n\n'

        # Recurse through child nodes if any
        if node.hasChildNodes():
            for e in node.childNodes:
                rst+=self.element_to_rst(e,indent)

        return rst

    def get_xml(self):
        """Return DOM tree as well formed XML."""

        return self.dom.toprettyxml(indent='   ', newl='\n')

    def get_rst(self):
        """Return ReStructuredText for entire DOM tree."""

        rst=""

        # Get root element
        rst+=self.element_to_rst(self.dom)

        return rst

def help2ReStructuredText(file):
    """Converts IRAF help file to ReStructuredText format."""

    parser=IrafHelpParser(file)

    return parser.get_rst()

def help2Xml(file):
    """Converts IRAF help file to XML format."""

    parser=IrafHelpParser(file)

    return parser.get_xml()

if __name__ == "__main__":
    import sys
    print help2ReStructuredText(str(sys.argv[1]))
