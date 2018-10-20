#pylint: disable=missing-docstring
#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html
#pylint: enable=missing-docstring
import cgi
import anytree
from base import NodeBase, Property

class Tag(NodeBase):
    """
    A node representing an HTML tag (e.g., h1, strong).
    """
    def __init__(self, parent=None, name=None, **kwargs):
        kwargs.setdefault('close', True)
        kwargs.setdefault('string', None)
        kwargs['class'] = kwargs.pop('class_', u'')
        kwargs['style'] = kwargs.pop('style_', u'')
        kwargs['id'] = kwargs.pop('id_', u'')
        super(Tag, self).__init__(name=name, parent=parent, **kwargs)

        string = self.get('string', None)
        if string is not None:
            String(self, content=string)

    def addStyle(self, style):
        """
        Add to the existing style settings.
        """
        s = self.get('style', '').split(';')
        s += style.split(';')
        self['style'] = ';'.join(s)

    def addClass(self, *args):
        """
        Add to the existing class list.
        """
        c = self.get('class', '').split(' ')
        c += args
        self.set('class', ' '.join(c))

    def write(self):
        """Write the HTML as a string, e.g., <foo>...</foo>."""
        out = ''
        attr_list = []
        for key, value in self.iteritems():
            if value:# and (key != 'class'):
                attr_list.append('{}="{}"'.format(key, str(value).strip()))

        attr = ' '.join(attr_list)
        if attr:
            out += '<{} {}>'.format(self.name, attr)
        else:
            out += '<{}>'.format(self.name)

        for child in self.children:
            out += child.write()
        if self.get('close'): #pylint: disable=no-member
            out += '</{}>'.format(self.name)
        return out

    def text(self):
        """
        Convert String objects into a single string.
        """
        strings = []
        for node in anytree.PreOrderIter(self):
            if node.name == 'String':
                strings.append(node['content'])
        return u' '.join(strings)

class String(NodeBase):
    """
    A node for containing string content.
    """
    def __init__(self, parent=None, **kwargs):
        kwargs.setdefault('content', u'')
        kwargs.setdefault('escape', u'')
        super(String, self).__init__('String', parent, **kwargs)

    def write(self):
        if self.get('escape'):
            return cgi.escape(self.get('content'), quote=True)
        else:
            return self.get('content')
