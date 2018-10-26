#pylint: disable=missing-docstring
import uuid
import collections
import anytree
import MooseDocs
from MooseDocs.base import components
from MooseDocs.extensions import core
from MooseDocs.tree import tokens, html
from MooseDocs.tree.base import Property

def make_extension(**kwargs):
    return FloatExtension(**kwargs)

def create_float(parent, extension, reader, page, settings, **kwargs):
    """Helper for optionally creating a float based on the existence of caption and/or id."""
    cap = add_caption(None, extension, reader, page, settings)
    if cap:
        flt = Float(parent, **kwargs)
        cap.parent = flt
        return flt
    return parent

def caption_settings():
    """Return settings necessary for captions."""
    settings = dict()
    settings['caption'] = (None, "The caption text for the float object.")
    settings['prefix'] = (None, "The numbered caption label to include prior to the caption text.")
    return settings

def add_caption(parent, extension, reader, page, settings):
    """Helper for adding captions to float tokens."""
    cap = settings['caption']
    key = settings['id']
    caption = None
    if key:
        if settings['prefix'] is not None:
            prefix = settings['prefix']
        else:
            prefix = extension.get('prefix', None)
        caption = Caption(parent, key=key, prefix=prefix)
        if cap:
            reader.tokenize(caption, cap, page, MooseDocs.INLINE)
    elif cap:
        caption = Caption(parent)
        reader.tokenize(caption, cap, page, MooseDocs.INLINE)

    #key = caption.get('key', None) if caption else None
    #if key:
    #    tokens.Shortcut(parent, caption=key, link=u'#{}'.format(key),
    #                    string=u'{} {}'.format(caption['prefix'].title(), caption['number']))

    return caption

def create_modal(parent, title=None, content=None, **kwargs):
    """
    Create the necessary Modal tokens for creating modal windows with materialize.
    """
    modal = ModalLink(parent, **kwargs)
    if isinstance(title, unicode):
        ModalLinkTitle(modal, string=title)
    elif isinstance(title, tokens.Token):
        title.parent = ModalLinkTitle(modal)

    if isinstance(content, unicode):
        ModalLinkContent(modal, string=content)
    elif isinstance(content, tokens.Token):
        content.parent = ModalLinkContent(modal)

    return parent

def create_modal_link(parent, title=None, content=None, string=None, **kwargs):
    """
    Create the necessary tokens to create a link to a modal window with materialize.
    """
    kwargs.setdefault('bookmark', unicode(uuid.uuid4()))
    link = tokens.Link(parent,
                       url=u'#{}'.format(kwargs['bookmark']),
                       class_='modal-trigger',
                       string=string)
    create_modal(parent, title, content, **kwargs)
    return link

Float = tokens.newToken('Float', img=False)
Caption = tokens.newToken('Caption', key=u'', prefix=u'', number=1)
ModalLink = tokens.newToken('ModalLink', bookmark=True, bottom=False, close=True)
ModalLinkTitle = tokens.newToken('ModalLinkTitle')
ModalLinkContent = tokens.newToken('ModalLinkContent')

class FloatExtension(components.Extension):
    """
    Provides ability to add caption float elements (e.g., figures, table, etc.). This is only a
    base extension. It does not provide tables for example, just the tools to make floats
    in a uniform manner.
    """
    COUNTS = collections.defaultdict(int)
    def extend(self, reader, renderer):
        renderer.add('Float', RenderFloat())
        renderer.add('Caption', RenderCaption())
        renderer.add('ModalLink', RenderModalLink())
        renderer.add('ModalLinkTitle', RenderModalLinkTitle())
        renderer.add('ModalLinkContent', RenderModalLinkContent())

    def preTokenize(self, ast, page):
        """Reset float counters."""
        FloatExtension.COUNTS.clear()

    def postTokenize(self, ast, page):
        """Set float number for each counter."""
        for node in anytree.PreOrderIter(ast, filter_=lambda n: n.name == 'CountToken'):
            if node.prefix is not None:
                FloatExtension.COUNTS[node.prefix] += 1
                node.set('number', FloatExtension.COUNTS[node.prefix])

class RenderFloat(components.RenderComponent):
    def createHTML(self, parent, token, page): #pylint: disable=no-self-use
        div = html.Tag(parent, 'div', token)
        div.addClass('moose-float-div')
        return div

    def createMaterialize(self, parent, token, page): #pylint: disable=no-self-use
        div = html.Tag(parent, 'div', token)
        div.addClass('card')
        content = html.Tag(div, 'div')
        if token['img']:
            content.addClass('card-image')
        else:
            content.addClass('card-content')

        return content

    def createLatex(self, parent, token, page):
        pass

class RenderCaption(components.RenderComponent):
    def createHTML(self, parent, token, page): #pylint: disable=no-self-use
        caption = html.Tag(parent, 'p', class_="moose-caption")
        prefix = token.get('prefix', None)
        if prefix:
            heading = html.Tag(caption, 'span', class_="moose-caption-heading")
            html.String(heading, content=u"{} {}: ".format(prefix, token['number']))

        text = html.Tag(caption, 'span', class_="moose-caption-text")
        return text

    def createLatex(self, parent, token, page):
        pass

class RenderModalLink(core.RenderLink):
    def createMaterialize(self, parent, token, page):

        cls = "modal bottom-sheet" if token.bottom else "modal"
        modal = html.Tag(parent, 'div', class_=cls, id_=token.bookmark)
        modal_content = html.Tag(modal, 'div', class_="modal-content")

        if token.close:
            footer = html.Tag(modal, 'div', class_='modal-footer')
            html.Tag(footer, 'a', class_='modal-close btn-flat', string=u'Close')
        return modal_content

class RenderModalLinkTitle(components.RenderComponent):

    def createMaterialize(self, parent, token, page):
        return html.Tag(parent, 'h4')

class RenderModalLinkContent(components.RenderComponent):
    def createMaterialize(self, parent, token, page):
        return parent
