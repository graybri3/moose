"""
Contains base classes intended to be used internal to this module.
"""
import uuid
import copy
import MooseDocs
from MooseDocs import common
from MooseDocs.common import exceptions

#: A value for allowing ConfigObject.get method to work with a default of None
UNSET = uuid.uuid4()

class ConfigObject(object):
    """
    Base class for objects that contain configure options.
    """
    @staticmethod
    def defaultConfig():
        """
        Return the default configuration for this object. The configuration should be a dict of
        tuples, where the tuple contains the default value and the description.
        """
        return dict()

    def __init__(self, **kwargs):
        self.__config = self.defaultConfig()
        if not isinstance(self.__config, dict):
            msg = "The return type from 'defaultConfig' must be a 'dict', but a {} was provided."
            raise exceptions.MooseDocsException(msg, type(self.__config))
        self.update(**kwargs)

        # Stores the configuration that was established upon object creation, this allows the
        # config extension to mess with configuration items during execute but then restore
        # prior to processing the next page.
        self.__initial_config = copy.copy(self.__config)

    def update(self, **kwargs):
        """
        Update the configuration with the supplied key-value pairs.
        """
        unknown = []
        for key, value in kwargs.iteritems():

            if key not in self.__config:
                unknown.append(key)
            else:
                self.__config[key] = (value, self.__config[key][1]) #TODO: type check???

        if unknown:
            msg = "The following config options were not found in the default config options for " \
                  "the {} object:"
            for key in unknown:
                msg += '\n{}{}'.format(' '*4, key)
            raise exceptions.MooseDocsException(msg.format(type(self)))

    def resetConfig(self):
        """
        Reset configuration to original state.
        """
        self.__config = copy.copy(self.__initial_config)

    def keys(self):
        """
        Return the available configuration items.
        """
        return self.__config.keys()

    def __getitem__(self, name):
        """
        Return a configuration value by name using the [] operator.
        """
        return self.get(name)

    def get(self, name, default=UNSET):
        """
        Return a configuration value by name, with an optional default.
        """
        if (default is not UNSET) and (name not in self.__config):
            return default
        else:
            return self.__config[name][0]

class ReaderObject(object):
    """
    Basic functions for objects that have a Reader object.
    """
    def __init__(self):
        self.__reader = None

    def init(self, reader):
        """Initialize the class with the Reader object."""
        if self.initialized():
            msg = "The {} object has already been initialized, this method should not " \
                  "be called twice."
            raise MooseDocs.common.exceptions.MooseDocsException(msg, type(self))

        common.check_type('reader', reader, MooseDocs.base.readers.Reader)
        self.__reader = reader

    def initialized(self):
        """Returns True if the init method was called."""
        return self.__reader is not None

    @property
    def reader(self):
        """Return the Reader object."""
        if self.__reader is None:
            msg = "The init() method of the {} object must be called prior to accessing the " \
                  "reader property."
            raise MooseDocs.common.exceptions.MooseDocsException(msg, type(self))

        return self.__reader

class RendererObject(object):
    """
    Basic functions for objects that have a Renderer object.
    """
    def __init__(self):
        self.__renderer = None

    def init(self, renderer):
        """Initialize the class with the Renderer object."""
        if self.initialized():
            msg = "The {} object has already been initialized, this method should not " \
                  "be called twice."
            raise MooseDocs.common.exceptions.MooseDocsException(msg, type(self))

        common.check_type('renderer', renderer, MooseDocs.base.renderers.Renderer)
        self.__renderer = renderer

    def initialized(self):
        """Returns True if the init method was called."""
        return self.__renderer is not None

    @property
    def renderer(self):
        """Return the Renderer object."""
        if self.__renderer is None:
            msg = "The init() method of the {} object must be called prior to accessing this " \
                  "property."
            raise MooseDocs.common.exceptions.MooseDocsException(msg, type(self))

        return self.__renderer

class ComponentObject(object):
    """
    Class for objects that require a list of components (e.g., Reader and Renderers).
    """
    def __init__(self):
        self.__components = []

    @property
    def components(self):
        """
        Return the list of Component objects.
        """
        return self.__components

    def addComponent(self, comp):
        """
        Add a Component object.
        """
        common.check_type("component", comp, MooseDocs.base.components.Component)
        self.__components.append(comp)
