from __future__ import print_function

import numpy as np
import glob      # Imported for handling command line arguments.
import sys       # Imported for file directory listing capabilities.
import os        # Imported for file path manipulation capabilities.

from addict import Dict as aDict
from IPython.core.magic import Magics, magics_class, line_magic
from IPython.display import display, HTML

DEBUG = False

#===============================================================================
# Load all plugins
#===============================================================================
def _loadPlugins():
    """
    load existing plugins
    """
    # find the path of the lcslib module
    import masai
    directory = os.path.join(os.path.dirname(masai.__file__),'plugins')

    plugins = glob.glob(os.path.join(directory,"*","*_plugin.py"))

    dic = aDict()

    plugins.sort()

    for file in plugins:

        # Extract just the first part of the .py file name.
        path, basename = os.path.split(file)
        name = basename.split('_plugin')[0]

        if DEBUG:
            print ("Loading plugin: " + name)

        # Dynamically set the PYTHONPATH so the user doesn't have to. It assumes
        # the plugins are contained in subdirectories where the main file lives.
        sys.path.append(path)
        # Import the plugin module temporarily long enough to instantiate an object
        # which is stored in a globally accessible dictionary.

        tempModule = __import__(name+"_plugin")
        attrs = dir(tempModule)
        if '__plugin_main__' in attrs:
            attr = getattr(tempModule,'__plugin_main__')
            cls = getattr(tempModule, attr)
            if cls:
                dic[name.capitalize()] = cls()

    return dic

#===============================================================================
# List available plugins
#===============================================================================
def _availablePlugins(dictionary, *types):
    """
    list existing plugins and print a description of the available plugins
    """
    def _get(_type=None):
        stype = ""
        if _type:
            stype = " of type %s"%_type
        text = "<h6>Available plugins%s:</h6>"%stype
        text += "<ul>"
        keys = dictionary.keys()
        keys.sort()
        for i in keys:
            plugin = dictionary[i]
            if not _type or plugin.type in _type:
                text += "<li>"+plugin.description
        text += "</ul>"
        return text

    if types:
        text = ""
        for item in types:
            text += _get(item)
    else:
        text = _get()

    return HTML(text)


################################################################################
# Magic functions
################################################################################
@magics_class
class NMRMagics(Magics):
    
    plugs = {}
    
    def __init__(self, *args, **kargs):
        ip = kargs.get('shell')
        self.plugins = _loadPlugins()
        ip.user_ns.update(self.plugins)
        #print("Variables in the user namespace:", list(ip.user_ns.keys()))
        display(HTML("<p>ipython nmr tools loaded</p>"))
        super(NMRMagics, self).__init__(*args, **kargs)
        
    @line_magic
    def list_plugins(self, line):
        """
        list plugins
        usage:
            %list_plugins             # to list all plugins
            %list_plugins -i <type>   # to list only plugin of type <type>
            %list_plugins -i (<type1>,<type2>) # to list several type of plugins

        To date, available types are [sources, plots, processes, analysis, simulations]
        """
        ip = self.shell
        opts, argsl = self.parse_options(line, 'i:', mode='string')
        cat = opts.get('i')
        if not isinstance(cat, list):
            cat = [cat,]
        display(_availablePlugins(self.plugins, *cat))

def load_ipython_extension(ipython):
    ipython.register_magics(NMRMagics)


if __name__ == '__main__':
    import IPython
    IPython.start_ipython(argv=[])

