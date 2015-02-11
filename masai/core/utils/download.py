# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.utils.download
# ===============================================================================
# Copyright (C) 2012-2015 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
#===============================================================================
"""

Various methods to download distant file

"""
from __future__ import print_function

import os
import sys
import tarfile

from traits.api import (HasTraits, Directory, Instance, Str, Dict)
from traits.etsconfig.api import ETSConfig

from .ssh import SSHSession

from .preferences import Preferences
preferences = Preferences(filename=os.path.join(os.path.expanduser('~'), '.masai','masai.ini'))

class Spectrometer(HasTraits):

    distant = Directory()
    local = Directory()

    session = Instance(SSHSession)
    user = Str("nmr")
    password = Str("")

    excludefile = Str('exclude.txt')
    listingscript = Str('script_listing.py')

    catalog = Dict

    def __init__(self, host=None, user="nmr", password="", local=None):

        self.user = user
        self.password = password

        if local:
            self.local = local
        else:
            #look in the environment variables
            self.local = os.environ.get('NMRDATA',"")
            if not self.local:
                # look into the preferences
                self.local = preferences.get('env.nmrdata',"")
        if not self.local:
            print('Please, set the local directory!. Usage: Spectrometer(host, [user=?, ][password=? , ][local=?]')
        else:
            preferences.set('env.nmrdata',self.local)
            preferences.flush()

        self.session = SSHSession(user, host, password, verbose=True)


    def listing(self, user=None, name=None):

        self.distant = os.path.join('data')

        if user:
            self.distant = os.path.join(self.distant,user,'nmr')
            if name:
                self.distant = os.path.join(self.distant,user,'nmr')
        else:
            print("A 'user' is required in the argument of this function.  Usage: list_dataset('xxx', [name='yyyy'])  ")
            return

        self.session.scp(self.listingscript, 'to', self.listingscript)
        self.session.scp(self.excludefile, 'to', self.excludefile)

        _catalog = self.session.ssh('python script_listing.py %s'%self.distant)

        self.catalog = eval(_catalog)

        for k, v in self.catalog.iteritems():
            print(k, v, sep=': ')

        return

    def get(self, user=None, name=None, expno=None):
        """
        Get and donwnload distant data

        Parameters
        ==========
        user: Str
            user name
        name: Str
            name of the dataset
        expno: Int, Str or List of Str, optional
            expno(s) to download, if not present all expnos are taken

        """
        if not name or not user:
            print()
            print("ERROR:","We need to specify a *user* and a *name* for the dataset", file=sys.stderr)
            return

        path = os.path.join("data", user, "nmr", name)

        if not expno:
            print("We download all expno present in the *%s* dataset"%name)
            expno = self.catalog.get(path)

        if isinstance(expno, str):
            expno = [expno, ]

        expnos = tuple(expno)

        if not self.local:
            print('no local directory was specified', file=sys.stderr)
            return

        for item in expnos:
            p = os.path.join(self.local, path,item)
            if not os.path.exists(p):
                f= os.path.join(path, item)
                print("%s does not exist -> download."%item)
                self.session.ssh("rm tempo.tar.gz; tar cvzf tempo.tar.gz -X exclude.txt '%s'"%(f.strip()))
                self.session.scp(self.local,"from","tempo.tar.gz")
                try:
                    tar = tarfile.open(os.path.join(self.local,"tempo.tar.gz"))
                    tar.extractall(path=self.local)
                    tar.close()
                except:
                    print('Error: bypassed', file=sys.stderr)
            else:
                print("%s already exists. Bypass"%p)

if __name__ == '__main__':

    d = Spectrometer(host="192.168.223.175", local=os.path.expanduser(os.path.join('~','Dropbox','storage','bruker')))
    d.listing("Eddy")
    d.get("Eddy", "150114","7")

