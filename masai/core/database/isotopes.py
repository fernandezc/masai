# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.database.isotopes
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
# Laboratoire Catalyse et Spectrochimie, Caen, France.
# christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================

"""This module essentially define a class :class:`Isotopes` which handle all 
necessary features of NMR nuclei, such as their spin, larmor frequency and so on.

"""
#===============================================================================
# python imports.
#===============================================================================

#===============================================================================
# Enthought imports.
#===============================================================================
from traits.api import (HasTraits, Property, Str, Float, Enum, Instance)

#===============================================================================
# Numpy 
#===============================================================================
from numpy import nan

#: Dictionary containing the isotopes properties
dict_isotopes = {
    '1H': {'A': 1, 'abundance': 99.9885, 'name': 'hydrogen', 'symbol': 'H', 'stability': 'stable', 'quadrupole': 0,
           'gn': 5.5856912, 'Z': 1, 'spin': 0.5},
    '2H': {'A': 2, 'abundance': 0.0115, 'name': 'hydrogen', 'symbol': 'H', 'stability': 'stable', 'quadrupole': 0.00286,
           'gn': 0.8574376, 'Z': 1, 'spin': 1.0},
    '3H': {'A': 3, 'abundance': 0.0, 'name': 'hydrogen', 'symbol': 'H', 'stability': 'radioactive', 'quadrupole': 0,
           'gn': 5.95792, 'Z': 1, 'spin': 0.5},
    '3He': {'A': 3, 'abundance': 0.000137, 'name': 'helium', 'symbol': 'He', 'stability': 'stable', 'quadrupole': 0,
            'gn': -4.255248, 'Z': 2, 'spin': 0.5},
    '4He': {'A': 4, 'abundance': 99.999863, 'name': 'helium', 'symbol': 'He', 'stability': 'stable', 'quadrupole': 0,
            'gn': 0.0, 'Z': 2, 'spin': 0.0},
    '6Li': {'A': 6, 'abundance': 7.59, 'name': 'lithium', 'symbol': 'Li', 'stability': 'stable',
            'quadrupole': -0.000808,
            'gn': 0.8220514, 'Z': 3, 'spin': 1.0},
    '7Li': {'A': 7, 'abundance': 92.41, 'name': 'lithium', 'symbol': 'Li', 'stability': 'stable', 'quadrupole': -0.0401,
            'gn': 2.170961, 'Z': 3, 'spin': 1.5},
    '9Be': {'A': 9, 'abundance': 100.0, 'name': 'beryllium', 'symbol': 'Be', 'stability': 'stable',
            'quadrupole': 0.05288,
            'gn': -0.785, 'Z': 4, 'spin': 1.5},
    '10B': {'A': 10, 'abundance': 19.9, 'name': 'boron', 'symbol': 'B', 'stability': 'stable', 'quadrupole': 0.08459,
            'gn': 0.600216, 'Z': 5, 'spin': 3.0},
    '11B': {'A': 11, 'abundance': 80.1, 'name': 'boron', 'symbol': 'B', 'stability': 'stable', 'quadrupole': 0.04059,
            'gn': 1.792424, 'Z': 5, 'spin': 1.5},
    '12C': {'A': 12, 'abundance': 98.93, 'name': 'carbon', 'symbol': 'C', 'stability': 'stable', 'quadrupole': 0,
            'gn': 0.0, 'Z': 6, 'spin': 0.0},
    '13C': {'A': 13, 'abundance': 1.07, 'name': 'carbon', 'symbol': 'C', 'stability': 'stable', 'quadrupole': 0,
            'gn': 1.40482, 'Z': 6, 'spin': 0.5},
    '14C': {'A': 14, 'abundance': 0.0, 'name': 'carbon', 'symbol': 'C', 'stability': 'radioactive', 'quadrupole': nan,
            'gn': 0.273, 'Z': 6, 'spin': 3.0},
    '14N': {'A': 14, 'abundance': 99.632, 'name': 'nitrogen', 'symbol': 'N', 'stability': 'stable',
            'quadrupole': 0.02044,
            'gn': 0.4037607, 'Z': 7, 'spin': 1.0},
    '15N': {'A': 15, 'abundance': 0.368, 'name': 'nitrogen', 'symbol': 'N', 'stability': 'stable', 'quadrupole': 0,
            'gn': -0.5663784, 'Z': 7, 'spin': 0.5},
    '16O': {'A': 16, 'abundance': 99.757, 'name': 'oxygen', 'symbol': 'O', 'stability': 'stable', 'quadrupole': 0,
            'gn': 0.0, 'Z': 8, 'spin': 0.0},
    '17O': {'A': 17, 'abundance': 0.038, 'name': 'oxygen', 'symbol': 'O', 'stability': 'stable', 'quadrupole': -0.02558,
            'gn': -0.757516, 'Z': 8, 'spin': 2.5},
    '18O': {'A': 18, 'abundance': 0.205, 'name': 'oxygen', 'symbol': 'O', 'stability': 'radioactive', 'quadrupole': 0,
            'gn': 0.0, 'Z': 8, 'spin': 0.0},
    '19F': {'A': 19, 'abundance': 100.0, 'name': 'fluorine', 'symbol': 'F', 'stability': 'stable', 'quadrupole': 0,
            'gn': 5.257732, 'Z': 9, 'spin': 0.5},
    '20Ne': {'A': 20, 'abundance': 90.48, 'name': 'neon', 'symbol': 'Ne', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 10, 'spin': 0.0},
    '21Ne': {'A': 21, 'abundance': 0.27, 'name': 'neon', 'symbol': 'Ne', 'stability': 'stable', 'quadrupole': 0.10155,
             'gn': -0.441197, 'Z': 10, 'spin': 1.5},
    '22Ne': {'A': 22, 'abundance': 9.25, 'name': 'neon', 'symbol': 'Ne', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 10, 'spin': 0.0},
    '22Na': {'A': 22, 'abundance': 0.0, 'name': 'sodium', 'symbol': 'Na', 'stability': 'radioactive', 'quadrupole': nan,
             'gn': 0.582, 'Z': 11, 'spin': 3.0},
    '23Na': {'A': 23, 'abundance': 100.0, 'name': 'sodium', 'symbol': 'Na', 'stability': 'stable', 'quadrupole': 0.104,
             'gn': 1.478391, 'Z': 11, 'spin': 1.5},
    '24Mg': {'A': 24, 'abundance': 78.99, 'name': 'magnesium', 'symbol': 'Mg', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 12, 'spin': 0.0},
    '25Mg': {'A': 25, 'abundance': 10.0, 'name': 'magnesium', 'symbol': 'Mg', 'stability': 'stable',
             'quadrupole': 0.1994,
             'gn': -0.34218, 'Z': 12, 'spin': 2.5},
    '26Mg': {'A': 26, 'abundance': 11.01, 'name': 'magnesium', 'symbol': 'Mg', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 12, 'spin': 0.0},
    '27Al': {'A': 27, 'abundance': 100.0, 'name': 'aluminium', 'symbol': 'Al', 'stability': 'stable',
             'quadrupole': 0.1466,
             'gn': 1.456601, 'Z': 13, 'spin': 2.5},
    '28Si': {'A': 28, 'abundance': 92.2297, 'name': 'silicon', 'symbol': 'Si', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 14, 'spin': 0.0},
    '29Si': {'A': 29, 'abundance': 4.6832, 'name': 'silicon', 'symbol': 'Si', 'stability': 'stable', 'quadrupole': 0,
             'gn': -1.1106, 'Z': 14, 'spin': 0.5},
    '30Si': {'A': 30, 'abundance': 3.0872, 'name': 'silicon', 'symbol': 'Si', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 14, 'spin': 0.0},
    '31P': {'A': 31, 'abundance': 100.0, 'name': 'phosphorus', 'symbol': 'P', 'stability': 'stable', 'quadrupole': 0,
            'gn': 2.2632, 'Z': 15, 'spin': 0.5},
    '32S': {'A': 32, 'abundance': 94.93, 'name': 'sulfur', 'symbol': 'S', 'stability': 'stable', 'quadrupole': 0,
            'gn': 0.0, 'Z': 16, 'spin': 0.0},
    '33S': {'A': 33, 'abundance': 0.76, 'name': 'sulfur', 'symbol': 'S', 'stability': 'stable', 'quadrupole': -0.0678,
            'gn': 0.42911, 'Z': 16, 'spin': 1.5},
    '34S': {'A': 34, 'abundance': 4.29, 'name': 'sulfur', 'symbol': 'S', 'stability': 'stable', 'quadrupole': 0,
            'gn': 0.0, 'Z': 16, 'spin': 0.0},
    '36S': {'A': 36, 'abundance': 0.02, 'name': 'sulfur', 'symbol': 'S', 'stability': 'stable', 'quadrupole': 0,
            'gn': 0.0, 'Z': 16, 'spin': 0.0},
    '35Cl': {'A': 35, 'abundance': 75.78, 'name': 'chlorine', 'symbol': 'Cl', 'stability': 'stable',
             'quadrupole': -0.08165,
             'gn': 0.5479157, 'Z': 17, 'spin': 1.5},
    '36Cl': {'A': 36, 'abundance': 0.0, 'name': 'chlorine', 'symbol': 'Cl', 'stability': 'radioactive',
             'quadrupole': -0.018,
             'gn': 0.642735, 'Z': 17, 'spin': 2.0},
    '37Cl': {'A': 37, 'abundance': 24.22, 'name': 'chlorine', 'symbol': 'Cl', 'stability': 'stable',
             'quadrupole': -0.06435,
             'gn': 0.456082, 'Z': 17, 'spin': 1.5},
    '36Ar': {'A': 36, 'abundance': 0.3365, 'name': 'argon', 'symbol': 'Ar', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 18, 'spin': 0.0},
    '38Ar': {'A': 38, 'abundance': 0.0632, 'name': 'argon', 'symbol': 'Ar', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 18, 'spin': 0.0},
    '39Ar': {'A': 39, 'abundance': 0.0, 'name': 'argon', 'symbol': 'Ar', 'stability': 'radioactive',
             'quadrupole': -0.12,
             'gn': -0.37, 'Z': 18, 'spin': 3.5},
    '40Ar': {'A': 40, 'abundance': 99.6003, 'name': 'argon', 'symbol': 'Ar', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 18, 'spin': 0.0},
    '39K': {'A': 39, 'abundance': 93.2581, 'name': 'potassium', 'symbol': 'K', 'stability': 'stable',
            'quadrupole': 0.0585,
            'gn': 0.2609909, 'Z': 19, 'spin': 1.5},
    '40K': {'A': 40, 'abundance': 0.0117, 'name': 'potassium', 'symbol': 'K', 'stability': 'stable',
            'quadrupole': -0.073,
            'gn': -0.32453, 'Z': 19, 'spin': 4.0},
    '41K': {'A': 41, 'abundance': 6.7302, 'name': 'potassium', 'symbol': 'K', 'stability': 'stable',
            'quadrupole': 0.0711,
            'gn': 0.1432542, 'Z': 19, 'spin': 1.5},
    '40Ca': {'A': 40, 'abundance': 96.941, 'name': 'calcium', 'symbol': 'Ca', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 20, 'spin': 0.0},
    '41Ca': {'A': 41, 'abundance': 0.0, 'name': 'calcium', 'symbol': 'Ca', 'stability': 'radioactive',
             'quadrupole': -0.0665,
             'gn': -0.4556514, 'Z': 20, 'spin': 3.5},
    '42Ca': {'A': 42, 'abundance': 0.647, 'name': 'calcium', 'symbol': 'Ca', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 20, 'spin': 0.0},
    '43Ca': {'A': 43, 'abundance': 0.135, 'name': 'calcium', 'symbol': 'Ca', 'stability': 'stable',
             'quadrupole': -0.0408,
             'gn': -0.376414, 'Z': 20, 'spin': 3.5},
    '44Ca': {'A': 44, 'abundance': 2.086, 'name': 'calcium', 'symbol': 'Ca', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 20, 'spin': 0.0},
    '46Ca': {'A': 46, 'abundance': 0.004, 'name': 'calcium', 'symbol': 'Ca', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 20, 'spin': 0.0},
    '48Ca': {'A': 48, 'abundance': 0.187, 'name': 'calcium', 'symbol': 'Ca', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 20, 'spin': 0.0},
    '45Sc': {'A': 45, 'abundance': 100.0, 'name': 'scandium', 'symbol': 'Sc', 'stability': 'stable',
             'quadrupole': -0.22,
             'gn': 1.35906, 'Z': 21, 'spin': 3.5},
    '46Ti': {'A': 46, 'abundance': 8.25, 'name': 'titanium', 'symbol': 'Ti', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 22, 'spin': 0.0},
    '47Ti': {'A': 47, 'abundance': 7.44, 'name': 'titanium', 'symbol': 'Ti', 'stability': 'stable', 'quadrupole': 0.302,
             'gn': -0.31539, 'Z': 22, 'spin': 2.5},
    '48Ti': {'A': 48, 'abundance': 73.72, 'name': 'titanium', 'symbol': 'Ti', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 22, 'spin': 0.0},
    '49Ti': {'A': 49, 'abundance': 5.41, 'name': 'titanium', 'symbol': 'Ti', 'stability': 'stable', 'quadrupole': 0.247,
             'gn': -0.315477, 'Z': 22, 'spin': 3.5},
    '50Ti': {'A': 50, 'abundance': 5.18, 'name': 'titanium', 'symbol': 'Ti', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 22, 'spin': 0.0},
    '50V': {'A': 50, 'abundance': 0.25, 'name': 'vanadium', 'symbol': 'V', 'stability': 'stable', 'quadrupole': 0.21,
            'gn': 0.556593, 'Z': 23, 'spin': 6.0},
    '51V': {'A': 51, 'abundance': 99.75, 'name': 'vanadium', 'symbol': 'V', 'stability': 'stable', 'quadrupole': -0.052,
            'gn': 1.46836, 'Z': 23, 'spin': 3.5},
    '50Cr': {'A': 50, 'abundance': 4.345, 'name': 'chromium', 'symbol': 'Cr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 24, 'spin': 0.0},
    '52Cr': {'A': 52, 'abundance': 83.789, 'name': 'chromium', 'symbol': 'Cr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 24, 'spin': 0.0},
    '53Cr': {'A': 53, 'abundance': 9.501, 'name': 'chromium', 'symbol': 'Cr', 'stability': 'stable',
             'quadrupole': -0.15,
             'gn': -0.3147, 'Z': 24, 'spin': 1.5},
    '54Cr': {'A': 54, 'abundance': 2.365, 'name': 'chromium', 'symbol': 'Cr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 24, 'spin': 0.0},
    '53Mn': {'A': 53, 'abundance': 0.0, 'name': 'manganese', 'symbol': 'Mn', 'stability': 'radioactive',
             'quadrupole': nan,
             'gn': 1.435, 'Z': 25, 'spin': 3.5},
    '55Mn': {'A': 55, 'abundance': 100.0, 'name': 'manganese', 'symbol': 'Mn', 'stability': 'stable',
             'quadrupole': 0.33,
             'gn': 1.3819, 'Z': 25, 'spin': 2.5},
    '54Fe': {'A': 54, 'abundance': 5.845, 'name': 'iron', 'symbol': 'Fe', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 26, 'spin': 0.0},
    '56Fe': {'A': 56, 'abundance': 91.754, 'name': 'iron', 'symbol': 'Fe', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 26, 'spin': 0.0},
    '57Fe': {'A': 57, 'abundance': 2.119, 'name': 'iron', 'symbol': 'Fe', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.1806, 'Z': 26, 'spin': 0.5},
    '58Fe': {'A': 58, 'abundance': 0.282, 'name': 'iron', 'symbol': 'Fe', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 26, 'spin': 0.0},
    '59Co': {'A': 59, 'abundance': 100.0, 'name': 'cobalt', 'symbol': 'Co', 'stability': 'stable', 'quadrupole': 0.42,
             'gn': 1.318, 'Z': 27, 'spin': 3.5},
    '60Co': {'A': 60, 'abundance': 0.0, 'name': 'cobalt', 'symbol': 'Co', 'stability': 'radioactive',
             'quadrupole': 0.44,
             'gn': 0.7589, 'Z': 27, 'spin': 5.0},
    '58Ni': {'A': 58, 'abundance': 68.0769, 'name': 'nickel', 'symbol': 'Ni', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 28, 'spin': 0.0},
    '60Ni': {'A': 60, 'abundance': 26.2231, 'name': 'nickel', 'symbol': 'Ni', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 28, 'spin': 0.0},
    '61Ni': {'A': 61, 'abundance': 1.1399, 'name': 'nickel', 'symbol': 'Ni', 'stability': 'stable', 'quadrupole': 0.162,
             'gn': -0.50001, 'Z': 28, 'spin': 1.5},
    '62Ni': {'A': 62, 'abundance': 3.6345, 'name': 'nickel', 'symbol': 'Ni', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 28, 'spin': 0.0},
    '64Ni': {'A': 64, 'abundance': 0.9256, 'name': 'nickel', 'symbol': 'Ni', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 28, 'spin': 0.0},
    '63Cu': {'A': 63, 'abundance': 69.17, 'name': 'copper', 'symbol': 'Cu', 'stability': 'stable', 'quadrupole': -0.22,
             'gn': 1.484, 'Z': 29, 'spin': 1.5},
    '65Cu': {'A': 65, 'abundance': 30.83, 'name': 'copper', 'symbol': 'Cu', 'stability': 'stable', 'quadrupole': -0.204,
             'gn': 1.588, 'Z': 29, 'spin': 1.5},
    '64Zn': {'A': 64, 'abundance': 48.63, 'name': 'zinc', 'symbol': 'Zn', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 30, 'spin': 0.0},
    '66Zn': {'A': 66, 'abundance': 27.9, 'name': 'zinc', 'symbol': 'Zn', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 30, 'spin': 0.0},
    '67Zn': {'A': 67, 'abundance': 4.1, 'name': 'zinc', 'symbol': 'Zn', 'stability': 'stable', 'quadrupole': 0.15,
             'gn': 0.350312, 'Z': 30, 'spin': 2.5},
    '68Zn': {'A': 68, 'abundance': 18.75, 'name': 'zinc', 'symbol': 'Zn', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 30, 'spin': 0.0},
    '70Zn': {'A': 70, 'abundance': 0.62, 'name': 'zinc', 'symbol': 'Zn', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 30, 'spin': 0.0},
    '69Ga': {'A': 69, 'abundance': 60.108, 'name': 'gallium', 'symbol': 'Ga', 'stability': 'stable',
             'quadrupole': 0.171,
             'gn': 1.34439, 'Z': 31, 'spin': 1.5},
    '71Ga': {'A': 71, 'abundance': 39.892, 'name': 'gallium', 'symbol': 'Ga', 'stability': 'stable',
             'quadrupole': 0.107,
             'gn': 1.70818, 'Z': 31, 'spin': 1.5},
    '70Ge': {'A': 70, 'abundance': 20.84, 'name': 'germanium', 'symbol': 'Ge', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 32, 'spin': 0.0},
    '72Ge': {'A': 72, 'abundance': 27.54, 'name': 'germanium', 'symbol': 'Ge', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 32, 'spin': 0.0},
    '73Ge': {'A': 73, 'abundance': 7.73, 'name': 'germanium', 'symbol': 'Ge', 'stability': 'stable',
             'quadrupole': -0.196,
             'gn': -0.1954371, 'Z': 32, 'spin': 4.5},
    '74Ge': {'A': 74, 'abundance': 36.28, 'name': 'germanium', 'symbol': 'Ge', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 32, 'spin': 0.0},
    '76Ge': {'A': 76, 'abundance': 7.61, 'name': 'germanium', 'symbol': 'Ge', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 32, 'spin': 0.0},
    '75As': {'A': 75, 'abundance': 100.0, 'name': 'arsenic', 'symbol': 'As', 'stability': 'stable', 'quadrupole': 0.314,
             'gn': 0.959647, 'Z': 33, 'spin': 1.5},
    '74Se': {'A': 74, 'abundance': 0.89, 'name': 'selenium', 'symbol': 'Se', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 34, 'spin': 0.0},
    '76Se': {'A': 76, 'abundance': 9.37, 'name': 'selenium', 'symbol': 'Se', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 34, 'spin': 0.0},
    '77Se': {'A': 77, 'abundance': 7.63, 'name': 'selenium', 'symbol': 'Se', 'stability': 'stable', 'quadrupole': 0,
             'gn': 1.0693, 'Z': 34, 'spin': 0.5},
    '78Se': {'A': 78, 'abundance': 23.77, 'name': 'selenium', 'symbol': 'Se', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 34, 'spin': 0.0},
    '79Se': {'A': 79, 'abundance': 0.0, 'name': 'selenium', 'symbol': 'Se', 'stability': 'radioactive',
             'quadrupole': 0.8,
             'gn': -0.2908, 'Z': 34, 'spin': 3.5},
    '80Se': {'A': 80, 'abundance': 49.61, 'name': 'selenium', 'symbol': 'Se', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 34, 'spin': 0.0},
    '82Se': {'A': 82, 'abundance': 8.73, 'name': 'selenium', 'symbol': 'Se', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 34, 'spin': 0.0},
    '79Br': {'A': 79, 'abundance': 50.69, 'name': 'bromine', 'symbol': 'Br', 'stability': 'stable', 'quadrupole': 0.313,
             'gn': 1.404266, 'Z': 35, 'spin': 1.5},
    '81Br': {'A': 81, 'abundance': 49.31, 'name': 'bromine', 'symbol': 'Br', 'stability': 'stable', 'quadrupole': 0.262,
             'gn': 1.513706, 'Z': 35, 'spin': 1.5},
    '78Kr': {'A': 78, 'abundance': 0.35, 'name': 'krypton', 'symbol': 'Kr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 36, 'spin': 0.0},
    '80Kr': {'A': 80, 'abundance': 2.28, 'name': 'krypton', 'symbol': 'Kr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 36, 'spin': 0.0},
    '82Kr': {'A': 82, 'abundance': 11.58, 'name': 'krypton', 'symbol': 'Kr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 36, 'spin': 0.0},
    '83Kr': {'A': 83, 'abundance': 11.49, 'name': 'krypton', 'symbol': 'Kr', 'stability': 'stable', 'quadrupole': 0.259,
             'gn': -0.215704, 'Z': 36, 'spin': 4.5},
    '84Kr': {'A': 84, 'abundance': 57.0, 'name': 'krypton', 'symbol': 'Kr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 36, 'spin': 0.0},
    '85Kr': {'A': 85, 'abundance': 0.0, 'name': 'krypton', 'symbol': 'Kr', 'stability': 'radioactive',
             'quadrupole': 0.44,
             'gn': 0.2233, 'Z': 36, 'spin': 4.5},
    '86Kr': {'A': 86, 'abundance': 17.3, 'name': 'krypton', 'symbol': 'Kr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 36, 'spin': 0.0},
    '85Rb': {'A': 85, 'abundance': 72.17, 'name': 'rubidium', 'symbol': 'Rb', 'stability': 'stable',
             'quadrupole': 0.276,
             'gn': 0.541253, 'Z': 37, 'spin': 2.5},
    '87Rb': {'A': 87, 'abundance': 27.83, 'name': 'rubidium', 'symbol': 'Rb', 'stability': 'stable',
             'quadrupole': 0.1335,
             'gn': 1.83427, 'Z': 37, 'spin': 1.5},
    '84Sr': {'A': 84, 'abundance': 0.56, 'name': 'strontium', 'symbol': 'Sr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 38, 'spin': 0.0},
    '86Sr': {'A': 86, 'abundance': 9.86, 'name': 'strontium', 'symbol': 'Sr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 38, 'spin': 0.0},
    '87Sr': {'A': 87, 'abundance': 7.0, 'name': 'strontium', 'symbol': 'Sr', 'stability': 'stable', 'quadrupole': 0.335,
             'gn': -0.24291, 'Z': 38, 'spin': 4.5},
    '88Sr': {'A': 88, 'abundance': 82.58, 'name': 'strontium', 'symbol': 'Sr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 38, 'spin': 0.0},
    '89Y': {'A': 89, 'abundance': 100.0, 'name': 'yttrium', 'symbol': 'Y', 'stability': 'stable', 'quadrupole': 0,
            'gn': -0.2748361, 'Z': 39, 'spin': 0.5},
    '90Zr': {'A': 90, 'abundance': 51.45, 'name': 'zirconium', 'symbol': 'Zr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 40, 'spin': 0.0},
    '91Zr': {'A': 91, 'abundance': 11.22, 'name': 'zirconium', 'symbol': 'Zr', 'stability': 'stable',
             'quadrupole': -0.176,
             'gn': -0.521448, 'Z': 40, 'spin': 2.5},
    '92Zr': {'A': 92, 'abundance': 17.15, 'name': 'zirconium', 'symbol': 'Zr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 40, 'spin': 0.0},
    '94Zr': {'A': 94, 'abundance': 17.38, 'name': 'zirconium', 'symbol': 'Zr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 40, 'spin': 0.0},
    '96Zr': {'A': 96, 'abundance': 2.8, 'name': 'zirconium', 'symbol': 'Zr', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 40, 'spin': 0.0},
    '93Nb': {'A': 93, 'abundance': 100.0, 'name': 'niobium', 'symbol': 'Nb', 'stability': 'stable', 'quadrupole': -0.32,
             'gn': 1.3712, 'Z': 41, 'spin': 4.5},
    '92Mo': {'A': 92, 'abundance': 14.84, 'name': 'molybdenum', 'symbol': 'Mo', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 42, 'spin': 0.0},
    '94Mo': {'A': 94, 'abundance': 9.25, 'name': 'molybdenum', 'symbol': 'Mo', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 42, 'spin': 0.0},
    '95Mo': {'A': 95, 'abundance': 15.92, 'name': 'molybdenum', 'symbol': 'Mo', 'stability': 'stable',
             'quadrupole': -0.022,
             'gn': -0.3656, 'Z': 42, 'spin': 2.5},
    '96Mo': {'A': 96, 'abundance': 16.68, 'name': 'molybdenum', 'symbol': 'Mo', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 42, 'spin': 0.0},
    '97Mo': {'A': 97, 'abundance': 9.55, 'name': 'molybdenum', 'symbol': 'Mo', 'stability': 'stable',
             'quadrupole': 0.255,
             'gn': -0.3734, 'Z': 42, 'spin': 2.5},
    '98Mo': {'A': 98, 'abundance': 24.13, 'name': 'molybdenum', 'symbol': 'Mo', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 42, 'spin': 0.0},
    '100Mo': {'A': 100, 'abundance': 9.63, 'name': 'molybdenum', 'symbol': 'Mo', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 42, 'spin': 0.0},
    '99Tc': {'A': 99, 'abundance': 0.0, 'name': 'technetium', 'symbol': 'Tc', 'stability': 'radioactive',
             'quadrupole': -0.129,
             'gn': 1.2632, 'Z': 43, 'spin': 4.5},
    '96Ru': {'A': 96, 'abundance': 5.54, 'name': 'ruthenium', 'symbol': 'Ru', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 44, 'spin': 0.0},
    '98Ru': {'A': 98, 'abundance': 1.87, 'name': 'ruthenium', 'symbol': 'Ru', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 44, 'spin': 0.0},
    '99Ru': {'A': 99, 'abundance': 12.76, 'name': 'ruthenium', 'symbol': 'Ru', 'stability': 'stable',
             'quadrupole': 0.079,
             'gn': -0.249, 'Z': 44, 'spin': 2.5},
    '100Ru': {'A': 100, 'abundance': 12.6, 'name': 'ruthenium', 'symbol': 'Ru', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 44, 'spin': 0.0},
    '101Ru': {'A': 101, 'abundance': 17.06, 'name': 'ruthenium', 'symbol': 'Ru', 'stability': 'stable',
              'quadrupole': 0.457,
              'gn': -0.279, 'Z': 44, 'spin': 2.5},
    '102Ru': {'A': 102, 'abundance': 31.55, 'name': 'ruthenium', 'symbol': 'Ru', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 44, 'spin': 0.0},
    '104Ru': {'A': 104, 'abundance': 18.62, 'name': 'ruthenium', 'symbol': 'Ru', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 44, 'spin': 0.0},
    '102Rh': {'A': 102, 'abundance': 0.0, 'name': 'rhodium', 'symbol': 'Rh', 'stability': 'radioactive',
              'quadrupole': nan,
              'gn': 0.685, 'Z': 45, 'spin': 6.0},
    '103Rh': {'A': 103, 'abundance': 100.0, 'name': 'rhodium', 'symbol': 'Rh', 'stability': 'stable', 'quadrupole': 0,
              'gn': -0.1768, 'Z': 45, 'spin': 0.5},
    '102Pd': {'A': 102, 'abundance': 1.02, 'name': 'palladium', 'symbol': 'Pd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 46, 'spin': 0.0},
    '104Pd': {'A': 104, 'abundance': 11.14, 'name': 'palladium', 'symbol': 'Pd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 46, 'spin': 0.0},
    '105Pd': {'A': 105, 'abundance': 22.33, 'name': 'palladium', 'symbol': 'Pd', 'stability': 'stable',
              'quadrupole': 0.66,
              'gn': -0.256, 'Z': 46, 'spin': 2.5},
    '106Pd': {'A': 106, 'abundance': 27.33, 'name': 'palladium', 'symbol': 'Pd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 46, 'spin': 0.0},
    '108Pd': {'A': 108, 'abundance': 26.46, 'name': 'palladium', 'symbol': 'Pd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 46, 'spin': 0.0},
    '110Pd': {'A': 110, 'abundance': 11.72, 'name': 'palladium', 'symbol': 'Pd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 46, 'spin': 0.0},
    '107Ag': {'A': 107, 'abundance': 51.839, 'name': 'silver', 'symbol': 'Ag', 'stability': 'stable', 'quadrupole': 0,
              'gn': -0.227249, 'Z': 47, 'spin': 0.5},
    '109Ag': {'A': 109, 'abundance': 48.161, 'name': 'silver', 'symbol': 'Ag', 'stability': 'stable', 'quadrupole': 0,
              'gn': -0.261743, 'Z': 47, 'spin': 0.5},
    '106Cd': {'A': 106, 'abundance': 1.25, 'name': 'cadmium', 'symbol': 'Cd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 48, 'spin': 0.0},
    '108Cd': {'A': 108, 'abundance': 0.89, 'name': 'cadmium', 'symbol': 'Cd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 48, 'spin': 0.0},
    '110Cd': {'A': 110, 'abundance': 12.49, 'name': 'cadmium', 'symbol': 'Cd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 48, 'spin': 0.0},
    '111Cd': {'A': 111, 'abundance': 12.8, 'name': 'cadmium', 'symbol': 'Cd', 'stability': 'stable', 'quadrupole': 0,
              'gn': -1.19043, 'Z': 48, 'spin': 0.5},
    '112Cd': {'A': 112, 'abundance': 24.13, 'name': 'cadmium', 'symbol': 'Cd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 48, 'spin': 0.0},
    '113Cd': {'A': 113, 'abundance': 12.22, 'name': 'cadmium', 'symbol': 'Cd', 'stability': 'stable', 'quadrupole': 0,
              'gn': -1.2454, 'Z': 48, 'spin': 0.5},
    '114Cd': {'A': 114, 'abundance': 28.73, 'name': 'cadmium', 'symbol': 'Cd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 48, 'spin': 0.0},
    '116Cd': {'A': 116, 'abundance': 7.49, 'name': 'cadmium', 'symbol': 'Cd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 48, 'spin': 0.0},
    '113In': {'A': 113, 'abundance': 4.29, 'name': 'indium', 'symbol': 'In', 'stability': 'stable', 'quadrupole': 0.799,
              'gn': 1.22864, 'Z': 49, 'spin': 4.5},
    '115In': {'A': 115, 'abundance': 95.71, 'name': 'indium', 'symbol': 'In', 'stability': 'stable', 'quadrupole': 0.81,
              'gn': 1.23129, 'Z': 49, 'spin': 4.5},
    '112Sn': {'A': 112, 'abundance': 0.97, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 50, 'spin': 0.0},
    '114Sn': {'A': 114, 'abundance': 0.66, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 50, 'spin': 0.0},
    '115Sn': {'A': 115, 'abundance': 0.34, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': -1.83766, 'Z': 50, 'spin': 0.5},
    '116Sn': {'A': 116, 'abundance': 14.54, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 50, 'spin': 0.0},
    '117Sn': {'A': 117, 'abundance': 7.68, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': -2.00208, 'Z': 50, 'spin': 0.5},
    '118Sn': {'A': 118, 'abundance': 24.22, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 50, 'spin': 0.0},
    '119Sn': {'A': 119, 'abundance': 8.59, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': -2.09456, 'Z': 50, 'spin': 0.5},
    '120Sn': {'A': 120, 'abundance': 32.58, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 50, 'spin': 0.0},
    '122Sn': {'A': 122, 'abundance': 4.63, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 50, 'spin': 0.0},
    '124Sn': {'A': 124, 'abundance': 5.79, 'name': 'tin', 'symbol': 'Sn', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 50, 'spin': 0.0},
    '121Sb': {'A': 121, 'abundance': 57.21, 'name': 'antimony', 'symbol': 'Sb', 'stability': 'stable',
              'quadrupole': -0.36,
              'gn': 1.3455, 'Z': 51, 'spin': 2.5},
    '123Sb': {'A': 123, 'abundance': 42.79, 'name': 'antimony', 'symbol': 'Sb', 'stability': 'stable',
              'quadrupole': -0.49,
              'gn': 0.72876, 'Z': 51, 'spin': 3.5},
    '125Sb': {'A': 125, 'abundance': 0.0, 'name': 'antimony', 'symbol': 'Sb', 'stability': 'radioactive',
              'quadrupole': nan,
              'gn': 0.7514, 'Z': 51, 'spin': 3.5},
    '120Te': {'A': 120, 'abundance': 0.09, 'name': 'tellurium', 'symbol': 'Te', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 52, 'spin': 0.0},
    '122Te': {'A': 122, 'abundance': 2.55, 'name': 'tellurium', 'symbol': 'Te', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 52, 'spin': 0.0},
    '123Te': {'A': 123, 'abundance': 0.89, 'name': 'tellurium', 'symbol': 'Te', 'stability': 'stable', 'quadrupole': 0,
              'gn': -1.4736, 'Z': 52, 'spin': 0.5},
    '124Te': {'A': 124, 'abundance': 4.74, 'name': 'tellurium', 'symbol': 'Te', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 52, 'spin': 0.0},
    '125Te': {'A': 125, 'abundance': 7.07, 'name': 'tellurium', 'symbol': 'Te', 'stability': 'stable', 'quadrupole': 0,
              'gn': -1.7766, 'Z': 52, 'spin': 0.5},
    '126Te': {'A': 126, 'abundance': 18.84, 'name': 'tellurium', 'symbol': 'Te', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 52, 'spin': 0.0},
    '128Te': {'A': 128, 'abundance': 31.74, 'name': 'tellurium', 'symbol': 'Te', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 52, 'spin': 0.0},
    '130Te': {'A': 130, 'abundance': 34.08, 'name': 'tellurium', 'symbol': 'Te', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 52, 'spin': 0.0},
    '127I': {'A': 127, 'abundance': 100.0, 'name': 'iodine', 'symbol': 'I', 'stability': 'stable', 'quadrupole': -0.71,
             'gn': 1.1253, 'Z': 53, 'spin': 2.5},
    '129I': {'A': 129, 'abundance': 0.0, 'name': 'iodine', 'symbol': 'I', 'stability': 'radioactive',
             'quadrupole': -0.482,
             'gn': 0.74886, 'Z': 53, 'spin': 3.5},
    '124Xe': {'A': 124, 'abundance': 0.09, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 54, 'spin': 0.0},
    '126Xe': {'A': 126, 'abundance': 0.09, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 54, 'spin': 0.0},
    '128Xe': {'A': 128, 'abundance': 1.92, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 54, 'spin': 0.0},
    '129Xe': {'A': 129, 'abundance': 26.44, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable', 'quadrupole': 0,
              'gn': -1.55595, 'Z': 54, 'spin': 0.5},
    '130Xe': {'A': 130, 'abundance': 4.08, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 54, 'spin': 0.0},
    '131Xe': {'A': 131, 'abundance': 21.18, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable',
              'quadrupole': -0.114,
              'gn': 0.46124, 'Z': 54, 'spin': 1.5},
    '132Xe': {'A': 132, 'abundance': 26.89, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 54, 'spin': 0.0},
    '134Xe': {'A': 134, 'abundance': 10.44, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 54, 'spin': 0.0},
    '136Xe': {'A': 136, 'abundance': 8.87, 'name': 'xenon', 'symbol': 'Xe', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 54, 'spin': 0.0},
    '133Cs': {'A': 133, 'abundance': 100.0, 'name': 'caesium', 'symbol': 'Cs', 'stability': 'stable',
              'quadrupole': -0.00343,
              'gn': 0.7378377, 'Z': 55, 'spin': 3.5},
    '134Cs': {'A': 134, 'abundance': 0.0, 'name': 'caesium', 'symbol': 'Cs', 'stability': 'radioactive',
              'quadrupole': 0.389,
              'gn': 0.74842, 'Z': 55, 'spin': 4.0},
    '135Cs': {'A': 135, 'abundance': 0.0, 'name': 'caesium', 'symbol': 'Cs', 'stability': 'radioactive',
              'quadrupole': 0.05,
              'gn': 0.78069, 'Z': 55, 'spin': 3.5},
    '137Cs': {'A': 137, 'abundance': 0.0, 'name': 'caesium', 'symbol': 'Cs', 'stability': 'radioactive',
              'quadrupole': 0.051,
              'gn': 0.8118, 'Z': 55, 'spin': 3.5},
    '130Ba': {'A': 130, 'abundance': 0.106, 'name': 'barium', 'symbol': 'Ba', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 56, 'spin': 0.0},
    '132Ba': {'A': 132, 'abundance': 0.101, 'name': 'barium', 'symbol': 'Ba', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 56, 'spin': 0.0},
    '133Ba': {'A': 133, 'abundance': 0.0, 'name': 'barium', 'symbol': 'Ba', 'stability': 'radioactive', 'quadrupole': 0,
              'gn': -1.54, 'Z': 56, 'spin': 0.5},
    '134Ba': {'A': 134, 'abundance': 2.417, 'name': 'barium', 'symbol': 'Ba', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 56, 'spin': 0.0},
    '135Ba': {'A': 135, 'abundance': 6.592, 'name': 'barium', 'symbol': 'Ba', 'stability': 'stable', 'quadrupole': 0.16,
              'gn': 0.55884, 'Z': 56, 'spin': 1.5},
    '136Ba': {'A': 136, 'abundance': 7.854, 'name': 'barium', 'symbol': 'Ba', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 56, 'spin': 0.0},
    '137Ba': {'A': 137, 'abundance': 11.232, 'name': 'barium', 'symbol': 'Ba', 'stability': 'stable',
              'quadrupole': 0.245,
              'gn': 0.62515, 'Z': 56, 'spin': 1.5},
    '138Ba': {'A': 138, 'abundance': 71.698, 'name': 'barium', 'symbol': 'Ba', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 56, 'spin': 0.0},
    '137La': {'A': 137, 'abundance': 0.0, 'name': 'lanthanum', 'symbol': 'La', 'stability': 'radioactive',
              'quadrupole': 0.26,
              'gn': 0.77, 'Z': 57, 'spin': 3.5},
    '138La': {'A': 138, 'abundance': 0.09, 'name': 'lanthanum', 'symbol': 'La', 'stability': 'stable',
              'quadrupole': 0.45,
              'gn': 0.74278, 'Z': 57, 'spin': 5.0},
    '139La': {'A': 139, 'abundance': 99.91, 'name': 'lanthanum', 'symbol': 'La', 'stability': 'stable',
              'quadrupole': 0.2,
              'gn': 0.7952, 'Z': 57, 'spin': 3.5},
    '136Ce': {'A': 136, 'abundance': 0.185, 'name': 'cerium', 'symbol': 'Ce', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 58, 'spin': 0.0},
    '138Ce': {'A': 138, 'abundance': 0.251, 'name': 'cerium', 'symbol': 'Ce', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 58, 'spin': 0.0},
    '140Ce': {'A': 140, 'abundance': 88.45, 'name': 'cerium', 'symbol': 'Ce', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 58, 'spin': 0.0},
    '142Ce': {'A': 142, 'abundance': 11.114, 'name': 'cerium', 'symbol': 'Ce', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 58, 'spin': 0.0},
    '141Pr': {'A': 141, 'abundance': 100.0, 'name': 'praesodymium', 'symbol': 'Pr', 'stability': 'stable',
              'quadrupole': -0.0589,
              'gn': 1.6, 'Z': 59, 'spin': 2.5},
    '142Nd': {'A': 142, 'abundance': 27.2, 'name': 'neodymium', 'symbol': 'Nd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 60, 'spin': 0.0},
    '143Nd': {'A': 143, 'abundance': 12.2, 'name': 'neodymium', 'symbol': 'Nd', 'stability': 'stable',
              'quadrupole': -0.63,
              'gn': -0.3076, 'Z': 60, 'spin': 3.5},
    '144Nd': {'A': 144, 'abundance': 23.8, 'name': 'neodymium', 'symbol': 'Nd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 60, 'spin': 0.0},
    '145Nd': {'A': 145, 'abundance': 8.3, 'name': 'neodymium', 'symbol': 'Nd', 'stability': 'stable',
              'quadrupole': -0.33,
              'gn': -0.19, 'Z': 60, 'spin': 3.5},
    '146Nd': {'A': 146, 'abundance': 17.2, 'name': 'neodymium', 'symbol': 'Nd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 60, 'spin': 0.0},
    '148Nd': {'A': 148, 'abundance': 5.7, 'name': 'neodymium', 'symbol': 'Nd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 60, 'spin': 0.0},
    '150Nd': {'A': 150, 'abundance': 5.6, 'name': 'neodymium', 'symbol': 'Nd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 60, 'spin': 0.0},
    '147Pm': {'A': 147, 'abundance': 0.0, 'name': 'promethium', 'symbol': 'Pm', 'stability': 'radioactive',
              'quadrupole': 0.74,
              'gn': 0.752, 'Z': 61, 'spin': 3.5},
    '144Sm': {'A': 144, 'abundance': 3.07, 'name': 'samarium', 'symbol': 'Sm', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 62, 'spin': 0.0},
    '147Sm': {'A': 147, 'abundance': 14.99, 'name': 'samarium', 'symbol': 'Sm', 'stability': 'stable',
              'quadrupole': -0.259,
              'gn': -0.2322, 'Z': 62, 'spin': 3.5},
    '148Sm': {'A': 148, 'abundance': 11.24, 'name': 'samarium', 'symbol': 'Sm', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 62, 'spin': 0.0},
    '149Sm': {'A': 149, 'abundance': 13.82, 'name': 'samarium', 'symbol': 'Sm', 'stability': 'stable',
              'quadrupole': 0.074,
              'gn': 0.1915, 'Z': 62, 'spin': 3.5},
    '150Sm': {'A': 150, 'abundance': 7.38, 'name': 'samarium', 'symbol': 'Sm', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 62, 'spin': 0.0},
    '151Sm': {'A': 151, 'abundance': 0.0, 'name': 'samarium', 'symbol': 'Sm', 'stability': 'radioactive',
              'quadrupole': 0.67,
              'gn': 0.142, 'Z': 62, 'spin': 2.5},
    '152Sm': {'A': 152, 'abundance': 26.75, 'name': 'samarium', 'symbol': 'Sm', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 62, 'spin': 0.0},
    '154Sm': {'A': 154, 'abundance': 22.75, 'name': 'samarium', 'symbol': 'Sm', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 62, 'spin': 0.0},
    '151Eu': {'A': 151, 'abundance': 47.81, 'name': 'europium', 'symbol': 'Eu', 'stability': 'stable',
              'quadrupole': 0.903,
              'gn': 1.389, 'Z': 63, 'spin': 2.5},
    '152Eu': {'A': 152, 'abundance': 0.0, 'name': 'europium', 'symbol': 'Eu', 'stability': 'radioactive',
              'quadrupole': 2.5,
              'gn': 0.64713, 'Z': 63, 'spin': 3.0},
    '153Eu': {'A': 153, 'abundance': 52.19, 'name': 'europium', 'symbol': 'Eu', 'stability': 'stable',
              'quadrupole': 2.412,
              'gn': 0.6134, 'Z': 63, 'spin': 2.5},
    '154Eu': {'A': 154, 'abundance': 0.0, 'name': 'europium', 'symbol': 'Eu', 'stability': 'radioactive',
              'quadrupole': 2.84,
              'gn': 0.6683, 'Z': 63, 'spin': 3.0},
    '155Eu': {'A': 155, 'abundance': 0.0, 'name': 'europium', 'symbol': 'Eu', 'stability': 'radioactive',
              'quadrupole': 2.3,
              'gn': 0.772, 'Z': 63, 'spin': 2.5},
    '152Gd': {'A': 152, 'abundance': 0.2, 'name': 'gadolinium', 'symbol': 'Gd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 64, 'spin': 0.0},
    '154Gd': {'A': 154, 'abundance': 2.18, 'name': 'gadolinium', 'symbol': 'Gd', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 64, 'spin': 0.0},
    '155Gd': {'A': 155, 'abundance': 14.8, 'name': 'gadolinium', 'symbol': 'Gd', 'stability': 'stable',
              'quadrupole': 1.27,
              'gn': -0.1723, 'Z': 64, 'spin': 1.5},
    '156Gd': {'A': 156, 'abundance': 20.47, 'name': 'gadolinium', 'symbol': 'Gd', 'stability': 'stable',
              'quadrupole': 0,
              'gn': 0.0, 'Z': 64, 'spin': 0.0},
    '157Gd': {'A': 157, 'abundance': 15.65, 'name': 'gadolinium', 'symbol': 'Gd', 'stability': 'stable',
              'quadrupole': 1.35,
              'gn': -0.2253, 'Z': 64, 'spin': 1.5},
    '158Gd': {'A': 158, 'abundance': 24.84, 'name': 'gadolinium', 'symbol': 'Gd', 'stability': 'stable',
              'quadrupole': 0,
              'gn': 0.0, 'Z': 64, 'spin': 0.0},
    '160Gd': {'A': 160, 'abundance': 21.86, 'name': 'gadolinium', 'symbol': 'Gd', 'stability': 'stable',
              'quadrupole': 0,
              'gn': 0.0, 'Z': 64, 'spin': 0.0},
    '157Tb': {'A': 157, 'abundance': 0.0, 'name': 'terbium', 'symbol': 'Tb', 'stability': 'radioactive',
              'quadrupole': 1.4,
              'gn': 1.3, 'Z': 65, 'spin': 1.5},
    '159Tb': {'A': 159, 'abundance': 100.0, 'name': 'terbium', 'symbol': 'Tb', 'stability': 'stable',
              'quadrupole': 1.432,
              'gn': 1.342, 'Z': 65, 'spin': 1.5},
    '160Tb': {'A': 160, 'abundance': 0.0, 'name': 'terbium', 'symbol': 'Tb', 'stability': 'radioactive',
              'quadrupole': 3.85,
              'gn': 0.597, 'Z': 65, 'spin': 3.0},
    '156Dy': {'A': 156, 'abundance': 0.06, 'name': 'dysprosium', 'symbol': 'Dy', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 66, 'spin': 0.0},
    '158Dy': {'A': 158, 'abundance': 0.1, 'name': 'dysprosium', 'symbol': 'Dy', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 66, 'spin': 0.0},
    '160Dy': {'A': 160, 'abundance': 2.34, 'name': 'dysprosium', 'symbol': 'Dy', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 66, 'spin': 0.0},
    '161Dy': {'A': 161, 'abundance': 18.91, 'name': 'dysprosium', 'symbol': 'Dy', 'stability': 'stable',
              'quadrupole': 2.507,
              'gn': -0.189, 'Z': 66, 'spin': 2.5},
    '162Dy': {'A': 162, 'abundance': 25.51, 'name': 'dysprosium', 'symbol': 'Dy', 'stability': 'stable',
              'quadrupole': 0,
              'gn': 0.0, 'Z': 66, 'spin': 0.0},
    '163Dy': {'A': 163, 'abundance': 24.9, 'name': 'dysprosium', 'symbol': 'Dy', 'stability': 'stable',
              'quadrupole': 2.648,
              'gn': 0.266, 'Z': 66, 'spin': 2.5},
    '164Dy': {'A': 164, 'abundance': 28.18, 'name': 'dysprosium', 'symbol': 'Dy', 'stability': 'stable',
              'quadrupole': 0,
              'gn': 0.0, 'Z': 66, 'spin': 0.0},
    '165Ho': {'A': 165, 'abundance': 100.0, 'name': 'holmium', 'symbol': 'Ho', 'stability': 'stable',
              'quadrupole': 3.58,
              'gn': 1.192, 'Z': 67, 'spin': 3.5},
    '162Er': {'A': 162, 'abundance': 0.14, 'name': 'erbium', 'symbol': 'Er', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 68, 'spin': 0.0},
    '164Er': {'A': 164, 'abundance': 1.61, 'name': 'erbium', 'symbol': 'Er', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 68, 'spin': 0.0},
    '166Er': {'A': 166, 'abundance': 33.61, 'name': 'erbium', 'symbol': 'Er', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 68, 'spin': 0.0},
    '167Er': {'A': 167, 'abundance': 22.93, 'name': 'erbium', 'symbol': 'Er', 'stability': 'stable',
              'quadrupole': 3.565,
              'gn': -0.1618, 'Z': 68, 'spin': 3.5},
    '168Er': {'A': 168, 'abundance': 26.78, 'name': 'erbium', 'symbol': 'Er', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 68, 'spin': 0.0},
    '170Er': {'A': 170, 'abundance': 14.93, 'name': 'erbium', 'symbol': 'Er', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 68, 'spin': 0.0},
    '169Tm': {'A': 169, 'abundance': 100.0, 'name': 'thulium', 'symbol': 'Tm', 'stability': 'stable', 'quadrupole': 0,
              'gn': -0.466, 'Z': 69, 'spin': 0.5},
    '171Tm': {'A': 171, 'abundance': 0.0, 'name': 'thulium', 'symbol': 'Tm', 'stability': 'radioactive',
              'quadrupole': 0,
              'gn': -0.4606, 'Z': 69, 'spin': 0.5},
    '168Yb': {'A': 168, 'abundance': 0.13, 'name': 'ytterbium', 'symbol': 'Yb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 70, 'spin': 0.0},
    '170Yb': {'A': 170, 'abundance': 3.04, 'name': 'ytterbium', 'symbol': 'Yb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 70, 'spin': 0.0},
    '171Yb': {'A': 171, 'abundance': 14.28, 'name': 'ytterbium', 'symbol': 'Yb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.9885, 'Z': 70, 'spin': 0.5},
    '172Yb': {'A': 172, 'abundance': 21.83, 'name': 'ytterbium', 'symbol': 'Yb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 70, 'spin': 0.0},
    '173Yb': {'A': 173, 'abundance': 16.13, 'name': 'ytterbium', 'symbol': 'Yb', 'stability': 'stable',
              'quadrupole': 2.8,
              'gn': -0.27195, 'Z': 70, 'spin': 2.5},
    '174Yb': {'A': 174, 'abundance': 31.83, 'name': 'ytterbium', 'symbol': 'Yb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 70, 'spin': 0.0},
    '176Yb': {'A': 176, 'abundance': 12.76, 'name': 'ytterbium', 'symbol': 'Yb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 70, 'spin': 0.0},
    '173Lu': {'A': 173, 'abundance': 0.0, 'name': 'lutetium', 'symbol': 'Lu', 'stability': 'radioactive',
              'quadrupole': 3.56,
              'gn': 0.669, 'Z': 71, 'spin': 3.5},
    '174Lu': {'A': 174, 'abundance': 0.0, 'name': 'lutetium', 'symbol': 'Lu', 'stability': 'radioactive',
              'quadrupole': nan,
              'gn': 1.94, 'Z': 71, 'spin': 1.0},
    '175Lu': {'A': 175, 'abundance': 97.41, 'name': 'lutetium', 'symbol': 'Lu', 'stability': 'stable',
              'quadrupole': 3.49,
              'gn': 0.63943, 'Z': 71, 'spin': 3.5},
    '176Lu': {'A': 176, 'abundance': 2.59, 'name': 'lutetium', 'symbol': 'Lu', 'stability': 'stable',
              'quadrupole': 4.97,
              'gn': 0.452, 'Z': 71, 'spin': 7.0},
    '174Hf': {'A': 174, 'abundance': 0.16, 'name': 'hafnium', 'symbol': 'Hf', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 72, 'spin': 0.0},
    '176Hf': {'A': 176, 'abundance': 5.26, 'name': 'hafnium', 'symbol': 'Hf', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 72, 'spin': 0.0},
    '177Hf': {'A': 177, 'abundance': 18.6, 'name': 'hafnium', 'symbol': 'Hf', 'stability': 'stable',
              'quadrupole': 3.365,
              'gn': 0.2267, 'Z': 72, 'spin': 3.5},
    '178Hf': {'A': 178, 'abundance': 27.28, 'name': 'hafnium', 'symbol': 'Hf', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 72, 'spin': 0.0},
    '179Hf': {'A': 179, 'abundance': 13.62, 'name': 'hafnium', 'symbol': 'Hf', 'stability': 'stable',
              'quadrupole': 3.793,
              'gn': -0.1424, 'Z': 72, 'spin': 4.5},
    '180Hf': {'A': 180, 'abundance': 35.08, 'name': 'hafnium', 'symbol': 'Hf', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 72, 'spin': 0.0},
    '180Ta': {'A': 180, 'abundance': 0.012, 'name': 'tantalum', 'symbol': 'Ta', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 73, 'spin': 0.0},
    '181Ta': {'A': 181, 'abundance': 99.988, 'name': 'tantalum', 'symbol': 'Ta', 'stability': 'stable',
              'quadrupole': 3.17,
              'gn': 0.67729, 'Z': 73, 'spin': 3.5},
    '180W': {'A': 180, 'abundance': 0.12, 'name': 'tungsten', 'symbol': 'W', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 74, 'spin': 0.0},
    '182W': {'A': 182, 'abundance': 26.5, 'name': 'tungsten', 'symbol': 'W', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 74, 'spin': 0.0},
    '183W': {'A': 183, 'abundance': 14.31, 'name': 'tungsten', 'symbol': 'W', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.2355694, 'Z': 74, 'spin': 0.5},
    '184W': {'A': 184, 'abundance': 30.64, 'name': 'tungsten', 'symbol': 'W', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 74, 'spin': 0.0},
    '186W': {'A': 186, 'abundance': 28.43, 'name': 'tungsten', 'symbol': 'W', 'stability': 'stable', 'quadrupole': 0,
             'gn': 0.0, 'Z': 74, 'spin': 0.0},
    '185Re': {'A': 185, 'abundance': 37.4, 'name': 'rhenium', 'symbol': 'Re', 'stability': 'stable', 'quadrupole': 2.18,
              'gn': 1.2748, 'Z': 75, 'spin': 2.5},
    '187Re': {'A': 187, 'abundance': 62.6, 'name': 'rhenium', 'symbol': 'Re', 'stability': 'stable', 'quadrupole': 2.07,
              'gn': 1.2878, 'Z': 75, 'spin': 2.5},
    '184Os': {'A': 184, 'abundance': 0.02, 'name': 'osmium', 'symbol': 'Os', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 76, 'spin': 0.0},
    '186Os': {'A': 186, 'abundance': 1.59, 'name': 'osmium', 'symbol': 'Os', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 76, 'spin': 0.0},
    '187Os': {'A': 187, 'abundance': 1.96, 'name': 'osmium', 'symbol': 'Os', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.1311, 'Z': 76, 'spin': 0.5},
    '188Os': {'A': 188, 'abundance': 13.24, 'name': 'osmium', 'symbol': 'Os', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 76, 'spin': 0.0},
    '189Os': {'A': 189, 'abundance': 16.15, 'name': 'osmium', 'symbol': 'Os', 'stability': 'stable',
              'quadrupole': 0.856,
              'gn': 0.488, 'Z': 76, 'spin': 1.5},
    '190Os': {'A': 190, 'abundance': 26.26, 'name': 'osmium', 'symbol': 'Os', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 76, 'spin': 0.0},
    '192Os': {'A': 192, 'abundance': 40.78, 'name': 'osmium', 'symbol': 'Os', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 76, 'spin': 0.0},
    '191Ir': {'A': 191, 'abundance': 37.3, 'name': 'iridium', 'symbol': 'Ir', 'stability': 'stable',
              'quadrupole': 0.816,
              'gn': 0.097, 'Z': 77, 'spin': 1.5},
    '193Ir': {'A': 193, 'abundance': 62.7, 'name': 'iridium', 'symbol': 'Ir', 'stability': 'stable',
              'quadrupole': 0.751,
              'gn': 0.107, 'Z': 77, 'spin': 1.5},
    '190Pt': {'A': 190, 'abundance': 0.014, 'name': 'platinum', 'symbol': 'Pt', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 78, 'spin': 0.0},
    '192Pt': {'A': 192, 'abundance': 0.784, 'name': 'platinum', 'symbol': 'Pt', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 78, 'spin': 0.0},
    '194Pt': {'A': 194, 'abundance': 32.967, 'name': 'platinum', 'symbol': 'Pt', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 78, 'spin': 0.0},
    '195Pt': {'A': 195, 'abundance': 33.832, 'name': 'platinum', 'symbol': 'Pt', 'stability': 'stable', 'quadrupole': 0,
              'gn': 1.219, 'Z': 78, 'spin': 0.5},
    '196Pt': {'A': 196, 'abundance': 25.242, 'name': 'platinum', 'symbol': 'Pt', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 78, 'spin': 0.0},
    '198Pt': {'A': 198, 'abundance': 7.163, 'name': 'platinum', 'symbol': 'Pt', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 78, 'spin': 0.0},
    '197Au': {'A': 197, 'abundance': 100.0, 'name': 'gold', 'symbol': 'Au', 'stability': 'stable', 'quadrupole': 0.547,
              'gn': 0.097968, 'Z': 79, 'spin': 1.5},
    '196Hg': {'A': 196, 'abundance': 0.15, 'name': 'mercury', 'symbol': 'Hg', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 80, 'spin': 0.0},
    '198Hg': {'A': 198, 'abundance': 9.97, 'name': 'mercury', 'symbol': 'Hg', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 80, 'spin': 0.0},
    '199Hg': {'A': 199, 'abundance': 16.87, 'name': 'mercury', 'symbol': 'Hg', 'stability': 'stable', 'quadrupole': 0,
              'gn': 1.01177, 'Z': 80, 'spin': 0.5},
    '200Hg': {'A': 200, 'abundance': 23.1, 'name': 'mercury', 'symbol': 'Hg', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 80, 'spin': 0.0},
    '201Hg': {'A': 201, 'abundance': 13.18, 'name': 'mercury', 'symbol': 'Hg', 'stability': 'stable',
              'quadrupole': 0.386,
              'gn': -0.373483, 'Z': 80, 'spin': 1.5},
    '202Hg': {'A': 202, 'abundance': 29.86, 'name': 'mercury', 'symbol': 'Hg', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 80, 'spin': 0.0},
    '204Hg': {'A': 204, 'abundance': 6.87, 'name': 'mercury', 'symbol': 'Hg', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 80, 'spin': 0.0},
    '203Tl': {'A': 203, 'abundance': 29.524, 'name': 'thallium', 'symbol': 'Tl', 'stability': 'stable', 'quadrupole': 0,
              'gn': 3.244514, 'Z': 81, 'spin': 0.5},
    '204Tl': {'A': 204, 'abundance': 0.0, 'name': 'thallium', 'symbol': 'Tl', 'stability': 'radioactive',
              'quadrupole': nan,
              'gn': 0.0454, 'Z': 81, 'spin': 2.0},
    '205Tl': {'A': 205, 'abundance': 70.476, 'name': 'thallium', 'symbol': 'Tl', 'stability': 'stable', 'quadrupole': 0,
              'gn': 3.2754, 'Z': 81, 'spin': 0.5},
    '204Pb': {'A': 204, 'abundance': 1.4, 'name': 'lead', 'symbol': 'Pb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 82, 'spin': 0.0},
    '206Pb': {'A': 206, 'abundance': 24.1, 'name': 'lead', 'symbol': 'Pb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 82, 'spin': 0.0},
    '207Pb': {'A': 207, 'abundance': 22.1, 'name': 'lead', 'symbol': 'Pb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 1.1748, 'Z': 82, 'spin': 0.5},
    '208Pb': {'A': 208, 'abundance': 52.4, 'name': 'lead', 'symbol': 'Pb', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 82, 'spin': 0.0},
    '207Bi': {'A': 207, 'abundance': 0.0, 'name': 'bismuth', 'symbol': 'Bi', 'stability': 'radioactive',
              'quadrupole': nan,
              'gn': 0.97, 'Z': 83, 'spin': 4.5},
    '209Bi': {'A': 209, 'abundance': 100.0, 'name': 'bismuth', 'symbol': 'Bi', 'stability': 'stable',
              'quadrupole': -0.516,
              'gn': 0.938, 'Z': 83, 'spin': 4.5},
    '209Po': {'A': 209, 'abundance': 0.0, 'name': 'polonium', 'symbol': 'Po', 'stability': 'radioactive',
              'quadrupole': 0,
              'gn': 1.5, 'Z': 84, 'spin': 0.5},
    '227Ac': {'A': 227, 'abundance': 0.0, 'name': 'actinium', 'symbol': 'Ac', 'stability': 'radioactive',
              'quadrupole': 1.7,
              'gn': 0.73, 'Z': 89, 'spin': 1.5},
    '229Th': {'A': 229, 'abundance': 0.0, 'name': 'thorium', 'symbol': 'Th', 'stability': 'radioactive',
              'quadrupole': 4.3,
              'gn': 0.16, 'Z': 90, 'spin': 2.5},
    '232Th': {'A': 232, 'abundance': 100.0, 'name': 'thorium', 'symbol': 'Th', 'stability': 'stable', 'quadrupole': 0,
              'gn': 0.0, 'Z': 90, 'spin': 0.0},
    '234U': {'A': 234, 'abundance': 0.0055, 'name': 'uranium', 'symbol': 'U', 'stability': 'radioactive',
             'quadrupole': 0,
             'gn': 0.0, 'Z': 92, 'spin': 0.0},
    '235U': {'A': 235, 'abundance': 0.72, 'name': 'uranium', 'symbol': 'U', 'stability': 'radioactive',
             'quadrupole': 4.936,
             'gn': -1.09, 'Z': 92, 'spin': 3.5},
    '238U': {'A': 238, 'abundance': 99.2745, 'name': 'uranium', 'symbol': 'U', 'stability': 'radioactive',
             'quadrupole': 0,
             'gn': 0.0, 'Z': 92, 'spin': 0.0},
    'off': {'A': 0, 'abundance': 0, 'name': 'undefined', 'symbol': '_', 'stability': 'stable', 'quadrupole': 0,
            'gn': 0.0, 'Z': 0, 'spin': 0.0},
}


def isotope_list():
    """Helper function returning an ordered list of the isotope names
    
    """
    d = sorted(dict_isotopes.items(), key=lambda (k, v): (v, k))
    k, _v = zip(*d)
    return k


#===============================================================================
# Isotopes
#===============================================================================
class Isotopes(HasTraits):
    """
    This class defines usefull properties of nuclei

    Parameters
    ----------
    nucleus : String, optional, default='1H'
        In the AX form where A is the atomic mass and X the atom symbol

    Attributes
    ----------
    name : string
        Long name of the nucleus
    symbol : string
        Symbol of the nucleus
    A : int
        Atomic mass / Number of nucleons
    Z : int
        Atomic number / number of protons
    spin : float
        The qpin quantum number
    abundance : Float
        The natural abundance in %
    gamma : float
        Gyromagnetic ratio divided by nuclear magneton
    Q : Float
        The electric quadrupole moment in barns [1]
    stability : string
        Stability of the nucleus (radioactive or not)

    Examples
    --------

    How to use it? 

    >>> isotope = Isotopes() 
    >>> isotope.spin
    0.5
    >>> isotope.nucleus = '27Al'
    >>> isotope.name
    'aluminium'
    >>> isotope.symbol
    'Al'
    >>> isotope = Isotopes('off')

    References
    ----------   
    .. [1] Nuclear magnetic moments are taken from Stone, Table of Nuclear
        Magnetic Dipole and Electric Quadrupole Moments, Atomic Data
        and Nuclear Data Tables 90, 75-176 (2005). 
        Nuclear quadrupole moments are taken from P.Pyykko, Mol.Phys.
        99, 1617-1629 (2001) and the 2002 edition of the CRC Handbook
        of Physics and Chemistry (which took it from Pyykko and others).

    """

    nucleus = Str(Desc='Isotope symbol')

    #===========================================================================
    # _get_spin
    #===========================================================================
    spin = Property(Float, Desc="Spin quantum number", depends_on='nucleus')

    def _get_spin(self):
        return dict_isotopes[self.nucleus]['spin']

    #===========================================================================
    # _get_Z
    #===========================================================================
    Z = Property(Float, Desc='Atomic number', depends_on='nucleus')

    def _get_Z(self):
        return dict_isotopes[self.nucleus]['Z']

    #===========================================================================
    # _get_A
    #===========================================================================
    A = Property(Float, Desc='Atomic mass', depends_on='nucleus')

    def _get_A(self):
        return dict_isotopes[self.nucleus]['A']

    #===========================================================================
    # _get_name
    #===========================================================================
    name = Property(Float, Desc='Current name', depends_on='nucleus')

    def _get_name(self):
        return dict_isotopes[self.nucleus]['name']

    #===========================================================================
    # _get_gamma
    #===========================================================================
    gamma = Property(Float,
                     Desc='gyromagnetic ratio divided by nuclear magneton',
                     depends_on='nucleus')

    def _get_gamma(self):
        return dict_isotopes[self.nucleus]['gn']

    #===========================================================================
    # _get_abundance
    #===========================================================================
    abundance = Property(Float, Desc='natural abundance in percent',
                         depends_on='nucleus')

    def _get_abundance(self):
        return dict_isotopes[self.nucleus]['abundance']

    #===========================================================================
    # _get_Q
    #===========================================================================
    Q = Property(Float, Desc='Electric quadrupole moment in barn',
                 depends_on='nucleus')

    def _get_Q(self):
        return dict_isotopes[self.nucleus]['quadrupole']

    #===========================================================================
    # _get_symbol
    #===========================================================================
    symbol = Property(Float, Desc="Symbol", depends_on='nucleus')

    def _get_symbol(self):
        return dict_isotopes[self.nucleus]['symbol']

    def __str__(self):
        return self.nucleus

    #===========================================================================
    # __init__
    #===========================================================================
    def __init__(self, nucleus='1H'):
        """
        Class constructor
        """
        self.nucleus = nucleus


class NmrProperties(HasTraits):
    """
    This class displays the properties of a given nucleus
    
    """
    #: This property return a list of isotopes

    isotope = Enum(isotope_list())
    current = Property(Instance(Isotopes), depends_on='isotope')
    _current = Instance(Isotopes, ())

    def _get_current(self):
        self._current = Isotopes(self.isotope)
        return self._current

if __name__ == '__main__':
    nmr = NmrProperties(isotope='27Al')
    