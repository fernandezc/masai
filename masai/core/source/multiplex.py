"""
****************************************************************
MSM (Multiplex Spam MQMAS) Au program
version: 20141109 (Python version)
****************************************************************
Authors:
C. Fernandez, A. Vicente
LCS (Laboratoire Catalyse et Spectrochimie)
UMR6506 CNRS/ENSICAEN,
Universite de Caen-Basse Normandie
6 bd du Marechal Juin,
14050 Caen,
France
mailto: Christian.fernandez@ensicaen.fr

References:
-----------
(1) N.Ivchenko et al., J.Magn.Reson. 160 (2003) 52-58
(2) Z.Gan, H.T.Kwak, J.Magn.Reson. 168 (2004) 346-351.
(3) N.Malicki, et al., Solid State NMR, 2005
"""

import os

import numpy as np


# local 
import nmrglue as ng


def read_multiplex(path, expno=None, procno=1, coherence=3, acqus_files=None, debug=False):
    """
    coherence: coherence selected during the evolution period
    """
    # ----------------------------------------------------------------
    # Reading 3D SER file
    # ----------------------------------------------------------------
    if expno is not None:
        # this mean that the path is stopped to name
        path = os.path.join(path, str(expno))

    dic, ser = ng.bruker.read(path, acqus_files=acqus_files)  # read_prog

    # Check if is a 3D spectra
    parmode = dic['acqus']['PARMODE']
    if ( parmode != 2 ):
        raise Exception("Not a 3D spectrum!")

    # Read size of the acquisition (F3) dimension
    td3 = dic['acqus']['TD']
    if 256 * (td3 / 256) != td3:
        raise Exception(
            "Wrong accumulated data. Please use a TD in the F3 acquisition dimension which is a multiple of 256")

    # Read size of the phase (F1) dimension
    td1 = dic['acqu3s']['TD']
    if debug: print("td1", td1)

    # Read size of the MQ (F2) dimension
    td2 = dic['acqu2s']['TD']
    if debug: print("td2", td2)

    # define more explicit variables
    nphases = td1 / 2
    nphasespam = 2
    npoints = td3 / 2
    nrows = ser.size / nphasespam / nphases / npoints

    # reshape the data array
    ser = ser.reshape(nrows, nphases, nphasespam, npoints)  # Acquisition order for multiplex 3-1-2
    #print "We use this data shape (nrows, nphases, nphasespam, npoints)=", data.shape
    td3 = dic['acqus']['TD']
    if debug: print("td3", td3)
    npoints = td3 / 2

    #----------------------------------------------------------------
    # Processing of the 3D SER file
    #----------------------------------------------------------------

    data = np.zeros((nrows * 2, npoints), dtype=np.complex128)

    for row in range(nrows):

        # read (1-3) planes
        # => number of separated phases (TD1) * size of FID in char (td3)

        # Process the n phases to separate and recombine the various pathways
        # This lead to the creation of the Cosine and Sine Part (SRH)
        # of the resulting SER file
        C0 = np.zeros(npoints, dtype=np.complex128)
        S0 = np.zeros(npoints, dtype=np.complex128)
        C1 = np.zeros(npoints, dtype=np.complex128)
        S1 = np.zeros(npoints, dtype=np.complex128)

        for ip in range(nphases):
            # Compute the numerical phase phi
            phi = 2.0 * np.pi * float(coherence * ip) / float(nphases)
            c = np.cos(phi)
            s = np.sin(phi)

            # Contribution from phi3 = 0 
            tc = c * ser[row, ip, 0]
            ts = s * ser[row, ip, 0]
            C0 = C0 + tc
            S0 = S0 + ts
            C1 = C1 - ts
            S1 = S1 + tc

            # Contribution from phi3 = 180 
            tc = c * ser[row, ip, 1]
            ts = s * ser[row, ip, 1]
            C0 = C0 - tc
            S0 = S0 - ts
            C1 = C1 - ts
            S1 = S1 + tc

        # End of loop over phases ip

        # Cosine Part (SRH)
        data[2 * row] = (-C1.imag + C0.real) + (C1.real + C0.imag) * 1.0j
        # Sine Part (SRH)
        data[2 * row + 1] = (-S1.imag + S0.real) + (S1.real + S0.imag) * 1.0j

    #reverse
    #data = -data[::-1]
    # End of loop over (1-3) planes (rows)
    dic['acqus']['PARMODE'] = 1

    return dic, data  