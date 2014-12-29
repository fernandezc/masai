#===============================================================================
# pulsar makefile
#===============================================================================
# Copyright (C) 2012-2013 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.  
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law 
# and abiding by the rules of distribution of free software.  
# You can  use, modify and/ or redistribute the software under 
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
#===============================================================================
# Makefile for compiling pulsar
#-------------------------------

# You can set these variables from the command line.
DEBUG =

help:
	@echo "Please use  'make <target>' where <target> is one of"
	@echo "  compile        to compile the forpulsar module with f2py (need the forpulsar.pyf file)"
	@echo "  compwtime      to compile the forpulsar module with f2py (need the forpulsar.pyf file) and reporting time information (do not work on macosx)"
	@echo "  verbcompile    to compile (verbose) the forpulsar module with f2py (need the forpulsar.pyf file)"
#	@echo "  pyf        to create the forpulsar.pyf signature file"

compile:
	f2py -c f90/forpulsar.pyf  f90/forpulsar.f90
	@echo "compiling pulsar is finished\n"
	-rm -rf *.mod
	@echo "cleaning mod temporary finished\n"

compile64:
	f2py --f90flags="-m64" -c f90/forpulsar.pyf  f90/forpulsar.f90
	@echo "compiling pulsar is finished\n"
	-rm -rf *.mod
	@echo "cleaning mod temporary finished\n"

compwtime:
	f2py -c -DF2PY_REPORT_ATEXIT f90/forpulsar.pyf  f90/forpulsar.f90
	@echo "\n\ncompiling of forpulsar is finished"
	
verbcompile:
	f2py --verbose -c f90/forpulsar.pyf  f90/forpulsar.f90
	@echo "\n\ncompiling of forpulsar is finished"

#do not use this command as we have manually edited the *.pyf file
#pyf:
#	f2py --verbose -m forpulsar -h f90/forpulsar.pyf f90/forpulsar.f90 --overwrite-signature
#	@echo "forpulsar is ready to be revised"
	
