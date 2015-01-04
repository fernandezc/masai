######
MASAÏ 
######

|masai|, a framework for processing, analysing and modelling 
of solid state NMR spectra

:Version: 0.2.alpha.dev

Warning
=======

	|masai| is still experimental and under active development. 
	It is not mature nor stable. 
	Its current design is subject to major changes, reorganizations, bugs and crashes!!!.
	
	**Should I say? be patient as I expect a fully working version with the upcoming 0.2 stable release.**
	
	
What is |masai|?
==================

|masai| is a software written in python and partially in fortran.

It (will) provides a general framework for the processing, analysis and modelling 
of Solid State NMR spectra.

A *GUI* application allows users to make basic processing and fitting of NMR spectra
(for now limited to the Bruker format).

The *API* should also be accessible to be included in other programs.

Can't wait? Here is how you can get it ready in three steps
============================================================

1. Install python distribution from `Enthougth Python Distribution (EPD) <http://www.enthought.com/products/epd.php>`_ 
	
Be sure that you have a Fortran compiler on your machine or get it from `gfortran <http://gcc.gnu.org/wiki/GFortran>`_
	 
2. Then in a terminal, issue this command::

	$ pip install masai

assuming that you have pip installed, or::

	$ easy_install pip
	$ pip install masai
	
To update a previous installed version::

	$ pip install -U masai
	
3. To run the program, just type in a terminal::

	$ masai
	
and you should have access to the GUI interface.
	
Look at the `documentation <http://fernandezc.github.com/masai>`_ for more information about installation.


Documentation
===============

For more information, have a look to the (under work!) web pages at 
`http://fernandezc.github.com/masai <http://fernandezc.github.com/masai>`_

Citing |masai|
===============

When using Masaï for your own work, and before a paper is published 
on this software, you are kindly requested to cite the program this way::

	Christian Fernandez, 
	Masaï: a framework for processing, analysing and modelling of solid state NMR spectra, 
	http://fernandezc.github.com/masai, version 0.2, ENSICAEN/Université de Caen/CNRS, 2012

.. |masai| replace:: **Masaï**