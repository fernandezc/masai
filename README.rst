######
MASAÏ 
######

|masai|, processing, analysing and modelling of solid state NMR spectra using IPython notebook.

:Version: 0.2.1

Warning
=======

	|masai| is still experimental and under active development. 
	It is not mature nor stable. 
	Its current design is subject to major changes, reorganizations, bugs and crashes!!!.
	
What is |masai|?
=================

|masai| is a software written in python (and partially in fortran.)

It (will) provides a general *API* for the processing, analysis and modelling
of Solid State NMR spectra.

The *API* is designed to work with IPython but can also be included in other programs.

Can't wait? Here is how you can get it ready in three steps
============================================================

1. Install python distribution from `Enthougth Python Distribution (EPD) <http://www.enthought.com/products/epd.php>`_ 
	
Be sure that you have a Fortran compiler on your machine or get it from `gfortran <http://gcc.gnu.org/wiki/GFortran>`_

2. Then in a terminal, issue this command::

    $ cd masai
    $ python setup.py develop


3. To run the program, open a terminal, cd to the directory where masai/Doc is installed, and then issue the following command::

    $ ipython notebook &


This should open your default Browser (*e.g*. Safari), and open a page `<http://localhost:8888>`_.

Now you can create a notebook or start with the tutorial in the doc folder.

if the notebook was open as above suggested, you can find it here: 
`Tutorial <http://localhost:8888/notebooks/masai/doc/tutorial/Index.ipynb>`_

Look at the `documentation <http://www.numare.fr/masai>`_ for more information about installation.


Documentation
===============

For more information, have a look to the (under work!) web pages at `<http://www.numare.fr/masai>`_

Citing |masai|
===============

When using Masaï for your own work, you are kindly requested to cite the program this way::

	Christian Fernandez, 
	Masaï: a framework for processing, analysing and modelling of solid state NMR spectra, 
	http://www.numare.fr/masai, version 0.2, ENSICAEN/Université de Caen/CNRS, 2015

.. |masai| replace:: **Masaï**