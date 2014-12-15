#!/usr/bin/env python
"""
SPECVIEW is a model to display spectroscopic data


Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0        8 Nov 2009

TODO
----


LIMITATIONS
-----------

"""
from pylab import *


def specview(wave, flux, xtext='$\lambda$', ytext='Flux',
             ls='-', lw=2, color='#FF0000'):
    """Create the display a spectrum

       returns a figure
    """
    f = figure(figsize=(8, 8), dpi=72)
    ax = f.axes([0.1, 0.1, 0.8, 0.8])
    ax.plot(wave, flux, ls=ls, lw=lw, color=color)

    ax.ylabel(ytext)
    ax.xlabel(xtext)
    return f
