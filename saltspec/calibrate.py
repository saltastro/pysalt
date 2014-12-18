import sys
import numpy as np
import pyfits


class Spectra:

    def __init__(self, warr, farr):
        self.warr = warr
        self.farr = farr

    def interp(self, warr):
        self.farr = np.interp(warr, self.warr, self.farr)
        self.warr = warr


def calibrate(obs_spectra, std_spectra, ext_spectra, airmass, exptime):

    # re-interpt the std_spectra over the same wavelength
    std_spectra.interp(obs_spectra.warr)

    # re-interp the ext_spetra over the sam ewavelength
    ext_spectra.interp(obs_spectra.warr)

    # create the calibration spectra
    cal_spectra = Spectra(obs_spectra.warr, obs_spectra.farr.copy())

    # set up the bandpass
    bandpass = np.diff(obs_spectra.warr).mean()

    # correct for extinction
    cal_spectra.farr = cal_spectra.farr / \
        10 ** (-0.4 * airmass * ext_spectra.farr)

    # correct for the exposure time and calculation the sensitivity curve
    cal_spectra.farr = cal_spectra.farr / exptime / bandpass / std_spectra.farr

    return cal_spectra


def magtoflux(marr, fzero):
    return fzero * 10 ** (-0.4 * marr)


def fnutofwave(warr, farr):
    """Converts farr in ergs/s/cm2/Hz to ergs/s/cm2/A"""
    c = 2.99792458e18  # spped of light in Angstroms/s
    return farr * c / warr ** 2

if __name__ == '__main__':
    obsfile = sys.argv[1]
    stdfile = sys.argv[2]
    extfile = sys.argv[3]
    airmass = float(sys.argv[4])
    exptime = float(sys.argv[5])

    # read in the obsfile
    warr, farr = np.loadtxt(obsfile, usecols=(0, 1), unpack=True)
    obs_spectra = Spectra(warr, farr)

    # read in the std file
    warr, farr = np.loadtxt(stdfile, usecols=(0, 1), unpack=True)
    farr = np.convolve(farr, np.ones(25), mode='same')
    std_spectra = Spectra(warr, farr)

    # read in the extfile
    warr, ext = np.loadtxt(extfile, usecols=(0, 1), unpack=True)
    ext_spectra = Spectra(warr, ext)

    # produce the spectra
    cal_spectra = calibrate(
        obs_spectra,
        std_spectra,
        ext_spectra,
        airmass,
        exptime)

    # plot the results
    from pylab import *
    figure()
    plot(cal_spectra.warr, cal_spectra.farr)
    show()

    # write it out
    outfile = obsfile[:-4] + ".flux"
    fout = open(outfile, 'w')
    for i in range(len(cal_spectra.warr)):
        fout.write(
            '%7.3f %6.5e \n' %
            (cal_spectra.warr[i], cal_spectra.farr[i]))
