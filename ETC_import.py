
#TODO: remove asserts

#import numpy as np
from numpy import array, arange, pi, hstack, vstack
from pickle import load as pload  # Only needed in point source case
from synphot import SourceSpectrum, SpectralElement  #TODO <<---- SLOW IMPORT!
from synphot.models import Empirical1D, Gaussian1D, Box1D
import astropy.units as u
import synphot.units as uu  # used in main code, not this file

from ETC_config import *

# Check config file inputs are valid and make some derived parameters

def lists_are_same(list_1, list_2):
    """ Check if lists are same regardless of order """
    ''' use dict.keys() to check if dicts have same keys'''
    if len(list_1) != len(list_2):
        return False
    return sorted(list_1) == sorted(list_2)

for d in channel_dicts:
    assert lists_are_same(channels ,d.keys()), "Mismatched channel key names"

plate_scale={} # Unit equivalence
for k in channels: plate_scale[k]=u.pixel_scale(platescale[k])

dispersion_scale_nobin={}  # Unit equivalence
for k in channels:
    lambdamin,lambdamax = channelRange[k]
    dispersion_scale_nobin[k]=u.pixel_scale( (lambdamax-lambdamin)/(Npix_dispers[k]*u.pix) )


# Min and max wavelength of all channels together
totalRange = array([channelRange[k].to('nm') for k in channels]) #np.array() removes units
totalRange = [totalRange.min(),totalRange.max()]*u.nm

telescope_Area = (1. - Obscuration**2)*pi*(telescope_D/2.)**2 #Collecting area 


# Function definitions

def rangeQ(q0, q1, dq=None):
    '''Construct an evenly spaced array using unitful quantities'''
    
    assert isinstance(q0, u.Quantity) and isinstance(q1, u.Quantity), "Inputs must be Quantities"
    
    unit = q0.unit
    if dq is None: dq = 1.*unit
    else: assert isinstance(dq, u.Quantity), "Inputs must be Quantities"

    # Make sure all input units are same dimension
    assert (unit.physical_type == q1.unit.physical_type) and (unit.physical_type == dq.unit.physical_type), \
        "All inputs must have same physical dimensions (different units are OK)"

    v0 = q0.value
    v1 = q1.to_value(unit)    
    dv = dq.to_value(unit)
        
    return arange(v0,v1,dv)*unit


def LoadCSVSpec(filename ,CSVdir=CSVdir):
    '''Helper to load throughput/spectrum element from CSV file'''
    return SpectralElement.from_file(CSVdir+filename ,wave_unit=default_waveunit)

# Something to try...
#import functools
#@functools.lru_cache(maxsize=128)
def seeingLambda(w ,FWHM ,pivot=500.*u.nm):
    '''Seeing law scaled to wavelength'''
    assert u.get_physical_type(w) == 'length', "w must have units of length"
    #pivot = 500.*u.nm
    return FWHM*(w/pivot)**0.2

def Extinction_atm(airmass):
    '''Compute Transmission vs. wavelength modulated by airmass.  Returns bandpass object'''
    '''Is there another parameter we can use to correct atm model using nightly data?'''
    bandpass = LoadCSVSpec(throughputFile_atm)
    bandpass.model.lookup_table = bandpass.model.lookup_table**airmass

    return bandpass

def convolveLSF_OLD(spectrum, LSFsigma ,kernel_upsample=5. ,kernel_range=4.):
    '''Placeholder until we have LSF data. Convolve spectrum at focal plane with LSF'''
    
    assert isinstance(LSFsigma,u.Quantity), "LSF sigma needs units"
    assert isinstance(spectrum ,(SourceSpectrum,SpectralElement)), \
        "Input spectrum must be SourceSpectrum or SpectralElement class"
    
    from scipy.signal import convolve

    # LSF placeholder model is just a gaussian
    LSF = Gaussian1D(amplitude=1. ,mean=0. ,stddev=LSFsigma)
    dlamda = LSFsigma/kernel_upsample
    
    ## KERNEL RANGE MUST HAVE 0 IN THE EXACT CENTER (need odd array length)
    xk = rangeQ(0.*LSFsigma,kernel_range*LSFsigma+dlamda,dlamda)
    xk = hstack((-xk[::-1],xk[1:]))  # Symmetrizes range, e.g. [-2,-1,0,1,2]

    lamdamin, lamdamax = spectrum.waverange
    x = rangeQ(lamdamin ,lamdamax ,dlamda)

    # Dimensionless arrays for use with convolve; BEWARE of boundary behavior
    # convolve(mode=same) matches size of 1st argument
    kernel = LSF(xk).value
    spec2 = convolve(spectrum(x).value, kernel, mode='same', method='auto')/kernel.sum()

    # Convert array to spectrum class; Model will be Empirical1D even if input was compound model
    if isinstance(spectrum ,SourceSpectrum):
        newspec = SourceSpectrum(Empirical1D, points=x, lookup_table=spec2*spectrum(1).unit, keep_neg=True)
    elif isinstance(spectrum ,SpectralElement):
        newspec = SpectralElement(Empirical1D, points=x, lookup_table=spec2, keep_neg=True)
    else:
        raise Exception("Unsupported input spectrum class")
        
    return newspec

def convolveLSF(spectrum ,slit_w ,seeing ,ch ,kernel_upsample=5. ,kernel_range_factor=4. ,pivot=500*u.nm):
    '''Placeholder until we have LSF data. Convolve spectrum at focal plane with LSF'''
    '''
    Approximates seeing for each channel as Gaussian with scale at channel center wavelength
    Final kernel is INSTRUMENT * (SLIT X SEEING)  (*=convolve)
    TODO: vary LSF in dispersion and/or spatial directions

    spectrum: synphot Spectrum object
    slit_w: slit width (angular size on sky)
    seeing: FWHM of seeing profile at pivot
    ch: channel name
    kernel_upsample: factor by which to upsample kernel scale
    kernel_range_factor: factor by which to extend range of kernel sampling
    pivot: wavelength where seeing FWHM is defined

    RETURNS: Spectrum object (Emprical1D) after LSF convolution
    '''
    
    assert isinstance(slit_w,u.Quantity), "slit_w needs units"
    assert isinstance(spectrum ,(SourceSpectrum,SpectralElement)), \
        "Input spectrum must be SourceSpectrum or SpectralElement class"
    
    from scipy.signal import convolve

    # Set seeing to that of center wavelength of channel
    # TODO: let seeing vary across channel?
    midlam = channelRange[ch].mean()
    sigma_seeing = seeingLambda(midlam ,seeing ,pivot=pivot)/2.35

    # Convert seeing and slit scales to wavelength
    # LSFsigma is already in wavelength units
    sigma_seeing_lam = sigma_seeing.to('pix', equivalencies=plate_scale[ch]) \
                                    .to(u.AA ,equivalencies=dispersion_scale_nobin[ch])
    slit_w_lam = slit_w.to('pix', equivalencies=plate_scale[ch]) \
                                    .to(u.AA ,equivalencies=dispersion_scale_nobin[ch])

    # Sample kernel out to slit width plus extra for instrument PSF
    # Not ideal if slit width >> seeing but should work
    kernel_range = slit_w_lam/2 + LSFsigma[ch]*kernel_range_factor

    # MIN of seeing and slit sets the slit scale; MAX of slit and instument sets total scale
    dlambda = max(min(sigma_seeing_lam, slit_w_lam/2.) ,LSFsigma[ch])/kernel_upsample

    # Make wavelength array for sampling kernel
    ## KERNEL RANGE MUST HAVE 0 IN THE EXACT CENTER (need odd array length)
    xk = rangeQ(0.*dlambda,kernel_range+dlambda,dlambda)
    xk = hstack((-xk[::-1],xk[1:]))  # Symmetrizes range, e.g. [-2,-1,0,1,2]

    # Make wavelength array for sampling spectrum; same spacing, larger range
    x = rangeQ(channelRange[ch][0] ,channelRange[ch][1] ,dlambda)

    # Unitless arrays for use with convolve; BEWARE of boundary behavior
    # convolve(mode=same) matches size of 1st argument
    slitLSF = Gaussian1D(stddev=sigma_seeing_lam)*Box1D(width=slit_w_lam) # Multiply seeing profile with slit "tophat"
    instLSF = Gaussian1D(stddev=LSFsigma[ch])
    kernel1 = slitLSF(xk).value
    kernel2 = instLSF(xk).value

    # Total kernel
    kernel = convolve(kernel1, kernel2, mode='same', method='auto')

    # Convolved spectrum as unitless array
    spec2 = convolve(spectrum(x).value, kernel, mode='same', method='auto')/kernel.sum()

    # Convert array to spectrum class; Model will be Empirical1D even if input was compound model
    if isinstance(spectrum ,SourceSpectrum):
        newspec = SourceSpectrum(Empirical1D, points=x, lookup_table=spec2*spectrum(1).unit, keep_neg=True)
    elif isinstance(spectrum ,SpectralElement):
        newspec = SpectralElement(Empirical1D, points=x, lookup_table=spec2, keep_neg=True)
    else:
        raise Exception("Unsupported input spectrum class")
        
    return newspec

def Moffat_scalefree(x,y ,beta=moffat_beta):
    '''Moffat PSF profile; x and y are dimensionless; not normalized here - we do that after tabulating'''
    return (1.+x**2+y**2)**(-beta)

# TABULATED IN PSF-profile-scratch.ipynb
PSFsum2D = pload(open(PSFsum2DFile ,'rb'))
from scipy.interpolate import dfitpack

def evaluate2Dinterp(f, x, y):
    '''Trick for quickly evaluating 2D interpolation on (x,y) pairs, not a 2D grid'''
    # https://stackoverflow.com/questions/47087109/evaluate-the-output-from-scipy-2d-interpolation-along-a-curve
    return dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], x, y)[0]

def slitFractions(lam, w ,h ,FWHM ,pivot=500.*u.nm):
    '''Compute fraction of PSF passing through slit and side slices. Assumes Moffat PSF'''
    
    ts = moffat_theta_factor * seeingLambda(lam ,FWHM ,pivot=pivot)  #specific to Moffat PSF
    #if ts.isscalar: ts=[ts]
        
    centerFrac = evaluate2Dinterp(PSFsum2D, w/ts/2., h/ts/2.)
    totalFrac = evaluate2Dinterp(PSFsum2D, 3.*w/ts/2., h/ts/2.)

    sideFrac = (totalFrac-centerFrac)/2.
    return {'total':totalFrac, 'center':centerFrac, 'side':sideFrac}

def profileOnDetector(channel ,slit_w ,seeing ,pivot ,lams ,spatial_range=None ,bin_spatial=1):
    ''' USE TABULATED MOFFAT INTEGRAL TO BACK OUT PIXELIZED SPATIAL PROFILES '''
    '''

    channel_plate_scale: Astropy unit equivalence between arcsec/pixels.  See plate_scale
    lams: wavelengths over which to sample the spatial profile
    seeing:  FWHM of seeing profile at wavelength pivot
    pivot:  pivot wavelength at which seeing is defined
    slit_w:  width of slit in arcsec (user selected)
    spatial_range:  how far to compute profile in the spatial direction
    bin_spatial:  step size (number of spatial pixels) to sample spatial profile; simulates CCD binning

    returns: profile_slit{'center' , side' } - for each key, a numpy array of the pixelized slit profile with shape (Npix,Nlambda)
             Npix is the number of pixels in the spatial direction counting half the profile from the center out
             Output is normalized to total flux entering slit section for each lambda;
             spatial direction should sum to 0.5 since we return only 1/2 the symmetric profile
    '''

    # # Define wavelengths over which to sample the spatial profile
    # # Profile is a slow function of wavelength; lambda_step_max is a lower bound for large ranges
    # # For ranges < lambda_step_max, use half of range as step
    # lambda_step_max = 10*u.nm  #maximum sampling of wavelength range
    # lambda_step = min(lambda_step_max, (lambda_range.max()-lambda_range.min())/2.)
    # lams = rangeQ(lambda_range[0],lambda_range[1]+lambda_step/2.,lambda_step)  # pad range max to ensure max is included

    if spatial_range is None: spatial_range = 3*seeing
    #lams = binCenters[channel]

    # Step through integer pixel widths from 0 (profile center) to some max spatial height
    assert bin_spatial - int(bin_spatial) == 0, "binning must be an integer"
    dx = (1*u.pix).to(u.arcsec ,equivalencies=plate_scale[channel]) * bin_spatial
    x=rangeQ(0*u.arcsec ,spatial_range ,dx)

    # For each wavelength, integrate over the slit width and from origin to each pixel boundary
    # slitFractions() halves w and h in its integrals so double the pixel height argument
    # This returns Npix+1 dicts containing arrays of length Nlambda
    PSFsums=[slitFractions(lams, slit_w, 2*xi ,seeing) for xi in x]  #SLOW?

    # Normalize to sum over full slit
    # This sums only from 0 to slit height, so double it below to get full slit value
    profile_norm=slitFractions(lams, slit_w, slit_h ,seeing)  #shape is (Nlambda)

    # Subtract sums at adjacent pixel boundaries to get sum in each pixel
    profile_slit={}
    for k in ['center','side']:

        #profile_slit[k]=array([(pfs1-pfs0)/profile_norm[k] for pfs0, pfs1 in zip(PSFsums[:-1][k],PSFsums[1:][k]) ]) 
        profile_slit[k]=array([(PSFsums[i+1][k]-PSFsums[i][k])/(2*profile_norm[k]) for i in range(len(x)-1)])
        # NB doubled normalization cf. above note

    # shapes are (Npix, Nlambda)

    return profile_slit

def plotAllChannels(spec ,lambda_range=None ,binned=False ,spec_allchan=None ,binCenters=None):
    '''
    spec: dictionary with a spectrum for each channel
    lambda_range: shade a wavelength range of interest
    binned:  if true, calculate flux in counts; otherwise use flux density
    spec_allchan: multiply spectra in all channels by this spectrum (bandpass)
    x: wavelength values to use for plot (good to use binCenters)
    '''

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12,4))
    for k in spec.keys():
        if binCenters: x = binCenters[k]
        else: x = spec[k].waveset

        y = spec[k]
        if spec_allchan is not None: y *= spec_allchan

        if binned:
            y = y(x ,flux_unit='count', area=telescope_Area)
            ax.plot(x ,y ,drawstyle='steps-mid', label=k ,color=channelColor[k])
        else:
            y = y(x)
            ax.plot(x, y ,label=k ,color=channelColor[k])

    ax.set_ylabel('Flux (%s)' % y[0].unit ) #evaluate the spectrum once to get units
    ax.set_xlim(totalRange)
    ax.legend()    
    
    if lambda_range is not None:
        ax.axvspan(lambda_range[0], lambda_range[1], alpha=0.2, color='grey') # shade user range
