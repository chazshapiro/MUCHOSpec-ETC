
#TODO: remove asserts
#TODO: rewite LSF convolution using astropy/specutils

#import numpy as np
from numpy import array, arange, pi, hstack, vstack
from pickle import load as pload  # Only needed in point source case
from synphot import SourceSpectrum, SpectralElement  #TODO <<---- SLOW IMPORT!
from synphot.models import Empirical1D, Gaussian1D, Box1D
import astropy.units as u
import synphot.units as uu  # used in main code, not this file
from functools import cache
from copy import deepcopy

from ETC.ETC_config import *

import ETC.path as p
from os import path
ETCdir = path.dirname(p.__file__)
sourcesdir = ETCdir+'/sources/'
#CSVdir = path[0]+'/CSV/'
PSFsum2DFile = ETCdir+'/PSFsum2D.pkl' #pre-tabulated integral of PSF over slit and side slices

# Check config file inputs are valid and make some derived parameters

def lists_are_same(list_1, list_2):
    """ Check if lists are same regardless of order """
    ''' use dict.keys() to check if dicts have same keys'''
    if len(list_1) != len(list_2):
        return False
    return sorted(list_1) == sorted(list_2)

for d in channel_dicts:
    assert lists_are_same(channels ,d.keys()), "Mismatched channel key names"

# Unit equivalence
plate_scale = { k : u.pixel_scale(platescale[k]) for k in channels }

dispersion_scale_nobin={}  # Unit equivalence
for k in channels:
    lambdamin,lambdamax = channelRange[k]
    dispersion_scale_nobin[k]=u.pixel_scale( (lambdamax-lambdamin)/(Npix_dispers[k]*u.pix) )

# Min and max wavelength of all channels together
totalRange = array([channelRange[k].to('nm') for k in channels]) #np.array() removes units
totalRange = [totalRange.min(),totalRange.max()]*u.nm

telescope_Area = (1. - Obscuration**2)*pi*(telescope_D/2.)**2 #Collecting area 


# Function definitions

@cache  # clever, but saves only ms
def makeBinCenters(binspec, chanlist=channels):
    '''Compute the center wavelength of each pixel along the spectral direction, accounting for binning'''
    binCenters={}
    for k in chanlist:
        lambdamin,lambdamax = channelRange[k]
        binCenters[k] = rangeQ(lambdamin,lambdamax, binspec*(lambdamax-lambdamin)/Npix_dispers[k])
        # dispersion_scale[k]=u.pixel_scale( args.binning[0]*(lambdamax-lambdamin)/(Npix_dispers[k]*u.pix) )
    return binCenters

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

def seeingLambda(w ,FWHM ,pivot=500.*u.nm):
    '''Seeing law scaled to wavelength'''
    assert u.get_physical_type(w) == 'length', "w must have units of length"
    #pivot = 500.*u.nm
    return FWHM*(w/pivot)**0.2

def makeSource(args):
    ''' Load the source model, mix with astrophysics, normalize .
    Returns a spectrum object for the source at the top of the atmosphere
    '''
    # Load source spectrum model
    if args.srcmodel[0].lower() == 'template':
        label = args.srctemp.split('/')[-1]  #used for plot labels
        sourceSpectrum = SourceSpectrum.from_file(args.srctemp)
    elif args.srcmodel[0].lower() == 'blackbody':
        from synphot.models import BlackBodyNorm1D
        label = 'blackbody '+str(args.tempK)
        sourceSpectrum = SourceSpectrum(BlackBodyNorm1D, temperature=args.tempK)
    elif args.srcmodel[0].lower() == 'constant':
        from astropy.modeling.models import Const1D
        label = 'constant'
        sourceSpectrum = SourceSpectrum(Const1D)

    else: raise Exception('Invalid source model')

    # Redshift it
    sourceSpectrum.z = args.z

    # Galactic extinction
    # https://pysynphot.readthedocs.io/en/latest/spectrum.html#pysynphot-extinction
    if args.E_BV != 0:
        from synphot import ReddeningLaw
        redlaw = ReddeningLaw.from_extinction_model(args.extmodel)
        sourceSpectrum *= SpectralElement(redlaw.extinction_curve(args.E_BV))

    # Load bandpass for normalization
    if args.magref[1].lower() == 'user':
        # Use the wavelength range from command line
        from synphot.models import Box1D
        norm_band = SpectralElement(Box1D, amplitude=1, x_0=args.wrange.mean(), 
                                    width=(args.wrange[1]-args.wrange[0]) )
    else:
        norm_band = SpectralElement.from_filter('johnson_'+args.magref[1].lower())

    # Normalize source - this is done after all astrophysical adjustments and before the atmosphere
    # So we are fixing the magnitude "at the top of the atmosphere"
    if args.magref[0].upper() == 'AB':
        sourceSpectrum = sourceSpectrum.normalize(args.mag*u.ABmag ,band=norm_band )
                                                  #,wavelengths=sourceSpectrum.waveset)
    elif args.magref[0].upper() == 'VEGA':
        sourceSpectrum = sourceSpectrum.normalize(args.mag*uu.VEGAMAG ,band=norm_band 
                                                  ,vegaspec=SourceSpectrum.from_vega())

    return sourceSpectrum


def Extinction_atm(airmass):
    '''Compute Transmission vs. wavelength modulated by airmass.  Returns bandpass object'''
    '''Is there another parameter we can use to correct atm model using nightly data?'''
    bandpass = LoadCSVSpec(throughputFile_atm)
    bandpass.model.lookup_table = bandpass.model.lookup_table**airmass

    return bandpass

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

from scipy.special import gamma
def sharpess_scalefree(theta, x=0., beta=moffat_beta):
    ''' 1D sharpness of seeing-limited PSF (Moffat) along spatial direction
    x = horizontal offset; changes y projection for e.g. side image slice
    '''

    # Normalization of PSF
    integral1 = 3.14159**.5 * theta**(2*beta) * (theta**2+x**2)**(0.5-beta) * gamma(beta-0.5) / gamma(beta)

    # Integral over PSF squared (sharpness)
    integral2 = 3.14159**.5 * theta**(4*beta) * (theta**2+x**2)**(0.5-2*beta) * gamma(2*beta-0.5) / gamma(2*beta)

    return integral2/integral1**2


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

def slitEfficiency(w ,h ,FWHM ,pivot=500.*u.nm ,optics=None):
    '''Compute fraction of PSF passing through slit and side slices, assuming Moffat PSF
    Return as bandpass objects. 

    optics = Bandpass object for slicer side optics'''

    # Slow function of wavelength so choose 10nm sampling
    lams = rangeQ(totalRange[0],totalRange[1],10*u.nm)

    ts = moffat_theta_factor * seeingLambda(lams ,FWHM ,pivot=pivot)  #specific to Moffat PSF
    #if ts.isscalar: ts=[ts]
        
    centerFrac = evaluate2Dinterp(PSFsum2D, w/ts/2., h/ts/2.)
    totalFrac = evaluate2Dinterp(PSFsum2D, 3.*w/ts/2., h/ts/2.)

    sideFrac = (totalFrac-centerFrac)/2.

    if optics is not None:
        sideFrac *= optics(lams).value  # Bandpass object is dimensionless Quantity
        totalFrac = centerFrac + sideFrac*2

    throughput_slicer = {}
    throughput_slicer['center'] = SpectralElement(Empirical1D, points=lams, lookup_table=centerFrac)
    throughput_slicer['side'] = SpectralElement(Empirical1D, points=lams, lookup_table=sideFrac)
    throughput_slicer['total'] = SpectralElement(Empirical1D, points=lams, lookup_table=totalFrac)

    return throughput_slicer

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

def applySlit(slitw, source_at_slit, sky_at_slit, throughput_slicerOptics, args ,chanlist):
    '''Convert source and sky spectra at slit entrance to spectra at focal plane.
    Applies throughput of slit and slicer, convolves with LSF, computes sharpness parameters
    1/sharpness is effective spatial extent of profile in (possibly binned) pixels 
    Assumes seeing-limited PSF

    INPUTS
    slitw:  slit width (unitful)
    source_at_slit: dict of source spectra for each channel, including effects of atm and all throughputs except slit/slicer
    sky_at_slit:  dict of ky spectra for each channel, including all throughputs except slit/slicer
    throughput_slicerOptics: throughput of slicer optics (bandpass object)
    args:  argparse object with all the command line arguments
    '''

    binCenters = makeBinCenters(args.binning[0], chanlist=chanlist)
    if args.noslicer or args.fastSNR:   slicer_paths = ['center']
    else:                               slicer_paths = ['center','side']

    # Combine slit fractions arrays with Optics to make throughput elements
    throughput_slicer = slitEfficiency(slitw ,slit_h ,args.seeing[0] ,pivot=args.seeing[1] ,optics=throughput_slicerOptics)

    # Compute pixelized spatial profiles for a flat spectrum
    # Multiplying spectra by these profiles "distributes" counts over pixels in spatial direction
    # THIS IS ONLY HALF THE (symmetric) PROFILE, so it is normalized to 0.5
    # profile_slit[k][lightpath] sums to 0.5 in each w bin and each path individually
    # profile_slit[k][lightpath] shape is (Nspatial, Nspectral)

    profile_slit = { k: profileOnDetector(k ,slitw ,args.seeing[0] ,args.seeing[1] ,binCenters[k]
                                            ,spatial_range=None ,bin_spatial=args.binning[1])
                    for k in chanlist }

    sharpness = { k : { s: 
        array([0.5]*len(binCenters[k])) if args.fastSNR  #1/sharpness = [2, 2, 2, 2...]
        else 2*(profile_slit[k][s]**2).sum(0)
        for s in slicer_paths }  for k in chanlist }

    # Multiply source spectrum by all throughputs, atmosphere, slit loss, and convolve with LSF
    # These are the flux densities at the focal plane array (FPA)
    # Side slice throughputs are for a SINGLE side slice

    sourceSpectrumFPA={}  #Total flux summed over all spatial pixels and used slices
    skySpectrumFPA={}     #Flux PER spatial pixel, depends on slice bc of slicer optics

    # background flux/pixel ~ slit_width * spatial_pixel_height; we'll scale sky flux by this later
    for k in chanlist:
        # area of sky projected onto 1 pixel in units of arcsec^2
        bg_pix_area = slitw * (1*u.pix).to('arcsec' ,equivalencies=plate_scale[k])
        bg_pix_area = bg_pix_area.to('arcsec2').value

        sourceSpectrumFPA[k]={}
        skySpectrumFPA[k]={}

        for s in slicer_paths:
            spec = source_at_slit[k] * throughput_slicer[s]
            sourceSpectrumFPA[k][s] = convolveLSF(spec, slitw ,args.seeing[0] ,k ,pivot=args.seeing[1])

            if args.fastSNR:
                # scale signal down to 2 center pixels
                sourceSpectrumFPA[k][s] *= SpectralElement(Empirical1D, points=binCenters[k], lookup_table=2*profile_slit[k][s][0])

            # Doesn't include atmosphere or slitloss for sky flux
            # Scale sky flux by effective area: slit_width*pixel_height
            spec = sky_at_slit[k] * bg_pix_area
            if s == 'side': spec *= throughput_slicerOptics
            skySpectrumFPA[k][s] = convolveLSF(spec, slitw ,args.seeing[0] ,k ,pivot=args.seeing[1])

    return sourceSpectrumFPA, skySpectrumFPA, sharpness

def computeSNR(exptime, slitw, args, SSSfocalplane, allChans=False):
        '''Compute SNR from inputs'''
        '''
        exptime: exposure time (unitful)
        slitw: slit width (unitful)
        args: argparse object containing user inputs
        SSSfocalplane: function that returns source spectra, sky spectra, and sharpness at the FPA
        allChans: boolean
        
        RETURNS:
            allChans --> dictionary of SNR/wavelength bin for all channels
            not allChans --> single SNR value            
        '''
    
        # Only loop over channels we need
        if allChans: chanlist = channels # master list from config
        else: chanlist = (args.channel)  # user's input channel; use tuple not list to allow function caching

        binCenters = makeBinCenters(args.binning[0], chanlist=chanlist)
        if args.noslicer or args.fastSNR:   slicer_paths = ['center']
        else:                               slicer_paths = ['center','side']

        sourceSpectrumFPA, skySpectrumFPA, sharpness = SSSfocalplane(slitw, chanlist)  #Cached

        # Convert signal and background to counts (total per wavelength), including noise
        # flux_unit='count' actually makes Spectrum object return counts/second,
        #   so adjust units so we can scale by unitfull exptime

        signal={}  #Total fluence (counts) summed over all spatial pixels
        bgvar={}   #Variance (counts) PER spatial pixel read
        SNR2={}
        SNR={}

        for k in chanlist:
            signal[k] = {}
            bgvar[k] = {}
            SNR2[k] = {}

            for s in slicer_paths:
                signal[k][s] = sourceSpectrumFPA[k][s](binCenters[k] ,flux_unit='count' ,area=telescope_Area) \
                                /u.s*exptime 

                bgvar[k][s] = skySpectrumFPA[k][s](binCenters[k] ,flux_unit='count' ,area=telescope_Area) \
                                /u.s*exptime
                bgvar[k][s] += darkcurrent[k]*exptime*u.pix 
                bgvar[k][s] *= args.binning[1]  #account for more background per pixel if binning
                bgvar[k][s] += (readnoise[k]*u.pix)**2/u.ct  #read noise per read pixel is unchanged
                                
                SIGNAL = signal[k][s]
                NOISE2 = SIGNAL + bgvar[k][s]/sharpness[k][s]
                SNR2[k][s] = (SIGNAL**2/NOISE2 /u.ct).value  #**.5/u.ct**.5

            SNR[k] = deepcopy(SNR2[k]['center'])
            if not (args.noslicer or args.fastSNR): SNR[k] += 2*SNR2[k]['side']
            SNR[k] = SNR[k]**0.5

        # Return SNR data for each channel
        if allChans: return SNR

        # Return 1 number - the SNR caclulated from user input
        else:
            # Find the bins where the target wavelengths live
            ch = args.channel
            closest_bin_i = [abs(binCenters[ch]-wr).argmin() for wr in args.wrange]
            return (SNR[ch][closest_bin_i[0]:closest_bin_i[1]+1]).mean()

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
