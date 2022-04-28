# Fixed parameters for this instrument

#TODO: hide sourcesdir and PSFsum2DFile with imports (?)
#TODO: asymmetric plate scale?

import astropy.units as u

slit_w_range=[0.2,10.]*u.arcsec
slit_h=60.*u.arcsec

telescope_D = 508.*u.cm  #Diameter
Obscuration = 0.3 # As a ratio of M1

moffat_beta=4.765
moffat_theta_factor = 0.5/(2**(1./moffat_beta) - 1.)**.5  # theta = factor*seeingFWHM

# Paths to data for this instrument

CSVdir='/home/developer/Software/ETC/CSV/'
default_waveunit=u.nm  #assume units for all CSV files

sourcesdir='/home/developer/Software/ETC/sources/'

skybackground_file = 'Gemini_skybg_50_10.txt'  #placeholder sky model

throughputFile_atm = 'atm-extinction-Palomar.csv'  #dimensionless T  (Flux/Flux_above_atmosphere)
throughputFile_telescope = 'throughput-Palomar-200inch.csv'
throughputFile_slicer = 'throughput_slicersides_temp.csv'

PSFsum2DFile = '/home/developer/Software/ETC/PSFsum2D.pkl' #pre-tabulated integral of PSF over slit and side slices

# Channel-wise fixed parameters for this instrument; all channel keys must match

channels=['U','G','R','I']
#channel_dicts = []

# Colors for plotting
channelColor={
    'U':'blue',
    'G':'green',
    'R':'red',
    'I':'magenta'
}

throughputFile_spectrograph={
    'U':'throughput-NGPS-spectrograph-U.csv',
    'G':'throughput-NGPS-spectrograph-G.csv',
    'R':'throughput-NGPS-spectrograph-R.csv',
    'I':'throughput-NGPS-spectrograph-I.csv'
}

QEFile={
    'U':'QE-LBNL-CCD-blue.csv',
    'G':'QE-LBNL-CCD-red.csv',
    'R':'QE-LBNL-CCD-red.csv',
    'I':'QE-LBNL-CCD-red.csv'
}

#LSFFile={}  # Wait for data

#1/2 of FWHM requirement --> sigma
LSFsigma={
    'U': 1.0/2/2.35*u.AA,
    'G': 1.4/2/2.35*u.AA,
    'R': 1.85/2/2.35*u.AA,
    'I': 2.25/2/2.35*u.AA,
}

channelRange={
    'U': [310.,436.]*u.nm,
    'G': [417.,590.]*u.nm,
    'R': [561.,794.]*u.nm,
    'I': [756.,1040.]*u.nm
}

# Width of detector (px) in the dispersion direction
Npix_dispers={
    'U':4096,
    'G':4096,
    'R':4096,
    'I':4096
}

platescale={
    'U':0.191*u.arcsec/u.pix,
    'G':0.191*u.arcsec/u.pix,
    'R':0.191*u.arcsec/u.pix,
    'I':0.191*u.arcsec/u.pix
}

darkcurrent={
    'U':2.0*u.count/u.pix/u.hr,
    'G':2.0*u.count/u.pix/u.hr,
    'R':2.0*u.count/u.pix/u.hr,
    'I':2.0*u.count/u.pix/u.hr
}

readnoise={
    'U':4.0*u.count/u.pix,
    'G':4.0*u.count/u.pix,
    'R':4.0*u.count/u.pix,
    'I':4.0*u.count/u.pix
}

# List all the dictionaries here so we can check their channels match
channel_dicts = [
channelColor
,throughputFile_spectrograph
,QEFile
,LSFsigma
,channelRange
,Npix_dispers
,platescale
,darkcurrent
,readnoise
]