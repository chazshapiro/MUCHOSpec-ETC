# Fixed parameters for this instrument

#TODO: asymmetric plate scale?

import astropy.units as u
from astropy.table import QTable

slit_w_range=[0.2,10.]*u.arcsec
slit_h=60.*u.arcsec
slit_efficiency_max = 0.97  # max fraction of PSF to enclose when optimizing for SNR

telescope_D = 508.*u.cm  #Diameter
Obscuration = 0.3 # As a ratio of M1

moffat_beta=4.765
moffat_theta_factor = 0.5/(2**(1./moffat_beta) - 1.)**.5  # theta = factor*seeingFWHM

# Paths to data for this instrument

CSVdir='/home/developer/Software/ETC/CSV/'
default_waveunit=u.nm  #assume units for all CSV files

skybackground_file = 'Gemini_skybg_50_10.txt'  #placeholder sky model

throughputFile_atm = 'atm-extinction-Palomar.csv'  #dimensionless T  (Flux/Flux_above_atmosphere)
throughputFile_telescope = 'throughput-Palomar-200inch.csv'
throughputFile_slicer = 'throughput_slicersides_temp.csv'

# Channel-wise fixed parameters for this instrument; all channel keys must match

channels=('U','G','R','I')  # use tuple not list to allow function caching

chanConfig=QTable([channels], names=['channel'])
chanConfig.add_index('channel')  # Allows us to specify rows by channel

chanConfig['channelRange']=[[310.,436.], [417.,590.], [561.,794.], [756.,1040.]] * u.nm

# Width of detector (px) in the dispersion direction
chanConfig['Npix_dispers']=(4096, 4096, 4096, 4096)

chanConfig['platescale']=(0.191, 0.191, 0.191, 0.191)*u.arcsec/u.pix

chanConfig['darkcurrent']=(2.0, 2.0, 2.0, 2.0)*u.count/u.pix/u.hr

chanConfig['readnoise']=(4.0, 4.0, 4.0, 4.0)*u.count/u.pix

#LSFFile={}  # Wait for data

#1/2 of FWHM requirement --> sigma
chanConfig['LSFsigma']=(1.0, 1.4, 1.85, 2.25)*u.AA /2/2.35

# Colors for plotting
chanConfig['channelColor']=('blue','green','red','magenta')

chanConfig['throughputFile_spectrograph']=(
    'throughput-NGPS-spectrograph-U.csv',
    'throughput-NGPS-spectrograph-G.csv',
    'throughput-NGPS-spectrograph-R.csv',
    'throughput-NGPS-spectrograph-I.csv'
)

chanConfig['QEFile']=(
    'QE-LBNL-CCD-blue.csv',
    'QE-LBNL-CCD-red.csv',
    'QE-LBNL-CCD-red.csv',
    'QE-LBNL-CCD-red.csv'
)

# Make standalone dicts from the columns in the data table
chanConfigi=chanConfig.loc  # lets us use the index column
for col in chanConfig.colnames:
    globals()[col] = { ch : chanConfigi[ch][col] for ch in channels}

# Create dict from a list
def dictify(vals, keys=channels):
    return { k : v for k,v in zip(keys,vals)}

# class DictObj:
#     def __init__(self, chanConfig:QTable):
#         assert isinstance(chanConfig, QTable)
#         chanConfigi=chanConfig.loc  # lets us use the index column
#         self.channels = chanConfig['channel']
#         for key in chanConfig.colnames:
#             coldict = { ch : chanConfigi[ch][key] for ch in self.channels}
#             setattr(self, key, coldict)
#             #recursive stuff
#             # if isinstance(val, (list, tuple)):
#             #    setattr(self, key, [DictObj(x) if isinstance(x, dict) else x for x in val])
#             # else:
#                # setattr(self, key, DictObj(val) if isinstance(val, dict) else val)
 
# P = DictObj(chanConfig)