import argparse

help = 'Run the Exposure Time Calculator.  Outputs are SNR, EXPTIME, wavelength range, and optional plots. '
help += 'The model assumes that signals from 3 image slicer paths are summed for the SNR calculation.'
#parser = argparse.ArgumentParser(description=help)
epilog = 'Example minimum argument set: \n./main.py G 500 510 SNR 10 -slit .5 -seeing 1 500 -airmass 1 -skymag 21.4 -srcmodel blackbody -tempK 6000 -mag 18. -magref VEGA johnson_v'

parser = argparse.ArgumentParser(  # Make printed help text wider
  formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=40) ,description=help ,epilog=epilog)

SNRparam = parser.add_argument_group('SNR parameters')

from ETC_config import channels
help = 'Spectrograph channel used for SNR'
parser.add_argument('channel', type=str, choices=channels ,help=help)

help = 'Min and max wavelength (nm) for SNR avg, e.g. "500 510". Will be rounded up to a whole number of bins'
parser.add_argument('wrange', type=float ,nargs=2 ,help=help )

help = 'Fix SNR or EXPTIME and calculate the other'
parser.add_argument('ETCmode', type=str, choices=['SNR','EXPTIME'], help=help)

help = 'Value of the fixed parameter: SNR (dimensionless) or EXPTIME (s)'
parser.add_argument('ETCfixed', type=float, help=help)

help = 'On-chip binning along dispersion (D) and spatial (S) axes'
parser.add_argument('-binning', type=int, nargs=2, default=[1,1], metavar=('BIN_D','BIN_S') ,help=help)

help = 'Number of spatial pixels on either side of profile center to use for SNR. If none, fit spatial profile'
parser.add_argument('-SNR_pix', type=int, default=None, help=help)

help = 'Plot SNR vs. wavelength for the solution'
parser.add_argument('-plotSNR', action='store_true', help=help)

help = 'Make diagnostic plots'
parser.add_argument('-plotdiag', action='store_true', help=help)

obsparam = parser.add_argument_group('REQUIRED Observation conditions')

help = 'Slit width (arcsec)'
obsparam.add_argument('-slit', type=float, required=True, help=help)

help = 'Seeing FWHM (arcsec) defined at pivot wavelength (nm)'
obsparam.add_argument('-seeing', type=float, nargs=2, metavar=('SEEING','PIVOT'), required=True, help=help)

help = 'Airmass (dimensionless)'
obsparam.add_argument('-airmass', type=float, required=True, help=help)

help = 'Sky brightness magnitude per arcsec^2 (VEGA, johnson_v)'
obsparam.add_argument('-skymag', type=float, required=True, help=help)

sourceparam_req = parser.add_argument_group('REQUIRED Source parameters')

help = 'Source magnitude (observed)'
sourceparam_req.add_argument('-mag', type=float, required=True, help=help)

help = 'Reference system and filter for source magnitude, e.g. "VEGA johnson_v"'
sourceparam_req.add_argument('-magref', type=str, nargs=2, metavar=('SYSTEM','FILTER'), required=True, help=help)

help = 'Astronomical source model type'
sourceparam_req.add_argument('-srcmodel', type=str, required=True, choices=['blackbody', 'template'], help=help)

sourceparam_add = parser.add_argument_group('Additional source parameters')

help = 'Source template category and filename; REQUIRED with -srcmodel template'
sourceparam_add.add_argument('-srctemp', type=str, nargs=2, metavar=('CAT','FILE'), help=help)

help = 'Blackbody temperature (K); REQUIRED with -srcmodel blackbody'
sourceparam_add.add_argument('-tempK', type=float, help=help)

help = 'Redshift'
sourceparam_add.add_argument('-z', type=float, default=0., help=help)

help = 'Selective Extinction E(B-V); default=0'
sourceparam_add.add_argument('-E_BV', type=float, default=0., help=help)

help = 'Normalized Extinction R_V; default=3.1'
sourceparam_add.add_argument('-R_V', type=float, default=3.1, help=help)
