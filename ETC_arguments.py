import argparse
from sys import path
from numpy.ma import is_masked
from ETC.ETC_import import sourcesdir

# Hacks to avoid exiting the program when there's a parser error
class ArgumentParserError(Exception): pass
def raiseerror(self, message): raise ArgumentParserError(message)

help = 'Run the Exposure Time Calculator.  Outputs are SNR, EXPTIME, wavelength range, and optional plots. '
help += 'The model assumes that signals from 3 image slicer paths are summed for the SNR calculation.'
epilog = 'Example minimum argument set: \n./ETC_main.py G 500 510 SNR 10 -slit .5 -seeing 1 500 -airmass 1 -skymag 21.4 -srcmodel blackbody 6000 -mag 18. -magref AB user'

parser = argparse.ArgumentParser(  # Make printed help text wider
  formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=40) ,description=help ,epilog=epilog)

# Define some extra types to enforce positive values for arguments

def posfloat(value): # require > 0
    fvalue = float(value)
    if fvalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive float" % value)
    return fvalue

def nonegfloat(value): # require >= 0
    fvalue = float(value)
    if fvalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid non-negative float" % value)
    return fvalue

def posint(value): # require > 0
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int" % value)
    return ivalue

SNRparam = parser.add_argument_group('SNR parameters')

from ETC.ETC_config import channels
help = 'Spectrograph channel used for SNR'
parser.add_argument('channel', type=str, choices=channels ,help=help)

help = 'Min and max wavelength (nm) for SNR avg, e.g. "500 510". Will be rounded up to a whole number of bins'
parser.add_argument('wrange', type=posfloat ,nargs=2 ,help=help )

help = 'Fix SNR or EXPTIME and calculate the other'
parser.add_argument('ETCmode', type=str, choices=['SNR','EXPTIME'], help=help)

help = 'Value of the fixed parameter: SNR (dimensionless) or EXPTIME (s)'
parser.add_argument('ETCfixed', type=posfloat, help=help)

help = 'On-chip binning along dispersion (D) and spatial (S) axes'
parser.add_argument('-binning', type=posint, nargs=2, default=[1,1], metavar=('BIN_D','BIN_S') ,help=help)

help = 'Number of spatial pixels on either side of profile center to use for SNR. If none, fit spatial profile'
parser.add_argument('-SNR_pix', type=posint, default=None, help=help)

help = 'Only use flux from the center slit, not side slices'
parser.add_argument('-noslicer', action='store_true', help=help)

help = 'Plot SNR vs. wavelength for the solution'
parser.add_argument('-plotSNR', action='store_true', help=help)

help = 'Make diagnostic plots'
parser.add_argument('-plotdiag', action='store_true', help=help)

obsparam = parser.add_argument_group('REQUIRED Observation conditions')

help = 'Slit width (arcsec)'
obsparam.add_argument('-slit', '-slitwidth', type=posfloat, required=True, help=help)

help = 'Seeing FWHM (arcsec) defined at pivot wavelength (nm)'
obsparam.add_argument('-seeing', type=posfloat, nargs=2, metavar=('SEEING','PIVOT'), required=True, help=help)

help = 'Airmass (dimensionless)'
obsparam.add_argument('-airmass', type=float, required=True, help=help)

help = 'Sky brightness magnitude per arcsec^2 (VEGA, johnson_v)'
obsparam.add_argument('-skymag', type=float, required=True, help=help)

sourceparam_req = parser.add_argument_group('REQUIRED Source parameters')

help = 'Source magnitude (observed at top of atmosphere)'
sourceparam_req.add_argument('-mag', type=float, required=True, help=help)

help = '''Reference system (AB, VEGA) and Johnson filter (UBVRIJK) for source magnitude, e.g. "VEGA V".  
Use FILTER="user" to normalize to the WRANGE input'''
sourceparam_req.add_argument('-magref', type=str, nargs=2, metavar=('SYSTEM','FILTER'), required=True, help=help)

help = 'Astronomical source model.  Examples: "blackbody 5000", "template spiral_001"'
sourceparam_req.add_argument('-srcmodel', nargs='+', required=True, help=help)

sourceparam_add = parser.add_argument_group('Additional source parameters')

help = 'Redshift'
sourceparam_add.add_argument('-z', type=nonegfloat, default=0., help=help)

help = 'Selective Extinction E(B-V); default=0'
sourceparam_add.add_argument('-E_BV', type=float, default=0., help=help)

help = 'Extinction model; default="mwavg" (Diffuse Milky Way, R_V=3.1)'
sourceparam_add.add_argument('-extmodel', type=str, default='mwavg', help=help)

# ETC parameter summary for external modules
etc_args = ['channel', 'wrange', 'SNR'] # Order is important
etc_kwargs = ['slitwidth', 'airmass', 'skymag','seeing', 'mag', 'magref', 'srcmodel']
etc_optkwargs = ['binning', 'SNR_pix', 'z', 'E_BV', 'extmodel']  ### -noslicer takes no argument in ETC command

def formETCcommmand(row):  ### Maybe make command as list not a big string
	'''Form the ETC command line string from an astropy table row'''
	cmd = '%s %s SNR %s ' % tuple([row[k] for k in etc_args])
	cmd_kwargs = [ '-%s %s'%(k,row[k]) for k in etc_kwargs if not is_masked(row[k]) ]

	# These columns are optional, so first check if they exist
	cmd_optkwargs = []
	for k in etc_optkwargs:
		if k in row.colnames:
			if not is_masked(row[k]):
				cmd_optkwargs.append('-%s %s'%(k,row[k]))
	return cmd + ' '.join(cmd_kwargs+cmd_optkwargs)

# Check that inputs are valid and append units where applicable
def check_inputs_add_units(args):

	if len(args.srcmodel) < 2: parser.error('-srcmodel requires at least 2 arguments')

	choices = ['blackbody','template']
	if args.srcmodel[0].lower() not in choices: parser.error('-srcmodel first argument must be in '+str(choices))

	# Check for valid template if using a template model
	if args.srcmodel[0].lower()=='template':

		args.srctemp = args.srcmodel[1]  # copy template name to new attribute

		# Automatically go looking for the template in the sources directory
		from os import walk
		from os.path import join as joinpath
		foundTemplate = False
		for root, dirs, files in walk(sourcesdir):
			if not foundTemplate:
				for name in files:
					if name == args.srctemp+".fits":
						dir_and_name = '/'.join(joinpath(root, name).split('/')[-2:])
						#print("Found template: "+dir_and_name)
						args.srctemp = joinpath(root, name)
						foundTemplate = True
						continue

		if not foundTemplate: parser.error("Could not find source template: "+args.srctemp+'.fits')

	# Check for temperature if using blackbody model
	if args.srcmodel[0].lower()=='blackbody':
		try: args.tempK = posfloat(args.srcmodel[1])  # copy blackbody temperature to new attribute
		except: parser.error("-srcmodel blackbody TEMPK requires TEMPK to be a positive float")

	# Valid mag system
	choices = ['AB','VEGA']
	if args.magref[0].upper() not in choices: parser.error('-magref SYSTEM must be in '+str(choices))

	# Valid mag filter
	choices = 'UBVRIJK'
	if args.magref[1].upper() != 'USER':
		if args.magref[1].upper() not in choices:
			parser.error('-magref FILTER must be "USER" or in '+str(choices))

	# Valid extinction model
	choices = ('lmc30dor', 'lmcavg', 'mwavg', 'mwdense', 'mwrv21', 'mwrv40', 'smcbar', 'xgalsb')
	if args.extmodel not in choices: parser.error('-extmodel must be in '+str(choices))

	# Append units to inputs where applicable
	import astropy.units as u
	args.slit *= u.arcsec
	args.wrange *= u.nm
	if args.ETCmode == 'EXPTIME': args.ETCfixed *= u.s
	args.seeing[0] *= u.arcsec
	args.seeing[1] *= u.nm
	if hasattr(args, 'tempK'): args.tempK*=u.K

	# Check wavelength range is (min, max) and within specified channel
	from ETC.ETC_config import channelRange
	if args.wrange[0] >= args.wrange[1]: parser.error("Wavelength range must be in form [min, max]")
	if args.wrange[0] < channelRange[args.channel][0]: parser.error("Wavelength range not in channel %s"%args.channel)
	if args.wrange[1] > channelRange[args.channel][1]: parser.error("Wavelength range not in channel %s"%args.channel)

