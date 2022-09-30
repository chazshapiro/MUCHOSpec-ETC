#!/usr/bin/env python
# Author: Chaz Shapiro (2022)
#
# Notes on SNR estimation: https://www.stsci.edu/instruments/wfpc2/Wfpc2_hand/HTML/W2_61.html

# TODOS: K-CORRECTION
#        USER SED, EXTENDED SOURCES (mag/arcsec**2), HOST GALAXY (mag/area)
#        MOON PHASE/POSITION, ZODIACAL LIGHT
#        Unique FILENAMES for plots?
#
# Clean up global-ish variables?
# Error handling
# Should background variances be 2x to acount for sky subtraction?

from ETC.ETC_config import *
from ETC.ETC_arguments import *
from ETC.ETC_import import *
from numpy import array, arange, vstack, log
from scipy import optimize

''' LOAD DATA that doesn't depend on user input '''

# Load sky background; eventually will need to interpolate a larger table
# background units dimensions are not same as other flux units
# File units = u.photon/u.s/u.nm/u.arcsec**2/u.m**2 ~ phot/s/wavelength  VS  nm
skySpec0 = SourceSpectrum.from_file(CSVdir+skybackground_file ,wave_unit='nm') #HARDCODED UNIT
# assumes units = phot/s/cm^2/Angstrom 

# Don't need to match mag reference and filter to source; will use whatever data we have
skyFilter = SpectralElement.from_filter('johnson_v')

# Load telescope throughput
throughput_telescope = LoadCSVSpec(throughputFile_telescope)

# Load throughputs and detector QE for all spectrograph channels
throughput_spectrograph = { k : LoadCSVSpec(throughputFile_spectrograph[k]) for k in channels }
QE = { k : LoadCSVSpec(QEFile[k]) for k in channels }

# Combine spectra with all throughputs except for slit/slicer
TP = { k : throughput_spectrograph[k]*QE[k]*throughput_telescope for k in channels }

# Load throughput for slicer optics (either side of slit)
throughput_slicerOptics = LoadCSVSpec(throughputFile_slicer)

# tt.stop('Setup')

''' MAIN FUNCTION RETURNS DICT OF RESULTS '''
def main(args ,quiet=False ,ETCextras=False):
    if args.timer:
        import timer
        tt = timer.Timer()
        def stopstart_timer(msg):
            tt.stop(msg)
            tt.start()
        tt.start()

    from ETC.ETC_config import channels

    # # Check for valid inputs
    # # If running as command line script, exit gracefully from a problem, otherwise raise an exception
    # try:
    #     check_inputs_add_units(args)
    # except Exception as e:
    #     if __name__ == "__main__": parser.error(e)
    #     else: raise e

    # Only bother with channels being used for SNR;  Saves ~0.3s
    if not args.plotSNR and not args.plotdiag: channels = [args.channel]

    # Load source spectrum model and normalize
    sourceSpectrum = makeSource(args)

    # Normalize the sky; normalization is wrong but proportional to phot/s/wavelength, same as file
    skySpec = skySpec0.normalize(args.skymag*uu.VEGAMAG ,band=skyFilter ,vegaspec=SourceSpectrum.from_vega() ) ###18.4ms
    # new units = VEGAMAG/arcsec^2 since skymag is really mag/arcsec^2

    # Load "throughput" for atmosphere
    throughput_atm = Extinction_atm(args.airmass)

    # Source is modulated by atm, sky is not
    sourceSpec_wTP = { k : sourceSpectrum * throughput_atm * TP[k] for k in channels }
    skySpec_wTP = { k : skySpec * TP[k] for k in channels }

    # Print the bins where the target wavelengths live; won't exactly match input range
    binCenters = makeBinCenters(args.binspect ,chanlist=channels ,wrange=args.wrange)
    # closest_bin_i = [abs(binCenters[args.channel]-wr).argmin() for wr in args.wrange]
    if not quiet:
        print( args.wrange, "-->", binCenters[args.channel][[0,-1]].round(3) )

    if args.noslicer or args.fastSNR:   slicer_paths = ['center']
    else:                               slicer_paths = ['center','side']


    '''  ALL SLIT DEPENDENCE BELOW THIS LINE '''

    # Setup the slicer function; cache saves time in loops where slit width doesn't change
    @cache
    def SSSfocalplane(w, chanlist):

        sourceSpectrumFPA, skySpectrumFPA, sharpness = \
            applySlit(w, sourceSpec_wTP, skySpec_wTP, throughput_slicerOptics, args ,chanlist)

        # Convert signal and background to counts (total per wavelength), including noise
        # flux_unit='count' actually makes Spectrum object return counts/second,
        #   so adjust units so we can scale by unitfull exptime

        signal = {k: {s: sourceSpectrumFPA[k][s](binCenters[k] ,flux_unit='count' ,area=telescope_Area)/u.s
                        for s in slicer_paths} for k in chanlist }

        bgvar = {k: {s: skySpectrumFPA[k][s](binCenters[k] ,flux_unit='count' ,area=telescope_Area)/u.s
                        for s in slicer_paths} for k in chanlist }

        return signal, bgvar, sharpness

    # Solve this function to find slitwidth that gives 97% efficiency (3% slit loss)
    def efffunc(slitw_arcsec):
        eff = slitEfficiency(slitw_arcsec*u.arcsec ,slit_h ,args.seeing[0] ,pivot=args.seeing[1] ,optics=throughput_slicerOptics)
        if args.noslicer: eff = eff['center']
        else:             eff = eff['total']
        return eff(args.wrange).mean() - slit_efficiency_max

    '''  ALL EXPTIME DEPENDENCE BELOW THIS LINE '''

    # Unpack SNR or EXPTIME and set the other to None
    if args.ETCmode == 'SNR':       SNR_target, exptime = (args.ETCfixed, None)
    elif args.ETCmode == 'EXPTIME': SNR_target, exptime = (None, args.ETCfixed)
    else: raise Exception('Invalid ETC mode')

    if args.slitmode=='SET': args.ETCmode += '..SLIT'  # makes this section more readable

    if args.timer: stopstart_timer('Main setup')

    ''' Compute SNR or solve for EXPTIME or SLIT depending on what is fixed '''

    SNRfunc = None

    # EXPTIME and SLIT are fixed; just compute SNR
    if args.ETCmode == 'EXPTIME..SLIT':
        t = exptime
        slitw_result = args.slit
        SNR_result = computeSNR_test(t, args.slit, args, SSSfocalplane).astype('float16')

    # SNR and SLIT are fixed; solve for EXPTIME
    elif args.ETCmode == 'SNR..SLIT':

        slitw_result = args.slit
        SNR_result = SNR_target

        # Find root in log-log space where SNR(t) is ~linear; doesn't save much time but more likely to converge
        def SNRfunc(t_sec):
            return log( computeSNR_test(10**t_sec*u.s, args.slit, args, SSSfocalplane) / SNR_target)
            # return computeSNR(t_sec*u.s, args.slit, args, SSSfocalplane) - SNR_target

        #ans = optimize.root_scalar(SNRfunc ,x0=1 ,x1=100)  ### Very bright stars may not converge here; try log-log
        ans = optimize.root_scalar(SNRfunc ,x0=0 ,x1=3 ,xtol=.01)  ###52ms

        # Check for converged answer
        if ans.converged:
            # t = (ans.root).astype('float16')*u.s
            t = (10.**(ans.root).astype('float16'))*u.s
        else:
            raise RuntimeError('ETC calculation did not converge')

    # EXPTIME is fixed; find SLIT that maximizes SNR
    elif args.ETCmode == 'EXPTIME':

        t = exptime

        ans = optimize.root_scalar(efffunc ,bracket=tuple(slit_w_range.to('arcsec').value) ,x0=args.seeing[0].to('arcsec').value)

        if args.timer: stopstart_timer('Find 97')

        # Check for converged answer
        if ans.converged:
            slitw_max_arcsec = ans.root # unitless, in arcsec
        else:
            raise RuntimeError('Max slit efficiency calculation did not converge')

        # Solve for 'best' slitwidth (maximize SNR)
        def SNRfunc(slitw_arcsec):
            # print(slitw_arcsec)
            ret = computeSNR_test(t, slitw_arcsec*u.arcsec, args, SSSfocalplane)
            return -ret

        ans = optimize.minimize_scalar(SNRfunc ,bounds=slit_w_range.value 
                                        ,method='bounded' ,options={'xatol':0.02}) # absolute tolerance in slitw
        slitw_best_arcsec = ans.x

        if slitw_best_arcsec <= slitw_max_arcsec:
            SNR_result = -ans.fun
            slitw_result = slitw_best_arcsec*u.arcsec
        else:
            SNR_result = -SNRfunc(slitw_max_arcsec)
            slitw_result = slitw_max_arcsec*u.arcsec

        if not quiet:
            print('SLIT=%s  limit=%.2f  best=%.2f'%(slitw_result.round(2), slitw_max_arcsec, slitw_best_arcsec))

    # SNR is fixed; solve for EXPTIME; for each EXPTIME tried, optimize SLIT
    # elif args.ETCmode == 'SNR':
        
    if args.timer: stopstart_timer('Main calculation')

    # RETURN FROM MAIN()

    result = {
        'exptime':t,
        'SNR':SNR_result*u.Unit(1),
        'slitwidth':slitw_result,
        'slit efficiency': (efffunc(slitw_result.to('arcsec').value)+slit_efficiency_max)*u.Unit(1)
        }

    if not quiet:
        print(  '  '.join(['%s=%s' % (k,v.round(2)) for (k,v) in result.items()])  )

    result['wrange'] = args.wrange

    if args.plotSNR:
        result['plotSNR'] = computeSNR_test(t, slitw_result ,args, SSSfocalplane, allChans=True)

    if ETCextras: # return extra functions and data for plotting
        return result, efffunc, SNRfunc
    else:
        return result


def runETC(row ,check=False):
    '''Convert an astropy table row into an ETC command and run the ETC to get an exposure time'''
    cmd = formETCcommand(row)
    args = parser.parse_args(cmd.split())  ### Change this parser.error ? Should be OK since we check ETC cols at the beginning
    check_inputs_add_units(args) # Check for valid inputs; append units

    # If only checking input, we're done
    if check: return True

    result = main(args ,quiet=True)  # Unitful (s)
    t = result['exptime']
    assert isinstance(t,u.Quantity), "Got invalid result from ETC"
    return result

def plotSNR_vs_slit(args, plt):

    w_arcsec = arange(.2,4.5,.2)

    ### HARDCODED UNITS
    fwhm = [makeLSFkernel(wi ,args.seeing[0] ,args.channel ,pivot=args.seeing[1], kernel_upsample=10.)[1] for wi in w_arcsec*u.arcsec]
    fwhm = array(fwhm)*fwhm[0].unit
    R = array(args.wrange).mean()*u.nm/fwhm
    R = R.to(1).value

    slicer = [False, True]
    colors = ['orange','blue']
    labels = ['slit only','with slicer']

    fig, ax1 = plt.subplots(figsize=(10,6))

    for s, c, l in zip(slicer, colors, labels):

        # Compute SNR for all slit widths; find max SNR and where slit encloses 97% of PSF

        args.noslicer = not s
        result, efffunc, SNRfunc = main(args ,quiet=True, ETCextras=True)  #SNRfunc will be SNR vs slit in arcsec

        # Check that main() returned a function we can loop over for plotting
        if SNRfunc is None:
            from sys import exit
            exit('No SNR function for plotting generated in this mode')

        snrs = -array(list(map(SNRfunc,w_arcsec)))
        ans = optimize.root_scalar(efffunc ,bracket=(w_arcsec.min(),w_arcsec.max()) ,x0=args.seeing[0])
        w97 = ans.root

        #Plot SNR for each case with vertical markers
        ax1.plot(w_arcsec, snrs ,color=c ,label=l)
        if s: ll = '97%'
        else: ll = None
        ax1.axvline(w97 ,color=c ,ls=':' ,label=ll)
        if s: ll = 'max SNR'
        else: ll = None
        ax1.axvline(w_arcsec[snrs.argmax()] ,color=c ,ls='--',label=ll)

    # Add labels for SNR plots
    ax1.set_xlabel('Slit width (arcsec)')
    ax1.set_ylabel('SNR')
    ax1.set_title(str(result['wrange'])+'   mag=%s'%args.mag+'   exptime='+str(result['exptime']))
    ax1.legend(loc='center right')

    # Add plot and 2nd y-axis
    ax2 = ax1.twinx()
    ax2.plot(w_arcsec, R, label='R' ,color='green')
    ax2.set_ylabel('R', color='green')
    ax2.tick_params(axis ='y', labelcolor = 'green') 
    ax2.tick_params(axis='y', direction='in', length=6, width=2, colors='g')#, grid_color='r', grid_alpha=0.5)
    ax2.grid(False)

    plt.savefig('SNR_vs_slit.png')
    print('Wrote', 'SNR_vs_slit.png')

    #plt.show()


if __name__ == "__main__":

    args = parser.parse_args()

    # Check for valid inputs; append units
    try: check_inputs_add_units(args)
    except Exception as e: parser.error(e)  # Exit gracefully

    # Run the main program
    result = main(args ,quiet=False)
    print('')

    # Plot SNR vs. wavelength if requested
    if args.plotSNR or args.plotdiag:
        print('Plotting...')

        import matplotlib.pyplot as plt
        #matplotlib.rcParams.update({'font.size': 14})
        from astropy.visualization import astropy_mpl_style, quantity_support
        plt.style.use(astropy_mpl_style)
        quantity_support()

    if args.plotSNR:
        SNR1 = result['plotSNR']

        fig, ax = plt.subplots(figsize=(15,4))

        binCenters = makeBinCenters(args.binspect)
        for k in channels:
            ax.plot(binCenters[k], SNR1[k] ,color=channelColor[k] ,label=k)

        plt.ylabel('SNR / wavelength bin')
        plt.legend()
        ax.axvspan(args.wrange[0], args.wrange[1], alpha=0.2, color='grey') # shade user range
        plt.title('EXPTIME = '+str(result['exptime']))
        plt.savefig('plotSNR.png')
        print('Wrote', 'plotSNR.png')

    if args.plotdiag:
        # Plot SNR and resolution (R) vs. slit width for this target
        plotSNR_vs_slit(args, plt)
