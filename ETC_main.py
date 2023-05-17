#!/usr/bin/python3
# Author: Chaz Shapiro (2022)
#
# Notes on SNR estimation: https://www.stsci.edu/instruments/wfpc2/Wfpc2_hand/HTML/W2_61.html

# TODOS: 
### LSF is incorrect for side slices
### K-CORRECTION
### USER SED, HOST GALAXY (mag/area)
### Unique FILENAMES for plots?
#
# Clean up global-ish variables?
# Should background variances be 2x to acount for sky subtraction?

from ETC.ETC_config import *
from ETC.ETC_arguments import *
from ETC.ETC_import import *
from numpy import array, arange, vstack, log, where
from scipy import optimize

''' LOAD DATA that doesn't depend on user input '''

# Load sky default background; used if rubin_sim spectrum is not provided
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
def main(args ,quiet=False ,ETCextras=False ,plotSNR=False ,plotdiag=False, skyspec=None):
    '''
    skyspec = 2xN array tabulating a sky background spectrum; units = (nm , ergs/s/cm^2/Å )
              Generated in OTM using SkyModel.return_wave_spec() from rubin_sim package
              This currently overrides args.skymag which is still required by argparse. 
    '''

    from ETC.ETC_config import channels

    # Only bother with channels being used for SNR to save time
    if plotSNR or plotdiag:
        args.wrange_ = None
    else:
        channels = [args.channel]
        args.wrange_ = args.wrange

    binCenters = makeBinCenters(args.binspect ,chanlist=channels ,wrange=args.wrange_)
    args.binCenters_ = binCenters

    # Print the bins where the target wavelengths live; won't exactly match input range
    if not quiet:
        print( args.wrange, "-->", binCenters[args.channel][[0,-1]].round(3) )

    if args.noslicer or args.fastSNR:   slicer_paths = ['center']
    else:                               slicer_paths = ['center','side']

    # Load source spectrum model and normalize
    sourceSpectrum = makeSource(args)

    if skyspec: 
        skySpec = SourceSpectrum(Empirical1D, keep_neg=True, 
                                  points=skyspec[0]*u.nm, lookup_table=skyspec[1][0]*uu.FLAM) # (nm , ergs/s/cm^2/Å )

    else:
    # Normalize the sky; normalization is wrong but proportional to phot/s/wavelength, same as file
    # skySpec = skySpec0.normalize(args.skymag*uu.VEGAMAG ,band=skyFilter ,vegaspec=vegaspec ) 
        skySpec = skySpec0*1.575e-17 * 10.**((21.4-args.skymag)/2.5) ### HARDCODE NORMALIZED TO VEGA mag 21.4 JOHNSON V
    # new units = VEGAMAG/arcsec^2 since skymag is really mag/arcsec^2

    # Load "throughput" for atmosphere
    throughput_atm = Extinction_atm(args.airmass)

    # Source is modulated by atm, sky is not
    if args.hires:
        sourceSpec_wTP = { k : sourceSpectrum * throughput_atm * TP[k] for k in channels }
        skySpec_wTP = { k : skySpec * TP[k] for k in channels }
    else:
        # Replace analytic models with lookup tables  ### just as fast to make tables hi-res??
        sourceSpec_wTP = { k : make_empirical(sourceSpectrum*throughput_atm*TP[k], binCenters[k]) for k in channels }
        skySpec_wTP = { k : make_empirical(skySpec*TP[k], binCenters[k]) for k in channels }


    '''  ALL SLIT DEPENDENCE BELOW THIS LINE '''

    # Setup the slicer function; cache saves time in loops where slit width doesn't change
    @cache
    def SSSfocalplane(w):

        sourceSpectrumFPA, skySpectrumFPA, sharpness = \
            applySlit(w, sourceSpec_wTP, skySpec_wTP, throughput_slicerOptics, args ,channels)

        # Convert signal and background to counts (total per wavelength), including noise
        # flux_unit='count' actually makes Spectrum object return counts/second,
        #   so adjust units so we can scale by unitfull exptime

        if args.hires:
            # Integrate spectrum over each bin
            signal = {k: {s: sourceSpectrumFPA[k][s](binCenters[k] ,flux_unit='count' ,area=telescope_Area)/u.s
                           for s in slicer_paths} for k in channels }

            bgvar = {k: {s: skySpectrumFPA[k][s](binCenters[k] ,flux_unit='count' ,area=telescope_Area)/u.s
                           for s in slicer_paths} for k in channels }

        else:
            # Slightly faster approximation for smooth functions; assume spectrum is constant across bins
            # binspect is the width of the bin in pixels
            signal = {k: {s: (sourceSpectrumFPA[k][s](binCenters[k])*telescope_Area/u.ph*u.ct
                            *(args.binspect*u.pix).to('nm',  equivalencies=dispersion_scale_nobin[k])).decompose()
                            for s in slicer_paths} for k in channels }

            bgvar = {k: {s: (skySpectrumFPA[k][s](binCenters[k])*telescope_Area/u.ph*u.ct
                            *(args.binspect*u.pix).to('nm',  equivalencies=dispersion_scale_nobin[k])).decompose()
                            for s in slicer_paths} for k in channels }

        return signal, bgvar, sharpness

    # Solve this function to find slitwidth that gives 97% efficiency (3% slit loss)
    def efffunc(slitw_arcsec):
        eff = slitEfficiency(slitw_arcsec*u.arcsec ,slit_h ,args.seeing[0] ,pivot=args.seeing[1] ,optics=throughput_slicerOptics)
        if args.noslicer: eff = eff['center']
        else:             eff = eff['total']
        return eff(args.wrange).mean() - slit_efficiency_max

    '''  ALL EXPTIME DEPENDENCE BELOW THIS LINE '''

    # Unpack SNR or EXPTIME and set the other to None
    ETCmode = args.ETCmode
    if ETCmode == 'SNR':       SNR_target, exptime = (args.ETCfixed, None)
    elif ETCmode == 'EXPTIME': SNR_target, exptime = (None, args.ETCfixed)
    else: raise Exception('Invalid ETC mode')

    if args.slitmode=='SET': ETCmode += '..SLIT'  # makes this section more readable

    ''' Compute SNR or solve for EXPTIME or SLIT depending on what is fixed '''

    SNRfunc = None

    # EXPTIME and SLIT are fixed; just compute SNR
    if ETCmode == 'EXPTIME..SLIT':
        t = exptime
        slitw_result = args.slit
        SNR_result = computeSNR(t, args.slit, args, SSSfocalplane)

    # SNR and SLIT are fixed; solve for EXPTIME
    elif ETCmode == 'SNR..SLIT':

        slitw_result = args.slit
        SNR_result = SNR_target

        # Find root in log-log space where SNR(t) is ~linear; doesn't save much time but more likely to converge
        def SNRfunc(t_sec):
            return log( computeSNR(10**t_sec*u.s, args.slit, args, SSSfocalplane) / SNR_target)
            # return computeSNR(t_sec*u.s, args.slit, args, SSSfocalplane) - SNR_target

        #ans = optimize.root_scalar(SNRfunc ,x0=1 ,x1=100)  # Very bright stars may not converge here; try log-log
        ans = optimize.root_scalar(SNRfunc ,x0=0 ,x1=3 ,xtol=.01)

        # Check for converged answer
        if ans.converged:
            # t = (ans.root).astype('float16')*u.s
            t = (10.**(ans.root).astype('float16'))*u.s
        else:
            raise RuntimeError('ETC calculation did not converge')

    # EXPTIME is fixed; find SLIT that maximizes SNR
    elif ETCmode == 'EXPTIME':

        t = exptime

        # Find "max" slit width containing 97% PSF
        # ans = optimize.root_scalar(efffunc ,bracket=tuple(slit_w_range.to('arcsec').value) ,x0=args.seeing[0].to('arcsec').value)

        # Check for converged answer
        # if ans.converged:
        #     slitw_max_arcsec = ans.root # unitless, in arcsec
        # else:
        #     raise RuntimeError('Max slit efficiency calculation did not converge')

        # Solve for 'best' slitwidth (maximize SNR)
        def SNRfunc(slitw_arcsec):
            # print(slitw_arcsec)
            ret = computeSNR(t, slitw_arcsec*u.arcsec, args, SSSfocalplane)
            return -ret

        # When optimizing SNR, don't look beyond max slit width set by 97%
        slitbound = slit_w_range.to('arcsec').value
        # if slitw_max_arcsec < slitbound[1]: slitbound[1] = slitw_max_arcsec

        ans = optimize.minimize_scalar(SNRfunc ,bounds=slitbound 
                                        ,method='bounded' ,options={'xatol':0.05}) # absolute tolerance in slitw

        SNR_result = -ans.fun
        slitw_result = ans.x*u.arcsec
        slitw_max_arcsec = slitbound[-1]

        if not quiet:
            print('SLIT=%s  limit=%.2f  best=%.2f'%(slitw_result.round(2), slitw_max_arcsec, slitw_result.value))

    # SNR is fixed; solve for EXPTIME; for each EXPTIME tried, optimize SLIT
    # elif args.ETCmode == 'SNR':

    # Compute final resolution R = lambda/fwhm;  approximate fwhm at center of channel, lambda at center of wrange
    # For FWHM, use kernel result or spectral bin width, whichever is larger
    kernel, kernel_fwhm, kernel_dlambda = makeLSFkernel(slitw_result ,args.seeing[0] ,args.channel ,pivot=args.seeing[1])
    # kernel_fwhm=kernel_fwhm[0]
    binres = (args.binspect*u.pix).to('nm',  equivalencies=dispersion_scale_nobin[args.channel])
    Res = (args.wrange.mean()/max(binres, kernel_fwhm)).to(1)
        
    # RETURN FROM MAIN()

    result = {
        'exptime':t,
        'SNR':SNR_result*u.Unit(1),
        'slitwidth':slitw_result,
        'slit efficiency': (efffunc(slitw_result.to('arcsec').value)+slit_efficiency_max)*u.Unit(1),
        'resolution':Res
        }

    if not quiet:
        print(  '  '.join(['%s=%s' % (k,v.round(2)) for (k,v) in result.items()])  )

    result['wrange'] = args.wrange

    if args.plotSNR:
        result['plotSNR'] = computeSNR(t, slitw_result ,args, SSSfocalplane, allChans=True)

    if ETCextras: # return extra functions and data for plotting
        return result, efffunc, SNRfunc
    else:
        return result


def runETC(row ,check=False, skyspec=None):
    '''Convert an astropy table row into an ETC command and run the ETC to get an exposure time'''
    cmd = formETCcommand(row)
    args = parser.parse_args(cmd.split())  ### Change this parser.error ? Should be OK since we check ETC cols at the beginning
    check_inputs_add_units(args) # Check for valid inputs; append units

    # skyspec should be a tuple of arrays:  (Nlambda, (Nspec, Nlambda)) where Nspec > 1 for multiple filters
    if skyspec:
        try:
            assert len(skyspec) == 2
            for ss in skyspec[1]:
                assert len(skyspec[0]) == len(ss)
        except:
            raise Exception('skyspec must be a 2-item tuple; 0: wavelength samples, 1: N arrays of spectrum values ')

    # If only checking input, we're done
    if check: return True

    result = main(args ,quiet=True, skyspec=skyspec)  # Unitful (s)
    t = result['exptime']
    assert isinstance(t,u.Quantity), "Got invalid result from ETC"
    return result




if __name__ == "__main__":

    args = parser.parse_args()

    # Check for valid inputs; append units
    try: check_inputs_add_units(args)
    except Exception as e: parser.error(e)  # Exit gracefully

    if args.timer:
        from timer import Timer
        tt = Timer()
        tt.start()

    # Run the main program
    result = main(args ,quiet=False)
    print('')

    if args.timer: tt.stop()

    # Plot SNR vs. wavelength if requested
    if args.plotSNR or args.plotdiag:
        print('Plotting...')

        from ETC.ETC_plots import *
        import matplotlib.pyplot as plt
        #matplotlib.rcParams.update({'font.size': 14})
        from astropy.visualization import astropy_mpl_style, quantity_support
        plt.style.use(astropy_mpl_style)
        quantity_support()

    if args.plotSNR:
        args.ETCmode = 'EXPTIME'
        args.ETCfixed = result['exptime']
        args.slitmode = 'SET'
        args.slit = result['slitwidth']
        result_plot = main(args ,quiet=True ,plotSNR=True)
        SNR1 = result_plot['plotSNR']

        fig, ax = plt.subplots(figsize=(15,4))

        binCenters = makeBinCenters(args.binspect)
        for k in channels:
            ax.plot(binCenters[k], SNR1[k] ,color=channelColor[k] ,label=k)

        plt.ylabel('SNR / wavelength bin')
        plt.legend()
        ax.axvspan(args.wrange[0], args.wrange[1], alpha=0.2, color='grey') # shade user range
        plt.title('EXPTIME = '+str(result['exptime'].round(3)))
        plt.savefig('plotSNR.png')
        print('Wrote', 'plotSNR.png')

    if args.plotdiag:
        # Plot SNR and resolution (R) vs. slit width for this target
        plotSNR_vs_slit_mags(args, plt)
