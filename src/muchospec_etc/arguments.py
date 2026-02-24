import argparse
from numpy.ma import is_masked
from .config import slitmodes, channels


# Parser Error Handling
class ArgumentParserError(Exception):
    pass


def raiseerror(self, message):
    raise ArgumentParserError(message)


def noQuitETCparser():
    """Change behavior of ArgumentParser to raise exception instead of exit"""
    argparse.ArgumentParser.error = classmethod(raiseerror)


# Parser Definition
description = (
    "Run the Exposure Time Calculator. Outputs are SNR, EXPTIME, "
    "wavelength range, and optional plots. "
    "The model assumes signals from 3 image slicer paths are summed."
)

epilog = (
    "Example minimum argument set:\n"
    "./ETC_main.py G 500 510 SNR 10 "
    "-slit SET .5 -seeing 1 500 -airmass 1 "
    "-skymag 21.4 -mag 18 -magsystem AB -magfilter match"
)

parser = argparse.ArgumentParser(
    formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=40),
    description=description,
    epilog=epilog,
)


# Custom Types
def posfloat(value):
    f = float(value)
    if f <= 0:
        raise argparse.ArgumentTypeError(f"{value} must be > 0")
    return f


def nonegfloat(value):
    f = float(value)
    if f < 0:
        raise argparse.ArgumentTypeError(f"{value} must be >= 0")
    return f


def posint(value):
    i = int(value)
    if i <= 0:
        raise argparse.ArgumentTypeError(f"{value} must be > 0")
    return i


def slitfloat(value):
    from .config import slit_w_range

    f = float(value)
    slitmin, slitmax = slit_w_range.to("arcsec").value
    if f < slitmin or f > slitmax:
        raise argparse.ArgumentTypeError(
            f"slitwidth must be in range {slit_w_range}"
        )
    return f



# Positional Arguments
parser.add_argument(
    "channel",
    type=str,
    choices=channels,
    help="Spectrograph channel used for SNR",
)

parser.add_argument(
    "wrange",
    type=posfloat,
    nargs=2,
    help="Min and max wavelength (nm) for SNR avg",
)

parser.add_argument(
    "ETCmode",
    type=str,
    choices=["SNR", "EXPTIME", "SET"],
    help="Fix SNR or EXPTIME and calculate the other",
)

parser.add_argument(
    "ETCfixed",
    type=posfloat,
    help="Value of fixed parameter",
)


# Optional Arguments
parser.add_argument("-binspect", "-bindisp", type=posint, default=1)
parser.add_argument("-binspat", type=posint, default=1)
parser.add_argument("-noslicer", action="store_true")
parser.add_argument("-fastSNR", type=int)
parser.add_argument("-plotSNR", action="store_true")
parser.add_argument("-plotslit", action="store_true")
parser.add_argument("-timer", action="store_true")
parser.add_argument("-hires", action="store_true")
parser.add_argument("-hires_solve", action="store_true")


# Observation Parameters
obs = parser.add_argument_group("REQUIRED Observation conditions")

obs.add_argument(
    "-slit",
    "-slitwidth",
    required=True,
    nargs="*",
    metavar=("MODE", "VALUE"),
    help=f"Valid modes: {list(slitmodes.keys())}",
)

obs.add_argument(
    "-seeing",
    type=posfloat,
    nargs=2,
    metavar=("SEEING", "PIVOT"),
    required=True,
)

obs.add_argument("-airmass", type=float, required=True)
obs.add_argument("-skymag", type=float, required=True)


# Source Parameters
src_req = parser.add_argument_group("REQUIRED Source parameters")

src_req.add_argument("-mag", type=float, required=True)
src_req.add_argument("-magsystem", choices=["AB", "VEGA", "Vega"], required=True)
src_req.add_argument(
    "-magfilter",
    choices=list("UBVRIJK") + ["user", "USER", "User", "match", "MATCH", "Match"],
    required=True,
)

src_opt = parser.add_argument_group("Additional source parameters")

src_opt.add_argument("-model", nargs="+", default=["constant"])
src_opt.add_argument("-z", type=nonegfloat, default=0.0)
src_opt.add_argument("-E_BV", type=float, default=0.0)
src_opt.add_argument("-extmodel", default="mwavg")
src_opt.add_argument("-extended", type=posint)


# External ETC API Metadata
etc_args = ["channel", "wrange", "exptime"]
etc_kwargs = ["slitwidth", "airmass", "skymag", "seeing", "mag", "magsystem", "magfilter"]
etc_optargs = ["srcmodel"]
etc_optkwargs = ["binspect", "binspat"]



# Command Builder
def formETCcommand(row):
    """Form command string from astropy row"""

    cmd = "%s %s %s " % tuple(row[k] for k in etc_args)

    cmd_kwargs = [
        f"-{k} {row[k]}"
        for k in etc_kwargs
        if not is_masked(row[k])
    ]

    cols_exist = set(etc_optargs) & set(row.keys())
    cmd_optargs = [
        row[k] for k in cols_exist
        if not is_masked(row[k])
    ]

    cols_exist = set(etc_optkwargs) & set(row.keys())
    cmd_optkwargs = [
        f"-{k} {row[k]}"
        for k in cols_exist
        if not is_masked(row[k])
    ]

    return cmd + " ".join(cmd_kwargs + cmd_optkwargs + cmd_optargs)



# Input Validation / Unit Attachment
def check_inputs_add_units(args):
    import astropy.units as u
    from .config import channelRange

    model = args.model[0].lower()
    model_args = {"blackbody": 2, "template": 2, "constant": 1}

    if model not in model_args:
        parser.error(f"-model must be one of {model_args.keys()}")

    if len(args.model) != model_args[model]:
        parser.error(
            f'-model "{model}" requires {model_args[model]} args'
        )

    # Blackbody
    if model == "blackbody":
        args.tempK = posfloat(args.model[1]) * u.K

    # Slit validation
    args.slitmode = args.slit[0].upper()
    if args.slitmode not in slitmodes:
        parser.error(f"-slitwidth MODE must be in {list(slitmodes)}")

    if args.slitmode == "AUTO":
        args.slitmode = "SNR"
        args.slit = slitmodes["AUTO"][0]
    elif len(args.slit) > 1:
        args.slit = posfloat(args.slit[1])
    else:
        parser.error(f"-slitwidth {args.slitmode} requires parameter")

    if args.slitmode == "SET":
        args.slit = slitfloat(args.slit) * u.arcsec

    # Units
    args.wrange *= u.nm
    args.seeing[0] *= u.arcsec
    args.seeing[1] *= u.nm

    if args.ETCmode in ["EXPTIME", "SET"]:
        args.ETCfixed *= u.s
        args.ETCmode = "EXPTIME"

    # Seeing scaling
    args.seeing[0] *= args.airmass ** 0.6

    if args.extended:
        args.seeing[0] = 100 * u.arcsec

    # Channel range validation
    if args.wrange[0] >= args.wrange[1]:
        parser.error("Wavelength range must be [min,max]")

    if args.wrange[0] < channelRange[args.channel][0]:
        parser.error("Range outside channel")

    if args.wrange[1] > channelRange[args.channel][1]:
        parser.error("Range outside channel")
