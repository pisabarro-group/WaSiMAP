from __future__ import annotations

import argparse
import sys, traceback

from .wasimap import WaterMapper
from . import __version__

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="wasimap",
        description="""
        WaSiMap water-site mapper CLI. 
        
        - Use on toplevel directory containing MD trajectories and topologies.
        - More instructions at https://github.com/pisabarro-group/wasimap

        Example usage:

        wasimap --gui 
        (Executes analysis on all trajectories in the local dir and opens GUI after processing)

        wasimap --onlygui 
        (Don't process anything. Open GUI for existing results)

        wasimap sim1 
        (executes analysis on trajectory and topology called sim1)
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=0.35,
        help="H-bonding distance cutoff (default: 0.35 nm)",
    )
    parser.add_argument(
        "-r",
        "--persistence",
        type=float,
        default=5.0,
        help="Minimum occupancy to be shown in GUI (default: 5.0). Keep in mind this value is proportional to your total frames",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Launch GUI after scavenging.",
    )
    parser.add_argument(
        "--onlygui",
        action="store_true",
        help="Do not process anything. Show GUI for existing results",
    )

    parser.add_argument(
        "--testdata",
        action="store_true",
        help="Download light MD test data to local dir (approx. 320MB)",
    )
    
    parser.add_argument(
        "--testdata2",
        action="store_true",
        help="Download thorough MD test data to local dir (approx. 15 GB)",
    )

    parser.add_argument(
        "-m"
        "--min-residence-ps",
        type=float,
        default=100.0,
        dest="min_residence_ps",
        help=(
            "Minimum accumulated time (in picoseconds) a water must remain near "
            "a heavy-atom to be considered a resident (default: 100.0 ps). Only these"
            "will be taken into account. Shorter residence times will be discarded."
            "If the trajectory doesn't report the time step size, it defaults to use"
            "20 continuous frames" 
        ),
    )


    parser.add_argument(
        "inputs",
        nargs="*",
        help="Name of an individual simulation (name only, no extension. trajectory and topology files must have the same name).",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        wm = WaterMapper(
            distance_threshold=args.cutoff,
            persistence=args.persistence,
            gui=args.gui,
            onlygui=args.onlygui,
            inputs=args.inputs,
            testdata=args.testdata,
            testdata2=args.testdata2,
            min_residence_ps=args.min_residence_ps,
        )
        #Run main application
        wm.run()
        return 0

    except KeyboardInterrupt:
        return 130
    except Exception as exc:
        tb = traceback.extract_tb(exc.__traceback__)[-1]
        print(f"[wasimap] ERROR: {exc} ({tb.filename}:{tb.lineno})", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())