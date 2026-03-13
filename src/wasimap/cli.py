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

        wasimap --testdata
        (Download MD test data from Zenodo cloud)


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
        type=int,
        default=5,
        help="Minimum percentage of frames a water must reside to be resident (default: 5).",
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
        help="Download MD test data to local dir (approx. 320MB)",
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