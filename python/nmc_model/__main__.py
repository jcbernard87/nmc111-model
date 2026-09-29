"""Run the model described by an input file: python -m nmc_model input.nml"""
import sys

from .inputs import run_file


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: python -m nmc_model input.nml", file=sys.stderr)
        return 2
    r = run_file(argv[0])
    print(f"exit {r.exit_reason} after {r.steps} steps")
    return 0


if __name__ == "__main__":
    sys.exit(main())
