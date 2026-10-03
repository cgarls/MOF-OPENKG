"""Command-line entry point for the compact release package."""

import argparse
from .mapping import json_to_tsv


def main(argv=None):
    parser = argparse.ArgumentParser(prog="mofs-openkg")
    sub = parser.add_subparsers(dest="command", required=True)
    convert = sub.add_parser("json-to-tsv", help="convert synthesis JSON records to triples")
    convert.add_argument("input")
    convert.add_argument("output")
    args = parser.parse_args(argv)
    if args.command == "json-to-tsv":
        print(f"wrote {json_to_tsv(args.input, args.output)} triples")


if __name__ == "__main__":
    main()
