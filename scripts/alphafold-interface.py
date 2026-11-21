#!/usr/bin/env python3
"""Evaluate AlphaFold interfaces and write a tall patch TSV."""
import argparse
import csv
import sys

from sieve.alphafold_interface import run
from sieve.alphafold_patch_metrics import ERRORS
from sieve.alphafold_result_format import PatchInput, read_patch_inputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='fold locator, CIF, ZIP, directory, or - for locators/result TSV on stdin')
    parser.add_argument('--output', required=True, help='output TSV path')
    parser.add_argument('--chains', nargs=2, metavar=('CHAIN_A', 'CHAIN_B'))
    parser.add_argument('--domain', action='append', default=[], metavar='CHAIN:DOMAIN',
                        help='explicit domain assignment, repeatable')
    parser.add_argument('--domain-file', '--domains', help='TSV with Sequence Name and Domain Type columns')
    parser.add_argument('--manifest', help='manifest mapping fold locators and chains to sequence names')
    parser.add_argument('--full-data', help='matching confidence JSON for direct CIF input')
    parser.add_argument('--dssp', help='mkdssp executable')
    parser.add_argument('--motif-window', action='append', default=[], metavar='CHAIN:START-END',
                        help='TRAF-binding motif residues, inclusive, repeatable across chains')
    parser.add_argument('--traf-threshold', action='append', default=[], metavar='NAME=VALUE',
                        help='override mixed_area, packing_area, or paired_area')
    parser.add_argument('--math-window', action='append', default=[], metavar='CHAIN:START-END',
                        help='explicit TRAF MATH domain residues, inclusive')
    args = parser.parse_args(argv)
    try:
        requests = read_patch_inputs(sys.stdin) if args.input == '-' else [PatchInput(args.input)]
        return run(requests, args.output, manifest=args.manifest, domain_file=args.domain_file,
                   domains=args.domain, chains=args.chains, full_data=args.full_data, dssp=args.dssp,
                   motif_windows=args.motif_window, traf_thresholds=args.traf_threshold,
                   math_windows=args.math_window)
    except (*ERRORS, csv.Error) as error:
        parser.error(str(error))


if __name__ == '__main__':
    sys.exit(main())
