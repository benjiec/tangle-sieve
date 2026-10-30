#!/usr/bin/env python3
"""Classify death-family faces on shared geometric patches using explicit TSV rules."""
import argparse
import csv
import sys
from sieve.alphafold_death_faces import run, ERRORS
from sieve.alphafold_pdockq2_patches import Options
from sieve.alphafold_result_format import PatchInput, read_patch_inputs


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',help='fold locator, CIF, ZIP, directory, or - for locators/result rows on stdin')
    for flag in ('manifest','domains','rules','output'):
        parser.add_argument('--'+flag,required=True)
    parser.add_argument('--full-data')
    parser.add_argument('--dssp')
    parser.add_argument('--chains',nargs=2)
    parser.add_argument('--distance-cutoff',type=float,default=8.0)
    parser.add_argument('--max-bridge-residues',type=int,default=1)
    parser.add_argument('--min-residues-per-chain',type=int,default=3)
    parser.add_argument('--min-contacts',type=int,default=5)
    parser.add_argument('--max-separator-checks',type=int,default=1000000)
    args=parser.parse_args(argv)
    try:
        requests=list(read_patch_inputs(sys.stdin)) if args.input=='-' else [PatchInput(args.input)]
        return run(requests,args.manifest,args.domains,args.rules,args.output,
            options=Options(args.distance_cutoff,args.max_bridge_residues,args.min_residues_per_chain,
                            args.min_contacts,args.max_separator_checks),
            chains=args.chains,full_data=args.full_data,dssp=args.dssp)
    except (*ERRORS,csv.Error) as exc:
        parser.error(str(exc))


if __name__=='__main__':
    sys.exit(main())
