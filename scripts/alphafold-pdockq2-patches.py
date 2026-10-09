#!/usr/bin/env python3
"""Detect spatial contact patches and recursively split narrow residue bridges.

Algorithm, equations, input mapping, limitations, and output specification:
../docs/pDockQ2.md (relative to this script).
"""
import argparse
import io
from pathlib import Path
import sys
import zipfile

from sieve.alphafold_pdockq2_patches import Options, analyze, write_patches, write_summary
from sieve.patch_description_summary import PatchDescriptions, find_dssp
from sieve.interface_description import DescriptionOptions


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='CIF, ZIP, or directory of ZIPs')
    parser.add_argument('--full-data', help='matching confidence JSON, required for CIF')
    parser.add_argument('--output', help='optional patch TSV (recomputed; no cache); summary always printed')
    parser.add_argument('--chains', nargs=2, metavar=('CHAIN1', 'CHAIN2'))
    parser.add_argument('--distance-cutoff', type=float, default=8.0)
    parser.add_argument('--max-bridge-residues', type=int, default=1)
    parser.add_argument('--min-residues-per-chain', type=int, default=3)
    parser.add_argument('--min-contacts', type=int, default=5)
    parser.add_argument('--max-separator-checks', type=int, default=1000000)
    parser.add_argument('--sequence-flank', type=int, default=5)
    parser.add_argument('--spatial-radius', type=float, default=8)
    parser.add_argument('--sasa-points', type=int, default=240)
    dssp_group = parser.add_mutually_exclusive_group()
    dssp_group.add_argument('--dssp', help='mkdssp executable for secondary structure')
    dssp_group.add_argument('--no-dssp', action='store_true', help='omit secondary structure; retain other descriptors')
    args = parser.parse_args(argv)
    try:
        output = Path(args.output).resolve() if args.output else None
        if output is not None and (output == Path(args.input).resolve() or (args.full_data and output == Path(args.full_data).resolve())):
            raise ValueError('output must not overwrite an input')
        if output is not None and output.suffix.lower() in ('.zip', '.cif', '.json'):
            raise ValueError('output must not use an input-file extension')
        options = Options(args.distance_cutoff, args.max_bridge_residues,
                          args.min_residues_per_chain, args.min_contacts, args.max_separator_checks)
        targets = {}
        inputs = {}
        rows = analyze(args.input, args.full_data, options, args.chains, summary_targets=targets, summary_inputs=inputs)
        descriptions = PatchDescriptions(inputs, None if args.no_dssp else find_dssp(args.dssp),
            DescriptionOptions(cutoff=args.distance_cutoff, flank=args.sequence_flank,
                               radius=args.spatial_radius, sasa_points=args.sasa_points))
        summary = io.StringIO()
        write_summary(rows, targets, summary, descriptions)
        if output is not None:
            with output.open('w', newline='') as stream:
                write_patches(rows, stream)
        sys.stdout.write(summary.getvalue())
    except (OSError, ValueError, KeyError, TypeError, zipfile.BadZipFile) as error:
        parser.error(str(error))
    return 0


if __name__ == '__main__':
    sys.exit(main())
