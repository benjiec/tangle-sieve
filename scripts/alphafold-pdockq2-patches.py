#!/usr/bin/env python3
"""Detect spatial contact patches and recursively split narrow residue bridges.

Algorithm, equations, input mapping, limitations, and output specification:
../docs/pDockQ2.md (relative to this script).
"""
import argparse
import io
import sys
import zipfile

from sieve.alphafold_pdockq2_patches import Options, write_summary
from sieve.alphafold_patch_metrics import extract_patch_metrics
from sieve.alphafold_secondary import find_dssp
from sieve.interface_description import DescriptionOptions
from sieve.alphafold_result_format import PatchInput, read_patch_inputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='fold_locator, CIF, ZIP, directory, or - to read locators/result TSV from stdin')
    parser.add_argument('--full-data', help='matching confidence JSON, required for CIF')
    parser.add_argument('--chains', nargs=2, metavar=('CHAIN1', 'CHAIN2'))
    parser.add_argument('--distance-cutoff', type=float, default=8.0)
    parser.add_argument('--max-bridge-residues', type=int, default=1)
    parser.add_argument('--min-residues-per-chain', type=int, default=3)
    parser.add_argument('--min-contacts', type=int, default=5)
    parser.add_argument('--max-separator-checks', type=int, default=1000000)
    parser.add_argument('--structure', action='store_true', help='print the per-chain structure summary')
    parser.add_argument('--sequence-flank', type=int, default=5)
    parser.add_argument('--spatial-radius', type=float, default=8)
    parser.add_argument('--sasa-points', type=int, default=240)
    dssp_group = parser.add_mutually_exclusive_group()
    dssp_group.add_argument('--dssp', help='mkdssp executable; overrides SIEVE_DSSP, then PATH discovery')
    dssp_group.add_argument('--no-dssp', action='store_true', help='omit secondary structure; retain other descriptors')
    args = parser.parse_args(argv)
    try:
        requests = list(read_patch_inputs(sys.stdin)) if args.input == '-' else [PatchInput(args.input)]
        if not requests:
            return 0
        options = Options(args.distance_cutoff, args.max_bridge_residues,
                          args.min_residues_per_chain, args.min_contacts, args.max_separator_checks)
        summary = io.StringIO()
        dssp = None if args.no_dssp else find_dssp(args.dssp)
        for request in requests:
            try:
                metrics = extract_patch_metrics(
                    request.locator, options=options, chains=args.chains,
                    full_data=args.full_data, dssp=dssp, selected_model=request.model_id,
                    description_options=DescriptionOptions(
                        cutoff=args.distance_cutoff, flank=args.sequence_flank,
                        radius=args.spatial_radius, sasa_points=args.sasa_points,
                    ),
                )
                write_summary(
                    metrics.rows, metrics.targets, summary, metrics.descriptions,
                    structure=args.structure, model_scores=metrics.model_scores,
                    fold_ids=metrics.fold_ids, sequence_ids=metrics.sequence_ids, chains=args.chains,
                )
            except (OSError, ValueError, KeyError, TypeError, zipfile.BadZipFile) as error:
                if request.record_number is not None:
                    raise ValueError(f'stdin record {request.record_number}: {error}') from error
                raise
        sys.stdout.write(summary.getvalue())
    except (OSError, ValueError, KeyError, TypeError, zipfile.BadZipFile) as error:
        parser.error(str(error))
    return 0


if __name__ == '__main__':
    sys.exit(main())
