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
from sieve.alphafold_result_format import PatchInput, read_patch_inputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='fold_locator, CIF, ZIP, directory, or - to read locators/result TSV from stdin')
    parser.add_argument('--full-data', help='matching confidence JSON, required for CIF')
    parser.add_argument('--output', help='optional patch TSV (recomputed; no cache); summary printed unless --silent')
    parser.add_argument('--silent', action='store_true', help='suppress stdout summaries; errors and warnings remain on stderr')
    parser.add_argument('--dd', action='store_true', help='enable provisional six-helix DD face candidates')
    parser.add_argument('--dd-helix-presence-threshold', type=int, default=3, help='minimum contacting residues per required helix (default 3)')
    parser.add_argument('--dd-loop-presence-threshold', type=int, default=1, help='minimum contacting residues per required loop (default 1)')
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
        from sieve.alphafold_dd_faces import validate_threshold
        validate_threshold(args.dd_helix_presence_threshold)
        validate_threshold(args.dd_loop_presence_threshold, 'loop')
        requests = list(read_patch_inputs(sys.stdin)) if args.input == '-' else [PatchInput(args.input)]
        output = Path(args.output).resolve() if args.output else None
        if output is not None:
            sources = [Path(request.locator.split('::', 1)[0]).resolve() for request in requests]
            if args.full_data:
                sources.append(Path(args.full_data).resolve())
            if output in sources:
                raise ValueError('output must not overwrite an input')
            if output.suffix.lower() in ('.zip', '.cif', '.json'):
                raise ValueError('output must not use an input-file extension')
        if not requests:
            return 0
        options = Options(args.distance_cutoff, args.max_bridge_residues,
                          args.min_residues_per_chain, args.min_contacts, args.max_separator_checks)
        rows = []
        summary = io.StringIO()
        for request in requests:
            try:
                targets, inputs, model_scores, fold_ids, sequence_ids = {}, {}, {}, {}, {}
                current_rows = analyze(request.locator, args.full_data, options, args.chains,
                    summary_targets=targets, summary_inputs=inputs, model_scores=model_scores,
                    fold_ids=fold_ids, sequence_ids=sequence_ids, selected_model=request.model_id)
                descriptions = PatchDescriptions(inputs, None if args.no_dssp else find_dssp(args.dssp),
                    DescriptionOptions(cutoff=args.distance_cutoff, flank=args.sequence_flank,
                                       radius=args.spatial_radius, sasa_points=args.sasa_points),
                    dd=args.dd, dd_threshold=args.dd_helix_presence_threshold, dd_loop_threshold=args.dd_loop_presence_threshold)
                if not args.silent:
                    write_summary(current_rows, targets, summary, descriptions, structure=args.structure,
                                  model_scores=model_scores, fold_ids=fold_ids, sequence_ids=sequence_ids, chains=args.chains)
                if output is not None:
                    for row in current_rows:
                        if row['scope'] != 'patch' or row['status'] != 'leaf':
                            continue
                        row.update(descriptions.describe(row)['metrics'])
                        row.update(sequence_id=sequence_ids.get(row['source']), fold_id=fold_ids.get(row['source']))
                        row['dd enabled'] = args.dd
                        row['dd helix presence threshold'] = args.dd_helix_presence_threshold
                        row['dd loop presence threshold'] = args.dd_loop_presence_threshold
                rows.extend(current_rows)
            except (OSError, ValueError, KeyError, TypeError, zipfile.BadZipFile) as error:
                if request.record_number is not None:
                    raise ValueError(f'stdin record {request.record_number}: {error}') from error
                raise
        if output is not None:
            with output.open('w', newline='') as stream:
                write_patches(rows, stream)
        if not args.silent:
            sys.stdout.write(summary.getvalue())
    except (OSError, ValueError, KeyError, TypeError, zipfile.BadZipFile) as error:
        parser.error(str(error))
    return 0


if __name__ == '__main__':
    sys.exit(main())
