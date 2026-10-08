#!/usr/bin/env python3
"""Describe and compare two AF ensembles. See docs/interface-comparison.md."""
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path
import shutil
import sys
import Bio
import numpy

from sieve.interface_description import DescriptionOptions, describe_input
from sieve.interface_comparison import compare_ensembles


def write_tsv(path, rows):
    if not rows:
        path.write_text('')
        return
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    return 'NA' if value is None else f'{value:.3f}'


def report(reference, query, comparison):
    lines = ['# AF interface comparison', '',
             'Exploratory comparison of predicted structures; scores are not binding probabilities.',
             'Residue correspondence is provisional, especially for divergent or incomplete constructs.', '',
             '## Chain and residue correspondence', '']
    for m in comparison['primary']['mappings']:
        row = m['rows'][0] if m['rows'] else {}
        lines.append(f"- Reference {row.get('reference_chain','?')} → query {row.get('query_chain','?')}: {m['aligned_pairs']} pairs; identity {number(m['identity'])}; reference/query coverage {m['reference_coverage']:.3f}/{m['query_coverage']:.3f}; tied optimal alignments: {m['multiple_optimal_alignments']}.")
    lines += ['', '## Model comparisons', '', '| Ref | Query | Contact F1 | Precision | Recall | Unmapped contacts ref/query | Anchor RMSD Å | Partner RMSD after anchor fit Å | Partner independent RMSD Å |', '|---|---|---|---|---|---|---|---|---|']
    for row in comparison['primary']['comparisons']:
        lines.append(f"| {row['reference_model']} | {row['query_model']} | {number(row['f1'])} | {number(row['precision'])} | {number(row['recall'])} | {row['reference_unmapped_contacts']}/{row['query_unmapped_contacts']} | {number(row['anchor_rmsd'])} | {number(row['partner_rmsd_after_anchor_fit'])} | {number(row['partner_independent_rmsd'])} |")
    lines += ['', '## Mapping sensitivity', '']
    for v in [comparison['primary'], *comparison['mapping_sensitivity']]:
        scores = [r['f1'] for r in v['comparisons'] if r['f1'] is not None]
        lines.append(f"- Gap-open penalty {v['gap_open']}: contact F1 range {number(min(scores) if scores else None)}–{number(max(scores) if scores else None)}; mapped pairs {[m['aligned_pairs'] for m in v['mappings']]}.")
    for label, models in [('Reference',reference),('Query',query)]:
        lines += ['',f'## {label} ensemble', '']
        for m in models:
            full = next((p for p in m['patches'] if p['score']['patch id']=='full'), None)
            score = full['score']['pDockQ2 max'] if full else None
            lines += [f"### Model {m['model']}", '', f"Cα contacts: {len(m['contacts'])}; full pDockQ2 max: {number(score)}. DSSP: {m['dssp']['status']}.",
                      f"Buried surface area per chain (Å²): {', '.join(c+' '+number(v) for c,v in m['buried_area_by_chain'].items())}.", '']
            for chain in m['chains']:
                segments = [s for s in m['secondary_elements'] if s['chain']==chain and s['code'] in ('H','G','I','E','B','P')]
                lines.append(f"- {chain} elements: " + (', '.join(f"{s['label']} {s['start']}–{s['end']}" for s in segments) or 'none assigned'))
            for p in m['patches']:
                s = p['score']
                if s['patch id']!='full' and s['status']!='leaf':
                    continue
                lines += ['',f"**Patch {s['patch id']} ({s['status']})**, pDockQ2 max {s['pDockQ2 max']:.3f}"]
                if s['contact count']:
                    lines.append(f"PAE-only confidence: {s['chain 1']}→{s['chain 2']} {s['normalized PAE 1 to 2']:.3f}; {s['chain 2']}→{s['chain 1']} {s['normalized PAE 2 to 1']:.3f}.")
                    lines.append(f"Interface pLDDT: {s['chain 1']} {s['interface pLDDT 1']:.2f}; {s['chain 2']} {s['interface pLDDT 2']:.2f}.")
                for chain, d in p['description'].items():
                    part = d['patch']
                    elements = [s['label'] for s in p['secondary_elements'] if s['chain']==chain]
                    lines.append(f"- {chain}: residues {','.join(map(str,part['residues'])) or 'none'}; elements {', '.join(elements) or 'none'}; mean hydropathy {number(part['mean_hydropathy'])}; nominal charge {part['nominal_sidechain_charge']:+g}; sequence/spatial shell hydropathy {number(d['sequence_shell']['mean_hydropathy'])}/{number(d['spatial_shell']['mean_hydropathy'])}.")
                lines.append(f"- Typed residue-pair candidates: {p['typed_residue_pair_counts'] or 'none'}.")
    lines += ['', '## Interpretation limits and validation', '',
              'The weakest link is the residue mapping: inspect mapping.tsv and the independent domain fits before interpreting contact F1. Re-run with fuller constructs and a curated correspondence to test stability. Low overlap can reflect mapping errors, different folds, or a different predicted interface; this report does not distinguish them automatically.',
              'Five AF models are not independent biological replicates. Contact recurrence is descriptive. Compare the known-protein AF ensemble with an experimental complex before treating its interface as a validated structural reference.',
              'Salt bridges and nonpolar contacts are geometric candidates; nominal charge is not an electrostatic potential. No binding energy or affinity is computed.', '',
              f'Definitions and sources: [interface-comparison.md]({Path(__file__).resolve().parents[1] / "docs/interface-comparison.md"}). Full atom details, all patch nodes, confidence values, and alternate mappings are in the JSON files.', '']
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', help='known-protein AF ZIP')
    parser.add_argument('query', help='novel-protein AF ZIP')
    parser.add_argument('--reference-chains', nargs=2, required=True)
    parser.add_argument('--query-chains', nargs=2, required=True)
    parser.add_argument('--output-dir', required=True, help='new directory; existing directories are refused')
    parser.add_argument('--distance-cutoff', type=float, default=8)
    parser.add_argument('--sequence-flank', type=int, default=5)
    parser.add_argument('--spatial-radius', type=float, default=8)
    parser.add_argument('--sasa-points', type=int, default=240)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--dssp', help='mkdssp executable')
    group.add_argument('--no-dssp', action='store_true')
    args = parser.parse_args(argv)
    for path in (args.reference, args.query):
        if not Path(path).is_file() or Path(path).suffix.lower() != '.zip':
            parser.error('reference and query must be AF ZIP files')
    output = Path(args.output_dir)
    if output.exists():
        parser.error('output directory already exists; choose a new directory')
    local = Path(__file__).resolve().parents[1]/'tmp/dssp/bin/mkdssp'
    dssp = None if args.no_dssp else args.dssp or shutil.which('mkdssp') or (str(local) if local.is_file() else None)
    if not args.no_dssp and dssp is None:
        parser.error('mkdssp not found; supply --dssp or explicitly use --no-dssp')
    options = DescriptionOptions(cutoff=args.distance_cutoff,flank=args.sequence_flank,radius=args.spatial_radius,sasa_points=args.sasa_points)
    reference = describe_input(args.reference,args.reference_chains,options,dssp)
    query = describe_input(args.query,args.query_chains,options,dssp)
    comparison = compare_ensembles(reference,query)
    comparison['provenance'] = {
        'inputs': [{'path': str(Path(p).resolve()), 'sha256': hashlib.sha256(Path(p).read_bytes()).hexdigest()} for p in (args.reference,args.query)],
        'python':platform.python_version(), 'biopython':Bio.__version__, 'numpy':numpy.__version__,
        'source_sha256': {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__), Path(__file__).resolve().parents[1]/'sieve/interface_description.py',
             Path(__file__).resolve().parents[1]/'sieve/interface_comparison.py']}}
    output.mkdir(parents=True)
    for name, data in [('reference',reference),('query',query),('comparison',comparison)]:
        (output/f'{name}.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    write_tsv(output/'mapping.tsv',[r for m in comparison['primary']['mappings'] for r in m['rows']])
    write_tsv(output/'model-comparisons.tsv',[{k:v for k,v in r.items() if k!='typed_contacts'} for r in comparison['primary']['comparisons']])
    for label, models in [('reference',reference),('query',query)]:
        write_tsv(output/f'{label}-interactions.tsv',[{'model':m['model'],**x} for m in models for x in m['interactions']])
        write_tsv(output/f'{label}-recurrence.tsv',comparison[f'{label}_recurrence'])
    (output/'report.md').write_text(report(reference,query,comparison))
    print(f'Report: {(output/"report.md").resolve()}')
    return 0


if __name__=='__main__':
    sys.exit(main())
