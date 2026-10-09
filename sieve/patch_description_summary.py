"""Human-readable structural descriptors for the existing patch summary."""
from collections import Counter, defaultdict
import io
import math
import os
import shutil

from Bio.PDB import MMCIFParser

from sieve.interface_description import (
    DescriptionOptions, atom_interactions, neighborhood, protein_residues,
    residue_name, sasa_burial, secondary_structure,
)
from sieve.alphafold_pdockq2_patches import ranges


def find_dssp(explicit=None):
    configured = explicit if explicit is not None else os.environ.get('SIEVE_DSSP')
    if configured is not None:
        executable = shutil.which(os.path.expanduser(configured)) if configured else None
        if executable is None:
            source = '--dssp' if explicit is not None else 'SIEVE_DSSP'
            raise ValueError(f'{source} does not identify an executable: {configured!r}')
        return executable
    return shutil.which('mkdssp')


def number(value):
    return 'n/a' if value is None else f'{value:.2f}'


def residue_ranges(numbers):
    return ','.join(str(a) if a == b else f'{a}–{b}' for a, b in ranges(numbers))


class PatchDescriptions:
    def __init__(self, inputs, dssp, options=None):
        self.models = {}
        self.options = options or DescriptionOptions()
        self.dssp = dssp
        self.pairs = {}
        for key, cif in inputs.items():
            model = MMCIFParser(QUIET=True).get_structure('summary', io.StringIO(cif))[0]
            for chain in model:
                seen = set()
                for residue in protein_residues(chain):
                    if residue.id[2] != ' ' or residue.id[1] in seen:
                        raise ValueError('descriptors require unique residue numbers without insertion codes')
                    seen.add(residue.id[1])
                    residue_name(residue)
                    if any(not all(math.isfinite(float(v)) for v in atom.coord) for atom in residue):
                        raise ValueError('nonfinite atom coordinates')
            codes, segments, status = secondary_structure(model, cif, dssp)
            self.models[key] = (model, codes, segments, status)

    def write_chains(self, source, stream):
        print('Chain secondary structure (N→C; identical sequences grouped):', file=stream)
        if self.dssp is None:
            print('  DSSP unavailable/disabled: secondary structure not assigned.', file=stream)
        groups = defaultdict(list)
        for (src, model_id), (model, codes, segments, status) in self.models.items():
            if src != source:
                continue
            for chain in model:
                residues = protein_residues(chain)
                if not residues:
                    continue
                sequence = ''.join(residue_name(r) for r in residues)
                parts = [s for s in segments if s['chain'] == chain.id]
                constitution = ' → '.join(f"{s['label']}[{s['code']}]({s['start']}–{s['end']})" for s in parts) or 'unassigned'
                counts = Counter(codes.get((chain.id, r.id[1]), '?') for r in residues)
                # Same sequence does not imply the same fold: retain every distinct assignment.
                groups[sequence].append((model_id, chain.id, constitution, counts))
        for index, (sequence, entries) in enumerate(groups.items(), 1):
            chains = ', '.join(sorted({c for _, c, _, _ in entries}))
            print(f'  Sequence group {index} | chains {chains} | {len(sequence)} residues', file=stream)
            variants = defaultdict(list)
            for model_id, chain, constitution, counts in entries:
                variants[(constitution, tuple(sorted(counts.items())))].append(f"model {model_id if model_id != '' else 'unnumbered'} chain {chain}")
            for (constitution, counts), members in variants.items():
                print(f"    {'; '.join(members)}", file=stream)
                print(f'      {constitution}', file=stream)
                print('      Composition: ' + ', '.join(f'{code} {count}/{len(sequence)} ({100*count/len(sequence):.1f}%)' for code, count in counts), file=stream)
        print('  Codes: H alpha, G 3₁₀, I pi, E strand, B isolated beta bridge, P polyproline; T/S/– loop; ? unassigned.', file=stream)
        print(file=stream)

    def write_patch(self, row, stream):
        key = (row['source'], row['model'])
        model, codes, segments, _ = self.models[key]
        chains = (row['chain 1'], row['chain 2'])
        pair_key = (*key, *chains)
        if pair_key not in self.pairs:
            self.pairs[pair_key] = (atom_interactions(model, chains, self.options),
                                    sasa_burial(model, chains, self.options.sasa_points))
        interactions, burial = self.pairs[pair_key]
        selected = {c: {n for n in map(int, row[f'interface residues {i}'].split(','))} if row[f'interface residues {i}'] else set() for i, c in enumerate(chains, 1)}
        print('     Local pLDDT (unique contacting residues): ' + ' | '.join(f"{c} {row[f'interface pLDDT {i}']:.2f}" for i, c in enumerate(chains, 1)), file=stream)
        charges = {}
        for c in chains:
            context = neighborhood(model[c], selected[c], self.options.flank, self.options.radius)
            patch = context['patch']
            charges[c] = patch['nominal_sidechain_charge']
            parts = []
            for s in segments:
                hits = selected[c].intersection(s['residues']) if s['chain'] == c else set()
                if hits:
                    parts.append(f"{s['label']}[{s['code']}]:{residue_ranges(hits)}")
            missing = sorted(n for n in selected[c] if (c,n) not in codes)
            if missing:
                parts.append(f'unassigned:{residue_ranges(missing)}')
            print(f"     {c} patch sequence: {patch['sequence']} | secondary structure: {'; '.join(parts) or 'none'}", file=stream)
            for label in ('patch', 'upstream', 'downstream', 'sequence_shell', 'spatial_shell'):
                values = context[label]
                print(f"       {label}: n={values['count']} | mean hydropathy {number(values['mean_hydropathy'])} | nominal charge {values['nominal_sidechain_charge']:+g} | K/R {values['positive_count']}, D/E {values['negative_count']}, H {values['histidine_count']}", file=stream)
            area = sum(burial[(c,n)]['buried_area'] for n in selected[c])
            nonpolar = sum(burial[(c,n)]['buried_nonpolar_area'] for n in selected[c])
            print(f'       Buried area: {area:.1f} Å²; C/S contribution {nonpolar:.1f} Å²', file=stream)
        print(f'     Nominal patch charge difference ({chains[0]} − {chains[1]}): {charges[chains[0]]-charges[chains[1]]:+g}', file=stream)
        features = [x for x in interactions if x['residue_1'] in selected[chains[0]] and x['residue_2'] in selected[chains[1]]]
        counts = Counter((x['type'], x['residue_1'], x['residue_2']) for x in features)
        totals = Counter(kind for kind, _, _ in counts)
        print('     Contact chemistry (unique residue pairs): ' + ' | '.join(f'{kind} {totals[kind]}' for kind in ('nonpolar_contact','salt_bridge_candidate','steric_overlap_candidate')), file=stream)
        salts = {}
        for x in interactions:
            if x['type'] != 'salt_bridge_candidate':
                continue
            a, b = x['residue_1'], x['residue_2']
            if a not in selected[chains[0]] and b not in selected[chains[1]]:
                continue
            if (a,b) not in salts or x['distance'] < salts[(a,b)]['distance']:
                salts[(a,b)] = x
        for (a,b), x in sorted(salts.items()):
            scope = 'within patch' if a in selected[chains[0]] and b in selected[chains[1]] else 'one endpoint in patch'
            print(f"       Salt candidate ({scope}): {chains[0]}:{x['aa_1']}{a}/{x['atom_1']}–{chains[1]}:{x['aa_2']}{b}/{x['atom_2']} {x['distance']:.2f} Å", file=stream)
        print(f'     Context: ±{self.options.flank} sequence residues; {self.options.radius:g} Å same-chain heavy-atom shell; shells exclude patch.', file=stream)
        print('     Chemistry is geometric evidence; nominal charge is not electrostatic potential.', file=stream)
