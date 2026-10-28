"""Human-readable structural descriptors for the existing patch summary."""
from collections import Counter, defaultdict
import io
import math
import os
import shutil

from Bio.PDB import MMCIFParser
from Bio.Data.PDBData import protein_letters_3to1

from sieve.interface_description import (
    DescriptionOptions, atom_interactions, neighborhood, protein_residues,
    residue_name, sasa_burial, secondary_structure, PHOSPHORYLATED,
)
from sieve.alphafold_pdockq2_patches import ranges
from sieve.alphafold_dd_faces import candidate_faces, candidate_interface, validate_threshold


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


def charge_text(values):
    charge = values['nominal_sidechain_charge']
    if charge is not None:
        return f'{charge:+g}'
    if values.get('unknown_ptm_count'):
        return 'n/a (unsupported PTM)'
    return (f"{values['charge_phosphate_minus1']:+g} (phosphate −1) / "
            f"{values['charge_phosphate_minus2']:+g} (phosphate −2)")


def residue_ranges(numbers):
    return ','.join(str(a) if a == b else f'{a}–{b}' for a, b in ranges(numbers))


class PatchDescriptions:
    def __init__(self, inputs, dssp, options=None, *, dd=False, dd_threshold=3, dd_loop_threshold=1):
        validate_threshold(dd_threshold)
        validate_threshold(dd_loop_threshold, 'loop')
        self.dd_loop_threshold = dd_loop_threshold
        self.dd, self.dd_threshold = dd, dd_threshold
        self.details = {}
        self.models = {}
        self.options = options or DescriptionOptions()
        self.dssp = dssp
        self.pairs = {}
        self.modified = set()
        for key, cif in inputs.items():
            model = MMCIFParser(QUIET=True).get_structure('summary', io.StringIO(cif))[0]
            for chain in model:
                if any(r.resname not in protein_letters_3to1 for r in protein_residues(chain)):
                    self.modified.add(key)
                seen = set()
                for residue in protein_residues(chain):
                    if residue.id[2] != ' ' or residue.id[1] in seen:
                        raise ValueError('descriptors require unique residue numbers without insertion codes')
                    seen.add(residue.id[1])
                    if residue.resname in protein_letters_3to1:
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
            if (src, model_id) in self.modified:
                sites = [(c.id,r.id[1]) for c in model for r in protein_residues(c)
                         if r.resname not in protein_letters_3to1]
                print(f'  Model {model_id}: DSSP coverage at modified residues: '
                      f'{sum(site in codes for site in sites)}/{len(sites)} assigned.', file=stream)
            for chain in model:
                residues = protein_residues(chain)
                if not residues:
                    continue
                sequence = ''.join(protein_letters_3to1.get(r.resname, f'[{r.resname}]') for r in residues)
                parts = [s for s in segments if s['chain'] == chain.id]
                constitution = ' → '.join(f"{s['label']}[{s['code']}]({s['start']}–{s['end']})" for s in parts) or 'unassigned'
                counts = Counter(codes.get((chain.id, r.id[1]), '?') for r in residues)
                # Same sequence does not imply the same fold: retain every distinct assignment.
                groups[sequence].append((model_id, chain.id, constitution, counts))
        for index, (sequence, entries) in enumerate(groups.items(), 1):
            chains = ', '.join(sorted({c for _, c, _, _ in entries}))
            print(f'  Sequence group {index} | chains {chains} | {sum(entries[0][3].values())} residues', file=stream)
            variants = defaultdict(list)
            for model_id, chain, constitution, counts in entries:
                variants[(constitution, tuple(sorted(counts.items())))].append(f"model {model_id if model_id != '' else 'unnumbered'} chain {chain}")
            for (constitution, counts), members in variants.items():
                print(f"    {'; '.join(members)}", file=stream)
                print(f'      {constitution}', file=stream)
                nres = sum(count for _, count in counts)
                print('      Composition: ' + ', '.join(f'{code} {count}/{nres} ({100*count/nres:.1f}%)' for code, count in counts), file=stream)
        print('  Codes: H alpha, G 3₁₀, I pi, E strand, B isolated beta bridge, P polyproline; T/S/– loop; ? unassigned.', file=stream)
        print(file=stream)

    def describe(self, row):
        """One descriptor calculation shared by TSV and human-readable summaries."""
        key = (row['source'], row['model'])
        chains = (row['chain 1'], row['chain 2'])
        cache_key = (*key, *chains, row['interface residues 1'], row['interface residues 2'])
        if cache_key in self.details:
            return self.details[cache_key]
        model, codes, segments, _ = self.models[key]
        pair_key = (*key, *chains)
        if pair_key not in self.pairs:
            self.pairs[pair_key] = (atom_interactions(model, chains, self.options),
                                   sasa_burial(model, chains, self.options.sasa_points))
        interactions, burial = self.pairs[pair_key]
        selected = {c: {int(n) for n in row[f'interface residues {i}'].split(',') if n}
                    for i,c in enumerate(chains,1)}
        contexts = {c: neighborhood(model[c], selected[c], self.options.flank, self.options.radius) for c in chains}
        structures = {}
        for c in chains:
            parts = []
            for s in segments:
                hits = selected[c].intersection(s['residues']) if s['chain'] == c else set()
                if hits:
                    parts.append(f"{s['label']}[{s['code']}]:{residue_ranges(hits)}")
            missing = sorted(n for n in selected[c] if (c,n) not in codes)
            if missing:
                parts.append(f'unassigned:{residue_ranges(missing)}')
            structures[c] = '; '.join(parts) or None
        areas = {c: sum(burial[(c,n)]['buried_area'] for n in selected[c]) for c in chains}
        cs = {c: sum(burial[(c,n)]['buried_nonpolar_area'] for n in selected[c]) for c in chains}
        features = [x for x in interactions if x['residue_1'] in selected[chains[0]] and x['residue_2'] in selected[chains[1]]]
        totals = Counter(kind for kind,_,_ in {(x['type'],x['residue_1'],x['residue_2']) for x in features})
        salts = {}
        for x in features:
            if x['type'] == 'salt_bridge_candidate':
                pair = (x['residue_1'],x['residue_2'])
                if pair not in salts or x['distance'] < salts[pair]['distance']:
                    salts[pair] = x
        a,b = (contexts[c]['patch'] for c in chains)
        if a['nominal_sidechain_charge'] is not None and b['nominal_sidechain_charge'] is not None:
            delta = a['nominal_sidechain_charge']-b['nominal_sidechain_charge']
        elif not a.get('unknown_ptm_count',0) and not b.get('unknown_ptm_count',0):
            delta = '; '.join(f"{a.get(f'charge_phosphate_minus{n}',a['nominal_sidechain_charge']) - b.get(f'charge_phosphate_minus{n}',b['nominal_sidechain_charge']):+g} (phosphate −{n})" for n in (1,2))
        else:
            delta = None
        total_area = sum(areas.values())
        metrics = {'total buried surface area': total_area,
                   'total C/S contribution percent': 100*sum(cs.values())/total_area if total_area > 0 else None,
                   'nonpolar contacts': totals['nonpolar_contact'],
                   'salt bridge candidates': totals['salt_bridge_candidate'],
                   'steric overlap candidates': totals['steric_overlap_candidate'],
                   'phosphate/basic contacts': totals['phosphate_basic_contact'],
                   'salt bridge candidate heavy-atom distance min': min((x['distance'] for x in salts.values()),default=None),
                   'salt bridge residues': '; '.join(f"{chains[0]}:{x['aa_1']}{a}/{x['atom_1']}–{chains[1]}:{x['aa_2']}{b}/{x['atom_2']} ({x['distance']:.3f} Å)" for (a,b),x in sorted(salts.items())) or None,
                   'nominal patch charge difference': delta}
        for i,c in enumerate(chains,1):
            metrics.update({f'interface structure {i}': structures[c],
                            f'patch {i} hydropathy mean': contexts[c]['patch']['mean_hydropathy'],
                            f'candidate DD family face {i}': candidate_faces(c,selected[c],codes,segments,self.dd_threshold,self.dd_loop_threshold) if self.dd else None})
        metrics['candidate DD interface'] = candidate_interface(
            metrics['candidate DD family face 1'], metrics['candidate DD family face 2'])
        details = dict(metrics=metrics,contexts=contexts,structures=structures,areas=areas,cs=cs,
                       selected=selected,interactions=interactions,totals=totals)
        self.details[cache_key] = details
        return details

    def write_patch(self, row, stream):
        key = (row['source'], row['model'])
        model, codes, segments, _ = self.models[key]
        chains = (row['chain 1'], row['chain 2'])
        if key in self.modified:
            print('     PTM scoring: mean normalized token-pair PAE per residue pair; not a PTM-calibrated pDockQ2.', file=stream)
            print('     PTM atoms retained. Hydropathy is n/a for sets containing modified residues.', file=stream)
            print('     Phosphate charge scenarios: −1 or −2 per SEP/TPO/PTR; protonation is not inferred.', file=stream)
            print('     Phosphate/basic contacts are distance-only candidates, not assigned hydrogen bonds.', file=stream)
            print('     Nonpolar-contact classification covers standard residues only; PTM carbon/sulfur burial is retained.', file=stream)
            for c in chains:
                for r in protein_residues(model[c]):
                    if r.resname in PHOSPHORYLATED:
                        missing = sorted({'P','O1P','O2P','O3P'} - {a.name for a in r})
                        if missing:
                            print(f'     WARNING: incomplete phosphate {c}:{r.resname}{r.id[1]}: missing {", ".join(missing)}; contacts use observed atoms only.', file=stream)
        details = self.describe(row)
        interactions, selected = details['interactions'], details['selected']
        origin = 'Inherited source pLDDT' if row.get('confidence origin') else 'Local pLDDT'
        print(f'     {origin} (unique contacting residues): ' + ' | '.join(f"{c} {number(row[f'interface pLDDT {i}'])}" for i, c in enumerate(chains, 1)), file=stream)
        for c in chains:
            context = details['contexts'][c]
            patch = context['patch']
            print(f"     {c} patch sequence: {patch['sequence']} | secondary structure: {details['structures'][c] or 'none'}", file=stream)
            if self.dd:
                print(f"       Candidate DD family face: {details['metrics'][f'candidate DD family face {chains.index(c)+1}'] or 'NA'} (six-helix coverage heuristic)", file=stream)
            for label in ('patch', 'upstream', 'downstream', 'sequence_shell', 'spatial_shell'):
                values = context[label]
                print(f"       {label}: n={values['count']} | mean hydropathy {number(values['mean_hydropathy'])} | nominal charge {charge_text(values)} | K/R {values['positive_count']}, D/E {values['negative_count']}, H {values['histidine_count']}", file=stream)
            area = details['areas'][c]
            nonpolar = details['cs'][c]
            print(f'       Buried area: {area:.1f} Å²; C/S contribution {nonpolar:.1f} Å²', file=stream)
        delta = details['metrics']['nominal patch charge difference']
        delta = f'{delta:+g}' if isinstance(delta,(int,float)) else delta or 'n/a (unsupported PTM)'
        print(f'     Nominal patch charge difference ({chains[0]} − {chains[1]}): {delta}', file=stream)
        totals = details['totals']
        print('     Contact chemistry (unique residue pairs): ' + ' | '.join(f'{kind} {totals[kind]}' for kind in ('nonpolar_contact','salt_bridge_candidate','steric_overlap_candidate')), file=stream)
        if key in self.modified:
            print(f"     Phosphate/basic contacts (unique residue pairs within patch): {totals['phosphate_basic_contact']}", file=stream)
        salts = {}
        for x in interactions:
            if x['type'] not in ('salt_bridge_candidate', 'phosphate_basic_contact'):
                continue
            a, b = x['residue_1'], x['residue_2']
            if a not in selected[chains[0]] and b not in selected[chains[1]]:
                continue
            if (a,b) not in salts or x['distance'] < salts[(a,b)]['distance']:
                salts[(a,b)] = x
        for (a,b), x in sorted(salts.items()):
            scope = 'within patch' if a in selected[chains[0]] and b in selected[chains[1]] else 'one endpoint in patch'
            label = 'Phosphate/basic candidate' if x['type'] == 'phosphate_basic_contact' else 'Salt candidate'
            print(f"       {label} ({scope}): {chains[0]}:{x['aa_1']}{a}/{x['atom_1']}–{chains[1]}:{x['aa_2']}{b}/{x['atom_2']} {x['distance']:.2f} Å", file=stream)
        print(f'     Context: ±{self.options.flank} sequence residues; {self.options.radius:g} Å same-chain heavy-atom shell; shells exclude patch.', file=stream)
        print('     Chemistry is geometric evidence; nominal charge is not electrostatic potential.', file=stream)
