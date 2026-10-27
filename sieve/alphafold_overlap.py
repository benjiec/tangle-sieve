"""Rigid assembly through one uniquely shared protein sequence; no relaxation."""
import csv
import io
import itertools
import json
import math
import os
import shutil
import string
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFIO

from sieve.alphafold_files import model_number
from sieve.alphafold_manifest import sequence_id
from sieve.alphafold_structure import (ASSEMBLY_FORMAT, COMPLEX_FILE, assembly_id,
                                      digest, load_structure, open_structure_fold)
from sieve.interface_description import VDW

VERSION = '1.0'
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}


@dataclass(frozen=True)
class Options:
    clash_overlap: float = 0.4
    rmsd_warning: float = 2.0
    site_rmsd_warning: float = 1.5
    coverage_warning: float = 0.9
    site_distance: float = 8.0
    contact_distance: float = 4.5

    def __post_init__(self):
        for key, value in asdict(self).items():
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f'{key} must be finite and positive')
        if self.coverage_warning > 1:
            raise ValueError('coverage-warning must be <= 1')


def fit(fixed, mobile):
    """Row-vector proper rotation: mobile @ rotation + translation -> fixed."""
    if len(fixed) < 3:
        raise ValueError('alignment requires at least three matched C-alpha coordinates')
    f, m = fixed.mean(axis=0), mobile.mean(axis=0)
    if min(np.linalg.matrix_rank(fixed-f), np.linalg.matrix_rank(mobile-m)) < 2:
        raise ValueError('alignment coordinates are collinear or degenerate')
    u, _, vt = np.linalg.svd((mobile-m).T @ (fixed-f))
    correction = np.eye(3)
    correction[-1,-1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return rotation, f-m @ rotation


def coordinates(atoms):
    return np.array([a['coord'] for a in atoms], dtype=float).reshape(-1,3)


def heavy(atoms):
    selected = [a for a in atoms if a['element'] not in ('H','D')]
    if not selected:
        raise ValueError('chain contains no heavy-atom coordinates')
    unknown = {a['element'] for a in selected} - VDW.keys()
    if unknown:
        raise ValueError(f'unsupported van der Waals elements: {sorted(unknown)}')
    return selected


def transformed(atoms, rotation, translation):
    return [dict(a,coord=p.tolist()) for a,p in zip(atoms,coordinates(atoms) @ rotation + translation)]


def unique_shared(first, second):
    common = set(first.sequences.values()) & set(second.sequences.values())
    if len(common) != 1:
        raise ValueError(f'exactly one shared sequence is required; found {len(common)}')
    sequence = next(iter(common))
    left=[c for c,s in first.sequences.items() if s == sequence]
    right=[c for c,s in second.sequences.items() if s == sequence]
    if len(left) != 1 or len(right) != 1:
        raise ValueError('shared sequence must occur exactly once in each input (duplicated shared chain)')
    return left[0],right[0]


def pair_geometry(left, right, options):
    """Chunked all-heavy-atom distances; nonshared chains have no declared covalent bonds."""
    left,right=heavy(left),heavy(right)
    right_xyz=coordinates(right)
    right_radii=np.array([VDW[a['element']] for a in right])
    closest=math.inf;max_overlap=-math.inf;clashes=[];contacts=set()
    for start in range(0,len(left),256):
        block=left[start:start+256]
        distances=np.linalg.norm(coordinates(block)[:,None,:]-right_xyz[None,:,:],axis=2)
        overlap=np.array([VDW[a['element']] for a in block])[:,None]+right_radii[None,:]-distances
        closest=min(closest,float(distances.min()));max_overlap=max(max_overlap,float(overlap.max()))
        for i,j in zip(*np.where(distances <= options.contact_distance)):
            contacts.add((block[i]['number'],right[j]['number']))
        for i,j in zip(*np.where(overlap > options.clash_overlap)):
            a,b=block[i],right[j]
            clashes.append(dict(residue_a=a['number'],atom_a=a['name'],residue_b=b['number'],atom_b=b['name'],
                                distance_A=float(distances[i,j]),overlap_A=float(overlap[i,j]),
                                backbone_backbone=a['name'] in BACKBONE and b['name'] in BACKBONE))
    summary=dict(min_distance_A=closest,max_overlap_A=max_overlap,clash_atom_pairs=len(clashes),
                 clash_residue_pairs=len({(d['residue_a'],d['residue_b']) for d in clashes}),
                 backbone_backbone_clashes=sum(d['backbone_backbone'] for d in clashes),
                 contact_residue_pairs=len(contacts))
    return summary,sorted(clashes,key=lambda d:-d['overlap_A']),contacts


def nearby(ca, atoms, cutoff):
    xyz=coordinates(heavy(atoms))
    return np.array([bool(np.any(np.linalg.norm(xyz-p,axis=1)<=cutoff)) for p in ca])


def rmsd(deviations):
    return float(np.sqrt(np.mean(deviations**2))) if len(deviations) else None


def assessment_status(*, alignment=False, cross_clashes=False, replacement=False,
                      alternative_clashes=False, inherited=False):
    """Summarize warning categories without implying that other checks passed."""
    if replacement and not any((alignment, cross_clashes, alternative_clashes, inherited)):
        return 'WARNING — Partners do not clash; merged shared-chain interface requires refinement'
    issues = []
    if alignment:
        issues.append('shared-chain alignment requires review')
    if cross_clashes:
        issues.append('cross-partner clashes detected')
    if replacement:
        issues.append('merged shared-chain interface requires refinement')
    if alternative_clashes:
        issues.append('partner clashes detected with alternative alignments')
    if inherited:
        issues.append('input assembly warnings require review')
    return 'WARNING — ' + '; '.join(issues) if issues else 'No detected geometric conflict'


def compare(first, second, options):
    shared_a,shared_b=unique_shared(first,second)
    ca_a={a['number']:a for a in first.atoms[shared_a] if a['name']=='CA'}
    ca_b={a['number']:a for a in second.atoms[shared_b] if a['name']=='CA'}
    positions=sorted(ca_a.keys() & ca_b.keys())
    fixed=coordinates([ca_a[n] for n in positions]);mobile=coordinates([ca_b[n] for n in positions])
    rotation,translation=fit(fixed,mobile)
    deviations=np.linalg.norm(mobile @ rotation + translation-fixed,axis=1)
    moved={c:transformed(a,rotation,translation) for c,a in second.atoms.items()}
    length=len(first.sequences[shared_a]);coverage=len(positions)/length
    partner_a=[c for c in first.sequences if c!=shared_a]
    partner_b=[c for c in second.sequences if c!=shared_b]
    if not partner_a or not partner_b:
        raise ValueError('each input must contain at least one nonshared partner chain')
    sites=[];site_rows=[]
    for side,structure,shared,partners,ca in [('a',first,shared_a,partner_a,fixed),('b',second,shared_b,partner_b,mobile)]:
        full_ca={a['number']:a for a in structure.atoms[shared] if a['name']=='CA'}
        for partner in partners:
            mask=nearby(ca,structure.atoms[partner],options.site_distance)
            all_numbers=list(full_ca)
            all_mask=nearby(coordinates(list(full_ca.values())),structure.atoms[partner],options.site_distance)
            all_site={n for n,yes in zip(all_numbers,all_mask) if yes}
            site_rows.append(dict(input=side,partner_chain=partner,matched_site_CA=int(mask.sum()),
                                  observed_site_CA=len(all_site),site_coverage=len(set(positions)&all_site)/len(all_site) if all_site else None,
                                  site_RMSD_A=rmsd(deviations[mask])))
            sites.append((f'site_{side}_{partner}',mask))
    local=[dict(sequence_position=n,CA_deviation_A=float(d),
                source_plddt_a=ca_a[n]['bfactor'] if first.plddt_known[shared_a] else None,
                source_plddt_b=ca_b[n]['bfactor'] if second.plddt_known[shared_b] else None) for n,d in zip(positions,deviations)]
    cross=[];cross_clashes=[]
    for a,b in itertools.product(partner_a,partner_b):
        metrics,clashes,_=pair_geometry(first.atoms[a],moved[b],options)
        cross.append(dict(chain_a=a,chain_b=b,**metrics))
        cross_clashes.extend(dict(chain_a=a,chain_b=b,**c) for c in clashes)
    changes=[];replacement_clashes=[];contact_changes=[]
    key=lambda c:(c['residue_a'],c['atom_a'],c['residue_b'],c['atom_b'])
    for b in partner_b:
        before,old_clashes,old_contacts=pair_geometry(second.atoms[shared_b],second.atoms[b],options)
        after,new_clashes,new_contacts=pair_geometry(first.atoms[shared_a],moved[b],options)
        old={key(c):c for c in old_clashes};new={key(c):c for c in new_clashes}
        introduced=new.keys()-old.keys()
        worsened={k for k in new.keys()&old.keys() if new[k]['overlap_A']>old[k]['overlap_A']+1e-6}
        changes.append(dict(partner_chain_b=b,shared_chain_a=shared_a,shared_chain_b=shared_b,
                            original_clashes=before['clash_atom_pairs'],assembled_clashes=after['clash_atom_pairs'],
                            introduced_clashes=len(introduced),worsened_clashes=len(worsened),
                            resolved_clashes=len(old.keys()-new.keys()),
                            original_max_overlap_A=before['max_overlap_A'],assembled_max_overlap_A=after['max_overlap_A'],
                            contacts_before=len(old_contacts),contacts_after=len(new_contacts),
                            contacts_lost=len(old_contacts-new_contacts),contacts_gained=len(new_contacts-old_contacts)))
        for k,c in new.items():
            replacement_clashes.append(dict(partner_chain_b=b,change='introduced' if k in introduced else 'worsened' if k in worsened else 'preexisting',**c))
        for label,pairs in [('lost',old_contacts-new_contacts),('gained',new_contacts-old_contacts)]:
            contact_changes.extend(dict(partner_chain_b=b,change=label,shared_position=x,partner_position=y) for x,y in sorted(pairs))
    sensitivity=[]
    masks=[('whole_chain',np.ones(len(positions),dtype=bool)),
           ('exclude_terminal_three',np.array([3<n<=length-3 for n in positions]))]+sites
    if first.plddt_known[shared_a] and second.plddt_known[shared_b]:
        masks.append(('source_plddt_ge_80',np.array([ca_a[n]['bfactor']>=80 and ca_b[n]['bfactor']>=80 for n in positions])))
    for label,mask in masks:
        entry=dict(alignment=label,matched_CA=int(mask.sum()),fit_RMSD_A=None,min_partner_distance_A=None,clash_atom_pairs=None,status='insufficient or degenerate coordinates')
        try:
            r,t=fit(fixed[mask],mobile[mask])
        except ValueError:
            sensitivity.append(entry);continue
        metrics=[pair_geometry(first.atoms[a],transformed(second.atoms[b],r,t),options)[0] for a,b in itertools.product(partner_a,partner_b)]
        entry.update(fit_RMSD_A=rmsd(np.linalg.norm(mobile[mask]@r+t-fixed[mask],axis=1)),
                     min_partner_distance_A=min(x['min_distance_A'] for x in metrics),clash_atom_pairs=sum(x['clash_atom_pairs'] for x in metrics),status='computed')
        sensitivity.append(entry)
    warnings=[]
    if coverage<options.coverage_warning:
        warnings.append(f'alignment coverage {coverage:.1%} < {options.coverage_warning:.1%}')
    global_rmsd=rmsd(deviations)
    if global_rmsd>options.rmsd_warning:
        warnings.append(f'shared-chain RMSD {global_rmsd:.2f} > {options.rmsd_warning:g} Å')
    fraction=float(np.mean(deviations>3))
    if fraction>0.1:
        warnings.append(f'{fraction:.1%} of matched Cα atoms deviate >3 Å (heuristic warning at >10%)')
    for site in site_rows:
        label=f"input {site['input']} / partner {site['partner_chain']}"
        if site['site_coverage'] is None or site['site_coverage']<options.coverage_warning:
            warnings.append(f'{label}: missing or incomplete binding-site Cα coverage')
        if site['site_RMSD_A'] is not None and site['site_RMSD_A']>options.site_rmsd_warning:
            warnings.append(f"{label}: site deviation {site['site_RMSD_A']:.2f} > {options.site_rmsd_warning:g} Å")
    alignment_warning = bool(warnings)
    clashes=sum(r['clash_atom_pairs'] for r in cross)
    if clashes: warnings.append(f'{clashes} cross-partner heavy-atom overlaps exceed {options.clash_overlap:g} Å')
    introduced=sum(r['introduced_clashes'] for r in changes)
    worsened=sum(r['worsened_clashes'] for r in changes)
    if introduced or worsened:
        warnings.append(f'shared-chain replacement introduces {introduced} and worsens {worsened} overlapping atom pairs')
    alternative_clashes = any(r['clash_atom_pairs'] for r in sensitivity if r['status']=='computed') and not clashes
    if alternative_clashes:
        warnings.append('partner clashes appear with an alternative alignment subset')
    inherited_warning = False
    for label,source in [('A',first),('B',second)]:
        if source.metadata and source.metadata.get('assessment',{}).get('warnings'):
            inherited_warning = True
            warnings.append(f'input {label} carries warnings from an earlier assembly; see provenance history')
    summary=dict(model_a=first.number,model_b=second.number,shared_chain_a=shared_a,shared_chain_b=shared_b,
                 shared_length=length,matched_CA=len(positions),coverage=coverage,CA_RMSD_A=global_rmsd,
                 CA_max_deviation_A=float(deviations.max()),CA_max_deviation_position=positions[int(deviations.argmax())],
                 fraction_CA_deviation_gt_3_A=fraction,
                 min_partner_distance_A=min(r['min_distance_A'] for r in cross),cross_partner_clashes=clashes,
                 cross_partner_backbone_clashes=sum(r['backbone_backbone_clashes'] for r in cross),
                 cross_partner_max_overlap_A=max(r['max_overlap_A'] for r in cross),
                 shared_replacement_introduced_clashes=introduced,shared_replacement_worsened_clashes=worsened,
                 shared_replacement_max_overlap_A=max(r['assembled_max_overlap_A'] for r in changes),
                 status=assessment_status(alignment=alignment_warning, cross_clashes=bool(clashes),
                                          replacement=bool(introduced or worsened), alternative_clashes=alternative_clashes,
                                          inherited=inherited_warning),
                 warnings='; '.join(warnings))
    return dict(summary=summary,warnings=warnings,sites=site_rows,deviations=local,cross=cross,cross_clashes=cross_clashes,
                changes=changes,replacement_clashes=replacement_clashes,contact_changes=contact_changes,
                sensitivity=sensitivity,rotation=rotation,translation=translation,moved=moved)


def complex_cif(first, second, result):
    shared_a,shared_b=result['summary']['shared_chain_a'],result['summary']['shared_chain_b']
    retained=[('a',c) for c in first.sequences]+[('b',c) for c in second.sequences if c!=shared_b]
    alphabet=string.ascii_uppercase+string.ascii_lowercase+string.digits
    if len(retained)>len(alphabet):
        raise ValueError('assembled model exceeds 62 single-character chain IDs supported by Rosetta')
    used=set();mapping={};data={'data_':'sieve_assembled_complex'}
    def add(category,record):
        for key,value in record.items():
            data.setdefault(category+'.'+key,[]).append(str(value))
    entities={};serial=0
    for side,old in retained:
        source=first if side=='a' else second
        new=old if old in alphabet and old not in used else next(c for c in alphabet if c not in used)
        used.add(new)
        seq=source.sequences[old]
        if seq not in entities:
            e=str(len(entities)+1);entities[seq]=e
            add('_entity',dict(id=e,type='polymer'))
            add('_entity_poly',dict(entity_id=e,type='polypeptide(L)',pdbx_seq_one_letter_code_can=seq))
            for n,mon in sorted(source.monomers[old].items()):
                add('_entity_poly_seq',dict(entity_id=e,num=n,mon_id=mon,hetero='n'))
        e=entities[seq]
        add('_struct_asym',dict(id=new,entity_id=e))
        mapping[new]=dict(input=side,original_chain=old,sequence=seq,plddt_known=source.plddt_known[old],
                          residue_origins=source.residue_origins[old])
        atoms=source.atoms[old] if side=='a' else result['moved'][old]
        for a in atoms:
            serial+=1
            add('_atom_site',dict(group_PDB='ATOM',id=serial,type_symbol=a['element'],label_atom_id=a['name'],label_alt_id='.',
                                  label_comp_id=a['resname'],label_asym_id=new,label_entity_id=e,label_seq_id=a['number'],
                                  pdbx_PDB_ins_code='?',Cartn_x=f"{a['coord'][0]:.6f}",Cartn_y=f"{a['coord'][1]:.6f}",Cartn_z=f"{a['coord'][2]:.6f}",
                                  occupancy=a['occupancy'],B_iso_or_equiv=a['bfactor'],auth_seq_id=a['number'],auth_comp_id=a['resname'],
                                  auth_asym_id=new,auth_atom_id=a['name'],pdbx_PDB_model_num=1))
    writer=MMCIFIO();writer.set_dict(data);stream=io.StringIO();writer.save(stream)
    return stream.getvalue(),mapping


def quick_summary(result, options):
    s=result['summary']
    lines=[s['status'], f"Selected models: A={s['model_a']}, B={s['model_b']}",f"Shared-chain alignment: {s['matched_CA']}/{s['shared_length']} Cα ({s['coverage']:.1%}); RMSD {s['CA_RMSD_A']:.3f} Å (warning >{options.rmsd_warning:g} Å)",
           f"Largest Cα deviation: {s['CA_max_deviation_A']:.3f} Å at sequence position {s['CA_max_deviation_position']}"]
    for row in result['sites']:
        value='NA' if row['site_RMSD_A'] is None else f"{row['site_RMSD_A']:.3f} Å"
        lines.append(f"Binding site {row['input']}:{row['partner_chain']}: {value}, {row['matched_site_CA']}/{row['observed_site_CA']} observed site Cα matched (warning >{options.site_rmsd_warning:g} Å)")
    lines.extend([f"Cross-partner clashes: {s['cross_partner_clashes']} atom pairs; backbone–backbone {s['cross_partner_backbone_clashes']}; closest heavy atoms {s['min_partner_distance_A']:.3f} Å; maximum overlap {max(0,s['cross_partner_max_overlap_A']):.3f} Å",
                  f"Shared-chain replacement: {s['shared_replacement_introduced_clashes']} introduced / {s['shared_replacement_worsened_clashes']} worsened overlapping atom pairs; maximum overlap {max(0,s['shared_replacement_max_overlap_A']):.3f} Å",
                  f'Clash screen: radius overlap >{options.clash_overlap:g} Å; heavy atoms only, not MolProbity clashscore.',
                  'Alignment thresholds are screening heuristics, not biological acceptance criteria.'])
    lines.extend('WARNING: '+w for w in result['warnings'])
    lines.append('No relaxation performed. This assessment does not establish binding or full stereochemical validity.')
    return '\n'.join(lines)


def table(path,rows,columns):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in rows:
            writer.writerow({k:('NA' if row.get(k) is None else row[k]) for k in columns})


def model_numbers(locator):
    with open_structure_fold(locator) as fold:
        numbers=[model_number(n) if model_number(n) is not None else 0 for n in fold.models]
        if len(set(numbers)) != len(numbers):
            raise ValueError('duplicate model numbers')
        return sorted(numbers)


def run(locator_a,locator_b,output,*,model_a=0,model_b=0,all_models=False,options=None,stdout=None):
    import sys
    options=options or Options();stdout=stdout or sys.stdout
    output=Path(output).resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError('output must be a new or empty directory; existing results are never overwritten')
    for locator in (locator_a,locator_b):
        path=Path(str(locator).split('::',1)[0]).resolve()
        if path==output or output in path.parents:
            raise ValueError('output must not contain or replace an input')
    first=load_structure(locator_a,model_a);second=load_structure(locator_b,model_b)
    primary=compare(first,second,options)
    text,mapping=complex_cif(first,second,primary)
    # Validate full sequence metadata and hashing before publishing any directory.
    from sieve.alphafold_manifest import chain_sequences
    sequences=chain_sequences(text)
    if sequences != {c:v['sequence'] for c,v in mapping.items()}:
        raise ValueError('assembled CIF sequence round-trip failed')
    comparisons=[primary['summary']]
    if all_models:
        def load_models(locator, selected):
            models = {}
            for n in model_numbers(locator):
                try:
                    models[n] = selected if n == selected.number else load_structure(locator, n)
                except (ValueError, KeyError, TypeError, OSError) as exc:
                    models[n] = str(exc)
            return models
        left=load_models(locator_a,first)
        right=load_models(locator_b,second)
        comparisons=[]
        for a,b in itertools.product(left,right):
            if a==model_a and b==model_b:
                comparisons.append(primary['summary']);continue
            try:
                if isinstance(left[a],str) or isinstance(right[b],str):
                    raise ValueError(left[a] if isinstance(left[a],str) else right[b])
                comparisons.append(compare(left[a],right[b],options)['summary'])
            except ValueError as exc:
                comparisons.append(dict(model_a=a,model_b=b,status='ERROR — cannot assemble',warnings=str(exc)))
    identity=dict(algorithm='rigid-shared-sequence-v1',inputs=[dict(cif_sha256=digest(s.content),model=s.number) for s in (first,second)],
                  shared_chains=[primary['summary']['shared_chain_a'],primary['summary']['shared_chain_b']],
                  retained_shared_input='a',alignment='all-matched-sequence-position-CA',output_sha256=digest(text))
    provenance=dict(format=ASSEMBLY_FORMAT,source_type='assembled',tool_version=VERSION,
                    created_at=datetime.now(timezone.utc).isoformat(),fold_id=assembly_id(identity),sequence_id=sequence_id(sequences.values()),
                    identity=identity,cif_sha256=digest(text),model_id=0,chains=mapping,
                    inputs=[dict(locator=s.locator,model=s.number,cif_sha256=digest(s.content),identity=s.identity,history=s.metadata) for s in (first,second)],
                    rotation=primary['rotation'].tolist(),translation=primary['translation'].tolist(),
                    transform_convention='row vectors: coordinates @ rotation + translation',settings=asdict(options),vdw_radii_A=VDW,
                    assessment=dict(status=primary['summary']['status'],warnings=primary['warnings']),
                    confidence='No assembled PAE, ipTM, pTM or seed; known per-atom pLDDT values are inherited only.')
    output.parent.mkdir(parents=True,exist_ok=True)
    temporary=Path(tempfile.mkdtemp(prefix='.overlap-',dir=output.parent))
    try:
        (temporary/COMPLEX_FILE).write_text(text)
        (temporary/'provenance.json').write_text(json.dumps(provenance,indent=2,allow_nan=False)+'\n')
        table(temporary/'comparisons.tsv',comparisons,list(primary['summary']))
        specs=[('shared_chain_deviations', 'deviations',['sequence_position','CA_deviation_A','source_plddt_a','source_plddt_b']),
               ('shared_chain_sites','sites',['input','partner_chain','matched_site_CA','observed_site_CA','site_coverage','site_RMSD_A']),
               ('alignment_sensitivity','sensitivity',['alignment','matched_CA','fit_RMSD_A','min_partner_distance_A','clash_atom_pairs','status']),
               ('cross_partner_metrics','cross',['chain_a','chain_b','min_distance_A','max_overlap_A','clash_atom_pairs','clash_residue_pairs','backbone_backbone_clashes','contact_residue_pairs']),
               ('cross_partner_clashes','cross_clashes',['chain_a','chain_b','residue_a','atom_a','residue_b','atom_b','distance_A','overlap_A','backbone_backbone']),
               ('shared_interface_changes','changes',['partner_chain_b','shared_chain_a','shared_chain_b','original_clashes','assembled_clashes','introduced_clashes','worsened_clashes','resolved_clashes','original_max_overlap_A','assembled_max_overlap_A','contacts_before','contacts_after','contacts_lost','contacts_gained']),
               ('shared_interface_clashes','replacement_clashes',['partner_chain_b','change','residue_a','atom_a','residue_b','atom_b','distance_A','overlap_A','backbone_backbone']),
               ('shared_interface_contact_changes','contact_changes',['partner_chain_b','change','shared_position','partner_position'])]
        for name,key,cols in specs: table(temporary/(name+'.tsv'),primary[key],cols)
        input_summary = '\n'.join(
            f"Input {label}: Fold ID {source.identity.get('fold_id') or 'NA'} | "
            f"Sequence ID {source.identity.get('sequence_id') or 'NA'} | Model {source.number}"
            for label, source in [('A', first), ('B', second)])
        summary=input_summary + '\n' + quick_summary(primary,options)
        doc=Path(__file__).resolve().parents[1]/'docs/alphafold-overlap.md'
        readme=f"# Assembly {provenance['fold_id']}\n\nSequence ID: {provenance['sequence_id']}\n\n```text\n{summary}\n```\n\n"
        readme+='## Output chain mapping\n\n'+''.join(f"- {c}: input {v['input'].upper()}, original chain {v['original_chain']}, {len(v['sequence'])} residues.\n" for c,v in mapping.items())
        readme+='\n'+doc.read_text()
        (temporary/'README.md').write_text(readme)
        # Includes metadata, sequence, and CIF checksum validation on the final artifact.
        load_structure(str(temporary),0)
        os.replace(temporary,output)
    finally:
        if temporary.exists(): shutil.rmtree(temporary)
    print(f"Fold ID: {provenance['fold_id']}\nSequence ID: {provenance['sequence_id']}\nLocator: {output}\n"+summary,file=stdout)
    if all_models:
        print(f"Model combinations: {len(comparisons)} assessed; {sum(r.get('status','').startswith('WARNING') for r in comparisons)} combinations with warnings; {sum(r.get('status','').startswith('ERROR') for r in comparisons)} errors.",file=stdout)
        print(f"All combinations: {output / 'comparisons.tsv'}", file=stdout)
    return provenance
