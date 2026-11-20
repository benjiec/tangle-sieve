"""Standalone death-family patch annotation using shared geometric detection."""
import sys
from pathlib import Path

from sieve.alphafold_cache import write_cache
from sieve.alphafold_manifest import read_manifest_groups
from sieve.alphafold_patch_metrics import ERRORS, jobs_for
from sieve.alphafold_pdockq2_patches import input_models, analyze, Options
from sieve.alphafold_pdockq2 import parse_model
from sieve.alphafold_secondary import find_dssp, load_secondary
from sieve.alphafold_death_rules import read_domains, read_rules, classify, candidate_interface

from sieve.alphafold_helix_normalization import normalize_helices, validate_kink_settings

KINK_SETTINGS = ["merge kinked helices", "kink connector length", "kink angle"]
SETTINGS = ['distance cutoff','max bridge residues','min residues per chain','min contacts','max separator checks']
SIDE_COLUMNS = ['chain', 'sequence name', 'domain', 'candidate face', 'status']
COLUMNS = ['fold_id','sequence_id','fold_locator','model_id','patch_id','patch label',
           'sequence_index',*SIDE_COLUMNS,'candidate interface','rules_id',*SETTINGS,*KINK_SETTINGS]
def run(requests, manifest, domains_path, rules_path, output, *, options=None, chains=None,
        full_data=None, dssp=None, stderr=None, merge_kinked_helices=False,
        kink_connector_length=1, kink_angle=70):
    validate_kink_settings(kink_connector_length,kink_angle)
    stderr=sys.stderr if stderr is None else stderr
    options=options or Options()
    if chains is not None and (len(chains)!=2 or chains[0]==chains[1]):
        raise ValueError('--chains requires two distinct chains')
    requests=list(requests)
    target=Path(output).resolve()
    protected=[Path(p).resolve() for p in (manifest,domains_path,rules_path)]
    protected += [Path(r.locator.split('::',1)[0]).resolve() for r in requests]
    if full_data: protected.append(Path(full_data).resolve())
    if target in protected or target.suffix.lower()!='.tsv':
        raise ValueError('output must be a separate .tsv file')
    groups=read_manifest_groups(manifest)
    domains=read_domains(domains_path)
    rules,rules_id=read_rules(rules_path)
    executable=find_dssp(dssp)
    if executable is None:
        raise ValueError('DSSP is required; use --dssp or SIEVE_DSSP')
    failed=False
    def error(message):
        nonlocal failed
        failed=True
        print(f'ERROR: {message}',file=stderr)
    # Resolve and validate every selected model's mappings before annotating or replacing output.
    prepared=[]; seen=set()
    for request in requests:
        for locator,number in jobs_for(request,error):
            key=(locator,number)
            if key in seen: continue
            seen.add(key)
            # A direct CIF belongs to its parent manifest fold.
            lookup=str(Path(locator).parent.resolve()) if Path(locator).suffix.lower()=='.cif' and '::' not in locator else locator
            if lookup not in groups:
                raise ValueError(f'locator missing from manifest: {lookup}')
            rows=groups[lookup]
            assignment={r['input chain']:r['input description'] for r in rows}
            chosen=chains or list(assignment)
            for chain in chosen:
                if chain not in assignment:
                    raise ValueError(f'chain {chain} missing from manifest: {lookup}')
                name=assignment[chain]
                if name not in domains:
                    raise ValueError(f'sequence missing from domain TSV: {name}')
                if domains[name] not in rules:
                    raise ValueError(f'domain has no rules: {domains[name]} ({name})')
            seeds={}
            try:
                records=list(input_models(locator,full_data,selected_model=number,model_seeds=seeds))
                if len(records)!=1: raise ValueError('expected exactly one model')
                record=records[0]
                residues,_=parse_model(record[2],record[3])
            except ERRORS as exc:
                error(f'{locator} model {number}: {exc}');continue
            if set(residues)!=set(assignment):
                raise ValueError(f'model/manifest chain mismatch: {lookup} model {number}')
            # Identity checks use the same implementation as patch scoring.
            from sieve.alphafold_composition import chain_compositions
            from sieve.alphafold_manifest import sequence_id, fold_id
            try:
                compositions=chain_compositions(record[2],record[3].get('_job_request'))
                seqid=sequence_id(compositions.values())
                seed=seeds.get((record[0],'' if number is None else number))
                fid=record[3]['_assembly']['fold_id'] if '_assembly' in record[3] else (fold_id(seed,compositions.values()) if seed is not None else None)
            except ERRORS as exc:
                error(f'{locator} model {number}: {exc}');continue
            if seqid!=rows[0]['sequence_id'] or fid!=rows[0]['fold_id']:
                raise ValueError(f'model/manifest identity mismatch: {lookup} model {number}')
            prepared.append((lookup,number,record,seeds,assignment,rows[0]))
    output_rows=[]
    for locator,number,record,seeds,assignment,metadata in prepared:
        try:
            patches=analyze(locator,options=options,chains=chains,model_records=[record],model_seeds=seeds)
            model,codes,segments,status=load_secondary(record[2],executable)
            if status['status']!='assigned':
                raise ValueError(f'DSSP unavailable: {status}')
            normalized={}
            if merge_kinked_helices:
                for chain in (chains or assignment):
                    normalized[chain]=normalize_helices(model,chain,codes,segments,
                        connector_length=kink_connector_length,angle=kink_angle)
            current=[]
            for patch in patches:
                if patch['scope']!='patch' or patch['status']!='leaf': continue
                row=dict(fold_locator=locator,fold_id=metadata['fold_id'],sequence_id=metadata['sequence_id'],
                         model_id=number,patch_id=patch['patch_id'],rules_id=rules_id,
                         **{'patch label':patch['patch id']})
                row.update({k:patch[k] for k in SETTINGS})
                row.update(dict(zip(KINK_SETTINGS,(merge_kinked_helices,kink_connector_length,kink_angle))))
                sides=[]
                for side in (1,2):
                    chain=patch[f'chain {side}'];name=assignment[chain];domain=domains[name]
                    selected={int(n) for n in patch[f'interface residues {side}'].split(',') if n}
                    chain_segments,report=normalized.get(chain,(segments,""))
                    faces,_,eligibility=classify(chain,selected,codes,chain_segments,rules[domain])
                    if report: eligibility += "; " + report
                    sides.append({'sequence_index':side,'chain':chain,'sequence name':name,
                                  'domain':domain,'candidate face':faces,'status':eligibility})
                row['candidate interface']=candidate_interface(*(side['candidate face'] for side in sides))
                current.extend(dict(row, **side) for side in sides)
            output_rows.extend(current)
        except ERRORS as exc:
            error(f'{locator} model {number}: {exc}')
    write_cache(target,COLUMNS,({k:'NA' if v is None else v for k,v in row.items()} for row in output_rows))
    return int(failed)
