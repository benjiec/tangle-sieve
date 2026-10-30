"""Standalone death-family patch annotation using shared geometric detection."""
import sys
from pathlib import Path
import subprocess
import zipfile

from sieve.alphafold_cache import write_cache
from sieve.alphafold_files import discover, fold_locator, model_number
from sieve.alphafold_manifest import read_manifest_groups
from sieve.alphafold_structure import open_structure_fold
from sieve.alphafold_pdockq2_patches import input_models, analyze, Options
from sieve.alphafold_pdockq2 import parse_model
from sieve.alphafold_secondary import find_dssp, load_secondary
from sieve.alphafold_death_rules import read_domains, read_rules, classify, candidate_interface

SETTINGS = ['distance cutoff','max bridge residues','min residues per chain','min contacts','max separator checks']
COLUMNS = ['fold_locator','fold_id','sequence_id','model_id','patch_id','patch label','chain 1','chain 2',
           'sequence name 1','sequence name 2','domain 1','domain 2','patch residues 1','patch residues 2',
           'interface structure 1','interface structure 2','candidate face 1','candidate face 2',
           'candidate interface','status 1','status 2','rules_id',*SETTINGS]
ERRORS = (OSError,ValueError,KeyError,TypeError,RuntimeError,zipfile.BadZipFile,subprocess.SubprocessError)


def jobs_for(request, error):
    path=Path(request.locator)
    if '::' in request.locator or path.suffix.lower()=='.cif':
        locators=[request.locator]
    elif path.is_dir() and any(model_number(p.name) is not None for p in path.iterdir() if p.is_file()):
        locators=[str(path.resolve())]
    else:
        locators=[fold_locator(f.source_path,f.fold_path) for f in discover(path,error)]
    for locator in locators:
        try:
            with open_structure_fold(locator) as fold:
                models=[model_number(n) for n in fold.models]
            if len(models)!=len(set(models)):
                raise ValueError('duplicate model numbers')
            if request.model_id is not None:
                if request.model_id not in models:
                    raise ValueError(f'model {request.model_id} not found')
                models=[request.model_id]
            for model in sorted(models,key=lambda n:-1 if n is None else n):
                yield locator,model
        except ERRORS as exc:
            error(f'{locator}: {exc}')


def run(requests, manifest, domains_path, rules_path, output, *, options=None, chains=None,
        full_data=None, dssp=None, stderr=None):
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
            _,codes,segments,status=load_secondary(record[2],executable)
            if status['status']!='assigned':
                raise ValueError(f'DSSP unavailable: {status}')
            current=[]
            for patch in patches:
                if patch['scope']!='patch' or patch['status']!='leaf': continue
                row=dict(fold_locator=locator,fold_id=metadata['fold_id'],sequence_id=metadata['sequence_id'],
                         model_id=number,patch_id=patch['patch_id'],rules_id=rules_id,
                         **{'patch label':patch['patch id']})
                row.update({k:patch[k] for k in SETTINGS})
                for side in (1,2):
                    chain=patch[f'chain {side}'];name=assignment[chain];domain=domains[name]
                    selected={int(n) for n in patch[f'interface residues {side}'].split(',') if n}
                    faces,_,eligibility=classify(chain,selected,codes,segments,rules[domain])
                    structure='; '.join(f"{s['label']}[{s['code']}]:"+','.join(map(str,sorted(selected & set(s['residues']))))
                        for s in segments if s['chain']==chain and selected & set(s['residues']))
                    row.update({f'chain {side}':chain,f'sequence name {side}':name,f'domain {side}':domain,
                        f'patch residues {side}':','.join(map(str,sorted(selected))),
                        f'interface structure {side}':structure or None,f'candidate face {side}':faces,
                        f'status {side}':eligibility})
                row['candidate interface']=candidate_interface(row['candidate face 1'],row['candidate face 2'])
                current.append(row)
            output_rows.extend(current)
        except ERRORS as exc:
            error(f'{locator} model {number}: {exc}')
    write_cache(target,COLUMNS,({k:'NA' if v is None else v for k,v in row.items()} for row in output_rows))
    return int(failed)
