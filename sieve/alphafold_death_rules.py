"""User-specified death-domain-family coverage heuristics; not topology validation."""
import csv
import re
from sieve.alphafold_structure import canonical, digest
from sieve.alphafold_dd_faces import candidate_interface

FACES = ('Ia','Ib','IIa','IIb','IIIa','IIIb')
ELEMENTS = {f'alpha{i}' for i in range(1,7)} | {f'loop{i}-{i+1}' for i in range(1,6)}


def read_domains(path):
    with open(path, encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t', strict=True)
        columns = reader.fieldnames or []
        lowered = [c.lower() for c in columns]
        if len(set(lowered)) != len(lowered) or not {'sequence name','domain'} <= set(lowered):
            raise ValueError('domain TSV requires unique Sequence Name and Domain headers (case-insensitive)')
        name, domain = (columns[lowered.index(c)] for c in ('sequence name','domain'))
        result = {}
        for row in reader:
            if None in row or any(v is None for v in row.values()) or not row[name].strip() or not row[domain].strip():
                raise ValueError('incomplete domain assignment')
            if row[name] in result:
                raise ValueError(f'duplicate sequence assignment: {row[name]}')
            result[row[name]] = row[domain]
        return result


def read_rules(path):
    fields = ('domain','face','element','min_residues','requirement_type','notes')
    records = []
    rules = {}
    with open(path, encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t', strict=True)
        columns = reader.fieldnames or []
        if len(set(columns)) != len(columns) or not set(fields) <= set(columns):
            raise ValueError('rule TSV requires: '+', '.join(fields))
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError('incomplete rule row')
            domain,face,element,minimum,kind,notes = (row[c] for c in fields)
            if (not domain.strip() or face not in FACES or element not in ELEMENTS or
                not re.fullmatch('[0-9]+',minimum) or kind not in ('mandatory','alternative_or','optional')):
                raise ValueError(f'invalid rule: {row}')
            minimum = int(minimum)
            if kind != 'optional' and minimum < 1:
                raise ValueError('mandatory/alternative rules require positive minimum residues')
            group = rules.setdefault(domain, {}).setdefault(face, [])
            if any(r['element'] == element for r in group):
                raise ValueError(f'duplicate rule for {domain}/{face}/{element}')
            rule = dict(element=element, min_residues=minimum, requirement_type=kind, notes=notes)
            group.append(rule)
            records.append(dict(domain=domain,face=face,**rule))
    if not rules:
        raise ValueError('rule TSV is empty')
    for domain,faces in rules.items():
        for face,group in faces.items():
            if not any(r['requirement_type'] != 'optional' for r in group):
                raise ValueError(f'{domain}/{face} has only optional rules')
    records.sort(key=lambda r:(r['domain'],FACES.index(r['face']),r['element']))
    return rules, digest(canonical(records))


def elements_for_chain(chain, codes, segments):
    helices = sorted((s for s in segments if s['chain']==chain and s['code']=='H'),key=lambda s:s['start'])
    if len(helices)!=6:
        return None, f'ineligible: expected 6 H segments, found {len(helices)}'
    if any(codes.get((chain,n),'?')=='?' for n in range(helices[0]['start'],helices[-1]['end']+1)):
        return None, 'ineligible: missing DSSP assignments'
    elements = {f'alpha{i}':set(s['residues']) for i,s in enumerate(helices,1)}
    for i in range(1,6):
        interval=set(range(helices[i-1]['end']+1,helices[i]['start']))
        elements[f'loop{i}-{i+1}']=interval if all(codes[(chain,n)] in ('-','T','S','G') for n in interval) else set()
    return elements, 'eligible'


def classify(chain, selected, codes, segments, rules):
    elements,status=elements_for_chain(chain,codes,segments)
    if elements is None:
        return None, [], status
    ranked=[]; evidence=[]
    for face in FACES:
        if face not in rules:
            continue
        mandatory=[]; alternatives=[]; active=[]
        for rule in rules[face]:
            positions=elements[rule['element']]
            count=len(selected & positions)
            item=dict(rule,contacting_residues=count,element_residues=len(positions),
                      passed=count>=rule['min_residues'])
            evidence.append(dict(face=face,**item))
            coverage=count/len(positions) if positions else 0
            if rule['requirement_type']=='optional':
                continue
            active.append(count)
            (mandatory if rule['requirement_type']=='mandatory' else alternatives).append((item['passed'],coverage))
        if not any(active):
            continue
        full=all(p for p,_ in mandatory) and (not alternatives or any(p for p,_ in alternatives))
        coverage=[c for _,c in mandatory]+([max(c for _,c in alternatives)] if alternatives else [])
        ranked.append((not full,-min(coverage),-sum(coverage)/len(coverage),FACES.index(face),face if full else face+' (partial)'))
    return '; '.join(r[-1] for r in sorted(ranked)) or 'none', evidence, status
