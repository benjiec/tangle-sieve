"""Shared DSSP discovery and structure loading."""
import io
import os
import shutil
from Bio.PDB import MMCIFParser
from sieve.interface_description import secondary_structure


def find_dssp(explicit=None):
    configured = explicit if explicit is not None else os.environ.get('SIEVE_DSSP')
    if configured is not None:
        executable = shutil.which(os.path.expanduser(configured)) if configured else None
        if executable is None:
            source = '--dssp' if explicit is not None else 'SIEVE_DSSP'
            raise ValueError(f'{source} does not identify an executable: {configured!r}')
        return executable
    return shutil.which('mkdssp')



def load_secondary(cif, executable):
    model = MMCIFParser(QUIET=True).get_structure("faces", io.StringIO(cif))[0]
    for chain in model:
        seen = set()
        for residue in chain:
            if "CA" not in residue:
                continue
            if residue.id[2] != " " or residue.id[1] in seen:
                raise ValueError("secondary structure requires unique residue numbers without insertion codes")
            seen.add(residue.id[1])
    codes, segments, status = secondary_structure(model, cif, executable)
    return model, codes, segments, status
