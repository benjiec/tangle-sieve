#!/usr/bin/env python3
"""Assemble complexes through exactly one uniquely shared protein sequence."""
import argparse
import sys
import zipfile
from sieve.alphafold_overlap import Options, run


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('locator_a',help='AF locator, assembled directory, or sequence-annotated CIF')
    parser.add_argument('locator_b',help='AF locator, assembled directory, or sequence-annotated CIF')
    parser.add_argument('--output',required=True,help='new or empty output directory')
    parser.add_argument('--model-a',type=int,default=0)
    parser.add_argument('--model-b',type=int,default=0)
    parser.add_argument('--all-models',action='store_true',help='also assess all combinations; save only the selected pair')
    for name,default in [('clash-overlap',.4),('rmsd-warning',2),('site-rmsd-warning',1.5),('coverage-warning',.9),('site-distance',8),('contact-distance',4.5)]:
        parser.add_argument('--'+name,type=float,default=default)
    args=parser.parse_args(argv)
    try:
        options=Options(**{name:getattr(args,name) for name in Options.__dataclass_fields__})
        run(args.locator_a,args.locator_b,args.output,model_a=args.model_a,model_b=args.model_b,all_models=args.all_models,options=options)
    except (OSError,ValueError,KeyError,TypeError,zipfile.BadZipFile) as exc:
        print(f'ERROR — cannot assemble: {exc}',file=sys.stderr)
        return 1
    return 0


if __name__=='__main__':
    sys.exit(main())
