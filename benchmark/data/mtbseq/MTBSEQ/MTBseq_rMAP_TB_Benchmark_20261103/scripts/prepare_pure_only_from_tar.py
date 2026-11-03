#!/usr/bin/env python3
"""Extract only two pure-control MTBseq Called and Position_Tables files from
existing full MTBseq results.tar.gz; do not extract BAMs or FASTQ files.

Usage:
  python3 prepare_pure_only_from_tar.py \
    --archive ~/Downloads/MTBSEQ/mtbseq_out/mtbseq_outputs.tar.gz \
    --output ~/Downloads/MTBseq_pure_only
"""
import argparse, tarfile, shutil
from pathlib import Path

PARENTS=('PURE-ERR038255-100x','PURE-ERR13260175-100x')

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--archive',type=Path,required=True)
 p.add_argument('--output',type=Path,required=True)
 args=p.parse_args()
 arc=args.archive.expanduser().resolve(); out=args.output.expanduser().resolve()
 if not arc.is_file():raise SystemExit(f'Archive not found: {arc}')
 out.mkdir(parents=True,exist_ok=True)
 seen=[]
 with tarfile.open(arc,'r|gz') as tar:
  for member in tar:
   if not member.isfile():continue
   parts=Path(member.name).parts
   matching_dir=next((x for x in ('Called','Position_Tables') if x in parts),None)
   if matching_dir is None:continue
   name=Path(member.name).name
   if not name.endswith('.tab') or not any(name.startswith(x+'_lib1') for x in PARENTS):continue
   if matching_dir=='Called' and 'gatk_position_variants' not in name:continue
   if matching_dir=='Position_Tables' and 'gatk_position_table' not in name:continue
   path=out/matching_dir/name
   path.parent.mkdir(parents=True,exist_ok=True)
   source=tar.extractfile(member)
   if source is None:continue
   if path.exists():raise SystemExit(f'Will not overwrite existing file: {path}')
   with source,open(path,'wb') as dest:shutil.copyfileobj(source,dest)
   seen.append((matching_dir,name,path.stat().st_size))
   print(f'Extracted {matching_dir}/{name} ({path.stat().st_size:,} bytes)',flush=True)
 with open(out/'samples.tsv','w') as f:
  for n in PARENTS:f.write(f'{n}\tlib1\n')
 called=[x for x in seen if x[0]=='Called'];ptable=[x for x in seen if x[0]=='Position_Tables']
 print(f'Found Called={len(called)} and Position_Tables={len(ptable)}')
 if len(called)!=2 or len(ptable)!=2:
  raise SystemExit('Incomplete input extraction. Inspect archive names; expected 2 Called and 2 Position_Tables. Do not run TBjoin yet.')
 print(f'Pure-only input staging successful: {out}')

if __name__=='__main__':main()
