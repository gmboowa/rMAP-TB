#!/usr/bin/env python3
"""Reproducible original-output audit for 18-library rMAP-TB / MTBseq comparison.
Uses Python standard library only. Reads original uploaded tables; does NOT rerun MTBseq.
"""
import csv, re, json, hashlib, statistics, math, shutil, zipfile, collections, datetime, os
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
INPUT=Path(os.environ.get('BENCHMARK_INPUT_DIR', str(BASE/'baseline_original_outputs')))
OUT=Path(os.environ.get('BENCHMARK_OUTPUT_DIR', str(BASE/'reproduced_audit')))
RAW=OUT/'baseline_original_outputs'; ANAL=OUT/'analysis'; SCRIPTS=OUT/'scripts'; REPO=OUT/'reports'; INPUTS=OUT/'inputs'
for d in (RAW,ANAL,SCRIPTS,REPO,INPUTS): d.mkdir(parents=True,exist_ok=True)
originals=['MTBseq_Amend_Audit.zip','MTBseq_Extra_Audit.zip','mtbseq_classification.tsv','mtbseq_groups.tsv','mtbseq_distance_matrix.tsv','mtbseq_statistics.tsv','mtbseq_resistance.tsv','collection.log','mtbseq_run(1).log','sample_map.tsv','status(1).txt','tb_surveillance_metadata(1).tsv','clustering_eligibility.tsv','pairwise_snp_distance_pairs.tsv','pairwise_snp_distance_matrix.tsv','resistance_profile_summary.tsv','run_metadata(1).txt','script(20261103-041350)']
for name in originals:
    if not (INPUT/name).is_file(): print('MISSING',name);continue
    shutil.copy2(INPUT/name,RAW/name)
rawfiles=[]
for name in sorted(originals):
 src=RAW/name
 if src.is_file():
  h=hashlib.sha256(src.read_bytes()).hexdigest()
  rawfiles.append({'filename':src.name,'size_bytes':src.stat().st_size,'sha256':h})
(RAW/'SHA256SUMS.txt').write_text(''.join(f'{x["sha256"]}  {x["filename"]}\n' for x in rawfiles))
# parse files from upstream ZIP and Amend ZIP
def readzip_table(zfile,suffix):
 with zipfile.ZipFile(INPUT/zfile) as z:
  matches=[n for n in z.namelist() if n.endswith(suffix) and not n.startswith('__MACOSX/')]
  if len(matches)!=1:raise ValueError(f'{zfile} {suffix} matches {matches}')
  import io
  return list(csv.reader(io.StringIO(z.read(matches[0]).decode('utf-8-sig')),delimiter='\t'))
amended=readzip_table('MTBseq_Extra_Audit.zip','_amended.tab')
phylo=readzip_table('MTBseq_Amend_Audit.zip','_amended_u95_phylo.tab')
final=readzip_table('MTBseq_Amend_Audit.zip','_amended_u95_phylo_w12.tab')
removed=readzip_table('MTBseq_Amend_Audit.zip','_amended_u95_phylo_w12_removed.tab')
assert len(amended[1])==200 and len(amended)==3150
assert len(phylo)==248 and len(final)==239 and len(removed)==11
samples=[amended[0][18+8*i].strip().removesuffix('_lib1') for i in range(18)]
canon=lambda s:re.sub(r'^\d{4}_','',s.strip().strip("'").removesuffix('_lib1').replace('-','_'))
names=[canon(s) for s in samples]
assert len(set(names))==18
# Load tsv by stdlib

def rdict(file):
 with open(INPUT/file,newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
meta={canon(r['sample']):r for r in rdict('tb_surveillance_metadata(1).tsv')}
flag={canon(r['sample']):r for r in rdict('clustering_eligibility.tsv')}
classif={canon(r['SampleID']):r for r in rdict('mtbseq_classification.tsv') if r.get('SampleID')}
stat={canon(r['SampleID']):r for r in rdict('mtbseq_statistics.tsv') if r.get('SampleID')}
assert set(names)==set(meta)==set(flag)==set(classif)==set(stat)
original_mat=INPUT/'mtbseq_distance_matrix.tsv'
rows=list(csv.reader(open(original_mat),delimiter='\t'))
assert len(rows)==18, len(rows)
mtb_dist={}
for i,row in enumerate(rows):
 assert canon(row[0])==names[i], (i,row[0],names[i])
 vals=[x.strip() for x in row[1:] if x.strip()]
 assert len(vals)==i,(i,len(vals))
 for j,v in enumerate(vals):mtb_dist[tuple(sorted((names[j],names[i])))] = int(v)
assert len(mtb_dist)==153
rmap_pairs={}
for r in rdict('pairwise_snp_distance_pairs.tsv'):
 a,b=canon(r['sample1']),canon(r['sample2'])
 rmap_pairs[tuple(sorted((a,b)))]=r
assert set(mtb_dist)==set(rmap_pairs)

def each_site(rows):
 for r in rows[2:]:
  yield r

def key(r):return (r[0].strip(),r[1].strip())
physet={key(r) for r in each_site(phylo)}
finset={key(r) for r in each_site(final)}
rmset={key(r) for r in each_site(removed)}
assert len(physet)==246 and len(finset)==237 and len(rmset)==9 and finset|rmset==physet
FULL18=0
for r in each_site(amended):
 if r[7].strip()=='18':FULL18+=1
assert FULL18==291,FULL18

def to_float(t):
 try:return float(t)
 except (TypeError,ValueError):return 0.0

def call(r,i):
 offset=18+8*i
 allele=r[offset+1].strip().upper()
 typ=r[offset].strip().upper()
 forward,reverse,qual,freq,depth=map(to_float,r[offset+2:offset+7])
 return dict(allele=allele,type=typ,forward=forward,reverse=reverse,qual20=qual,freq=freq,depth=depth,
   pass_qc=(allele in 'ACGT' and len(allele)==1 and typ in ('NONE','SNP') and forward>=4 and reverse>=4 and qual>=4 and freq>=75))

pure_idx=[names.index('PURE_ERR038255_100x'),names.index('PURE_ERR13260175_100x')]
assert pure_idx==[16,17]
site_rows=[]; pure_sites=[]; pure_marker_keys=set(); marker_rmap_keys=set()
for r in each_site(amended):
 k=key(r);c1,c2=call(r,16),call(r,17)
 pair_pass=c1['pass_qc'] and c2['pass_qc'] and r[1].strip()=='0'
 diff=pair_pass and c1['allele']!=c2['allele']
 if diff: pure_marker_keys.add(k)
 site_rows.append({'position':int(r[0]),'insindex':r[1].strip(),'reference_base':r[2].strip(),
  'gene':r[4].strip(),'category':r[14].strip(),'resistance_annotation':r[15].strip(),
  'count_unamb_all18':int(r[7]),'pct_unamb_all18':float(r[8]),'all18_fully_callable':int(r[7])==18,
  'in_mtbseq_u95_phylo':k in physet,'in_mtbseq_final_w12':k in finset,
  'removed_window12':k in rmset,'pure_both_pass_call_QC':pair_pass,
  'pure_control_allele_difference':diff})
 if pair_pass:
  pure_sites.append({'position':int(r[0]),'insindex':r[1].strip(),'reference_base':r[2].strip(),
  'allele_ERR038255':c1['allele'],'allele_ERR13260175':c2['allele'],
  'freq_ERR038255':round(c1['freq'],2),'freq_ERR13260175':round(c2['freq'],2),
  'covf_ERR038255':int(c1['forward']),'covr_ERR038255':int(c1['reverse']),
  'covf_ERR13260175':int(c2['forward']),'covr_ERR13260175':int(c2['reverse']),
  'different':int(diff),'mtbseq_final_w12_included':int(k in finset),
  'mtbseq_u95_phylo_included':int(k in physet),'category':r[14].strip(),'gene':r[4].strip(),
  'note':'candidate site from original 18-sample joint table; pair-only post hoc QC, NOT new TBjoin'})
assert len(site_rows)==3148
assert len(pure_sites)==2201,len(pure_sites)
assert sum(r['different'] for r in pure_sites)==1278
assert not pure_marker_keys & physet # none of pair's informative markers survives

def write_tsv(file, rows, fields=None):
 if not rows and not fields:raise ValueError(file)
 fields=fields or list(rows[0].keys())
 with open(ANAL/file,'w',newline='',encoding='utf-8') as fh:
  w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
write_tsv('mtbseq_site_mask_audit.tsv',site_rows)
write_tsv('pure_control_pair_only_candidate_sites.tsv',pure_sites)
filter_summary=[
 {'stage':'All candidate rows in 18-sample amended joint table','retained_positions':3148,'removed_since_prior':0,'note':'Candidate variable/variant positions, not genome-wide coordinates'},
 {'stage':'18/18 unambiguous calls (u95 rounding)','retained_positions':291,'removed_since_prior':2857,'note':'With 18 samples, 17/18 =94.44% is below 95%'},
 {'stage':'MTBseq phylogenetic eligibility','retained_positions':246,'removed_since_prior':45,'note':'Official TBamend _u95_phylo stage; additional region, allele-type and other exclusions'},
 {'stage':'MTBseq 12-base window','retained_positions':237,'removed_since_prior':9,'note':'Official _u95_phylo_w12 stage; nine removed sites'},
]
write_tsv('filtering_stage_counts.tsv',filter_summary)
# validate final FASTA and the full 153-pair matrix
with zipfile.ZipFile(INPUT/'MTBseq_Amend_Audit.zip') as z:
 seq_entry=next(n for n in z.namelist() if n.endswith('_u95_phylo_w12.plainIDs.fasta'))
 fasta=z.read(seq_entry).decode()
seq={}
for record in fasta.split('>'):
 if not record.strip():continue
 h,*lines=record.strip().splitlines()
 seq[canon(h)]=''.join(lines).upper()
assert set(seq)==set(names)
assert set(map(len,seq.values()))=={237}
reconstruct=[]
for i,a in enumerate(names):
 for b in names[i+1:]:
  dist=sum(x!=y for x,y in zip(seq[a],seq[b]))
  expected=mtb_dist[tuple(sorted((a,b)))]
  assert dist==expected,(a,b,dist,expected)
  p=rmap_pairs[tuple(sorted((a,b)))]
  reconstruct.append({'sample_a':a,'sample_b':b,'mtbseq_filtered_237_site_distance':expected,
  'mtbseq_reconstructed_237_site_distance':dist,'mtbseq_reproduced_exactly':'YES',
  'rmap_full_core_snp_distance':int(p['snp_distance']), 'rmap_comparable_sites':int(p['comparable_sites']),
  'rmap_cluster_class':p['cluster_class'],'site_mask_aligned':'NO — distinct denominators',
  'reporting_restriction':'Do not calculate diagnostic accuracy from these unmatched masks'})
write_tsv('all153_pairs_DESCRIPTIVE_UNMATCHED_MASKS.tsv',reconstruct)
# summaries of diagnostic-marker allele proportions (only two available parental controls)
marker_records=[r for r in each_site(amended) if key(r) in pure_marker_keys]
assert len(marker_records)==1278
mixture=[];med_errors=[]
for i,s in enumerate(samples):
 name=names[i]
 m=meta[name]; c=classif[name]; cf=flag[name];st=stat[name]
 fraction=re.search(r'_(\d+)pct([AB])?_',name)
 nom=int(fraction.group(1)) if fraction else None
 is_mix=name.startswith('MIX_'); is_same=name.startswith('CTRL_SAMELIN_');is_pure=name.startswith('PURE_')
 applies=('ERR038255_ERR13260175' in name) or is_pure
 # Parent-diagnostic metrics only when both pure controls represent the specific mixture's contributors.
 passing=[];all_nonmissing=[]
 if applies:
  for r in marker_records:
   cc=call(r,i)
   if 0<=cc['freq']<=100:all_nonmissing.append(cc['freq'])
   if cc['depth']>=20 and cc['forward']>=4 and cc['reverse']>=4 and cc['qual20']>=4 and 0<=cc['freq']<=100:
    passing.append(cc['freq'])
 measured=round(statistics.median([100-v for v in passing]),2) if passing else None
 mid=round(sum(25 <= v <=75 for v in passing)/len(passing),4) if passing else None
 full_lineage=m['lineage_group']
 row={'sample':name,'study_truth':('mixed_same_lineage' if is_same else 'mixed_different_lineages' if is_mix else 'pure'),
      'known_mixture_nominal_percent':nom,'orientation_from_A_B_label':'NOT_CONFIRMED',
      'reported_rmap_tbprofiler_mixed_flag':cf['tbprofiler_mixed_infection_status'],
      'rmap_lineage_group':full_lineage,'mtbseq_single_lineage_coll_easy':c['Coll lineage (easy)'].strip("' "),
      'mtbseq_mixed_infection_classifier':'NOT AVAILABLE in supplied MTBseq classification output',
      'mtbseq_mapped_reads_pct':st['% Mapped Reads'].strip("' "),
      'rmap_mean_depth':m['mean_depth'],
      'diagnostic_parent_panel_applicable':'YES' if applies else 'NO — parents not in pure-control panel',
      'parent_diagnostic_positions_all':1278 if applies else None,
      'parent_diagnostic_positions_depthQC20':len(passing) if applies else None,
      'median_minor_allele_proxy_pct':measured,'fraction_markers_with_major_freq_25to75':mid,
      'site_proxy_note':'100 minus MTBseq-reported major allele frequency; not direct minority variant count' if applies else 'Not calculated'}
 mixture.append(row)
 if is_mix and applies:
  med_errors.append((name,nom,measured,abs(measured-nom)))
write_tsv('mixture_detection_and_marker_frequencies.tsv',mixture)
# Additional summary
same_pure=[r for r in reconstruct if {r['sample_a'],r['sample_b']}==set(('PURE_ERR038255_100x','PURE_ERR13260175_100x'))][0]
assert same_pure['rmap_full_core_snp_distance']==1148
true_mixes=[r for r in mixture if r['study_truth']=='mixed_different_lineages']
assert len(true_mixes)==15
flagged=sum('Possible mixed' in r['reported_rmap_tbprofiler_mixed_flag'] for r in true_mixes)
assert flagged==15
pure_neg=sum('No mixed-infection signal' in r['reported_rmap_tbprofiler_mixed_flag'] for r in mixture if r['study_truth']=='pure')
assert pure_neg==2
same_lineage=[r for r in mixture if r['study_truth']=='mixed_same_lineage']
assert len(same_lineage)==1
strict_variation=sum(len(set(s[i] for s in seq.values()))>1 for i in range(237))
assert strict_variation==4
numeric_summary={
 'n_samples':18,'n_different_lineage_mixtures_labeled_MIX':15,'n_same_lineage_mixture_controls':1,'n_pure_strain_controls':2,
 'candidate_joint_positions':3148,'all18_unambiguous_positions':291,'phylo_positions':246,'final_window12_positions':237,
 'final_informative_positions':4,'original_mtbseq_distance_values':dict(collections.Counter(mtb_dist.values())),
 'mtbseq_pairwise_matrix_reconstructed_pairs':len(reconstruct),'pure_pair_candidate_positions_callable_both':2201,
 'pure_pair_candidate_different_positions':1278,'pure_pair_candidate_differences_surviving_u95':0,
 'rmap_pure_controls_full_core_snp_distance':1148,
 'rmap_pure_controls_comparable_sites':same_pure['rmap_comparable_sites'],
 'rmap_tbprofiler_labeled_mix_flagged':flagged,'rmap_tbprofiler_pure_controls_without_flag':pure_neg,
 'different_lineage_mixtures_with_multilineage_lineage_group':sum(r['rmap_lineage_group'].startswith('Mixed') for r in true_mixes),
 'source_parents_12_mixture_med_minor_allele_proxy_mae':round(statistics.mean(v[3] for v in med_errors),3),
}
(ANAL/'checked_summary.json').write_text(json.dumps(numeric_summary,indent=2)+'\n')
manifest={'project':'rMAP-TB versus MTBseq 18-library mixed-strain benchmarking audit','generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'baseline_workflow':'MTBSEQ_BENCHMARK','baseline_task':'MTBSEQ_FULL','mtbseq_version':'1.1.0',
 'original_submission':'d11eb528-1d56-496c-901b-3213d38af7de',
 'original_workflow':'e76bf8e8-bc07-46c5-a8d6-0a390c1cc569',
 'cloud_archive':'gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/mtbseq_out/mtbseq_outputs.tar.gz',
 'full_21_6GB_archive_included':False,'raw_files':rawfiles,'parameters':{'step':'TBfull','mincovf':4,'mincovr':4,'minphred20':4,'minfreq':75,'unambig':95,'window':12,'distance':12,'threads':8,
 'reference':'M._tuberculosis_H37Rv_2015-11-13'},'rmap_18_workflow':'06425efe-8104-42ef-b8f4-7f4b69f4e3b5','rmap_18_submission':'3d10dc72-6796-4f83-85ac-6ebbbf7a816'}
(RAW/'BASELINE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('OUTPUT',OUT)
print(json.dumps(numeric_summary,indent=2))
print('mix summary:',*(f"{x[0]} nominal={x[1]} measured={x[2]} err={x[3]}" for x in med_errors),sep='\n')
