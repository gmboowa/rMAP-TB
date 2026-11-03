#!/usr/bin/env python3
"""Read-only validation of shipped analysis tables. No external dependencies."""
import csv, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
A=ROOT/'analysis'
s=json.loads((A/'checked_summary.json').read_text())
assert s['n_samples']==18
assert [s['candidate_joint_positions'],s['all18_unambiguous_positions'],s['phylo_positions'],s['final_window12_positions']]==[3148,291,246,237]
assert s['mtbseq_pairwise_matrix_reconstructed_pairs']==153
assert s['pure_pair_candidate_positions_callable_both']==2201
assert s['pure_pair_candidate_different_positions']==1278
assert s['rmap_pure_controls_full_core_snp_distance']==1148
for name,n in [('all153_pairs_DESCRIPTIVE_UNMATCHED_MASKS.tsv',153),('pure_control_pair_only_candidate_sites.tsv',2201),('mixture_detection_and_marker_frequencies.tsv',18),('mtbseq_site_mask_audit.tsv',3148)]:
 with open(A/name,newline='') as f:
  rows=list(csv.DictReader(f,delimiter='\t'))
 assert len(rows)==n,(name,len(rows),n)
 print(name,n,'rows: PASS')
print('Original 18-library audit: all tests passed.')
