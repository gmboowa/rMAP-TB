#!/usr/bin/env python3
"""Read-only diagnostic for MTBseq vs rMAP-TB reference-coordinate mismatches.

Does not suppress a mismatch, edit source files, or calculate cross-pipeline
SNP distance. Outputs a full mismatch table and counts by reference base.
"""
import argparse
import csv
import io
import json
import re
import zipfile
from collections import Counter
from pathlib import Path

CANON = set('ACGT')

def name_id(value):
    return re.sub(r'^\d{4}_', '', value.strip().strip("'").removesuffix('_lib1').replace('-', '_'))

def read_table(archive, suffix):
    with zipfile.ZipFile(archive) as handle:
        candidates = [x for x in handle.namelist() if x.endswith(suffix) and not x.startswith('__MACOSX/')]
        if len(candidates) != 1:
            raise SystemExit(f'Expected exactly one {suffix} in {archive}; found: {candidates}')
        return list(csv.reader(io.StringIO(handle.read(candidates[0]).decode('utf-8-sig')), delimiter='\t'))

def fasta_sequences(path):
    sequences = {}
    current = None
    chunks = []
    with path.open('r', encoding='utf-8') as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current is not None:
                    if current in sequences:
                        raise SystemExit(f'Duplicate FASTA ID: {current}')
                    sequences[current] = ''.join(chunks).upper()
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        if current in sequences:
            raise SystemExit(f'Duplicate FASTA ID: {current}')
        sequences[current] = ''.join(chunks).upper()
    return sequences

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rmap-fasta', type=Path, required=True)
    parser.add_argument('--audit-zip', type=Path, required=True)
    parser.add_argument('--amend-zip', type=Path, default=None,
                        help='Optional ZIP containing the final _amended_u95_phylo_w12.tab')
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    seqs_raw = fasta_sequences(args.rmap_fasta.expanduser())
    seqs = {name_id(key): value for key, value in seqs_raw.items()}
    rows = read_table(args.audit_zip.expanduser(), '_amended.tab')
    sample_ids = [name_id(rows[0][18 + 8 * idx]) for idx in range(18)]
    missing = sorted(set(sample_ids) - set(seqs))
    extra = sorted(set(seqs) - set(sample_ids))
    if missing or len(extra) != 1:
        raise SystemExit(f'Sample/reference ID mismatch: missing={missing}; extra={extra}')
    ref_label = extra[0]
    ref = seqs[ref_label]
    lengths = sorted({len(v) for v in seqs.values()})
    if len(lengths) != 1:
        raise SystemExit(f'Not all sequences have the same length: {lengths}')
    filtered_keys = None
    if args.amend_zip:
        final = read_table(args.amend_zip.expanduser(), '_amended_u95_phylo_w12.tab')
        filtered_keys = {(r[0].strip(), r[1].strip()) for r in final[2:]}

    categories = Counter()
    final_categories = Counter()
    ref_char_counts = Counter()
    offset_matches = Counter()
    disagreements = []
    checked = 0
    skipped_nonzero_ins = 0
    skipped_noncanonical_ref = 0
    skipped_outside = 0
    for row in rows[2:]:
        if len(row) < 3:
            continue
        position, insert_index, mtb_ref = row[0].strip(), row[1].strip(), row[2].strip().upper()
        if insert_index != '0':
            skipped_nonzero_ins += 1
            continue
        if mtb_ref not in CANON:
            skipped_noncanonical_ref += 1
            continue
        try:
            pos = int(position)
        except ValueError:
            continue
        if not (1 <= pos <= len(ref)):
            skipped_outside += 1
            continue
        base = ref[pos - 1]
        checked += 1
        ref_char_counts[base] += 1
        if base == mtb_ref:
            classification = 'MATCH'
        elif base not in CANON:
            classification = 'MASKED_OR_GAP'
        else:
            classification = 'DIFFERENT_ACGT_REFERENCE_BASE'
        categories[classification] += 1
        is_final = filtered_keys is not None and (position, insert_index) in filtered_keys
        if is_final:
            final_categories[classification] += 1
        if classification != 'MATCH':
            disagreements.append({
                'position_1_based': pos,
                'mtbseq_reference': mtb_ref,
                'rmap_reference_character': base,
                'category': classification,
                'in_final_MTBseq_237_site_table': 'YES' if is_final else 'NO',
                'rmap_reference_context': ref[max(0, pos - 6):min(len(ref), pos + 5)],
            })
        for shift in (-2, -1, 0, 1, 2):
            k = pos - 1 + shift
            if 0 <= k < len(ref) and ref[k] == mtb_ref:
                offset_matches[str(shift)] += 1

    out = args.outdir.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'reference_discrepancies.tsv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=['position_1_based', 'mtbseq_reference',
            'rmap_reference_character', 'category', 'in_final_MTBseq_237_site_table',
            'rmap_reference_context'], delimiter='\t')
        writer.writeheader()
        writer.writerows(disagreements)
    summary = {
        'rmap_fasta': str(args.rmap_fasta.expanduser()),
        'mtbseq_joint_amended_zip': str(args.audit_zip.expanduser()),
        'rmap_reference_name': ref_label,
        'alignment_length': len(ref),
        'n_fasta_sequences': len(seqs),
        'n_matching_sample_ids': len(sample_ids),
        'checked_noninsertion_canonical_reference_positions': checked,
        'reference_comparison': dict(categories),
        'rmap_reference_character_counts_at_checked_positions': dict(ref_char_counts),
        'final_filtered_table_reference_comparison': dict(final_categories),
        'positions_skipped_with_nonzero_insertion_index': skipped_nonzero_ins,
        'positions_skipped_with_noncanonical_mtbseq_ref': skipped_noncanonical_ref,
        'positions_skipped_out_of_alignment_bounds': skipped_outside,
        'reference_matches_with_nearby_offsets_DONT_USE_TO_CORRECT': dict(offset_matches),
        'conclusion': 'Diagnostic only. Do not compute matched-site distances until coordinate assumptions are confirmed.'
    }
    (out / 'reference_diagnostic_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))
    print('\nWrote:', out / 'reference_discrepancies.tsv')
    print('Wrote:', out / 'reference_diagnostic_summary.json')

if __name__ == '__main__':
    main()
