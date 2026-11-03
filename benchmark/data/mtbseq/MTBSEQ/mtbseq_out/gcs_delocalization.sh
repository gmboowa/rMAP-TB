source '/mnt/disks/cromwell_root/gcs_transfer.sh'

timestamped_message 'Delocalization script execution started...'

# fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1
delocalize_092b8e1d719cdfa14815e505d7e73347=(
  "terra-dd7d30c6"       # project
  "3"   # max attempts
  "0" # parallel composite upload threshold, will not be used for directory types
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/memory_retry_rc"
  "/mnt/disks/cromwell_root/memory_retry_rc"
  "optional"
  "text/plain; charset=UTF-8"
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/rc"
  "/mnt/disks/cromwell_root/rc"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/stdout"
  "/mnt/disks/cromwell_root/stdout"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/stderr"
  "/mnt/disks/cromwell_root/stderr"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/mtbseq_out/samples.tsv"
  "/mnt/disks/cromwell_root/mtbseq_out/samples.tsv"
  "required"
  ""
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/mtbseq_out/mtbseq_outputs.tar.gz"
  "/mnt/disks/cromwell_root/mtbseq_out/mtbseq_outputs.tar.gz"
  "required"
  ""
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/mtbseq_out/mtbseq_run.log"
  "/mnt/disks/cromwell_root/mtbseq_out/mtbseq_run.log"
  "required"
  ""
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/mtbseq_out/sample_map.tsv"
  "/mnt/disks/cromwell_root/mtbseq_out/sample_map.tsv"
  "required"
  ""
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/mtbseq_out/status.txt"
  "/mnt/disks/cromwell_root/mtbseq_out/status.txt"
  "required"
  ""
)

delocalize "${delocalize_092b8e1d719cdfa14815e505d7e73347[@]}"
      
timestamped_message 'Delocalization script execution complete.'