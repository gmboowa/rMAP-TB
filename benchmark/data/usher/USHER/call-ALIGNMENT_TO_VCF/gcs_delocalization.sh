source '/mnt/disks/cromwell_root/gcs_transfer.sh'

timestamped_message 'Delocalization script execution started...'

# fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1
delocalize_092b8e1d719cdfa14815e505d7e73347=(
  "terra-dd7d30c6"       # project
  "3"   # max attempts
  "0" # parallel composite upload threshold, will not be used for directory types
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/memory_retry_rc"
  "/mnt/disks/cromwell_root/memory_retry_rc"
  "optional"
  "text/plain; charset=UTF-8"
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/rc"
  "/mnt/disks/cromwell_root/rc"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/stdout"
  "/mnt/disks/cromwell_root/stdout"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/stderr"
  "/mnt/disks/cromwell_root/stderr"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/usher_prep/alignment.vcf"
  "/mnt/disks/cromwell_root/usher_prep/alignment.vcf"
  "required"
  ""
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/usher_prep/conversion.log"
  "/mnt/disks/cromwell_root/usher_prep/conversion.log"
  "required"
  ""
  "file"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/usher_prep/alignment.ordered.fa"
  "/mnt/disks/cromwell_root/usher_prep/alignment.ordered.fa"
  "required"
  ""
)

delocalize "${delocalize_092b8e1d719cdfa14815e505d7e73347[@]}"
      
timestamped_message 'Delocalization script execution complete.'