source '/mnt/disks/cromwell_root/gcs_transfer.sh'

timestamped_message 'Localization script execution started...'


# No reference disks mounted since not requested in workflow options.





# Localize files from source bucket 'fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1' to container parent directory '/mnt/disks/cromwell_root/fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/3d10dc72-6796-4f83-85ac-6bebbbf7a816/rMAP_TB/06425efe-8104-42ef-b8f4-7f4b69f4e3b5/call-SNIPPY_CORE_MTBC/cacheCopy/snippy_core'.
files_to_localize_88dcc499028a45dc173fd107fe81be45=(
  "terra-dd7d30c6"   # project to use if requester pays
  "3" # max transfer attempts
  "/mnt/disks/cromwell_root/fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/3d10dc72-6796-4f83-85ac-6bebbbf7a816/rMAP_TB/06425efe-8104-42ef-b8f4-7f4b69f4e3b5/call-SNIPPY_CORE_MTBC/cacheCopy/snippy_core/" # container parent directory
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/3d10dc72-6796-4f83-85ac-6bebbbf7a816/rMAP_TB/06425efe-8104-42ef-b8f4-7f4b69f4e3b5/call-SNIPPY_CORE_MTBC/cacheCopy/snippy_core/core.vcf"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/3d10dc72-6796-4f83-85ac-6bebbbf7a816/rMAP_TB/06425efe-8104-42ef-b8f4-7f4b69f4e3b5/call-SNIPPY_CORE_MTBC/cacheCopy/snippy_core/core.full.aln"
)

localize_files "${files_to_localize_88dcc499028a45dc173fd107fe81be45[@]}"
       


# Localize files from source bucket 'fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1' to container parent directory '/mnt/disks/cromwell_root'.
files_to_localize_113423a69d1a0f36416ed1ad8ddda8b1=(
  "terra-dd7d30c6"   # project to use if requester pays
  "3" # max transfer attempts
  "/mnt/disks/cromwell_root/" # container parent directory
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/script"
)

localize_files "${files_to_localize_113423a69d1a0f36416ed1ad8ddda8b1[@]}"
       

timestamped_message 'Localization script execution complete.'