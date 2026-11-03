source '/mnt/disks/cromwell_root/gcs_transfer.sh'

timestamped_message 'Localization script execution started...'


# No reference disks mounted since not requested in workflow options.





# Localize files from source bucket 'fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1' to container parent directory '/mnt/disks/cromwell_root'.
files_to_localize_113423a69d1a0f36416ed1ad8ddda8b1=(
  "terra-dd7d30c6"   # project to use if requester pays
  "3" # max transfer attempts
  "/mnt/disks/cromwell_root/" # container parent directory
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-USHER_PLACE/script"
)

localize_files "${files_to_localize_113423a69d1a0f36416ed1ad8ddda8b1[@]}"
       


# Localize files from source bucket 'fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1' to container parent directory '/mnt/disks/cromwell_root/fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/usher_prep'.
files_to_localize_303565ade85b1aa438d162a6b03638d0=(
  "terra-dd7d30c6"   # project to use if requester pays
  "3" # max transfer attempts
  "/mnt/disks/cromwell_root/fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/usher_prep/" # container parent directory
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-ALIGNMENT_TO_VCF/usher_prep/alignment.vcf"
)

localize_files "${files_to_localize_303565ade85b1aa438d162a6b03638d0[@]}"
       


# Localize files from source bucket 'fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1' to container parent directory '/mnt/disks/cromwell_root/fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-BUILD_STARTING_TREE/starting_tree'.
files_to_localize_13a7d22cb18395b5436e59be5e76ff93=(
  "terra-dd7d30c6"   # project to use if requester pays
  "3" # max transfer attempts
  "/mnt/disks/cromwell_root/fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-BUILD_STARTING_TREE/starting_tree/" # container parent directory
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/c798405b-0542-4aab-afdb-4a3126a7e223/USHER_TB_BENCHMARK/0a3b74e0-43fd-429a-851b-6a4afc6828af/call-BUILD_STARTING_TREE/starting_tree/starting_tree.nwk"
)

localize_files "${files_to_localize_13a7d22cb18395b5436e59be5e76ff93[@]}"
       

timestamped_message 'Localization script execution complete.'