source '/mnt/disks/cromwell_root/gcs_transfer.sh'

timestamped_message 'Localization script execution started...'


# No reference disks mounted since not requested in workflow options.





# Localize files from source bucket 'fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1' to container parent directory '/mnt/disks/cromwell_root/fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures'.
files_to_localize_9a27fa94ed342bfa0185d9066f471e3e=(
  "terra-dd7d30c6"   # project to use if requester pays
  "3" # max transfer attempts
  "/mnt/disks/cromwell_root/fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/" # container parent directory
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_30pctA_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/PURE_ERR13260175_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038266_ERR13260106_50pctB_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_50pctB_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/PURE_ERR13260175_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_05pctB_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctB_140x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_20pctB_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctA_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038266_ERR13260106_10pctB_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/CTRL_SAMELIN_ERR038254_ERR038255_50pct_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_20pctB_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_30pctA_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/PURE_ERR038255_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038266_ERR13260106_20pctB_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctB_70x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_20pctA_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_05pctB_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_05pctA_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctB_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_20pctA_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_30pctB_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctB_40x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_05pctA_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_50pctB_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/CTRL_SAMELIN_ERR038254_ERR038255_50pct_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_30pctB_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/PURE_ERR038255_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctB_140x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctB_40x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038266_ERR13260106_20pctB_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctB_70x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctB_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038255_ERR13260175_10pctA_100x_2.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038266_ERR13260106_10pctB_100x_1.fastq.gz"
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/uploads/TB_mixtures/MIX_ERR038266_ERR13260106_50pctB_100x_1.fastq.gz"
)

localize_files "${files_to_localize_9a27fa94ed342bfa0185d9066f471e3e[@]}"
       


# Localize files from source bucket 'fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1' to container parent directory '/mnt/disks/cromwell_root'.
files_to_localize_113423a69d1a0f36416ed1ad8ddda8b1=(
  "terra-dd7d30c6"   # project to use if requester pays
  "3" # max transfer attempts
  "/mnt/disks/cromwell_root/" # container parent directory
  "gs://fc-67b90df7-2f0b-4506-947f-aaf2f23ceaa1/submissions/d11eb528-1d56-496c-901b-3213d38af7de/MTBSEQ_BENCHMARK/e76bf8e8-bc07-46c5-a8d6-0a390c1cc569/call-MTBSEQ_FULL/script"
)

localize_files "${files_to_localize_113423a69d1a0f36416ed1ad8ddda8b1[@]}"
       

timestamped_message 'Localization script execution complete.'