#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna prep elastic -m $DIR/sra_batch_33.manifest --skip-bad-records --core-instance-type m3.xlarge --master-instance-type m3.xlarge -o s3://rail-sra-hg38/sra_prep_batch_33 -c 20 --core-instance-bid-price 0.25 --master-instance-bid-price 0.25 -f --max-task-attempts 6
