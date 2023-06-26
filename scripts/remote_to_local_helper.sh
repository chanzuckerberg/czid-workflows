set -euo pipefail

## A helper function that replaces paths run remotely with an s3 bucket path
## Usage ./scripts/remote_to_local_helper.sh inputs-copied-from-remote-run.json s3://bucket-where-outputs-really-are

INPUT_JSON=$1 
S3_BUCKET_PATH=$(echo $2 | sed 's/\/$//' | sed 's/\//\\\//g')
REPLACE_MATCH=$(echo /mnt | sed 's/\//\\\//g')


jq '.' $INPUT_JSON | sed -E "s/(.*)${REPLACE_MATCH}.*(\/[A-Za-z0-9]+)(\/?)/\1${S3_BUCKET_PATH}\2/g" \
	| sed -E "s/(.*\"docker_image_id\": \")(.*\/)([A-Za-z0-9_-]+)(:.*)(\",)/\1czid-\3\5/g"
