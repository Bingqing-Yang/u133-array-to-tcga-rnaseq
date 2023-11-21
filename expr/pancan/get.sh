#!/bin/bash
# Retrieve raw data
# modify Dr. David Shih get.sh script
# credit to Dr. David Shih

set -euo pipefail

mkdir gdc
# "gdc-client" can be replaced by the location of the gdc-client script
gdc-client download -m manifest/PanCan-General_Open_GDC-Manifest_2.txt -d gdc
