#!/usr/bin/env bash

if [[ -z "${ALLOC_NAME}" ]]; then
  echo "The variable ALLOC_NAME -- which should contain your allocation name -- is not set. See the README.md file."
  exit 1
fi

export CLUSTER_OPTIONS="--account=$ALLOC_NAME"
NXF_OPTS="-Xms500M -Xmx3G" ./nextflow $@ -profile cluster
