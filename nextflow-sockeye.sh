#!/usr/bin/env bash

if [[ -z "${ALLOC_NAME}" ]]; then
  echo "The variable ALLOC_NAME -- which should contain your allocation name -- is not set. See the README.md file."
  exit 1
fi

# explicit number of nodes needed since Nov-23
export CLUSTER_OPTIONS="--nodes=1 --account=$ALLOC_NAME"
echo "Using CLUSTER_OPTIONS=$CLUSTER_OPTIONS"
./nextflow $@ -profile cluster
