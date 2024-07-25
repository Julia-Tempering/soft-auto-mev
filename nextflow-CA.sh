#!/usr/bin/env bash

export CLUSTER_OPTIONS="--account=def-bouchar3"
NXF_OPTS="-Xms500M -Xmx3G" ./nextflow $@ -profile cluster
