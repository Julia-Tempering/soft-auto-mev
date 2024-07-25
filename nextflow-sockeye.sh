#!/bin/bash

# explicit number of nodes needed since Nov-23
cl_ops="--nodes=1 --account"

# decide alloc code based on user names
if [[ "$USER" == @("tgliu"|"mbironla") ]]; then
    cl_ops="$cl_ops=st-tdjc-1"
else
    cl_ops="$cl_ops=st-alexbou-1"
fi

export CLUSTER_OPTIONS="$cl_ops"
echo "Using CLUSTER_OPTIONS=$CLUSTER_OPTIONS"
./nextflow $@ -profile cluster

