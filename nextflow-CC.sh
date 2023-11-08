#!/bin/bash

NXF_OPTS="-Xms500M -Xmx3G" ./nextflow $@ -profile CC
