#!/bin/bash

CONDA=/tool_deps/_conda/bin/conda

${CONDA} create -q -y --name __siteseq@0.1.0

source /tool_deps/_conda/bin/activate __siteseq@0.1.0

${CONDA} install -q -y --override-channels --channel bioconda \
--channel conda-forge --channel defaults \
pysam=0.15.3=py27hda2845c_1

pip install pyfaidx

${CONDA} clean --all --yes
