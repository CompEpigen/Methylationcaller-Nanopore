#!/bin/bash
cluster_script="bsub {params.misc} -n {params.slots} -W {params.runtime} -M {params.memusage}  -R \"rusage[mem={params.memusage}]\" -o .snakemake/logs/{params.jobname}.log -e .snakemake/logs/{params.jobname}.err"
snakemake -k --cluster "$cluster_script" --jobs 64 --latency-wait 120 $@
