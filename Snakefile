configfile: "config.yaml"
sample = config["samples"]
method = config["methods"]

import os

def activate_environment():
    os.system("module load gcc/7.2.0")
    os.system("module load anaconda3/2019.07")
    os.system("source activate /home/mayerma/.conda/envs/deepsignalenv")

activate_environment()

rule all:
  input:
    "output/deepsignal_methrix.bedGraph",
    "output/tombo_methrix.bedGraph",
    "output/nanopolish_methrix.bedGraph",
    "output/bisulfite_methrix.bedGraph"


rule multi_to_single_fast5:
    input:
        path = "raw/{sample}/guppy"
    output:
        single = directory("output/{sample}/fast5s_guppy.single.al")
    params:
        jobname='multi_to_single',
        runtime='160:00',
        memusage='10000',
        slots='1',
        misc=''
    shell:
        "multi_to_single_fast5 --input_path {input.path} --save_path {output.single} --recursive"

rule tombo_resquiggle:
    input:
        fast5s = rules.multi_to_single_fast5.output.single,
        reference = "GCF_000001405.26_GRCh38_genomic.fna"
    output:
        dummy = "output/{sample}/resquiggled_done.txt"
    params:
        jobname='tombo_resquiggle',
        runtime='60:00',
        memusage='40000',
        slots='10',
        misc=''
    run:
        shell("tombo resquiggle {input.fast5s} {input.reference} --processes 10 --overwrite")
        shell("touch {output.dummy}")

rule tombo_modifications:
    input:
        fast5path = rules.multi_to_single_fast5.output.single,
        dummy = rules.tombo_resquiggle.output.dummy
    output:
        modifications = "output/{sample}/alt_model.CpG.tombo.stats",
        dummy = "output/{sample}/tombo_mods_done.txt"
    params:
        jobname='tombo_modifications',
        runtime='150:00',
        memusage='90000',
        slots='8',
        misc=''
    run:
        shell("tombo detect_modifications alternative_model --fast5-basedirs {input.fast5path} --alternate-bases CpG --statistics-file-basename output/{wildcards.sample}/alt_model")
        shell("touch output/{wildcards.sample}/tombo_mods_done.txt")

rule tombo_browserfiles:
    input:
        fast5path = rules.multi_to_single_fast5.output.single,
        modifications = rules.tombo_modifications.output.modifications,
        dummy = rules.tombo_modifications.output.dummy
    output:
        "output/{sample}/tombo_{sample}.fraction_modified_reads.minus.wig",
        "output/{sample}/tombo_{sample}.fraction_modified_reads.plus.wig",
        "output/{sample}/tombo_{sample}.valid_coverage.minus.wig",
        "output/{sample}/tombo_{sample}.valid_coverage.plus.wig"
    params:
        jobname='tombo_browserfiles',
        runtime='20:00',
        memusage='16000',
        slots='4',
        misc=''
    run:
        shell("tombo text_output browser_files --fast5-basedirs {input.fast5path} --statistics-filename {input.modifications} --file-types fraction valid_coverage --browser-file-basename output/{wildcards.sample}/tombo_{wildcards.sample} ")

rule gather_all_tombo:
    input:
        expand("output/{SAMPLE}/tombo_{SAMPLE}.fraction_modified_reads.minus.wig",SAMPLE = sample),
        expand("output/{SAMPLE}/tombo_{SAMPLE}.fraction_modified_reads.plus.wig",SAMPLE = sample)
    output:
        "output/all_tombo_frequency.tsv"
    params:
        jobname='gather_tombo',
        runtime='10:00',
        memusage='16000',
        slots='4',
        misc=''
    shell:
        "cat {input} >> {output}"

rule deepsignal:
    input:
        fast5s = rules.tombo_resquiggle.input.fast5s,
        resquiggled = rules.tombo_resquiggle.output.dummy,
        reference = "GCF_000001405.26_GRCh38_genomic.fna"
    output:
        calls = "output/{sample}/fast5s.al.CpG.signal_features.17bases.rawsignals_360.tsv",
        call_mods = "output/{sample}/fast5s.al.CpG.call_mods.tsv"
    params:
        jobname='deepsignal',
        runtime='180:00',
        memusage='80000',
        slots='10',
        misc=''
    run:
        shell("deepsignal extract --fast5_dir {input.fast5s} --reference_path {input.reference} --write_path {output.calls} --nproc 10")
        shell("deepsignal call_mods --input_path {output.calls} --model_path model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt --result_file {output.call_mods} --nproc 10 --is_gpu no")

rule deepsignal_call_mods:
    input:
        script = "call_modification_frequency.py",
        call_mods = rules.deepsignal.output.call_mods
    output:
        mods_frequency = "output/{sample}/deepsignal_frequency.tsv"
    params:
        jobname='deepsignal_call_mods',
        runtime='50:00',
        memusage='16000',
        slots='8',
        misc=''
    shell:
        "python {input.script} --input_path {input.call_mods} --result_file {output.mods_frequency} --prob_cf 0"

rule gather_all_deepsignal:
    input:
        expand("output/{SAMPLE}/deepsignal_frequency.tsv", SAMPLE = sample)
    output:
        "output/all_deepsignal_frequency.tsv"
    params:
        jobname='gather_deepsignal',
        runtime='10:00',
        memusage='8000',
        slots='4',
        misc=''
    shell:
        "cat {input} >> {output}"


rule frequency_to_methrix:
    input:
        "output/all_{method}_frequency.tsv"
    output:
        "output/{method}_methrix.bedGraph"
    params:
        jobname='frequency_to_methrix',
        runtime='30:00',
        memusage='20000',
        slots='4',
        misc=''
    run:
        shell("python3 calc_freq_for_methrix.py --input_file {input} --result_file {output}")
