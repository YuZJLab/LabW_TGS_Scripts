#!/usr/bin/env bash
# VERSION 1.0

set -e # Prevent blank variables and stop when error raised
REFERENCE_FASTA="/home/pub/References/GATK_hg38/hg38.fa"
THREAD=40
__conda_setup="$("${HOME}/conda/bin/conda" 'shell.bash' 'hook' 2> /dev/null)"
eval "${__conda_setup}"


######################## Simple error handlers ########################
# From YuZJLab LinuxMiniPrograms, <https://gitee.com/YuZJLab/LinuxMiniPrograms>
getuser() {
    if [ -n "${USER:-}" ]; then
        printf "${USER}"
        return
        elif [ -n "${USERNAME:-}" ]; then
        printf "${USERNAME}"
        return
    else
        return 1
    fi
}
getcorenumber() {
    if [[ "${myos}" == *"BSD"* ]] || [[ "${myos}" == *"Darwin"* ]]; then
        sysctl -a | grep hw.ncpu | awk '{print $2;}'
    else
        cat /proc/cpuinfo | grep '^processor\s: ' | wc -l | awk '{print $1}'
    fi
}
killtree() {
    kill -19 ${1} || return
    for CPID in $(ps -o pid --no-headers --ppid ${1}); do
        killtree ${CPID} ${2}
    done
    kill -${2} ${1} || return
}
trimstr() {
    : "${1#"${1%%[![:space:]]*}"}"
    : "${_%"${_##*[![:space:]]}"}"
    printf '%s\n' "$_"
}
errh() {
    echo -e "\033[31mERROR: ${*}\033[0m" >&2
    exit 1
}
warnh() {
    echo -e "\033[31mWARNING: ${*}\033[0m" >&2
}
infoh() {
    echo -e "\033[33mINFO: ${*}\033[0m" >&2
}
conda_activate(){
    conda env list | grep "${1}" 2> /dev/null # Determine whether this environment exists
    source activate "${1}"
}


######################## Scripts starting from here do the job ########################

preqc_nanoplot(){
    conda_activate nanoplot
    NanoPlot --fastq "${1}.fastq.gz" \
    -o "${1}.nanoplot" \
    -t 40 \
    --verbose \
    --plots kde hex dot pauvre
}

# This is not used because it rediculously dumped out SAM format.
align_flair(){
    conda_activate flair
    flair.py align -g "${REFERENCE_FASTA}" \
    -r "${1}.fastq.gz" \
    -t ${THREAD} \
    -o "${1}.aligned"
    conda deactivate
}

# This program uses the following steps to reduce disk space:
# firstly it asks miminap to dump SAM format
# Then uses samtools to convert it to compressed BAM
align_minimap(){
    minimap2 -ax splice \
    -t ${THREAD} \
    --secondary=no \
    "${REFERENCE_FASTA}" "${1}.fastq.gz" | \
    samtools view -h \
    -@ ${THREAD} > "${1}.aligned.unsorted.bam"
}


qc_alignqc(){
    conda_activate alignqc
    alignqc analyze "${1}.aligned.unsorted.bam" -g "${REFERENCE_FASTA}" \
    --no_transcriptome \
    --threads ${THREAD} \
    -o "${1}.alignqc.xhtml"
    conda deactivate
}

######################## Scripts starting from here build the skeleton ########################

align(){
    align_minimap "${@}"
}

qc(){
    qc_alignqc "${@}"
}

main(){
    [ -f data.conf ]
    cat data.conf | while read line;do
        # align "${line}"
        qc "${line}"
    done
}

main
