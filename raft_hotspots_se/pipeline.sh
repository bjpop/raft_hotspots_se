#!/bin/bash

program_name="pipeline.sh"
reference=""
fastq=""
prefix=""
blacklist=""
sites=""

# Help message for using the program.
function show_help {
cat << UsageMessage

${program_name}: execute the raft DSB site detection pipeline 

Usage:
    ${program_name} [-h] [-v] -r reference -f fastq -p prefix -b blacklist -s sites

-h shows this help message

-v verbose output

UsageMessage
}

# Parse the command line arguments and set the global variables language and new_project_name
function parse_args {
    local OPTIND opt

    while getopts "hvr:f:p:b:s:" opt; do
        case "${opt}" in
            h)
                show_help
                exit 0
                ;;
	    v)  verbose=true
		;;
            r)  reference="${OPTARG}"
                ;;
            f)  fastq="${OPTARG}"
                ;;
            p)  prefix="${OPTARG}"
                ;;
            b)  blacklist="${OPTARG}"
                ;;
            s)  sites="${OPTARG}"
                ;;
        esac
    done

    shift $((OPTIND-1))

    [ "$1" = "--" ] && shift

    if [[ -z ${reference} ]]; then
        echo "${program_name}: ERROR: missing command line argument: -r reference, use -h for help"
        exit 2
    fi

    if [[ -z ${fastq} ]]; then
        echo "${program_name}: ERROR: missing command line argument: -f fastq, use -h for help"
        exit 2
    fi

    if [[ -z ${prefix} ]]; then
        echo "${program_name}: ERROR: missing command line argument: -p prefix, use -h for help"
        exit 2
    fi

    if [[ -z ${blacklist} ]]; then
        echo "${program_name}: ERROR: missing command line argument: -b blacklist, use -h for help"
        exit 2
    fi

    if [[ -z ${sites} ]]; then
        echo "${program_name}: ERROR: missing command line argument: -s sites, use -h for help"
        exit 2
    fi
}

function verbose_message {
    if [ "${verbose}" = true ]; then
        echo "${program_name} $1"
    fi
}

function run_commands {

    step="Pre-processing fastq file, present DSBs at start"
    verbose_message "$step"
    raft_fastq_2sites_parse --fastq ${fastq} > "${prefix}.preprocess.fastq" || {
        echo "FAILED: $step"
	exit 1
    }

    step="Aligning reads to the reference"
    verbose_message "$step"
    bwa mem ${reference} "${prefix}.preprocess.fastq" > ${prefix}.sam || {
        echo $step 
	exit 1
    }

    step="Converting SAM file to sorted BAM"
    verbose_message "$step"
    samtools view -u  ${prefix}.sam | samtools sort -@ 4 -o ${prefix}.sort.bam || {
        echo "FAILED: $step"
	exit 1
    }

    step="Converting BAM file to BED"
    verbose_message "$step"
    bedtools bamtobed -cigar -i ${prefix}.sort.bam > ${prefix}.bed || {
        echo "FAILED: $step"
	exit 1
    }

    step="Remove blacklisted site-overlapping reads; yield bed file"
    verbose_message "$step"
    bedtools subtract -A -a ${prefix}.bed -b ${blacklist} > ${prefix}.blf.bed || {
        echo "FAILED: $step"
	exit 1
    }

    step="Counting DSB points; filter for no clipping and high MQ"
    verbose_message "$step"
    raft_bed_2sites_parse --bed ${prefix}.blf.bed > ${prefix}.counts.bed || {
        echo "FAILED: $step"
	exit 1
    }

    step="Filter out DSBs proximal to Sau3A1 sites"
    verbose_message "$step"
    bedtools subtract -A -a ${prefix}.counts.bed -b ${sites} > ${prefix}.counts.blf.sites.bed || {
        echo "FAILED: $step"
	exit 1
    }
}

parse_args $@
run_commands
