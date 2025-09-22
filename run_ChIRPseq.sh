#!/bin/bash

print_help() {
  echo "Usage: $0 --oddFq <odd_R1 [odd_R2]> --evenFq <even_R1 [even_R2]> --inputFq <input_R1 [input_R2]> --prefix <output_prefix> --outputdir <output_directory> --repeatIndex <repeats_index_directory> --genomeIndex <genome_index_directory> --adapterFa <illumina_adapter> --sif <singularity_environment> --threads <number_of_threads> --binSize <binSize_for_bw> --g <hs> --peaktype <broad_or_narrow> --pval <p_value> --qval <q_value> --broad_cutoff <cutoff> --llocal <local_region> --keepdup <dup_setting> --nomodel <on_or_off> --nolambda <on_or_off> --callsummits <on_or_off> --extsize_val <extension_size> --shift_val <shift_bp>"
  echo "--oddFq              : Path to the odd fastq"
  echo "--evenFq             : Path to the even fastq"
  echo "--inputFq            : (optinal) Path to the input control fastq"
  echo "--prefix             : Prefix to the output files"
  echo "--outputdir          : path to the directory where the output will be stored"
  echo "--repeatIndex        : path to the directory where bowtie repeats reference build with prefix"
  echo "--genomeIndex        : path to the directory where bowtie genome reference build with prefix"
  echo "--adapterFa          : path to the adapter fasta"
  echo "--sif                : path to the singularity environment file"
  echo "--threads            : (optinal) Number of threads to use, default 8"
  echo "--binSize            : (optinal) Number of binsize to use, default 10"
  echo "--g                  : (optional) provide species code accepted by MACS3 (e.g. hs, mm, ce, dm), (default: hs)"
  echo "  --peaktype         : (optional) For MACS3 use 'broad' or 'narrow'; for SEACR use 'relaxed' or 'stringent' (default: narrow)."
  echo "  --pval             : (optional) MACS3 p-value cutoff for narrow/strong peaks. If specified, overrides q-value."
  echo "  --qval             : (optional) MACS3 q-value (FDR) cutoff for narrow/strong peaks (default: 0.01)."
  echo "  --broad_cutoff     : (optional) MACS3 p- or q-value cutoff for broad peaks, used only when --peaktype=broad (default: 0.1)."
  echo "  --llocal           : (optional) MACS3 --llocal setting for samples without control (default: 100000)."
  echo "  --keepdup          : (optional) MACS3 --keep-dup setting (default: all)."
  echo "  --nomodel          : (optional) MACS3 disable model building: on/off (default: on)."
  echo "  --nolambda         : (optional) MACS3 disable dynamic lambda: on/off (default: on)."
  echo "  --callsummits      : (optional) MACS3 call summits: on/off (default: on)."
  echo "  --extsize_val      : (optional) MACS3 fragment/extension size (bp); effective only when --nomodel=on."
  echo "  --shift_val        : (optional) MACS3 read 5' end shift (bp); effective only when --nomodel=on. Positive shifts 5'→3'; negative shifts 3'→5'."
  
}

if [[ "$#" -eq 0 || "$1" == "--help" ]]; then
  print_help
  exit 0
fi

oddFq_R1=""; oddFq_R2=""
evenFq_R1=""; evenFq_R2=""
inputFq_R1=""; inputFq_R2=""
prefix=""
outputdir=""
repeatIndex=""
genomeIndex=""
adapterFa=""
sif=""
threads=8
binSize=10
nomodel=${nomodel:-on}
callsummits=${callsummits:-on}
nolambda=${nolambda:-on}
shift_val=${shift_val:-""}
extsize_val=${extsize_val:-""}
peaktype=${peaktype:-narrow}       # narrow / broad
pval=${pval:-""}
qval=${qval:-""}
broad_cutoff=${broad_cutoff:-0.1}
llocal=${llocal:-100000}
keepdup=${keepdup:-all}
g=${g:-hs}                        # hs, mm, ce, dm

while [[ $# -gt 0 ]]; do
  case $1 in
    --oddFq)
      shift
      [[ $# -ge 1 ]] || { echo "Error: --oddFq requires at least one file"; exit 1; }
      oddFq_R1="$1"
      [[ -f "$oddFq_R1" ]] || { echo "Error: $oddFq_R1 not found"; exit 1; }
      if [[ $# -ge 2 && ! "$2" =~ ^-- ]]; then
      oddFq_R2="$2"
      [[ -f "$oddFq_R2" ]] || { echo "Error: $oddFq_R2 not found"; exit 1; }
      shift
      fi
      shift
      ;;
    --evenFq)
      shift
      [[ $# -ge 1 ]] || { echo "Error: --evenFq requires at least one file"; exit 1; }
      evenFq_R1="$1"
      [[ -f "$evenFq_R1" ]] || { echo "Error: $evenFq_R1 not found"; exit 1; }
      if [[ $# -ge 2 && ! "$2" =~ ^-- ]]; then
       evenFq_R2="$2"
      [[ -f "$evenFq_R2" ]] || { echo "Error: $evenFq_R2 not found"; exit 1; }
      shift
      fi
      shift
      ;;
    --inputFq)
      shift
      if [[ $# -ge 1 && ! "$1" =~ ^-- ]]; then
      inputFq_R1="$1"
      [[ -f "$inputFq_R1" ]] || { echo "Error: $inputFq_R1 not found"; exit 1; }
      shift
      if [[ $# -ge 1 && ! "$1" =~ ^-- ]]; then
      inputFq_R2="$1"
      [[ -f "$inputFq_R2" ]] || { echo "Error: $inputFq_R2 not found"; exit 1; }
      shift
      fi
      fi
      ;;
    --prefix) 
      if [[ -z "$2" ]]; then
        echo "Error: output prefix $2 not found!"
        exit 1
      fi
      prefix="$2"; shift 2 ;;
    --outputdir) 
      if [[ -z "$2" ]]; then
        echo "Error: output directory $2 not found!"
        exit 1
      fi
      outputdir="$2"; shift 2 ;;
    --repeatIndex) 
      if [[ -z "$2" ]]; then
        echo "Error: reference directory $2 not found!"
        exit 1
      fi
      repeatIndex="$2"; shift 2 ;;
    --genomeIndex) 
      if [[ -z "$2" ]]; then
        echo "Error: reference directory $2 not found!"
        exit 1
      fi
      genomeIndex="$2"; shift 2 ;;
    --adapterFa) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: adapter fa file $2 not found!"
        exit 1
      fi
      adapterFa="$2"; shift 2 ;;
    --sif) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: singularity sif file $2 not found!"
        exit 1
      fi
      sif="$2"; shift 2 ;;
    --threads) 
      threads="$2"; shift 2 ;;
    --binSize) 
      binSize="$2"; shift 2 ;;
    --g) 
      g="$2"; shift 2 ;;
    --peaktype) 
      peaktype="$2"; shift 2 ;;
    --nolambda) 
      nolambda="$2"; shift 2 ;;
    --pval) 
      pval="$2"; shift 2 ;;
    --qval) 
      qval="$2"; shift 2 ;;
    --broad_cutoff) 
      broad_cutoff="$2"; shift 2 ;;
    --llocal) 
      llocal="$2"; shift 2 ;;
    --keepdup) 
      keepdup="$2"; shift 2 ;;
    --nomodel) 
      nomodel="$2"; shift 2 ;;
    --callsummits) 
      callsummits="$2"; shift 2 ;;
    --extsize_val) 
      extsize_val="$2"; shift 2 ;;
    --shift_val) 
      shift_val="$2"; shift 2 ;;
    *) 
      echo "Unknown option $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$oddFq_R1" || -z "$evenFq_R1" || -z "$prefix" || -z "$outputdir" || -z "$repeatIndex" || -z "$genomeIndex" || -z "$adapterFa" || -z "$sif" ]]; then
  echo "Error: Missing required parameters"
  print_help
  exit 1
fi

########### mkdir
SING_EXEC1="singularity exec --cleanenv $sif"
$SING_EXEC1 mkdir -p ${outputdir}
$SING_EXEC1 mkdir -p ${outputdir}/qc
FastQCdir=${outputdir}/qc
$SING_EXEC1 mkdir -p ${outputdir}/trim
trimedfadir=${outputdir}/trim
$SING_EXEC1 mkdir -p ${outputdir}/bam
bowtieoutdir=${outputdir}/bam
$SING_EXEC1 mkdir -p ${outputdir}/bw
bwoutdir=${outputdir}/bw
$SING_EXEC1 mkdir -p ${outputdir}/figure
FigureDir=${outputdir}/figure
$SING_EXEC1 mkdir -p ${outputdir}/peak
pcoutdir=${outputdir}/peak
$SING_EXEC1 mkdir -p ${outputdir}/multiqc
multiqcdir=${outputdir}/multiqc

########### singularity command

SING_EXEC="singularity exec --cleanenv"
[[ -n "$oddFq_R1" ]] && oddFq_R1=$(${SING_EXEC1} readlink -f "$oddFq_R1")
[[ -n "$evenFq_R1" ]] && evenFq_R1=$(${SING_EXEC1} readlink -f "$evenFq_R1")
[[ -n "$inputFq_R1" ]] && inputFq_R1=$(${SING_EXEC1} readlink -f "$inputFq_R1")
[[ -n "$oddFq_R2" ]] && oddFq_R2=$(${SING_EXEC1} readlink -f "$oddFq_R2")
[[ -n "$evenFq_R2" ]] && evenFq_R2=$(${SING_EXEC1} readlink -f "$evenFq_R2")
[[ -n "$inputFq_R2" ]] && inputFq_R2=$(${SING_EXEC1} readlink -f "$inputFq_R2")
[[ -n "$outputdir" ]] && outputdir=$(${SING_EXEC1} readlink -f "$outputdir")
[[ -n "$repeatIndex" ]] && repeatIndex=$(${SING_EXEC1} readlink -f "$repeatIndex")
[[ -n "$genomeIndex" ]] && genomeIndex=$(${SING_EXEC1} readlink -f "$genomeIndex")
[[ -n "$adapterFa" ]] && adapterFa=$(${SING_EXEC1} readlink -f "$adapterFa")

[[ -n "$sif" ]] && sif=$( readlink -f "$sif")

bind_dirs=()
[[ -n "$oddFq_R1" ]] && bind_dirs+=("$(dirname "$oddFq_R1")")
[[ -n "$evenFq_R1" ]] && bind_dirs+=("$(dirname "$evenFq_R1")")
[[ -n "$inputFq_R1" ]] && bind_dirs+=("$(dirname "$inputFq_R1")")
[[ -n "$outputdir" ]] && bind_dirs+=("$outputdir")
[[ -n "$repeatIndex" ]] && bind_dirs+=("$(dirname "$repeatIndex")")
[[ -n "$genomeIndex" ]] && bind_dirs+=("$(dirname "$genomeIndex")")
[[ -n "$adapterFa" ]] && bind_dirs+=("$(dirname "$adapterFa")")
bind_dirs_unique=($(${SING_EXEC1} printf "%s\n" "${bind_dirs[@]}" | ${SING_EXEC1} sort -u))
for dir in "${bind_dirs_unique[@]}"; do
    real_dir=$(readlink -f "$dir")
    [[ -n "$real_dir" ]] && SING_EXEC+=" -B $real_dir:$real_dir"
done
SING_EXEC+=" $sif"

########### gzip
compress_if_needed() {
    local f="$1"
    local outdir="${outputdir}/rawdata"
    [[ -z "$f" ]] && return

    $SING_EXEC1 mkdir -p "$outdir"

    if [[ -f "$f" ]]; then
        if [[ "$f" != *.gz ]]; then
            local basef=$(basename "$f")
            local target="$outdir/$basef.gz"
            echo "gzip $f -> $target ..." >&2
            $SING_EXEC gzip -c "$f" > "$target"
            echo "$target"
        else
            echo "$f"
        fi
    fi
}

oddFq_R1=$(compress_if_needed "$oddFq_R1")
oddFq_R2=$(compress_if_needed "$oddFq_R2")
evenFq_R1=$(compress_if_needed "$evenFq_R1")
evenFq_R2=$(compress_if_needed "$evenFq_R2")
inputFq_R1=$(compress_if_needed "$inputFq_R1")
inputFq_R2=$(compress_if_needed "$inputFq_R2")

########### processing
process_group(){
    local label="$1"
    local r1="$2"
    local r2="$3"

    [[ -z "$r1" ]] && { echo "Skip ${label}: no R1 provided." >&2; return; }
    
    local base=${label}
    echo "### Processing group: ${label} (base=${base})" >&2

    if [[ -n "$r2" ]]; then
        # ---- PAIRED-END workflow ----
        echo "Detected PAIRED-END for ${label}:" >&2
        echo "  R1: $r1" >&2
        echo "  R2: $r2" >&2

        # FastQC
        echo "Running FastQC for ${label}..." >&2
        $SING_EXEC fastqc -o "${FastQCdir}" -t "${threads}" -q "$r1" "$r2"
	
        # cutadapt paired
        local trimmed_r1="${trimedfadir}/${base}_val_1.fq.gz"
        local trimmed_r2="${trimedfadir}/${base}_val_2.fq.gz"

        # trim
        echo "Running trim_galore for ${label}..." >&2
        $SING_EXEC bash -c "trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 --stringency 3 --paired '$r1' '$r2' -a file:${adapterFa} -o '${trimedfadir}' --basename '${base}' -j ${threads}"
        
        # remove repeats
        echo "Running repeats mapping for ${label}..." >&2
        $SING_EXEC bowtie2 -p ${threads} -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --rg-id {base} -x $repeatIndex -1 ${trimmed_r1} -2 ${trimmed_r2} --very-sensitive --un-conc ${bowtieoutdir}/${base}.fastq 2> ${bowtieoutdir}/${base}.repeats.bowtie.stats | $SING_EXEC samtools view -bS - | $SING_EXEC samtools sort -@ ${threads} -o ${bowtieoutdir}/${base}_repeat_sorted.bam
        $SING_EXEC gzip -f ${bowtieoutdir}/${base}.1.fastq > ${bowtieoutdir}/${base}.1.fastq.gz
        $SING_EXEC gzip -f ${bowtieoutdir}/${base}.2.fastq > ${bowtieoutdir}/${base}.2.fastq.gz
        rm ${bowtieoutdir}/${base}.1.fastq
        rm ${bowtieoutdir}/${base}.2.fastq
        
        # Mapping
        echo "Running bowtie2 mapping for ${label}..." >&2
        $SING_EXEC bowtie2 -t -k 1 --end-to-end --sensitive -p ${threads} --fr --no-mixed --no-discordant -X 1000 -x ${genomeIndex} -1 ${bowtieoutdir}/${base}.1.fastq.gz -2 ${bowtieoutdir}/${base}.2.fastq.gz 2> ${bowtieoutdir}/${base}.bowtie.stats | $SING_EXEC samtools view -q 255 -bS - | $SING_EXEC samtools sort -n -@ ${threads} - -T $base -O BAM -o ${bowtieoutdir}/${base}.sort.bam
        rm ${bowtieoutdir}/${base}.1.fastq.gz
        rm ${bowtieoutdir}/${base}.2.fastq.gz
    else
        # ---- SINGLE-END workflow ----
        echo "Detected SINGLE-END for ${label}: R1=$r1" >&2

        # FastQC
        echo "Running FastQC for ${label}..." >&2
        $SING_EXEC fastqc -o "${FastQCdir}" -t "${threads}" -q "$r1"

        # trim single
        local trimmed_se="${trimedfadir}/${base}_trimmed.fq.gz"
        echo "Running trim_galore for ${label}..." >&2
        $SING_EXEC bash -c "trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 --stringency 3 '$r1' -a file:${adapterFa} -o '${trimedfadir}' --basename '${base}' -j ${threads}"

        # remove repeats
        echo "Running repeats mapping for ${label}..." >&2
        $SING_EXEC bowtie2 -p ${threads} -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --rg-id {base} -x $repeatIndex -U ${trimmed_se} --un ${bowtieoutdir}/${base}.fastq 2> ${bowtieoutdir}/${base}.repeats.bowtie.stats | $SING_EXEC samtools view -bS - | $SING_EXEC samtools sort -@ ${threads} -o ${bowtieoutdir}/${base}_repeat_sorted.bam
        $SING_EXEC gzip -f ${bowtieoutdir}/${base}.fastq > ${bowtieoutdir}/${base}.fastq.gz
        #$SING_EXEC rm ${bowtieoutdir}/${base}.fastq
				
        # Mapping
        echo "Running bowtie2 mapping for ${label}..." >&2
        $SING_EXEC bowtie2 -t -k 1 --end-to-end --sensitive -p ${threads} -x ${genomeIndex} -U ${bowtieoutdir}/${base}.fastq.gz 2> ${bowtieoutdir}/${base}.bowtie.stats | $SING_EXEC samtools view -q 255 -bS - | $SING_EXEC samtools sort -n -@ ${threads} - -T $base -O BAM -o ${bowtieoutdir}/${base}.sort.bam
       # $SING_EXEC rm ${bowtieoutdir}/${base}.S.fastq.gz
    fi

    # ---------- post ----------
    echo "Markdup / index / flagstat for ${base}" >&2
    $SING_EXEC bash -c "samtools fixmate -r -m '${bowtieoutdir}/${base}.sort.bam' - | samtools sort -@ '${threads}' - | samtools markdup -r -s - '${bowtieoutdir}/${base}.DeDup.bam' 2> '${bowtieoutdir}/${base}.markdup.log'"
    $SING_EXEC samtools flagstat "${bowtieoutdir}/${base}.DeDup.bam" > "${bowtieoutdir}/${base}.flagstat.txt"
    $SING_EXEC samtools index "${bowtieoutdir}/${base}.DeDup.bam"
    rm -f "${bowtieoutdir}/${base}.sort.bam"

    echo "Generate bigWig for ${base}" >&2
    $SING_EXEC bamCoverage -b "${bowtieoutdir}/${base}.DeDup.bam" \
        -o "${bwoutdir}/${base}.DeDup.bw" \
        -p "${threads}" --binSize "${binSize}" --normalizeUsing CPM

    echo "Done ${label} (${base})" >&2
}

# ---------- running ----------
echo "Starting processing for odd group..." >&2
process_group "${prefix}.odd"  "$oddFq_R1"  "$oddFq_R2"

echo "Starting processing for even group..." >&2
process_group "${prefix}.even" "$evenFq_R1" "$evenFq_R2"

echo "Starting processing for input group..." >&2
process_group "${prefix}.input" "$inputFq_R1" "$inputFq_R2"


########### correlation
bwfiles=("${bwoutdir}/${prefix}.odd.DeDup.bw" "${bwoutdir}/${prefix}.even.DeDup.bw")
bamfiles=("${bowtieoutdir}/${prefix}.odd.DeDup.bam" "${bowtieoutdir}/${prefix}.even.DeDup.bam")
if [[ -f "${bwoutdir}/${prefix}.input.DeDup.bw" ]]; then
    bwfiles+=("${bwoutdir}/${prefix}.input.DeDup.bw")
    bamfiles+=("${bowtieoutdir}/${prefix}.input.DeDup.bam")
fi
# labels
labels=()
for f in "${bwfiles[@]}"; do
    labels+=("$(basename -s .bw "$f")")
done

echo "Check the correlation"
$SING_EXEC multiBigwigSummary bins -b "${bwfiles[@]}" --labels "${labels[@]}" -out ${FigureDir}/BW_compare_PCA.npz -p ${threads}
$SING_EXEC plotPCA -in ${FigureDir}/BW_compare_PCA.npz -o ${FigureDir}/BW_compare_PCA.pdf
$SING_EXEC plotCorrelation -in ${FigureDir}/BW_compare_PCA.npz --corMethod pearson --skipZeros -o ${FigureDir}/BW_compare_cor.pdf --whatToPlot heatmap --colorMap RdYlBu --plotNumbers
$SING_EXEC rm ${FigureDir}/BW_compare_PCA.npz

echo "check the seq quality"
$SING_EXEC plotFingerprint -b "${bamfiles[@]}" --labels "${labels[@]}" --skipZeros --plotFile ${FigureDir}/fingerprints.pdf

########### merge bam
$SING_EXEC samtools merge "${bowtieoutdir}/${prefix}.merged.DeDup.bam" "${bowtieoutdir}/${prefix}.odd.DeDup.bam" "${bowtieoutdir}/${prefix}.even.DeDup.bam"
$SING_EXEC samtools index "${bowtieoutdir}/${prefix}.merged.DeDup.bam"
$SING_EXEC bamCoverage -b "${bowtieoutdir}/${prefix}.merged.DeDup.bam" -o "${bwoutdir}/${prefix}.merged.DeDup.bw" -p "${threads}" --binSize "${binSize}" --normalizeUsing CPM

########### peak compare
echo "peak calling"
# -------------------- set --------------------
treatment="${prefix}.merged"
control="${prefix}.input"

treat_bam="${bowtieoutdir}/${treatment}.DeDup.bam"
ctrl_bam="${bowtieoutdir}/${control}.DeDup.bam"

# pair-single
bam_format="BAM"
[[ -n "$oddFq_R2" ]] && bam_format="BAMPE"

# -------------------- macs3 --------------------
build_macs_opts() {
    local has_ctrl="$1"  # "yes" or "no"
    local opts=""

    [[ $nomodel == "on" ]] && opts+=" --nomodel"
    [[ $callsummits == "on" ]] && opts+=" --call-summits"
    [[ $nolambda == "on" ]] && opts+=" --nolambda"
    [[ -n "$shift_val" ]] && opts+=" --shift $shift_val"
    [[ -n "$extsize_val" ]] && opts+=" --extsize $extsize_val"

    if [[ "$peaktype" == "broad" ]]; then
        opts+=" --broad"
        [[ -n "$broad_cutoff" ]] && opts+=" --broad-cutoff $broad_cutoff"
    fi

    [[ -n "$pval" ]] && opts+=" -p $pval"
    [[ -n "$qval" ]] && opts+=" -q $qval"
    [[ -z "$pval" && -z "$qval" ]] && opts+=" -q 0.01"

    # no control, then add --llocal
    [[ "$has_ctrl" == "no" && -n "$llocal" ]] && opts+=" --llocal $llocal"

    echo "$opts"
}


# -------------------- enrichment --------------------
plot_peak_heatmap() {
    local treatment="$1"
    local control="$2"      # 可以为空，treat-only
    local peak_file="$3"
    local peaktype="$4"
    local out_prefix="$5"

    [[ ! -s "$peak_file" ]] && { echo "Skip: $peak_file is missing or empty."; return; }

    if [[ "$peaktype" == "broad" ]]; then
        if [[ -n "$control" ]]; then
            $SING_EXEC computeMatrix scale-regions -p ${threads} \
                -S ${bwoutdir}/${treatment}.DeDup.bw ${bwoutdir}/${control}.DeDup.bw \
                -R $peak_file -a 5000 -b 5000 --regionBodyLength 1000 \
                --missingDataAsZero -o ${FigureDir}/${out_prefix}.gz
        else
            $SING_EXEC computeMatrix scale-regions -p ${threads} \
                -S ${bwoutdir}/${treatment}.DeDup.bw \
                -R $peak_file -a 5000 -b 5000 --regionBodyLength 1000 \
                --missingDataAsZero -o ${FigureDir}/${out_prefix}.gz
        fi
        $SING_EXEC plotHeatmap -m ${FigureDir}/${out_prefix}.gz \
            -out ${FigureDir}/${out_prefix}.pdf \
            --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4 \
            --startLabel "Start" --endLabel "End"
    else
        if [[ -n "$control" ]]; then
            $SING_EXEC computeMatrix scale-regions -p ${threads} \
                -S ${bwoutdir}/${treatment}.DeDup.bw ${bwoutdir}/${control}.DeDup.bw \
                -R $peak_file -a 5000 -b 5000 --regionBodyLength 0 \
                --missingDataAsZero -o ${FigureDir}/${out_prefix}.gz
        else
            $SING_EXEC computeMatrix scale-regions -p ${threads} \
                -S ${bwoutdir}/${treatment}.DeDup.bw \
                -R $peak_file -a 5000 -b 5000 --regionBodyLength 0 \
                --missingDataAsZero -o ${FigureDir}/${out_prefix}.gz
        fi
        $SING_EXEC plotHeatmap -m ${FigureDir}/${out_prefix}.gz \
            -out ${FigureDir}/${out_prefix}.pdf \
            --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4 \
            --refPointLabel "Peak"
    fi
}

# -------------------- running MACS3 --------------------
extra_opts=$(build_macs_opts)

if [[ -f "$treat_bam" && -f "$ctrl_bam" ]]; then
    echo "Calling peaks for $treatment vs $control ..."
    extra_opts=$(build_macs_opts "yes")
    $SING_EXEC macs3 callpeak \
        -t "$treat_bam" -c "$ctrl_bam" \
        -f "$bam_format" -g "$g" \
        --outdir "$pcoutdir" -n "${treatment}.vs.${control}" \
        $extra_opts --keep-dup "$keepdup" \
        2> >(tee -a "${pcoutdir}/${treatment}.vs.${control}.macs3.stats" >&2) \
        1> "${pcoutdir}/${treatment}.vs.${control}.macs3.stdout"
    rm "${pcoutdir}/${treatment}.vs.${control}.macs3.stdout"

    peak_file="${pcoutdir}/${treatment}.vs.${control}_peaks.${peaktype}Peak"
    plot_peak_heatmap "$treatment" "$control" "$peak_file" "$peaktype" "${treatment}.vs.${control}.peak"

elif [[ -f "$treat_bam" && ! -f "$ctrl_bam" ]]; then
    echo "Control BAM not found. Running MACS3 in treatment-only mode ..."
    extra_opts=$(build_macs_opts "no")
    $SING_EXEC macs3 callpeak \
        -t "$treat_bam" \
        -f "$bam_format" -g "$g" \
        --outdir "$pcoutdir" -n "$treatment" \
        $extra_opts --keep-dup "$keepdup" \
        2> >(tee -a "${pcoutdir}/${treatment}.macs3.stats" >&2) \
        1> "${pcoutdir}/${treatment}.macs3.stdout"
    rm "${pcoutdir}/${treatment}.macs3.stdout"

    peak_file="${pcoutdir}/${treatment}_peaks.${peaktype}Peak"
    plot_peak_heatmap "$treatment" "" "$peak_file" "$peaktype" "${treatment}.peak"

else
    echo "Treatment BAM not found. Nothing to do."
fi

########### multiqc
$SING_EXEC multiqc ${FastQCdir}/* ${trimedfadir}/* ${bowtieoutdir}/*.bowtie.stats ${bowtieoutdir}/*flagstat.txt -o ${multiqcdir} --force

########### remove
$SING_EXEC rm -r ${FastQCdir}
$SING_EXEC rm -r ${trimedfadir}
$SING_EXEC rm -f ${FigureDir}/*.peak.gz



