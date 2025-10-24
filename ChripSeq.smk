from pathlib import Path

# Load configuration file
configfile: "config.yaml"

# Define output directories (created automatically by Snakemake)
OUTPUT_DIR = Path(config["outputdir"]).resolve()
QC_DIR = OUTPUT_DIR / "rawdata.qc"
TRIM_DIR = OUTPUT_DIR / "trim"
BAM_DIR = OUTPUT_DIR / "bam"
BW_DIR = OUTPUT_DIR / "bw"
FIGURE_DIR = OUTPUT_DIR / "figure"
PEAK_DIR = OUTPUT_DIR / "peak"
MULTIQC_DIR = OUTPUT_DIR / "multiqc"

rule all:
    """Final target to ensure all pipeline outputs are generated"""
    input:
        # Post-processing outputs
        BW_DIR / f"{config['prefix']}.odd.DeDup.bw",
        BW_DIR / f"{config['prefix']}.even.DeDup.bw",
        
        # Merged outputs
        BW_DIR / f"{config['prefix']}.merged.DeDup.bw",
        
        # Figures
        FIGURE_DIR / "BW_compare_PCA.pdf",
        FIGURE_DIR / "BW_compare_cor.pdf",
        FIGURE_DIR / "fingerprints.pdf",
        
        # Peak results
        PEAK_DIR / "peakcalling.log",
        FIGURE_DIR / f"{config['prefix']}.peak.pdf",

        MULTIQC_DIR / "multiqc_report.html"


# --------------------------
# QC for raw data (odd group)
# --------------------------
rule fastqc_odd:
    """Run FastQC on raw odd-group reads"""
    input:
        r1 = config["oddFq"]["R1"],
    output:
        log = temp(QC_DIR / f"{config['prefix']}.odd.fastqc.log")
    params:
        threads = config["threads"],
        out_dir = QC_DIR,
        r2 = config["oddFq"].get("R2", "")
    container: config["sif"]
    shell:
        """
        if [ -n "{params.r2}" ]; then
            fastqc -o {params.out_dir} -t {params.threads} -q {input.r1} {params.r2} 2> {output.log}
        else
            fastqc -o {params.out_dir} -t {params.threads} -q {input.r1} 2> {output.log}
        fi
        """


# ---------------------------
# QC for raw data (even group)
# ---------------------------
rule fastqc_even:
    """Run FastQC on raw even-group reads"""
    input:
        r1 = config["evenFq"]["R1"],
    output:
        log = temp(QC_DIR / f"{config['prefix']}.even.fastqc.log")
    params:
        threads = config["threads"],
        out_dir = QC_DIR,
        r2 = config["evenFq"].get("R2", "")
    container: config["sif"]
    shell:
        """
        if [ -n "{params.r2}" ]; then
            fastqc -o {params.out_dir} -t {params.threads} -q {input.r1} {params.r2} 2> {output.log}
        else
            fastqc -o {params.out_dir} -t {params.threads} -q {input.r1} 2> {output.log}
        fi
        """


# --------------------------------
# QC for raw data (input group, optional)
# --------------------------------
rule fastqc_input:
    """Run FastQC on raw input-group reads (if provided)"""
    input:
        oddr1 = config["oddFq"]["R1"],
    output:
        log = temp(QC_DIR / f"{config['prefix']}.input.fastqc.log")
    params:
        threads = config["threads"],
        out_dir = QC_DIR,
        r1 = config["inputFq"].get("R1", ""),
        r2 = config["inputFq"].get("R2", "")
    container: config["sif"]
    shell:
        """
        if [ -z "{params.r1}" ]; then
            echo "No input provided, skipping fastqc" > {output.log}
        else
            if [ -n "{params.r2}" ]; then
                fastqc -o {params.out_dir} -t {params.threads} -q {params.r1} {params.r2} 2> {output.log}
            else
                fastqc -o {params.out_dir} -t {params.threads} -q {params.r1} 2> {output.log}
            fi
        fi
        """


# --------------------------
# Trim odd group reads
# --------------------------
rule trim_galore_odd:
    """Trim adapters/low-quality bases for odd group"""
    input:
        # Depend on FastQC log (QC-trim connection)
        r1 = config["oddFq"]["R1"]
    output:
        log = temp(TRIM_DIR / f"{config['prefix']}.odd.trim.log")
    params:
        threads = config["threads"],
        adapter = config["adapterFa"],
        out_dir = TRIM_DIR,
        basename = f"{config['prefix']}.odd",
        r2 = config["oddFq"].get("R2", "")
    container: config["sif"]
    shell:
        """
        if [ -n "{params.r2}" ]; then
            trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 \
                --stringency 3 --paired {input.r1} {params.r2} \
                -a file:{params.adapter} -o {params.out_dir} \
                --basename {params.basename} -j {params.threads} > {output.log} 2>&1
        else
            trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 \
                --stringency 3 {input.r1} \
                -a file:{params.adapter} -o {params.out_dir} \
                --basename {params.basename} -j {params.threads} > {output.log} 2>&1
        fi
        """


# ---------------------------
# Trim even group reads
# ---------------------------
rule trim_galore_even:
    """Trim adapters/low-quality bases for even group"""
    input:
        r1 = config["evenFq"]["R1"]
    output:
        log = temp(TRIM_DIR / f"{config['prefix']}.even.trim.log")
    params:
        threads = config["threads"],
        adapter = config["adapterFa"],
        out_dir = TRIM_DIR,
        basename = f"{config['prefix']}.even",
        r2 = config["evenFq"].get("R2", "")
    container: config["sif"]
    shell:
        """
        if [ -n "{params.r2}" ]; then
            trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 \
                --stringency 3 --paired {input.r1} {params.r2} \
                -a file:{params.adapter} -o {params.out_dir} \
                --basename {params.basename} -j {params.threads} > {output.log} 2>&1
        else
            trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 \
                --stringency 3 {input.r1} \
                -a file:{params.adapter} -o {params.out_dir} \
                --basename {params.basename} -j {params.threads} > {output.log} 2>&1
        fi
        """


# --------------------------------
# Trim input group reads (optional)
# --------------------------------
rule trim_galore_input:
    """Trim adapters/low-quality bases for input group (if provided)"""
    input:
        r1 = config["oddFq"]["R1"]
    output:
        log = temp(TRIM_DIR / f"{config['prefix']}.input.trim.log")
    params:
        threads = config["threads"],
        adapter = config["adapterFa"],
        out_dir = TRIM_DIR,
        basename = f"{config['prefix']}.input",
        r1 = config["inputFq"].get("R1", ""),
        r2 = config["inputFq"].get("R2", "")
    container: config["sif"]
    shell:
        """
        if [ -z "{params.r1}" ]; then
            echo "No input provided, skipping trim_galore." > {output.log}
        else
            if [ -n "{params.r2}" ]; then
                trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 \
                    --stringency 3 --paired {params.r1} {params.r2} \
                    -a file:{params.adapter} -o {params.out_dir} \
                    --basename {params.basename} -j {params.threads} > {output.log} 2>&1
            else
                trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 \
                    --stringency 3 {params.r1} \
                    -a file:{params.adapter} -o {params.out_dir} \
                    --basename {params.basename} -j {params.threads} > {output.log} 2>&1
            fi
        fi
        """


# --------------------------
# Map odd reads to repeats
# --------------------------
rule bowtie2_repeats_odd:
    """Map odd-group trimmed reads to repeats (filter repeats)"""
    input:
        trim_log = TRIM_DIR / f"{config['prefix']}.odd.trim.log"
    output:
        log = BAM_DIR / f"{config['prefix']}.odd.repeats.stats"
    params:
        threads = config["threads"],
        index = config["repeatIndex"],
        basename = f"{config['prefix']}.odd"
    container: config["sif"]
    shell:
        """
        if [ -e "{TRIM_DIR}/{params.basename}_val_2.fq.gz" ]; then
            bowtie2 -p {params.threads} -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
                --rg-id {params.basename} -x {params.index} \
                -1 {TRIM_DIR}/{params.basename}_val_1.fq.gz -2 {TRIM_DIR}/{params.basename}_val_2.fq.gz --very-sensitive \
                --un-conc {BAM_DIR}/{params.basename}.fastq \
                2> {output.log} | \
                samtools view -bS - | samtools sort -@ {params.threads} -o {BAM_DIR}/{params.basename}_repeat_sorted.bam
            gzip -f {BAM_DIR}/{params.basename}.1.fastq > ${BAM_DIR}/{params.basename}.1.fastq.gz
            gzip -f {BAM_DIR}/{params.basename}.2.fastq > ${BAM_DIR}/{params.basename}.2.fastq.gz
            rm {TRIM_DIR}/{params.basename}_val_1.fq.gz
            rm {TRIM_DIR}/{params.basename}_val_2.fq.gz
        else
            bowtie2 -p {params.threads} -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
                --rg-id {params.basename} -x {params.index} \
                -U {TRIM_DIR}/{params.basename}_trimmed.fq.gz --un {BAM_DIR}/{params.basename}.fastq \
                2> {output.log} | \
                samtools view -bS - | samtools sort -@ {params.threads} -o {BAM_DIR}/{params.basename}_repeat_sorted.bam
            gzip -f {BAM_DIR}/{params.basename}.fastq > {BAM_DIR}/{params.basename}.fastq.gz
            rm {TRIM_DIR}/{params.basename}_trimmed.fq.gz
        fi
        rm {BAM_DIR}/{params.basename}_repeat_sorted.bam
        """


# ---------------------------
# Map even reads to repeats
# ---------------------------
rule bowtie2_repeats_even:
    """Map even-group trimmed reads to repeats (filter repeats)"""
    input:
        trim_log = TRIM_DIR / f"{config['prefix']}.even.trim.log"
    output:
        log = BAM_DIR / f"{config['prefix']}.even.repeats.stats"
    params:
        threads = config["threads"],
        index = config["repeatIndex"],
        basename = f"{config['prefix']}.even"
    container: config["sif"]
    shell:
        """
        if [ -e "{TRIM_DIR}/{params.basename}_val_2.fq.gz" ]; then
            bowtie2 -p {params.threads} -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
                --rg-id {params.basename} -x {params.index} \
                -1 {TRIM_DIR}/{params.basename}_val_1.fq.gz -2 {TRIM_DIR}/{params.basename}_val_2.fq.gz --very-sensitive \
                --un-conc {BAM_DIR}/{params.basename}.fastq \
                2> {output.log} | \
                samtools view -bS - | samtools sort -@ {params.threads} -o {BAM_DIR}/{params.basename}_repeat_sorted.bam
            gzip -f ${BAM_DIR}/{params.basename}.1.fastq > ${BAM_DIR}/{params.basename}.1.fastq.gz
            gzip -f ${BAM_DIR}/{params.basename}.2.fastq > ${BAM_DIR}/{params.basename}.2.fastq.gz
            rm {TRIM_DIR}/{params.basename}_val_1.fq.gz
            rm {TRIM_DIR}/{params.basename}_val_2.fq.gz
        else
            bowtie2 -p {params.threads} -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
                --rg-id {params.basename} -x {params.index} \
                -U {TRIM_DIR}/{params.basename}_trimmed.fq.gz --un {BAM_DIR}/{params.basename}.fastq \
                2> {output.log} | \
                samtools view -bS - | samtools sort -@ {params.threads} -o {BAM_DIR}/{params.basename}_repeat_sorted.bam
            gzip -f {BAM_DIR}/{params.basename}.fastq > {BAM_DIR}/{params.basename}.fastq.gz
            rm {TRIM_DIR}/{params.basename}_trimmed.fq.gz
        fi
        rm {BAM_DIR}/{params.basename}_repeat_sorted.bam
        """


# --------------------------------
# Map input reads to repeats (optional)
# --------------------------------
rule bowtie2_repeats_input:
    """Map input-group trimmed reads to repeats (filter repeats)"""
    input:
        trim_log = TRIM_DIR / f"{config['prefix']}.input.trim.log"
    output:
        log = BAM_DIR / f"{config['prefix']}.input.repeats.stats"
    params:
        threads = config["threads"],
        index = config["repeatIndex"],
        basename = f"{config['prefix']}.input"
    container: config["sif"]
    shell:
        """
        if [ ! -e "{TRIM_DIR}/{params.basename}_trimmed.fq.gz" -a ! -e "{TRIM_DIR}/{params.basename}_val_1.fq.gz" ]; then
            echo "No input provided, skipping bowtie2.repeat." > {output.log}
        else
            if [ -e "{TRIM_DIR}/{params.basename}_val_2.fq.gz" ]; then
                bowtie2 -p {params.threads} -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
                    --rg-id {params.basename} -x {params.index} \
                    -1 {TRIM_DIR}/{params.basename}_val_1.fq.gz -2 {TRIM_DIR}/{params.basename}_val_2.fq.gz --very-sensitive \
                    --un-conc {BAM_DIR}/{params.basename}.fastq \
                    2> {output.log} | \
                    samtools view -bS - | samtools sort -@ {params.threads} -o {BAM_DIR}/{params.basename}_repeat_sorted.bam
                gzip -f {BAM_DIR}/{params.basename}.1.fastq > {BAM_DIR}/{params.basename}.1.fastq.gz
                gzip -f {BAM_DIR}/{params.basename}.2.fastq > {BAM_DIR}/{params.basename}.2.fastq.gz
                rm {TRIM_DIR}/{params.basename}_val_1.fq.gz
                rm {TRIM_DIR}/{params.basename}_val_2.fq.gz
            else
                bowtie2 -p {params.threads} -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
                    --rg-id {params.basename} -x {params.index} \
                    -U {TRIM_DIR}/{params.basename}_trimmed.fq.gz --un {BAM_DIR}/{params.basename}.fastq \
                    2> {output.log} | \
                    samtools view -bS - | samtools sort -@ {params.threads} -o {BAM_DIR}/{params.basename}_repeat_sorted.bam
                gzip -f {BAM_DIR}/{params.basename}.fastq > {BAM_DIR}/{params.basename}.fastq.gz
                rm {TRIM_DIR}/{params.basename}_trimmed.fq.gz
            fi
            rm {BAM_DIR}/{params.basename}_repeat_sorted.bam
        fi
        """


# --------------------------
# Map odd non-repeat reads to genome
# --------------------------
rule bowtie2_genome_odd:
    """Map odd-group non-repeat reads to genome"""
    input:
        log = BAM_DIR / f"{config['prefix']}.odd.repeats.stats"
    output:
        sorted_bam = BAM_DIR / f"{config['prefix']}.odd.sort.bam",
        log = BAM_DIR / f"{config['prefix']}.odd.genome.stats"
    params:
        threads = config["threads"],
        index = config["genomeIndex"],
        basename = f"{config['prefix']}.odd"
    container: config["sif"]
    shell:
        """
        if [ -e "{BAM_DIR}/{params.basename}.2.fastq.gz" ]; then
            bowtie2 -t -k 1 --end-to-end --sensitive -p {params.threads} \
                --fr --no-mixed --no-discordant -X 1000 \
                -x {params.index} -1 {BAM_DIR}/{params.basename}.1.fastq.gz -2 {BAM_DIR}/{params.basename}.2.fastq.gz \
                2> {BAM_DIR}/{params.basename}.genome.stats | \
                samtools view -q 255 -bS - | samtools sort -n -@ {params.threads} \
                - -T {params.basename} -O BAM -o {output.sorted_bam}
            rm {BAM_DIR}/{params.basename}.1.fastq.gz {BAM_DIR}/{params.basename}.2.fastq.gz
        else
            bowtie2 -t -k 1 --end-to-end --sensitive -p {params.threads} \
                -x {params.index} -U {BAM_DIR}/{params.basename}.fastq.gz \
                2> {BAM_DIR}/{params.basename}.genome.stats | \
                samtools view -q 255 -bS - | samtools sort -n -@ {params.threads} \
                - -T {params.basename} -O BAM -o {output.sorted_bam}
            rm {BAM_DIR}/{params.basename}.fastq.gz
        fi
        """


# ---------------------------
# Map even non-repeat reads to genome
# ---------------------------
rule bowtie2_genome_even:
    """Map even-group non-repeat reads to genome"""
    input:
        log = BAM_DIR / f"{config['prefix']}.even.repeats.stats"
    output:
        sorted_bam = BAM_DIR / f"{config['prefix']}.even.sort.bam",
        log = BAM_DIR / f"{config['prefix']}.even.genome.stats"
    params:
        threads = config["threads"],
        index = config["genomeIndex"],
        basename = f"{config['prefix']}.even"
    container: config["sif"]
    shell:
        """
        if [ -e "{BAM_DIR}/{params.basename}.2.fastq.gz" ]; then
            bowtie2 -t -k 1 --end-to-end --sensitive -p {params.threads} \
                --fr --no-mixed --no-discordant -X 1000 \
                -x {params.index} -1 {BAM_DIR}/{params.basename}.1.fastq.gz -2 {BAM_DIR}/{params.basename}.2.fastq.gz \
                2> {BAM_DIR}/{params.basename}.genome.stats | \
                samtools view -q 255 -bS - | samtools sort -n -@ {params.threads} \
                - -T {params.basename} -O BAM -o {output.sorted_bam}
            rm {BAM_DIR}/{params.basename}.1.fastq.gz {BAM_DIR}/{params.basename}.2.fastq.gz
        else
            bowtie2 -t -k 1 --end-to-end --sensitive -p {params.threads} \
                -x {params.index} -U {BAM_DIR}/{params.basename}.fastq.gz \
                2> {BAM_DIR}/{params.basename}.genome.stats | \
                samtools view -q 255 -bS - | samtools sort -n -@ {params.threads} \
                - -T {params.basename} -O BAM -o {output.sorted_bam}
            rm {BAM_DIR}/{params.basename}.fastq.gz
        fi
        """


# --------------------------------
# Map input non-repeat reads to genome (optional)
# --------------------------------
rule bowtie2_genome_input:
    """Map input-group non-repeat reads to genome"""
    input:
        log = BAM_DIR / f"{config['prefix']}.input.repeats.stats"
    output:
        log = BAM_DIR / f"{config['prefix']}.input.genome.stats"
    params:
        threads = config["threads"],
        index = config["genomeIndex"],
        basename = f"{config['prefix']}.input"
    container: config["sif"]
    shell:
        """
        if [ ! -e "{BAM_DIR}/{params.basename}.1.fastq.gz" -a ! -e "{BAM_DIR}/{params.basename}.fastq.gz" ]; then
            echo "No input provided, skipping bowtie2.genome." > {output.log}
        else
            if [ -e "{BAM_DIR}/{params.basename}.2.fastq.gz" ]; then
                bowtie2 -t -k 1 --end-to-end --sensitive -p {params.threads} \
                    --fr --no-mixed --no-discordant -X 1000 \
                    -x {params.index} -1 {BAM_DIR}/{params.basename}.1.fastq.gz -2 {BAM_DIR}/{params.basename}.2.fastq.gz \
                    2> {output.log} | \
                    samtools view -q 255 -bS - | samtools sort -n -@ {params.threads} \
                    - -T {params.basename} -O BAM -o {BAM_DIR}/{params.basename}.sort.bam
                rm {BAM_DIR}/{params.basename}.1.fastq.gz {BAM_DIR}/{params.basename}.2.fastq.gz
            else
                bowtie2 -t -k 1 --end-to-end --sensitive -p {params.threads} \
                    -x {params.index} -U {BAM_DIR}/{params.basename}.fastq.gz \
                    2> {output.log} | \
                    samtools view -q 255 -bS - | samtools sort -n -@ {params.threads} \
                    - -T {params.basename} -O BAM -o {BAM_DIR}/{params.basename}.sort.bam
                rm {BAM_DIR}/{params.basename}.fastq.gz
            fi
        fi
        """


# --------------------------
# Post-process odd BAM (dedup)
# --------------------------
rule postprocess_bam_odd:
    """Deduplicate and index odd-group genome BAM"""
    input:
        sorted_bam = BAM_DIR / f"{config['prefix']}.odd.sort.bam"
    output:
        dedup_bam = BAM_DIR / f"{config['prefix']}.odd.DeDup.bam",
        dedup_bai = BAM_DIR / f"{config['prefix']}.odd.DeDup.bam.bai",
        flagstat = BAM_DIR / f"{config['prefix']}.odd.flagstat.txt"
    params:
        threads = config["threads"],
        basename = f"{config['prefix']}.odd"
    container: config["sif"]
    shell:
        """
        samtools fixmate -r -m {input.sorted_bam} - | \
            samtools sort -@ {params.threads} - | \
            samtools markdup -r -s - {output.dedup_bam} 2> {BAM_DIR}/{params.basename}.markdup.log
        samtools flagstat {output.dedup_bam} > {output.flagstat}
        samtools index {output.dedup_bam}
        rm {input.sorted_bam}
        """


# ---------------------------
# Post-process even BAM (dedup)
# ---------------------------
rule postprocess_bam_even:
    """Deduplicate and index even-group genome BAM"""
    input:
        sorted_bam = BAM_DIR / f"{config['prefix']}.even.sort.bam"
    output:
        dedup_bam = BAM_DIR / f"{config['prefix']}.even.DeDup.bam",
        dedup_bai = BAM_DIR / f"{config['prefix']}.even.DeDup.bam.bai",
        flagstat = BAM_DIR / f"{config['prefix']}.even.flagstat.txt"
    params:
        threads = config["threads"],
        basename = f"{config['prefix']}.even"
    container: config["sif"]
    shell:
        """
        samtools fixmate -r -m {input.sorted_bam} - | \
            samtools sort -@ {params.threads} - | \
            samtools markdup -r -s - {output.dedup_bam} 2> {BAM_DIR}/{params.basename}.markdup.log
        samtools flagstat {output.dedup_bam} > {output.flagstat}
        samtools index {output.dedup_bam}
        rm {input.sorted_bam}
        """


# --------------------------------
# Post-process input BAM (dedup, optional)
# --------------------------------
rule postprocess_bam_input:
    """Deduplicate and index input-group genome BAM"""
    input:
        log = BAM_DIR / f"{config['prefix']}.input.genome.stats"
    output:
        log = BAM_DIR / f"{config['prefix']}.input.markdup.log"
    params:
        threads = config["threads"],
        basename = f"{config['prefix']}.input",
        sorted_bam = BAM_DIR / f"{config['prefix']}.input.sort.bam"
    container: config["sif"]
    shell:
        """
        if [ ! -e {params.sorted_bam} ]; then
            echo "No input provided, skipping samtools markdup." > {output.log}
        else
            samtools fixmate -r -m {params.sorted_bam} - | \
                samtools sort -@ {params.threads} - | \
                samtools markdup -r -s - {BAM_DIR}/{params.basename}.DeDup.bam 2> {BAM_DIR}/{params.basename}.markdup.log
            samtools flagstat {BAM_DIR}/{params.basename}.DeDup.bam > {BAM_DIR}/{params.basename}.flagstat.txt
            samtools index {BAM_DIR}/{params.basename}.DeDup.bam
            rm {params.sorted_bam}
        fi
        """


# --------------------------
# Generate odd BigWig
# --------------------------
rule bamcoverage_odd:
    """Generate normalized BigWig for odd group"""
    input:
        dedup_bam = BAM_DIR / f"{config['prefix']}.odd.DeDup.bam"
    output:
        bw = BW_DIR / f"{config['prefix']}.odd.DeDup.bw"
    params:
        threads = config["threads"],
        bin_size = config["binSize"]
    container: config["sif"]
    shell:
        """
        bamCoverage -b {input.dedup_bam} -o {output.bw} \
            -p {params.threads} --binSize {params.bin_size} --normalizeUsing CPM
        """


# ---------------------------
# Generate even BigWig
# ---------------------------
rule bamcoverage_even:
    """Generate normalized BigWig for even group"""
    input:
        dedup_bam = BAM_DIR / f"{config['prefix']}.even.DeDup.bam"
    output:
        bw = BW_DIR / f"{config['prefix']}.even.DeDup.bw"
    params:
        threads = config["threads"],
        bin_size = config["binSize"]
    container: config["sif"]
    shell:
        """
        bamCoverage -b {input.dedup_bam} -o {output.bw} \
            -p {params.threads} --binSize {params.bin_size} --normalizeUsing CPM
        """


# --------------------------------
# Generate input BigWig (optional)
# --------------------------------
rule bamcoverage_input:
    """Generate normalized BigWig for input group"""
    input:
        log = BAM_DIR / f"{config['prefix']}.input.markdup.log"
    output:
        log = temp(BW_DIR / f"{config['prefix']}.input.log")
    params:
        threads = config["threads"],
        bin_size = config["binSize"],
        dedup_bam = BAM_DIR / f"{config['prefix']}.input.DeDup.bam",
        basename = f"{config['prefix']}.input",
    container: config["sif"]
    shell:
        """
        if [ ! -e {params.dedup_bam} ]; then
            echo "No input provided, skipping" > {output.log}
        else
            bamCoverage -b {params.dedup_bam} -o {BW_DIR}/{params.basename}.DeDup.bw \
                -p {params.threads} --binSize {params.bin_size} --normalizeUsing CPM 2> {output.log}
        fi
        """


# --------------------------
# Correlation analysis
# --------------------------
rule multi_bw_analysis:
    """PCA and correlation analysis using deepTools"""
    input:
        odd_bw = BW_DIR / f"{config['prefix']}.odd.DeDup.bw",
        even_bw = BW_DIR / f"{config['prefix']}.even.DeDup.bw",
        log = BW_DIR / f"{config['prefix']}.input.log"
    output:
        pca_pdf = FIGURE_DIR / "BW_compare_PCA.pdf",
        cor_pdf = FIGURE_DIR / "BW_compare_cor.pdf"
    params:
        threads = config["threads"],
        prefix = config["prefix"],
        input_bw = BW_DIR / f"{config['prefix']}.input.DeDup.bw"
    container: config["sif"]
    shell:
        """
        # Collect BigWig files
        bw_files="{input.odd_bw} {input.even_bw}"
        [ -e "{params.input_bw}" ] && bw_files+=" {params.input_bw}"
        
        # Collect labels
        labels="{params.prefix}.odd {params.prefix}.even"
        [ -e "{params.input_bw}" ] && labels+=" {params.prefix}.input"
        
        # Run analysis
        multiBigwigSummary bins -b $bw_files --labels $labels \
            -out {FIGURE_DIR}/bw_summary.npz -p {params.threads}
        plotPCA -in {FIGURE_DIR}/bw_summary.npz -o {output.pca_pdf}
        plotCorrelation -in {FIGURE_DIR}/bw_summary.npz --corMethod pearson \
            --skipZeros -o {output.cor_pdf} --whatToPlot heatmap --colorMap RdYlBu --plotNumbers
        rm {FIGURE_DIR}/bw_summary.npz
        """


# --------------------------
# Fingerprint plot
# --------------------------
rule plot_fingerprint:
    """Generate fingerprint plot for library quality"""
    input:
        odd_bam = BAM_DIR / f"{config['prefix']}.odd.DeDup.bam",
        even_bam = BAM_DIR / f"{config['prefix']}.even.DeDup.bam",
        log = BAM_DIR / f"{config['prefix']}.input.markdup.log"
    output:
        pdf = FIGURE_DIR / "fingerprints.pdf"
    params:
        prefix = config["prefix"],
        input_bam = BAM_DIR / f"{config['prefix']}.input.DeDup.bam"
    container: config["sif"]
    shell:
        """
        # Collect BAM files
        bam_files="{input.odd_bam} {input.even_bam}"
        [ -e "{params.input_bam}" ] && bam_files+=" {params.input_bam}"
        
        # Collect labels
        labels="{params.prefix}.odd {params.prefix}.even"
        [ -e "{params.input_bam}" ] && labels+=" {params.prefix}.input"
        
        # Generate plot
        plotFingerprint -b $bam_files --labels $labels \
            --skipZeros --plotFile {output.pdf}
        """


# --------------------------
# Merge odd and even BAMs
# --------------------------
rule merge_bams:
    """Merge odd and even deduplicated BAMs"""
    input:
        odd_bam = BAM_DIR / f"{config['prefix']}.odd.DeDup.bam",
        even_bam = BAM_DIR / f"{config['prefix']}.even.DeDup.bam"
    output:
        merged_bam = BAM_DIR / f"{config['prefix']}.merged.DeDup.bam",
        merged_bai = BAM_DIR / f"{config['prefix']}.merged.DeDup.bam.bai"
    params:
        threads = config["threads"]
    container: config["sif"]
    shell:
        """
        samtools merge -@ {params.threads} {output.merged_bam} {input.odd_bam} {input.even_bam}
        samtools index {output.merged_bam}
        """


# --------------------------
# Generate merged BigWig
# --------------------------
rule bamcoverage_merged:
    """Generate normalized BigWig for merged BAM"""
    input:
        merged_bam = BAM_DIR / f"{config['prefix']}.merged.DeDup.bam"
    output:
        bw = BW_DIR / f"{config['prefix']}.merged.DeDup.bw"
    params:
        threads = config["threads"],
        bin_size = config["binSize"]
    container: config["sif"]
    shell:
        """
        bamCoverage -b {input.merged_bam} -o {output.bw} \
            -p {params.threads} --binSize {params.bin_size} --normalizeUsing CPM
        """


# --------------------------
# Call peaks with MACS3
# --------------------------
rule macs3_callpeak:
    """Call peaks using MACS3"""
    input:
        treat_bam = BAM_DIR / f"{config['prefix']}.merged.DeDup.bam",
        log = BAM_DIR / f"{config['prefix']}.input.markdup.log"
    output:
        log = PEAK_DIR / "peakcalling.log"
    params:
        genome = config["g"],
        peaktype = config["peaktype"],
        pval = config["pval"],
        qval = config["qval"],
        broad_cutoff = config["broad_cutoff"],
        llocal = config["llocal"],
        keepdup = config["keepdup"],
        nomodel = config["nomodel"],
        nolambda = config["nolambda"],
        callsummits = config["callsummits"],
        extsize = config["extsize_val"],
        shift = config["shift_val"],
        bam_format = "BAMPE" if (config["oddFq"].get("R2") or config["evenFq"].get("R2")) else "BAM",
        prefix = config["prefix"],
        ctrl_bam = BAM_DIR / f"{config['prefix']}.input.DeDup.bam"
    container: config["sif"]
    shell:
        """
        # Build MACS3 options
        opts=""
        [ "{params.nomodel}" = "on" ] && opts+=" --nomodel"
        [ "{params.callsummits}" = "on" ] && opts+=" --call-summits"
        [ "{params.nolambda}" = "on" ] && opts+=" --nolambda"
        [ -n "{params.shift}" ] && opts+=" --shift {params.shift}"
        [ -n "{params.extsize}" ] && opts+=" --extsize {params.extsize}"
        
        if [ "{params.peaktype}" = "broad" ]; then
            opts+=" --broad"
            [ -n "{params.broad_cutoff}" ] && opts+=" --broad-cutoff {params.broad_cutoff}"
        fi
        
        [ -n "{params.pval}" ] && opts+=" -p {params.pval}"
        [ -n "{params.qval}" ] && opts+=" -q {params.qval}"
        [ -z "{params.pval}" ] && [ -z "{params.qval}" ] && opts+=" -q 0.01"
        
        # Run peak calling
        if [ -e "{params.ctrl_bam}" ]; then
            macs3 callpeak -t {input.treat_bam} -c {params.ctrl_bam} \
                -f {params.bam_format} -g {params.genome} \
                --outdir {PEAK_DIR} -n "{params.prefix}.merged.vs.{params.prefix}.input" \
                $opts --keep-dup {params.keepdup} 2> {output.log}
        else
            macs3 callpeak -t {input.treat_bam} \
                -f {params.bam_format} -g {params.genome} \
                --outdir {PEAK_DIR} -n "{params.prefix}.merged" \
                $opts --keep-dup {params.keepdup} --llocal {params.llocal} 2> {output.log}
        fi
        """

rule plot_peak_heatmap:
    """plot peak heatmap"""
    input:
        peak_log = PEAK_DIR / "peakcalling.log",
        treat_bw = BW_DIR / f"{config['prefix']}.merged.DeDup.bw",
        log = BW_DIR / f"{config['prefix']}.input.log"
    output:
        heatmap_pdf = FIGURE_DIR / f"{config['prefix']}.peak.pdf"
    params:
        threads = config["threads"],
        peaktype = config["peaktype"],
        region_params = "--regionBodyLength 1000 -a 5000 -b 5000" if config["peaktype"] == "broad" else "--regionBodyLength 0 -a 5000 -b 5000",
        label_params = "--startLabel 'Start' --endLabel 'End'" if config["peaktype"] == "broad" else "--refPointLabel 'Peak'",
        prefix = config["prefix"],
        ctrl_bw = BW_DIR / f"{config['prefix']}.input.DeDup.bw"
    container: config["sif"]
    shell:
        """
        peak_file=$(find {PEAK_DIR} -maxdepth 1 -type f \( -name "*_peaks.broadPeak" -o -name "*_peaks.narrowPeak" \))
        
        if [ -z "$peak_file" ] || [ ! -f "$peak_file" ] || [ ! -s "$peak_file" ]; then
            echo "Warning: Peak file not found or is empty in {PEAK_DIR}. Skipping heatmap generation." >&2
            # Create empty output file to satisfy Snakemake's output requirement
            touch {output.heatmap_pdf}
            exit 0
        fi

        sig_files="{input.treat_bw}"
        [ -e "{params.ctrl_bw}" ] && sig_files+=" {params.ctrl_bw}"

        computeMatrix scale-regions -p {params.threads} \
            -S $sig_files \
            -R $peak_file \
            {params.region_params} \
            --missingDataAsZero \
            -o {FIGURE_DIR}/{params.prefix}.peak_matrix.gz
        
        plotHeatmap -m {FIGURE_DIR}/{params.prefix}.peak_matrix.gz \
            -out {output.heatmap_pdf} \
            --colorMap viridis \
            --missingDataColor white \
            --heatmapHeight 12 \
            --heatmapWidth 4 \
            {params.label_params}
        
        rm {FIGURE_DIR}/{params.prefix}.peak_matrix.gz
        """


rule multiqc_report:
    input:
        QC_DIR / f"{config['prefix']}.odd.fastqc.log",
        QC_DIR / f"{config['prefix']}.even.fastqc.log",
        QC_DIR / f"{config['prefix']}.input.fastqc.log",
        TRIM_DIR / f"{config['prefix']}.odd.trim.log",
        TRIM_DIR / f"{config['prefix']}.even.trim.log",
        TRIM_DIR / f"{config['prefix']}.input.trim.log",
        BAM_DIR / f"{config['prefix']}.odd.repeats.stats",
        BAM_DIR / f"{config['prefix']}.even.repeats.stats",
        BAM_DIR / f"{config['prefix']}.input.repeats.stats",
        BAM_DIR / f"{config['prefix']}.odd.genome.stats",
        BAM_DIR / f"{config['prefix']}.even.genome.stats",
        BAM_DIR / f"{config['prefix']}.input.genome.stats",
        BAM_DIR / f"{config['prefix']}.odd.flagstat.txt",
        BAM_DIR / f"{config['prefix']}.even.flagstat.txt",
        BAM_DIR / f"{config['prefix']}.input.markdup.log"
    output:
        MULTIQC_DIR / "multiqc_report.html"
    container: config["sif"]
    shell:
        """
        multiqc {TRIM_DIR}/* {BAM_DIR}/*.repeats.stats {BAM_DIR}/*.genome.stats \
            {BAM_DIR}/*flagstat.txt -o {MULTIQC_DIR} --force
        """

