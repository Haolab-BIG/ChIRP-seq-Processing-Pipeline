# ChIRP-seq-Processing-Pipeline
This pipeline processes raw ChIRP-seq FASTQ data through sequential steps including adapter trimming, quality control, genome alignment, and peak calling. The workflow supports three input types: odd, even, and optional input control. Repeat-masked sequences are removed prior to alignment. Odd and even datasets are merged and then compared to the input control if available.

# Part I Introduction
## i. Workflow
Here stands an throughout workflow of data analysis.
<img width="1016" height="320" alt="ChIRPseq" src="https://github.com/user-attachments/assets/bcbb4ede-e6a4-485d-aacf-1e024876711f" />

## ii. Features
The workflow is fully containerized using Singularity, including all required tools and dependencies. With a single command, raw FASTQ data can be processed through trimming, quality control, alignment, merging, repeat removal, and peak calling in a reproducible manner. Input control is optional; if not provided, peak calling is performed on the merged ChIRP-seq dataset alone.

# Part II Requirements
1.  **Recommended System Configuration**:

      * 8-core CPU
      * 24 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
			libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3.  **Download Basement Files**:

      * `run_ChIRPseq.sh`
      * `ChIRPseq.sif` (The Singularity container)
      * `illumina_adapter.fa`

4.  **Reference Data**: A directory containing Bowtie2 indexes. Below are detailed steps for the human **hg38 genome**. For other reference genomes, download the corresponding files and replace paths as needed.

      * 4.1 Repeat Index

        ```bash
        # Set base directory
        basementdir=/mnt1/2.NAS2024/wutan/9.pipe/3.chrip-seq/basement_data
        cd $basementdir

        # Download repeats from UCSC Table Browser and move to $basementdir (as shown in followed picture)
        singularity exec --cleanenv ChIRPseq.sif gunzip $basementdir/repeats.hg38.fa.gz

        # Remove spaces from headers to avoid duplicate names
        singularity exec --cleanenv ChIRPseq.sif awk '/^>/{gsub(/ /,"_"); print; next} {print}' repeats.hg38.fa > repeats.hg38.unique.fa
        singularity exec --cleanenv ChIRPseq.sif rm repeats.hg38.fa
        
        # Build Bowtie2 index for repeats
        singularity exec --cleanenv ChIRPseq.sif mkdir -p ${basementdir}/hg38_repeats
        singularity exec --cleanenv ChIRPseq.sif bowtie2-build --threads 8 -f ${basementdir}/repeats.hg38.unique.fa ${basementdir}/hg38_repeats/bowtie2_index
        ```
      
        <img width="806" height="589" alt="图片" src="https://github.com/user-attachments/assets/bd554377-e831-4532-9911-61773c4ed0fd" />

        * 4.2 genome index
        ```bash
        mkdir basement_data
        cd basement_data

        # Download Genome FASTA
        singularity exec --cleanenv ChIRPseq.sif wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz

        # Unzip the files
        singularity exec --cleanenv ChIRPseq.sif gunzip GRCh38.primary_assembly.genome.fa.gz
        singularity exec --cleanenv ChIRPseq.sif gunzip gencode.v46.primary_assembly.annotation.gtf.gz

        # Remove scafford
        singularity exec --cleanenv ChIRPseq.sif awk '/^>/ {p=0} /^>chr[0-9XYM]/ {p=1} p' GRCh38.primary_assembly.genome.fa > GRCh38.primary_assembly.genome.chr.fa

        # Build index
        singularity exec --cleanenv ChIRPseq.sif mkdir hg38
        singularity exec --cleanenv ChIRPseq.sif bowtie2-build --threads 8 -f GRCh38.primary_assembly.genome.chr.fa ./hg38/bowtie2_index

        # get chromatin size (only used when you select SEACR as the peak calling method)
        singularity exec --cleanenv ChIRPseq.sif samtools faidx GRCh38.primary_assembly.genome.chr.fa
        singularity exec --cleanenv ChIRPseq.sif cut -f1,2 GRCh38.primary_assembly.genome.chr.fa.fai > chromatin.size

        # Remove unnecessary files
        singularity exec --cleanenv ChIRPseq.sif rm GRCh38.primary_assembly.genome.chr.fa
        singularity exec --cleanenv ChIRPseq.sif rm GRCh38.primary_assembly.genome.fa
        ```

5.   **Required File Structure**
      ```bash
      basement_data/
      ├── ChIRPseq.sif
      ├── hg38/
            ├── bowtie2_index.1.bt2
            ├── bowtie2_index.2.bt2
            ├── bowtie2_index.3.bt2
            ├── bowtie2_index.4.bt2
            ├── bowtie2_index.rev.1.bt2
            └── bowtie2_index.rev.2.bt2
      ├── hg38_repeats/
            ├── bowtie2_index.1.bt2
            ├── bowtie2_index.2.bt2
            ├── bowtie2_index.3.bt2
            ├── bowtie2_index.4.bt2
            ├── bowtie2_index.rev.1.bt2
            └── bowtie2_index.rev.2.bt2
      ├── illumina_adapter.fa
      ├── run_ChIRPseq.sh
      └── repeats.hg38.unique.fa
      ```

# Part III Running

   * **Example code for ChIRP-seq**

      ```bash
      bash ./basement_data/run_ChIRPseq.sh \
                --oddFq /path_to_rawdata/HeLa_Terc_odd.fastq.gz \
                --evenFq /path_to_rawdata/HeLa_Terc_even.fastq.gz \
                --inputFq /path_to_rawdata/HeLa_Terc_input.fastq.gz \
                --prefix HeLa_Terc \
                --outputdir ./result \
                --repeatIndex ./basement_data/hg38_repeats/bowtie2_index \
                --genomeIndex ./basement_data/hg38/bowtie2_index \
                --adapterFa ./basement_data/illumina_adapter.fa \
                --sif ./basement_data/DNAProteinSeq.sif \
                --threads 8 \
                --binSize 10 --g hs
      ```
      
   * **Command Parameters**

      - `--evenFq`:             (required) The path(s) to the even FASTQ file(s)
      - `--oddFq`:              (required) The path(s) to the odd FASTQ file(s)
      - `--inputFq`:            (optinal) The path(s) to the input FASTQ file(s)
      - `--prefix`:             (optinal) A prefix name of output files()
      - `--outputdir`:          (required) Path to the directory where the output will be stored
      - `--repeatIndex`:        (required) Path to the directory where bowtie reference build with prefix
      - `--genomeIndex`:        (required) Path to the directory where bowtie reference build with prefix
      - `--adapterFa`:          (required) Path to the adapter fasta
      - `--sif`:                (required) Path to the singularity environment file
      - `--threads`:            (optional) Number of threads to use (default: 8)
      - `--binSize`:            (optional) Number of binsize to use (default: 10)
      - `--peakcalling`:        (optional) Number of binsize to use (default: 10)
      - `--g`:                  (optional) provide the species code accepted by MACS3, for example: `hs` (human), `mm` (mouse), `ce` (C. elegans), `dm` (Drosophila melanogaster), etc.
      - `--peaktype`:           (optional) For MACS3, set `broad` or `narrow` (default: `narrow`)
      - `--pval`:               (optional) For MACS3, P-value cutoff for peak calling to determine significance of narrow/strong peaks. If specified, MACS3 uses p-value instead of q-value
      - `--qval`:               (optional) For MACS3, Q-value (FDR) cutoff for narrow/strong peaks, controlling false discovery rate (default: 0.01)
      - `--broad_cutoff`:       (optional) For MACS3, P-value or q-value cutoff for broad/weak peaks, effective only when `--peaktype` is `broad` (default: 0.1)
      - `--llocal`:             (optional) For MACS3, `--llocal` value for samples without input control (default: 100000)
      - `--keepdup`:            (optional) For MACS3, `--keep-dup` setting (default: `all`)
      - `--nomodel`:            (optional) For MACS3, Disable model building, set `on` or `off` (default: `on`)
      - `--nolambda`:           (optional) For MACS3, Disable dynamic lambda, set `on` or `off` (default: `on`)
      - `--callsummits`:        (optional) For MACS3, Enable calling of peak summits within enriched regions, set `on` or `off` (default: `on`)
      - `--extsize_val`:        (optional) For MACS3, Set fragment/extension size for MACS3 in base pairs, effective only when `--nomodel` is `on`
      - `--shift_val`:          (optional) For MACS3, Set shift for read 5' ends in base pairs for MACS3, effective only when `--nomodel` is `on`; positive moves 5'->3', negative moves 3'->5'

# Part IV Output

   * **Output Structure with MACS3**
      ```bash
      result_chip_histone /
      ├── bam/
            ├── HeLa_Terc.even.bowtie.stats
            ├── HeLa_Terc.even.DeDup.bam
            ├── HeLa_Terc.even.DeDup.bam.bai
            ├── HeLa_Terc.even.flagstat.txt
            ├── HeLa_Terc.even.markdup.log
            ├── HeLa_Terc.even.repeats.bowtie.stats
            ├── HeLa_Terc.even_repeat.bam
            ├── HeLa_Terc.input.bowtie.stats
            ├── HeLa_Terc.input.DeDup.bam
            ├── HeLa_Terc.input.DeDup.bam.bai
            ├── HeLa_Terc.input.flagstat.txt
            ├── HeLa_Terc.input.markdup.log
            ├── HeLa_Terc.input.repeats.bowtie.stats
            ├── HeLa_Terc.input_repeat.bam
            ├── HeLa_Terc.odd.bowtie.stats
            ├── HeLa_Terc.odd.DeDup.bam
            ├── HeLa_Terc.odd.DeDup.bam.bai
            ├── HeLa_Terc.odd.flagstat.txt
            ├── HeLa_Terc.odd.markdup.log
            ├── HeLa_Terc.odd.repeats.bowtie.stats
            └── HeLa_Terc.odd_repeat.bam
      ├── bw/
            ├── HeLa_Terc.even.DeDup.bw
            ├── HeLa_Terc.input.DeDup.bw
            ├── HeLa_Terc.merged.DeDup.bw
            └── HeLa_Terc.odd.DeDup.bw
      ├── figure/
            ├── HeLa_Terc.peak.pdf
            ├── fingerprints.pdf
            ├── BW_compare_PCA.pdf
            └── BW_compare_cor.pdf
      ├── peak/
            ├── HeLa_Terc.merged.vs.HeLa_Terc.input.macs3.stats
            ├── HeLa_Terc.merged.vs.HeLa_Terc.input_peaks.narrowPeak
            ├── HeLa_Terc.merged.vs.HeLa_Terc.input_peaks.xls
            └── HeLa_Terc.merged.vs.HeLa_Terc.input_summits.bed
      ├── multiqc/
            ├── multiqc_data/
            └── multiqc_report.html
      ```

   * **Output Interpretation**

      - **`*.repeats.bowtie.stats`**
        - **Content**: Bowtie alignment summary when reads are mapped **to the repeat reference database**. Reports total reads processed, reads aligned to repeats (unique/multiple), unaligned reads, and overall alignment rate.

        - **Application**:  
          - Evaluates the **repeat-removal step**.  
          - **Desired outcome**: a **higher percentage of unmapped reads** indicates that most sequencing reads are *not* repeats and can proceed to genome alignment.  
          - Useful for checking repeat content and filtering efficiency; results were aggregated with **MultiQC**.

          <img width="1589" height="516" alt="图片" src="https://github.com/user-attachments/assets/654b8289-0b5a-4e10-a544-e0063dc40e91" />

      - **`*.bowtie.stats`**

        - **Content**: Bowtie alignment summary when reads are mapped **to the reference genome**. Includes total reads processed, uniquely mapped reads, multi-mapped reads, unaligned reads, and overall alignment rate.
        - **Application**:
          - Assesses **final genome alignment quality**.  
          - **Desired outcome**: a **high proportion of uniquely mapped reads** reflects good library quality and accurate genome mapping.  
          - Helps identify low-quality libraries or alignment issues and were visualized with **MultiQC**.
          
          <img width="1589" height="508" alt="图片" src="https://github.com/user-attachments/assets/0120f2e6-1560-4035-9682-372ba1600d60" />

      - **`*.DeDup.bam`**

        - **Content**: This is the main alignment file in Binary Alignment Map (BAM) format. It contains all the sequencing reads and their mapping coordinates on the reference genome. This version has had duplicate reads (PCR duplicates) removed. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bam.html.
        - **Application**: It's the primary evidence for read alignment and can be used for detailed inspection in genome browsers or for downstream analyses.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.bw`**

        - **Content**: A BigWig file that represents the End-seq signal coverage across the genome. It shows the read density (how many reads cover each position) in a compressed format. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bigWig.html
        - **Application**: Primarily used for visualization. You can load this file into a genome browser (e.g., IGV, UCSC Genome Browser) to see a "signal track" that shows gene expression levels visually across chromosomes.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.flagstat.txt`**

        - **Content**: Contains alignment statistics generated by samtools flagstat after removing duplicate reads. It reports the total number of reads, mapped reads, properly paired reads, singletons, and the number of duplicate reads removed, providing a summary of the final, deduplicated BAM file.
        - **Application**: Used to evaluate the quality of the deduplicated alignment, check library complexity, and ensure that downstream analyses (e.g., peak calling, coverage calculation) are based on high-quality, non-redundant reads. Results were aggregated with **MultiQC**.

          <img width="1590" height="870" alt="图片" src="https://github.com/user-attachments/assets/fc858aab-73d6-4af1-9f81-8a61ef3543b0" />

		  <img width="1593" height="871" alt="图片" src="https://github.com/user-attachments/assets/1fa28e0d-5775-49f3-8413-042a5782afcc" />

      - **`*.markdup.log`**

        - **Content**: Log file generated by `samtools markdup`, summarizing read duplication. It includes READ (total number of input reads), WRITTEN (reads retained after removing duplicates), EXCLUDED, EXAMINED, counts of PAIRED and SINGLE reads, as well as DUPLICATE SINGLE/PAIR and DUPLICATE TOTAL.
        - **Application**: Used to evaluate library complexity and duplication rate. A high WRITTEN/READ ratio indicates low duplication and good library complexity, while a low ratio suggests high PCR duplication or low-complexity sequencing.

      - **`multiqc_report`** : Open multiqc_report.html in a web browser to explore all sections interactively.

        - **General Statistics**: A combined table summarizing important metrics for each sample:
	  
          <img width="1592" height="512" alt="图片" src="https://github.com/user-attachments/assets/37078704-dd2a-4dc9-bb6a-7dc1d6dd2a58" />

        - **FastQC**: Quality-control metrics on raw and trimmed reads, including 'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores', 'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content', 'Sequence Length Distribution', 'Sequence Duplication Levels', 'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content':
    
          - Sequence counts for each sample. Estimate duplicate read counts:

            <img width="1595" height="557" alt="图片" src="https://github.com/user-attachments/assets/594acb6e-313c-4661-909b-bea2624e133b" />

          - Sequence Quality Histograms: The mean quality value across each base position in the read.
	  
            <img width="1591" height="619" alt="图片" src="https://github.com/user-attachments/assets/f21ef97e-74c1-45f5-8453-33756c3a6e43" />

          - Adapter Content: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.
	  
            <img width="1586" height="613" alt="图片" src="https://github.com/user-attachments/assets/7a660c29-6dd4-4b60-ae01-6ecbbd8c0ede" />

        - **Samtools**: This module parses the output from samtools flagstat to report the percentage of total, mapped, and properly paired reads, providing a summary of alignment quality. Helps evaluate the effectiveness of deduplication and ensures that downstream analyses (e.g., peak calling, coverage profiling) are based on unique, non-redundant reads. (as shown in **`*.flagstat.txt`**)
	  
        - **Bowtie**: Alignment statistics such as total reads, uniquely mapped reads, and multi-mapping rates (as shown in **`*.bowtie.stats`** and **`*.repeats.bowtie.stats`**)

      - **`*BW_compare_PCA.pdf`**

        - **Content**: PDF file showing the principal component analysis (PCA) of BigWig signal profiles across multiple samples. It visualizes sample-to-sample similarity and variance based on genome-wide coverage or signal intensities.
        - **Application**: Used to assess the overall relationship between samples, detect outliers, and evaluate batch effects or experimental reproducibility.

          <img width="707" height="705" alt="图片" src="https://github.com/user-attachments/assets/e1ed3bc8-d562-4dfd-bf8b-43d64fca4971" />

      - **`*BW_compare_cor.pdf`**

        - **Content**: PDF file showing a heatmap of pairwise correlations between samples based on BigWig signal profiles. It typically includes correlation values (Pearson) and visually represents sample similarity across the genome.
        - **Application**: Used to assess consistency and reproducibility between samples, identify outliers, and evaluate experimental quality in End-seq.

          <img width="783" height="787" alt="图片" src="https://github.com/user-attachments/assets/ce07afc6-544f-4573-a5cc-141023185cd1" />

      - **`*fingerprints.pdf`**

        - **Content**: PDF file generated by `plotFingerprint` showing the cumulative read coverage across the genome for each BAM file. It visualizes enrichment patterns and sequencing depth consistency among samples.
        - **Application**: In the fingerprint plot, a larger separation between treatment and control curves, indicates stronger enrichment and higher signal-to-noise ratio. Also, the fingerprint plot can help decide whether to call narrow peaks or broad peaks: Narrow peaks are appropriate when the signal is sharp and localized; Broad peaks are used when the signal spans wide genomic regions with diffuse enrichment, such as histone modifications.

          <img width="923" height="708" alt="图片" src="https://github.com/user-attachments/assets/5961fc2d-de94-4433-b669-b69fe109a651" />

      - **`*..macs3.stats`**

        - **Content**: Contains summary statistics from MACS3 peak calling, including number of input reads, effective genome size, estimated fragment size, number of peaks called, and other runtime information.
        - **Application**: Used to check if MACS3 ran successfully and to detect any errors or warnings during the peak calling process.

      - **`*_peaks.narrowPeak`**

        - **Content**: BED6+4 format file containing peak locations along with peak summit, p-value, and q-value.  
    Suitable for direct loading into the UCSC Genome Browser when the `--trackline` option is enabled.

        - **Columns**:

        | Column | Description |
        |--------|------------|
        | chrom  | Chromosome name |
        | start  | Start position of the peak (0-based) |
        | end    | End position of the peak (not inclusive) |
        | name   | Peak name or ID |
        | score  | Integer score for display, calculated as `int(-10*log10(pvalue))` or `int(-10*log10(qvalue))` depending on whether `-p` or `-q` was used as the cutoff. |
        | strand | Strand information (‘+’, ‘-’, or ‘.’ if not applicable) |
        | signalValue | Fold enrichment at the peak summit |
        | pValue | -log10 p-value at the peak summit |
        | qValue | -log10 q-value (FDR) at the peak summit |
        | summit | Relative summit position to the peak start |

        - **Application**: Used for quantitative peak analysis, filtering peaks by significance, fold enrichment, or integrating with downstream functional annotation, motif analysis, and visualization.

      - **`*_peaks.broadPeak`**

        - **Content**: BED6+3 format file (similar to narrowPeak, but without the 10th column for peak summits). Only available when `--broad` is enabled. In broad peak mode, the peak summit isn’t called, so the 5th, 7th–9th columns are the mean values across the peak region. Can be loaded directly into UCSC Genome Browser with `--trackline`.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chrom  | Chromosome name |
        | start  | Start position of the broad peak (0-based) |
        | end    | End position of the broad peak (not inclusive) |
        | name   | Peak name or ID |
        | score  | Mean score across the broad peak (similar to narrowPeak 5th column) |
        | strand | Strand information (‘+’, ‘-’, or ‘.’ if not applicable) |
        | signalValue | Mean enrichment signal across the peak |
        | pValue | Mean -log10 p-value across the peak |
        | qValue | Mean -log10 q-value (FDR) across the peak |

        - **Application**: Used for quantitative peak analysis, filtering peaks by significance, fold enrichment, or integrating with downstream functional annotation, motif analysis, and visualization.

      - **`*_peaks.gappedPeak`**

        - **Content**: BED12+3 format file containing broad regions and narrow peaks within them. Only available when `--broad` is enabled. Can be loaded into UCSC Genome Browser. Columns 5, 7–9 may need adjustment if integrating with narrowPeak conventions.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chrom       | Chromosome name |
        | start       | Start of the broad region (0-based) |
        | end         | End of the broad region (not inclusive) |
        | name        | Peak name or ID |
        | score       | Score for display in UCSC browser (grey levels, similar to narrowPeak 5th column) |
        | strand      | Strand information (‘+’, ‘-’, or ‘.’) |
        | thickStart  | Start of the first narrow peak within the broad region |
        | thickEnd    | End of the first narrow peak within the broad region |
        | itemRgb     | RGB color for UCSC browser (0 uses default color) |
        | blockCount  | Number of blocks (including 1bp at start and end of broad regions) |
        | blockSizes  | Comma-separated lengths of each block |
        | blockStarts | Comma-separated start positions of each block relative to `start` |
        | foldChange  | Fold-change of enrichment within the peak |
        | -log10(pvalue) | -log10 p-value for the peak |
        | -log10(qvalue) | -log10 q-value (FDR) for the peak |

        - **Application**: Used to analyze subpeak structure, study internal peak features, or visualize complex enrichment patterns in broad regions.

      - **`*_peaks.xls`**

        - **Content**: Tab-delimited summary of all peaks called by MACS3 with detailed metrics.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chr    | Chromosome name |
        | start  | Peak start position (0-based) |
        | end    | Peak end position (not inclusive) |
        | length | Peak length (end - start) |
        | pileup | Maximum pileup (number of overlapping tags) at the peak |
        | -log10(pvalue) | -log10 of p-value for peak significance |
        | fold_enrichment | Fold enrichment of the peak over background |
        | -log10(qvalue) | -log10 of q-value (FDR) for peak significance |
        | name   | Peak name or ID |

        - **Application**: Used for quantitative peak analysis, filtering peaks by significance, fold enrichment, or integrating with downstream functional annotation, motif analysis, and visualization.

      - **`*.peak.pdf`**

        - **Content**: Heatmap visualizing read enrichment over peaks. Generated using `plotHeatmap` from deepTools with a `*_peaks.broadPeak` and `*bw` input. The heatmap shows signal intensity (color-coded, viridis colormap) across all peaks, with missing data represented in white. The height and width of the heatmap are set for clear visualization of peak patterns.
        - **Application**: Used to assess global enrichment patterns across peaks. Peaks with strong enrichment appear as high-intensity bands or curves; if the signal is higher than control samples, it indicates that peak calling was successful and represents true biological enrichment.

          <img width="491" height="772" alt="图片" src="https://github.com/user-attachments/assets/af523f27-1758-43d1-8f0d-ee5f7924a6a3" />

      - **`*_peaks.broadPeak` / `*_peaks.narrowPeak`**

        - **Content**: Genomic interval files produced by peak-calling pipelines:

        - **Application**:  
          Each of these files can be loaded into **IGV** or the **UCSC Genome Browser** together with the matched **bigWig track file** for direct, interactive genome-wide visualization of peak regions and signal intensity.

          <img width="733" height="378" alt="图片" src="https://github.com/user-attachments/assets/964ebf87-53b1-49b7-8568-eb9e8b3c6cc0" />

		  
