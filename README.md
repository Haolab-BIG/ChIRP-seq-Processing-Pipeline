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

        # Download repeats from UCSC Table Browser (Repeats --> Genome)
        gunzip $basementdir/repeats.hg38.fa.gz

        # Remove spaces from headers to avoid duplicate names
        awk '/^>/{gsub(/ /,"_"); print; next} {print}' repeats.hg38.fa > repeats.hg38.unique.fa

        # Build Bowtie2 index for repeats
        mkdir -p ${basementdir}/hg38_repeats
        bowtie2-build --threads 8 -f ${basementdir}/repeats.hg38.unique.fa ${basementdir}/hg38_repeats/bowtie2_inde
        ```
      
        <img width="806" height="589" alt="图片" src="https://github.com/user-attachments/assets/bd554377-e831-4532-9911-61773c4ed0fd" />

        * 4.2 genome index
        ```bash
        mkdir basement_data
        cd basement_data

        # Download Genome FASTA
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz

        # Unzip the files
        gunzip GRCh38.primary_assembly.genome.fa.gz
        gunzip gencode.v46.primary_assembly.annotation.gtf.gz

        # Remove scafford
        awk '/^>/ {p=0} /^>chr[0-9XYM]/ {p=1} p' GRCh38.primary_assembly.genome.fa > GRCh38.primary_assembly.genome.chr.fa

        # Build index
        mkdir hg38
        singularity exec --cleanenv End-seq.sif bowtie2-build --threads 8 -f GRCh38.primary_assembly.genome.chr.fa ./hg38/bowtie2_index

        # get chromatin size (only used when you select SEACR as the peak calling method)
        samtools faidx GRCh38.primary_assembly.genome.chr.fa
        cut -f1,2 GRCh38.primary_assembly.genome.chr.fa.fai > chromatin.size

        # Remove unnecessary files
        rm GRCh38.primary_assembly.genome.chr.fa
        rm GRCh38.primary_assembly.genome.fa
        ```

5.   **Required File Structure**
      ```bash
      basement_data/
      ├── hg38/
            ├── bowtie2_index.1.bt2
            ├── bowtie2_index.2.bt2
            ├── bowtie2_index.3.bt2
            ├── bowtie2_index.4.bt2
            ├── bowtie2_index.rev.1.bt2
            └── bowtie2_index.rev.2.bt2
      ├── chromatin.size
      ├── Comparison.txt
      ├── DNAProteinSeq.sif
      ├── illumina_adapter.fa
      ├── run_DNAProteinSeq.sh
      └── SampleInfor.txt
      ```
