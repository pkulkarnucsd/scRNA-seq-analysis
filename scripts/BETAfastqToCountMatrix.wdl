version 1.0

# In the works Pipeline that converts fastq files into count matrix that can later be used for DEG Analysis

workflow scRNAseqPipeline {
    input {
        Array[File] fastq_files
    }

    # FastQC before trimming
    scatter (file in fastq_files) {
        call FastQC {
            input:
                input_fastq = file
        }
    }

    # Trim using Fastp
    scatter (file in fastq_files) {
        call Fastp {
            input:
                input_fastq = file
        }
    }

    # FastQC after trimming
    scatter (file in Fastp.trimmed_files) {
        call FastQC {
            input:
                input_fastq = file
        }
    }

    # HISAT2 Alignment (considering Human Genome)
    call DownloadReference {
        input:
            url = "https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz",
            output_name = "grch38_genome"
    }

    scatter (pair in Fastp.trimmed_pairs) {
        call HISAT2 {
            input:
                index_prefix = DownloadReference.downloaded,
                read1 = pair[0],
                read2 = pair[1]
        }
    }

    # Download GTF
    call DownloadGTF {
        input:
            url = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz",
            output_name = "Homo_sapiens.GRCh38.112.gtf.gz"
    }

    # FeatureCounts
    scatter (bam in HISAT2.bam_files) {
        call FeatureCounts {
            input:
                bam_file = bam,
                annotation_file = DownloadGTF.downloaded
        }
    }

    output {
        Array[File] fastqc_reports = FastQC.fastqc_reports
        Array[File] trimmed_fastqs = Fastp.trimmed_files
        Array[File] aligned_bams = HISAT2.bam_files
        Array[File] count_matrices = FeatureCounts.count_matrices
    }
}

task FastQC {
    input {
        File input_fastq
    }
    command {
        fastqc ~{input_fastq}
    }
    output {
        File fastqc_report = glob("*.html")
    }
}

task Fastp {
    input {
        File input_fastq
    }
    command {
        fastp -i ~{input_fastq} -o trimmed_~{basename(input_fastq)}
    }
    output {
        File trimmed_fastq = "trimmed_~{basename(input_fastq)}"
    }
}

task DownloadReference {
    input {
        String url
        String output_name
    }
    command {
        wget ~{url} -O ~{output_name}.tar.gz
        tar -xvf ~{output_name}.tar.gz
    }
    output {
        String downloaded = "~{output_name}"
    }
}

task HISAT2 {
    input {
        String index_prefix
        File read1
        File read2
    }
    command {
        hisat2 -x ~{index_prefix} -1 ~{read1} -2 ~{read2} -S aligned.sam
        samtools view -bS aligned.sam > aligned.bam
        samtools sort -o sorted.bam aligned.bam
        samtools index sorted.bam
    }
    output {
        File bam_file = "sorted.bam"
    }
}

task DownloadGTF {
    input {
        String url
        String output_name
    }
    command {
        wget ~{url} -O ~{output_name}
        gunzip ~{output_name}
    }
    output {
        File downloaded = "~{output_name%.gz}"
    }
}

task FeatureCounts {
    input {
        File bam_file
        File annotation_file
    }
    command {
        featureCounts -a ~{annotation_file} -o counts.txt ~{bam_file}
    }
    output {
        File count_matrix = "counts.txt"
    }
}
