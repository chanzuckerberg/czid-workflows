version 1.1

workflow index_generation {
    input {
        File input_fasta
        File? input_gtf
        String host_name
        File ercc_fasta
        File ercc_gtf
        String docker_image_id
    }

    call GenerateHostGenome {
        input:
        input_fasta = input_fasta,
        input_gtf = input_gtf,
        host_name = host_name,
        ercc_fasta = ercc_fasta,
        ercc_gtf = ercc_gtf,
        docker_image_id = docker_image_id
    }

    output {
        File original_input_fasta = GenerateHostGenome.original_input_fasta
        File? original_input_gtf = GenerateHostGenome.original_input_gtf
        File fasta_with_ercc_fa = GenerateHostGenome.fasta_with_ercc_fa
        File? gtf_with_ercc_gtf = GenerateHostGenome.gtf_with_ercc_gtf
        File star_genome_tar = GenerateHostGenome.star_genome_tar
        File bowtie_genome_tar = GenerateHostGenome.bowtie_genome_tar
        File minimap2_dna = GenerateHostGenome.minimap2_dna
        File minimap2_rna = GenerateHostGenome.minimap2_rna
    }
}

task GenerateHostGenome {
    input {
        File input_fasta
        File? input_gtf
        String host_name
        File ercc_fasta
        File ercc_gtf
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        #
        # Create fasta_with_ercc
        #

        INPUT_FASTA_PATH="~{input_fasta}"

        # Download input fa
        if [ ${INPUT_FASTA_PATH: -3} == ".gz" ]
        then
            gunzip -c $INPUT_FASTA_PATH > input.fa
            INPUT_FASTA_PATH=input.fa
        else
            cp $INPUT_FASTA_PATH input.fa
            INPUT_FASTA_PATH=input.fa
        fi

        # Concatenate ercc and input
        cat "~{ercc_fasta}" $INPUT_FASTA_PATH > fasta_with_ercc.fa

        #
        # Create gtf_with_ercc
        #

        INPUT_GTF_PATH="~{input_gtf}"
        GTF_PATH="~{ercc_gtf}"

        # Download input gtf, if provided
        if [[ -n "${INPUT_GTF_PATH}" ]] ; then
            if [ ${INPUT_GTF_PATH: -3} == ".gz" ]
            then
                gunzip -c $INPUT_GTF_PATH > input.gtf
                INPUT_GTF_PATH=input.gtf
            else
                cp $INPUT_GTF_PATH input.gtf
                INPUT_GTF_PATH=input.gtf
            fi
            # Concatenate ercc and input
            cat "~{ercc_gtf}" $INPUT_GTF_PATH > gtf_with_ercc.gtf
            GTF_PATH=gtf_with_ercc.gtf
        fi

        #
        # Generate STAR genome
        #

        # Make directory for STAR genome
        STAR_GENOME="~{host_name}_STAR_genome"
        # HACK: we used to support splitting star indexes into many parts, this made things slower
        # Here we generate the index as if it is in many parts, but there is only ever one part for
        # backwards compatibility
        mkdir -p "$STAR_GENOME/part-0"

        AVAILABLE_MEMORY=$(free --bytes | head -n 2 | tail -n 1 | sed "s/  */ /g" | cut -d' ' -f 7)

        STAR \
          --sjdbGTFfile $GTF_PATH \
          --runThreadN $(nproc) \
          --runMode genomeGenerate \
          --genomeFastaFiles fasta_with_ercc.fa \
          --limitGenomeGenerateRAM $AVAILABLE_MEMORY \
          --genomeDir "$STAR_GENOME/part-0"

        # create a parts.txt file for backwards compatibility
        echo 1 > "$STAR_GENOME/parts.txt"

        # tar STAR genome
        tar cvf "$STAR_GENOME.tar" -C $(pwd) $STAR_GENOME

        #
        # Generate bowtie2 genome
        #

        # Make directory for bowtie2 genome
        BOWTIE2_GENOME="~{host_name}_bowtie2_genome"
        mkdir $BOWTIE2_GENOME

        # Change into the directory to contain the output and generate bowtie2 genome
        cd $BOWTIE2_GENOME
        bowtie2-build ../fasta_with_ercc.fa "~{host_name}"
        cd ..

        # tar bowtie2 genome
        tar cvf "$BOWTIE2_GENOME.tar" -C $(pwd) $BOWTIE2_GENOME

        minimap2 -x map-ont -d "~{host_name}_minimap2_genome_dna.mmi" fasta_with_ercc.fa
        minimap2 -x splice -d "~{host_name}_minimap2_genome_rna.mmi" fasta_with_ercc.fa
    >>>

    output {
        File original_input_fasta = "input.fa"
        File? original_input_gtf = "input.gtf"
        File fasta_with_ercc_fa = "fasta_with_ercc.fa"
        File? gtf_with_ercc_gtf = "gtf_with_ercc.gtf"
        File star_genome_tar = "~{host_name}_STAR_genome.tar"
        File bowtie_genome_tar = "~{host_name}_bowtie2_genome.tar"
        File minimap2_dna = "~{host_name}_minimap2_genome_dna.mmi"
        File minimap2_rna = "~{host_name}_minimap2_genome_rna.mmi"
    }

    runtime {
        docker: docker_image_id
    }
}
