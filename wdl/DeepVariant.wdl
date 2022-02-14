# Use DeepVariant to generate VCF & gVCF for one sample
# ref: https://github.com/google/deepvariant/blob/r0.5/docs/deepvariant-gvcf-support.md
#
#              +----------------------------------------------------------------------------+
#              |                                                                            |
#              |  DeepVariant.wdl                                                           |
#              |                                                                            |
#              |  +-----------------+    +-----------------+    +------------------------+  |
# sample.bam   |  |                 |    |                 |    |                        |  |
#  genome.fa ----->  make_examples  |---->  call_variants  |---->  postprocess_variants  |-----> gVCF
#      range   |  |                 |    |                 |    |                        |  |
#              |  +-----------------+    +--------^--------+    +------------------------+  |
#              |                                  |                                         |
#              |                                  |                                         |
#              +----------------------------------|-----------------------------------------+
#                                                 |
#                                        DeepVariant Model

workflow DeepVariant {
    # reference genome
    File ref_fasta_gz

    # Genomic range(s) to call. Provide at most one of range or ranges_bed.
    # If neither is provided, calls the whole reference genome.
    String? range       # e.g. chr12:111766922-111817529
    File? ranges_bed

    # Read alignments - bam & bai (bai auto-generated if omitted)
    # The output vcf/gvcf filename is derived from the bam's.
    File bam
    File? bai

    # DeepVariant model type
    String model_type #WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA

    # Advanced setting for DeepVariant make_examples
    Int? gvcf_gq_binsize

    # DeepVariant docker image tag
    # c.f. https://github.com/google/deepvariant/blob/master/docs/deepvariant-docker.md
    String? deepvariant_docker
    String deepvariant_docker_default = select_first([deepvariant_docker,"gcr.io/deepvariant-docker/deepvariant:1.3.0"])

    # call make_examples { input:
    #     ref_fasta_gz = ref_fasta_gz,
    #     range = range,
    #     ranges_bed = ranges_bed,
    #     bam = bam,
    #     bai = bai,
    #     gvcf_gq_binsize = gvcf_gq_binsize,
    #     deepvariant_docker = deepvariant_docker_default
    # }

    call run_deepvariant { input:
        bam = bam,
        bai = bai,
        ref_fasta_gz = ref_fasta_gz,
        model_type = model_type,
        deepvariant_docker = deepvariant_docker_default
    }

    # call postprocess_variants { input:
    #     ref_fasta_gz = ref_fasta_gz,
    #     call_variants_output = call_variants.call_variants_output,
    #     gvcf_tfrecords_tar = make_examples.gvcf_tfrecords_tar,
    #     deepvariant_docker = deepvariant_docker_default
    # }

    output {
        File vcf_gz = run_deepvariant.call_variants_output
        #File gvcf_gz = postprocess_variants.gvcf_gz
    }
}

# DeepVariant make_examples
task make_examples {
    File ref_fasta_gz

    # bam & bai (bai auto-generated if omitted)
    File bam
    File? bai

    # Genomic range(s) to run on. Provide at most one of range or ranges_bed.
    # If neither is provided, calls the whole reference genome.
    String? range       # e.g. chr12:111766922-111817529
    File? ranges_bed

    Int? gvcf_gq_binsize
    
    String deepvariant_docker

    parameter_meta {
        ref_fasta_gz: "stream"
    }

    command <<<
        range_arg=""
        if [ -n "${range}" ]; then
            range_arg="--regions ${range}"
        fi
        if [ -n "${ranges_bed}" ]; then
            range_arg="--regions ${ranges_bed}"
        fi
        binsize_arg=""
        if [ -n "${gvcf_gq_binsize}" ]; then
            binsize_arg="--gvcf_gq-binsize ${gvcf_gq_binsize}"
        fi

        set -ex -o pipefail
        export SHELL=/bin/bash
        apt-get update -qq && apt-get install -y -qq samtools
        (zcat "${ref_fasta_gz}" > ref.fasta && samtools faidx ref.fasta) &
        output_name=$(basename "${bam}" .bam)
        mv "${bam}" input.bam
        if [ -n "${bai}" ]; then
            mv "${bai}" input.bam.bai
        else
            samtools index "input.bam"
        fi
        wait -n
        ls -lh

        mkdir examples/ gvcf/ logs/
        output_fn="examples/$output_name.tfrecord@$(nproc).gz"
        gvcf_fn="gvcf/$output_name.gvcf.tfrecord@$(nproc).gz"

        seq 0 $(( `nproc` - 1 )) | NO_GCE_CHECK=True parallel --halt 2 -t --results logs/ \
            "/opt/deepvariant/bin/make_examples --mode calling --ref ref.fasta --reads input.bam --examples '$output_fn' --gvcf '$gvcf_fn' --task {} $range_arg $binsize_arg 2>&1" > /dev/null

        mkdir tar/
        tar -cvzf "tar/$output_name.logs.tar.gz" -C logs/ .
        tar -cvf "tar/$output_name.examples.tar" -C examples/ .
        tar -cvf "tar/$output_name.gvcf_tfrecords.tar" -C gvcf/ .
    >>>

    runtime {
        docker: "${deepvariant_docker}"
        cpu: "32"
    }

    output {
        File examples_tar = glob("tar/*.examples.tar")[0]
        File gvcf_tfrecords_tar = glob("tar/*.gvcf_tfrecords.tar")[0]
        File logs_tar_gz = glob("tar/*.logs.tar.gz")[0]
    }
}

# DeepVariant call_variants (CPU)
task run_deepvariant {
    # DeepVariant model files (with no folder component)
    File ref_fasta_gz
    String model_type
    File bam
    File? bai
    String deepvariant_docker

    parameter_meta {
        ref_fasta_gz: "stream"
    }

    command <<<
        set -ex -o pipefail

        apt-get update -qq && apt-get install -y -qq samtools
        (zcat "${ref_fasta_gz}" > /input/ref.fasta && samtools faidx /input/ref.fasta)

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=${model_type} \
            --ref=/input/ref.fasta \
            --reads=${bam} \
            --output_vcf=/output/YOUR_OUTPUT_VCF \
            --output_gvcf=/output/YOUR_OUTPUT_GVCF \
            --num_shards=$(nproc) \
            --logging_dir=/output/logs \
            --dry_run=false

        NO_GCE_CHECK=True /opt/deepvariant/bin/run_deepvariant \
            --outfile "output/$output_name.call_variants.tfrecord.gz" \
            --examples "examples/$output_name.tfrecord@$n_examples.gz" \
            --checkpoint model.ckpt
    >>>

    runtime {
        docker: "${deepvariant_docker}"
        cpu: "32"
    }

    output {
        File call_variants_output = glob("output/*.gz")[0]
    }
}

# DeepVariant postprocess_variants
task postprocess_variants {
    File ref_fasta_gz
    File gvcf_tfrecords_tar
    File call_variants_output
    String deepvariant_docker

    parameter_meta {
        ref_fasta_gz: "stream"
        gvcf_tfrecords_tar: "stream"
    }

    command <<<
        set -ex -o pipefail
        apt-get update -qq && apt-get install -y -qq samtools wget
        # download a multithreaded version of the tabix bgzip utility
        wget --quiet -O bgzip "https://github.com/dnanexus-rnd/GLnexus/blob/master/cli/dxapplet/resources/usr/local/bin/bgzip?raw=true"
        chmod +x bgzip

        zcat "${ref_fasta_gz}" > ref.fasta
        samtools faidx ref.fasta

        mkdir gvcf output
        tar xvf "${gvcf_tfrecords_tar}" -C gvcf/
        n_gvcf_tfrecords=$(find gvcf/ -type f | wc -l)
        output_name=$(basename "${call_variants_output}" .call_variants.tfrecord.gz)
        NO_GCE_CHECK=True /opt/deepvariant/bin/postprocess_variants \
            --ref ref.fasta --infile "${call_variants_output}" \
            --nonvariant_site_tfrecord_path "gvcf/$output_name.gvcf.tfrecord@$n_gvcf_tfrecords.gz" \
            --outfile "output/$output_name.vcf" \
            --gvcf_outfile "output/$output_name.gvcf"

        ./bgzip -@ $(nproc) output/*.vcf &
        ./bgzip -@ $(nproc) output/*.gvcf
        wait -n
    >>>

    runtime {
        docker: "${deepvariant_docker}"
        cpu: "8"
    }

    output {
        File vcf_gz = glob("output/*.vcf.gz")[0]
        File gvcf_gz = glob("output/*.gvcf.gz")[0]
    }
}
