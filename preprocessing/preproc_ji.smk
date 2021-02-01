import os
import pandas as pd


output_dir = '/scratch/users/tbencomo/cs191w'
bam_dir = '/scratch/users/tbencomo/cs191w/data/ji'
metadata_fp = '/home/users/tbencomo/cs191w/data/ji_metadata.csv'
samples = pd.read_csv(metadata_fp)
samples['bamfile'] = bam_dir + '/' + samples['bamfile']
samples['sample_id'] = samples['patient'] + '_' + samples['condition'] + '_' + samples['replicate'].astype(str)

ref_dir = '/home/groups/carilee/refs/scRNAseq/refdata-gex-GRCh38-2020-A'

cellranger_mem = 46
cellranger_threads = 20


#samples = samples.head(1)

def get_bam(wildcards):
    return samples[samples['sample_id'] == wildcards.sample_id]['bamfile']

def get_input(wildcards):
    sample_id = wildcards.sample_id
    fqdir = os.listdir(os.path.join(bam_dir, "fqs", sample_id))
    fqdir = [os.path.join(bam_dir, "fqs", sample_id, f) for f in fqdir]
    nested_dir = [d for d in fqdir if os.path.isdir(d)][0]
    fqs = os.listdir(os.path.join(bam_dir, sample_id, nested_dir))
    fqs = [os.path.join(nested_dir, f) for f in fqs]
    print(fqs)
    return fqs  + [ref_dir]

rule targets:
    input:
        expand(os.path.join(bam_dir, "fqs","{sample_id}"), sample_id=samples['sample_id']),
        expand(os.path.join(output_dir, "cellranger", "ji", "{sample_id}"), sample_id = samples['sample_id'])

rule cellranger_count:
    input:
        os.path.join(bam_dir, "fqs", "{sample_id}")
        #get_input
    output:
        outdir=directory(os.path.join(output_dir, "cellranger" , "ji", "{sample_id}"))
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.outdir),
        #fqdir = lambda wildcards, input: os.path.dirname(input[0]),
        refdir = ref_dir,
        mem = cellranger_mem,
        nthreads = cellranger_threads
    shell:
        """
        cd {params.outdir}
        cellranger count --id={wildcards.sample_id} \
            --fastqs={input} \
            --transcriptome={params.refdir} \
            --localcores={params.nthreads} \
            --localmem={params.mem}
        """

rule bam2fastq:
    input:
        get_bam
    output:
        directory(os.path.join(bam_dir, "fqs", "{sample_id}"))
    threads: 8
    shell:
        """
        bamtofastq-1.3.2 --nthreads {threads} {input} {output}
        """

