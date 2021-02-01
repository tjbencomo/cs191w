import os
import pandas as pd


output_dir = '/scratch/users/tbencomo/cs191w/data/ji'
bam_dir = '/scratch/users/tbencomo/cs191w/data/ji'
metadata_fp = '/home/users/tbencomo/cs191w/data/ji_metadata.csv'
samples = pd.read_csv(metadata_fp)
samples['bamfile'] = bam_dir + '/' + samples['bamfile']
samples['sample_id'] = samples['patient'] + '_' + samples['condition'] + '_' + samples['replicate'].astype(str)

samples = samples.head(1)

def get_bam(wildcards):
    return samples[samples['sample_id'] == wildcards.sample_id]['bamfile']

rule targets:
    input:
        expand(os.path.join(output_dir, "fqs","{sample_id}"), sample_id=samples['sample_id'])

rule bam2fastq:
    input:
        get_bam
        #os.path.join(bam_dir, "{sample_id}.bam") 
    output:
        directory(os.path.join(output_dir, "fqs", "{sample_id}"))
    threads: 8
    shell:
        """
        bamtofastq-1.3.2 --nthreads {threads} {input} {output}
        """

