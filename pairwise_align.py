import subprocess, os
import begin
from pyfaidx import Fasta


@begin.subcommand
def pairwise_align(folder,work_dir):
    try:
        os.mkdir(work_dir)
    except:
        pass
    fastas = [folder+'/'+file for file in os.listdir(folder) if file.endswith('.fasta')]

    for count,fasta in enumerate(fastas):
        subprocess.call("rm temp.fa.fai && reformat.sh in=%s out=temp.fa addunderscore overwrite=true"%fasta,shell=True)
        f = Fasta('temp.fa')
        seqs = f.keys()
        subprocess.call("samtools faidx temp.fa %s > %s/main.fa"%(seqs[0],work_dir),shell=True)
        subprocess.call('samtools faidx temp.fa %s > %s/temp.fa && lastz --format=maf %s/main.fa %s/temp.fa > %s/%d.maf'%(' '.join(seqs[1:]),work_dir,work_dir,work_dir,work_dir,count),shell=True)






@begin.start
def main():
    pass