import os
import re
from Bio import SeqIO


def mRNAFilter(basedir, refdir, simdir, shortest_length):
    '''
    Filter mRNA from seqences

    Input:
    * {refdir}/transcriptome.fasta
    * {refdir}/annotation.gff
    * {refdir}/genome.fasta

    Output:
    * {basedir}/chromosome.fasta
    * {basedir}/mRNA.gtf
    * {simdir}/flux_simulator.gtf
    * {simdir}/chromosome/*.fasta

    Parameters:
        basedir: str, Path to species directory
        refdir: str, Path to reference directory
        simdir: str, Paht to simulation directory
        shortest_length: int, length threshold
    '''
    # 1. Get transcriptome name which len > threashold
    print('load transcriptome.fasta')
    transcript_pool = set()  # store all transcriptome name
    transcript_length = 0
    total_transcript_length = 0
    for seq in SeqIO.parse(refdir + '/transcriptome.fasta', 'fasta'):
        if len(seq.seq) >= shortest_length:
            transcript_name = seq.id.split('.')[0]
            transcript_pool.add(transcript_name)
        total_transcript_length += len(seq.seq)

    # 1. summary
    print('    - number of transcripts: {:>16d}'.format(len(transcript_pool)))
    print('    - total length of transcripts: {:>10d}'.format(total_transcript_length))

    # 2. Get mrna from rna and not from mitochordial
    #    and collect chromosome and gene along mrna
    print('extract mRNA in annotation.gff')
    chromosome_pool = set()  # store all chromosome name
    annotation_pool = list()
    mRNA_parent = dict()
    gene_pool = set()
    mRNA_pool = set()
    total_transcript_length = 0

    with open(refdir + '/annotation.gff', 'r') as annotation_file:
        for annotation_line in annotation_file:
            if annotation_line[0] == '#':
                continue

            # read each line
            annotation_line = annotation_line.rstrip()
            annotation_data = annotation_line.split('\t')
            chromosome = annotation_data[0]
            category   = annotation_data[2]
            attributes = annotation_data[8]

            # filter
            if chromosome in ['MT', 'Mito']:
                continue
            if category not in ['gene', 'mRNA', 'exon']:
                continue
            annotation_pool.append(annotation_data)

            # collect relationship of exon -> mRNA -> transcriptome
            # TOCHECK: I think mRNA will occur before exon data
            if category == 'mRNA':
                regex = re.match(r'ID=transcript:(\S+?);', attributes)
                transcript_name = regex.group(1)
                regex = re.search(r'Parent=gene:(\S+?);', attributes)
                gene_name = regex.group(1)
                mRNA_parent[transcript_name] = gene_name

            elif category == 'exon':
                regex = re.match(r'Parent=transcript:(\S+?);', attributes)
                transcript_name = regex.group(1)
                if transcript_name in transcript_pool and transcript_name in mRNA_parent:
                    total_transcript_length += int(annotation_data[4]) - int(annotation_data[3]) + 1
                    chromosome_pool.add(chromosome)
                    gene_pool.add(mRNA_parent[transcript_name])
                    mRNA_pool.add(transcript_name)

    # 2. summary
    print('    - number of chromosomes: {:>16d}'.format(len(chromosome_pool)))
    print('    - number of genes: {:>22d}'.format(len(gene_pool)))
    print('    - number of transcripts: {:>16d}'.format(len(mRNA_pool)))
    print('    - total length of transcripts: {:>10d}'.format(total_transcript_length))

    # 3. Regenerate the gtf from gff for mRNA and flux_simulator
    #    Select those gene mrna exon had appeared in the pool
    print('output flux_simulator.gtf and mRNA.gtf')
    mRNA_count = 0
    total_transcript_length = 0
    os.makedirs(simdir, exist_ok=True)
    annotation_out = open(basedir + '/mRNA.gtf', 'w')
    flux_annotation_out = open(simdir + '/flux_simulator.gtf', 'w')
    for annotation_data in annotation_pool:
        category   = annotation_data[2]
        attributes = annotation_data[8]

        if category == 'gene':
            regex = re.match(r'ID=gene:(\S+?);', attributes)
            gene_name = regex.group(1)
            if gene_name in gene_pool:
                new_attributes = 'gene_id "' + gene_name + '";'
                print('\t'.join(annotation_data[:-1]) + '\t' + new_attributes, file=annotation_out)

        elif category == 'mRNA':
            regex = re.match(r'ID=transcript:(\S+?);', attributes)
            transcript_name = regex.group(1)
            regex = re.search(r'Parent=gene:(\S+?);', attributes)
            gene_name = regex.group(1)
            if transcript_name in mRNA_pool and gene_name in gene_pool:
                mRNA_count += 1
                new_attributes = 'gene_id "' + gene_name + '"; transcript_id "' + transcript_name + '";'
                print('\t'.join(annotation_data[:-1]) + '\t' + new_attributes, file=annotation_out)

        elif category == 'exon':
            regex = re.match(r'Parent=transcript:(\S+?);', attributes)
            transcript_name = regex.group(1)
            if transcript_name in mRNA_pool and transcript_name in mRNA_parent:
                gene_name = mRNA_parent[transcript_name]
                if gene_name not in gene_pool:
                    continue
                total_transcript_length += int(annotation_data[4]) - int(annotation_data[3]) + 1
                new_attributes = 'gene_id "' + gene_name + '"; transcript_id "' + transcript_name + '";'
                print('\t'.join(annotation_data[:-1]) + '\t' + new_attributes, file=annotation_out)
                print('\t'.join(annotation_data[:-1]) + '\t' + new_attributes, file=flux_annotation_out)

    # 3. summary
    print('    - number of chromosomes: {:>16d}'.format(len(chromosome_pool)))
    print('    - number of genes: {:>22d}'.format(len(gene_pool)))
    print('    - number of transcripts: {:>16d}'.format(mRNA_count))
    print('    - total length of transcripts: {:>10d}'.format(total_transcript_length))
    annotation_out.close()
    flux_annotation_out.close()

    # 4. Filter those genome that appeared in the pool
    #    Separate chromosome to /chromosome/*.fasta
    print('output genome.fasta and chromosome.fasta')
    os.makedirs(simdir + '/chromosome', exist_ok=True)
    chromosome_count = 0
    with open(basedir + '/chromosome.fasta', 'w') as genome_out:
        for genome in SeqIO.parse(refdir + '/genome.fasta', 'fasta'):
            genome_name = genome.id
            if genome_name not in chromosome_pool:
                continue
            SeqIO.write([genome], open(simdir + '/chromosome/' + genome_name + '.fa', 'w'), 'fasta')
            chromosome_count += 1
            print(genome.format('fasta'), file=genome_out, end="")

    # 4. summary
    print('    - number of chromosomes: {:>16d}'.format(chromosome_count))

    # 5. Filter those mrna that appeared in the pool
    print('output mRNA.fasta')
    mRNA_count = 0
    total_transcript_length = 0
    retain_transcript = False
    with open(basedir + '/mRNA.fasta', 'w') as transcript_out:
        transcript_sequences = []
        for transcript in SeqIO.parse(refdir + '/transcriptome.fasta', 'fasta'):
            transcript_name = transcript.id.split('.')[0]
            if transcript_name in mRNA_pool:
                mRNA_count += 1
                total_transcript_length += len(transcript.seq)
                print(transcript.format('fasta'), file=transcript_out, end="")

    # 5. summary
    print('    - number of transcripts: {:>16d}'.format(mRNA_count))
    print('    - total length of transcripts: {:>10d}'.format(total_transcript_length))


def splitInterleavedReads(file_fastq, output_name):
    """
    Separte read1 and read2 from simulation fastq
    Input: xx.fastq
    Output: {output_name}_r1.fastq, {output_name}_r1.fastq
    """
    # init
    fragment_count = 0
    r1 = open(output_name + "_r1.fastq", "w")
    r2 = open(output_name + "_r2.fastq", "w")
    
    for i, seq in enumerate(SeqIO.parse(file_fastq, 'fastq')):
        # rename description
        if i % 2 == 0:
            fragment_count += 1
        desc = seq.id
        seq.description = ""
        seq.id = f"flux_simulator_{fragment_count}:{desc[:-4]}/{desc[-1]}"

        # write
        if desc[-1] == '1':
            print(seq.format('fastq'), file=r1, end="")
        else:
            print(seq.format('fastq'), file=r2, end="")

    # close
    r1.close()
    r2.close()


def fasta_length_filter(file_fasta, file_result, threshold):
    '''
    Select sequences that lenght is larger than threshold

    Input:
    * file_fasta, str
    * threshold, int

    Output:
    * file_result
    '''

    # 1. Get transcriptome name which len > threashold
    print('load transcriptome.fasta')
    with open(file_result, 'w') as result:
        for seq in SeqIO.parse(file_fasta, 'fasta'):
            if len(seq.seq) >= threshold:
                print(seq.format('fasta'), file=result, end="")
    

if __name__ == '__main__':
    data_folder = 'data/yeast'
    # mRNAFilter(data_folder, refdir=f"{data_folder}/download", simdir=f"{data_folder}/simulation", shortest_length=500)
    # splitInterleavedReads(f"{data_folder}/simulation/flux_simulator_yeast_low.fastq", f"{data_folder}/simlow")
    # fasta_length_filter("data/yeast/rnaspades/transcripts.fasta", "data/yeast/simlow.rnaspades.fasta", 500)
    sys.exit(0)
