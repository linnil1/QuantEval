import os
import re
import logging
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import numpy as np


def setupLogger():
    """setup logger"""
    logger = logging.getLogger('QuantEval')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


def mRNAFilter(basedir, refdir, simdir, shortest_length):
    """
    Filter mRNA from seqences

    Inputs
    ------
    * {refdir}/transcriptome.fasta
    * {refdir}/annotation.gff
    * {refdir}/genome.fasta

    Outputs
    -------
    The reference file that contains only mRNA
        * {basedir}/chromosome.fasta
        * {basedir}/mRNA.gtf
        * {basedir}/mRNA.fasta
    The files for simulation
        * {simdir}/flux_simulator.gtf
        * {simdir}/chromosome/*.fasta

    Parameters
    ----------
    basedir: str
        Path to species directory
    refdir: str
        Path to reference directory
    simdir: str
        Path to simulation directory
    shortest_length: int
        Length threshold of mrna
    """
    logger = logging.getLogger('QuantEval')

    # 1. Get transcriptome name which len > threashold
    logger.info("load transcriptome.fasta")
    transcript_pool = set()  # store all transcriptome name
    transcript_length = 0
    total_transcript_length = 0
    for seq in SeqIO.parse(refdir + "/transcriptome.fasta", "fasta"):
        if len(seq.seq) >= shortest_length:
            transcript_name = seq.id.split(".")[0]
            transcript_pool.add(transcript_name)
        total_transcript_length += len(seq.seq)

    # 1. summary
    print("    - number of transcripts: {:>16d}".format(len(transcript_pool)))
    print("    - total length of transcripts: {:>10d}".format(total_transcript_length))

    # 2. Get mrna from rna and not from mitochordial
    #    and collect chromosome and gene along mrna
    logger.info("extract mRNA in annotation.gff")
    chromosome_pool = set()  # store all chromosome name
    annotation_pool = list()
    mRNA_parent = dict()
    gene_pool = set()
    mRNA_pool = set()
    total_transcript_length = 0

    with open(refdir + "/annotation.gff", "r") as annotation_file:
        for annotation_line in annotation_file:
            if annotation_line[0] == "#":
                continue

            # read each line
            annotation_line = annotation_line.rstrip()
            annotation_data = annotation_line.split("\t")
            chromosome = annotation_data[0]
            category   = annotation_data[2]
            attributes = annotation_data[8]

            # filter
            if chromosome in ["MT", "Mito"]:
                continue
            if category not in ["gene", "mRNA", "exon"]:
                continue
            annotation_pool.append(annotation_data)

            # collect relationship of exon -> mRNA -> transcriptome
            # TOCHECK: I think mRNA will occur before exon data
            if category == "mRNA":
                regex = re.match(r"ID=transcript:(\S+?);", attributes)
                transcript_name = regex.group(1)
                regex = re.search(r"Parent=gene:(\S+?);", attributes)
                gene_name = regex.group(1)
                mRNA_parent[transcript_name] = gene_name

            elif category == "exon":
                regex = re.match(r"Parent=transcript:(\S+?);", attributes)
                transcript_name = regex.group(1)
                if transcript_name in transcript_pool and transcript_name in mRNA_parent:
                    total_transcript_length += int(annotation_data[4]) - int(annotation_data[3]) + 1
                    chromosome_pool.add(chromosome)
                    gene_pool.add(mRNA_parent[transcript_name])
                    mRNA_pool.add(transcript_name)

    # 2. summary
    print("    - number of chromosomes: {:>16d}".format(len(chromosome_pool)))
    print("    - number of genes: {:>22d}".format(len(gene_pool)))
    print("    - number of transcripts: {:>16d}".format(len(mRNA_pool)))
    print("    - total length of transcripts: {:>10d}".format(total_transcript_length))

    # 3. Regenerate the gtf from gff for mRNA and flux_simulator
    #    Select those gene mrna exon had appeared in the pool
    logger.info("output flux_simulator.gtf and mRNA.gtf")
    mRNA_count = 0
    total_transcript_length = 0
    os.makedirs(simdir, exist_ok=True)
    annotation_out = open(basedir + "/mRNA.gtf", "w")
    flux_annotation_out = open(simdir + "/flux_simulator.gtf", "w")
    for annotation_data in annotation_pool:
        category   = annotation_data[2]
        attributes = annotation_data[8]

        if category == "gene":
            regex = re.match(r"ID=gene:(\S+?);", attributes)
            gene_name = regex.group(1)
            if gene_name in gene_pool:
                new_attributes = 'gene_id "' + gene_name + '";'
                print("\t".join(annotation_data[:-1]) + "\t" + new_attributes, file=annotation_out)

        elif category == "mRNA":
            regex = re.match(r"ID=transcript:(\S+?);", attributes)
            transcript_name = regex.group(1)
            regex = re.search(r"Parent=gene:(\S+?);", attributes)
            gene_name = regex.group(1)
            if transcript_name in mRNA_pool and gene_name in gene_pool:
                mRNA_count += 1
                new_attributes = 'gene_id "' + gene_name + '"; transcript_id "' + transcript_name + '";'
                print("\t".join(annotation_data[:-1]) + "\t" + new_attributes, file=annotation_out)

        elif category == "exon":
            regex = re.match(r"Parent=transcript:(\S+?);", attributes)
            transcript_name = regex.group(1)
            if transcript_name in mRNA_pool and transcript_name in mRNA_parent:
                gene_name = mRNA_parent[transcript_name]
                if gene_name not in gene_pool:
                    continue
                total_transcript_length += int(annotation_data[4]) - int(annotation_data[3]) + 1
                new_attributes = 'gene_id "' + gene_name + '"; transcript_id "' + transcript_name + '";'
                print("\t".join(annotation_data[:-1]) + "\t" + new_attributes, file=annotation_out)
                print("\t".join(annotation_data[:-1]) + "\t" + new_attributes, file=flux_annotation_out)

    # 3. summary
    print("    - number of chromosomes: {:>16d}".format(len(chromosome_pool)))
    print("    - number of genes: {:>22d}".format(len(gene_pool)))
    print("    - number of transcripts: {:>16d}".format(mRNA_count))
    print("    - total length of transcripts: {:>10d}".format(total_transcript_length))
    annotation_out.close()
    flux_annotation_out.close()

    # 4. Filter those genome that appeared in the pool
    #    Separate chromosome to /chromosome/*.fasta
    logger.info("output genome.fasta and chromosome.fasta")
    os.makedirs(simdir + "/chromosome", exist_ok=True)
    chromosome_count = 0
    with open(basedir + "/chromosome.fasta", "w") as genome_out:
        for genome in SeqIO.parse(refdir + "/genome.fasta", "fasta"):
            genome_name = genome.id
            if genome_name not in chromosome_pool:
                continue
            SeqIO.write([genome], open(simdir + "/chromosome/" + genome_name + ".fa", "w"), "fasta")
            chromosome_count += 1
            print(genome.format("fasta"), file=genome_out, end="")

    # 4. summary
    print("    - number of chromosomes: {:>16d}".format(chromosome_count))

    # 5. Filter those mrna that appeared in the pool
    logger.info("output mRNA.fasta")
    mRNA_count = 0
    total_transcript_length = 0
    retain_transcript = False
    with open(basedir + "/mRNA.fasta", "w") as transcript_out:
        transcript_sequences = []
        for transcript in SeqIO.parse(refdir + "/transcriptome.fasta", "fasta"):
            transcript_name = transcript.id.split('.')[0]
            if transcript_name in mRNA_pool:
                mRNA_count += 1
                total_transcript_length += len(transcript.seq)
                print(transcript.format("fasta"), file=transcript_out, end="")

    # 5. summary
    print("    - number of transcripts: {:>16d}".format(mRNA_count))
    print("    - total length of transcripts: {:>10d}".format(total_transcript_length))


def splitInterleavedReads(file_fastq, output_name):
    """
    Separte read1 and read2 from simulation fastq

    Parameters
    ----------
    file_fastq: str
        Path to reads data from simulation
    output_name: str
        Path for splited reads files
        It will save to two files:
            {output_name}_r1.fastq and {output_name}_r1.fastq
    """
    # init
    logger = logging.getLogger('QuantEval')
    logger.info(f"Split {file_fastq} to {output_name}")
    fragment_count = 0
    r1 = open(output_name + "_r1.fastq", "w")
    r2 = open(output_name + "_r2.fastq", "w")

    for i, seq in enumerate(SeqIO.parse(file_fastq, "fastq")):
        # rename description
        if i % 2 == 0:
            fragment_count += 1
        desc = seq.id
        seq.description = ""
        seq.id = f"flux_simulator_{fragment_count}:{desc[:-4]}/{desc[-1]}"

        # write
        if desc[-1] == "1":
            print(seq.format("fastq"), file=r1, end="")
        else:
            print(seq.format("fastq"), file=r2, end="")

    # close
    r1.close()
    r2.close()


def fasta_length_filter(file_fasta, file_result, threshold):
    """
    Select sequences that lenght is larger than threshold

    Inputs
    ------
    file_fasta: str
        Path to fasta file
    file_result: str
        Path to filtered fasta file
    threshold: int
        The minimal threashold for transcriptome length
    """
    logger = logging.getLogger('QuantEval')
    logger.info(f"filter {file_fasta} to {file_result} with len > {threshold}")
    # Filter len > threashold
    with open(file_result, "w") as result:
        for seq in SeqIO.parse(file_fasta, "fasta"):
            if len(seq.seq) >= threshold:
                print(seq.format("fasta"), file=result, end="")


def readQuantTsv(method, path):
    """
    Read quantification result tsv file

    Parameters
    ----------
    method: {kallisto, rsem, salmon}
        The quantification tool
    path: str
        The path to tsv file

    Returns
    -------
    data: pandas.Dataframe
        The data frame has three column:
            'name', 'tpm', 'count'
    """

    # read kallisto
    if method == "kallisto":
        data = pd.read_table(path, sep="\t")
        data = data.rename(columns={'target_id': "name",
                                    'tpm': "tpm",
                                    'est_counts': "count"})
        return data.loc[:, ("name", "tpm", "count")]

    # read rsem
    elif method == "rsem":
        data = pd.read_table(path, sep="\t")
        data = data.rename(columns={'transcript_id': "name",
                                    'TPM': "tpm",
                                    'expected_count': "count"})
        return data.loc[:, ("name", "tpm", "count")]

    # read salmon
    elif method == "salmon":
        data = pd.read_table(path, sep="\t")
        data = data.rename(columns={'Name': "name",
                                    'TPM': "tpm",
                                    'NumReads': "count"})
        return data.loc[:, ("name", "tpm", "count")]

    # read answer
    elif method == "answer":
        data = pd.read_table(path, sep="\t")
        data = data.rename(columns={'name': "name",
                                    'answer_tpm': "tpm",
                                    'answer_count': "count"})
        return data.loc[:, ("name", "tpm", "count")]

    else:
        raise ValueError("Unknown Method")


def average_tpm(files_abundance, file_output):
    """
    Average the tmp created by kallisto, rsem and salmon

    Example usage:
        average_tpm({'kallisto': "/data/abundance.tsv"}, "merged.tsv")

    Parameters
    ----------
    files_abundance: dict
        key: quantification method
        value: path to abundance.tsv
    output: str
        The path you store the average value. Format: tsv.
    """
    logger = logging.getLogger('QuantEval')
    logger.info(f"Average {list(files_abundance.keys())} to {file_output}")
    # read all tables
    tables = [readQuantTsv(method, path) for method, path in files_abundance.items()]

    # Average the inner join
    # -> only reserve the transcripts that got by all three methods
    tmp = pd.concat(tables).groupby("name")
    expression_table = tmp.mean().loc[tmp.count()['tpm'] == len(tables), :]

    # save it
    expression_table.reset_index(level=0, inplace=True)
    expression_table = expression_table.rename(columns={'tpm': "answer_tpm",
                                                        'count': "answer_count"})
    expression_table = expression_table.round({'answer_tpm': 3, 'answer_count': 3})
    expression_table.to_csv(file_output, sep="\t", index=False)


def readsCount(files_reads):
    """
    Count the reads in input fastq files

    Parameters
    ----------
    files_read: list of str
        The list of path to the input file

    Returns
    -------
    table_count: pd.Dataframe
        The table that have two colnums: name and counts
    """
    # count sequence ID in fastq
    count = defaultdict(int)
    for f in files_reads:
        for seq in SeqIO.parse(f, 'fastq'):
            name = re.match('\S+?:\S+?:\S+?:(\S+?):', seq.id).group(1)
            count[name] += 1

    # store in table
    read_count = pd.DataFrame.from_dict(count, orient='index')
    read_count.columns = ['count']
    read_count['name'] = read_count.index.tolist()
    return read_count


def readLib(file_lib):
    """
    Collect sequences length distribution in simulation

    Parameters
    ----------
    files_lib: str
        The path to the flux_simulator.pro

    Returns
    -------
    count_dict: dict
        The distribution of length of sequences.
        key: length
        value: counts
    """
    count_dict = defaultdict(int)
    for lib_in in open(file_lib):
        data = lib_in.rstrip().split()
        start, end, count = int(data[0]), int(data[1]), int(data[-1])
        length = abs(start - end)
        count_dict[length] += count

    return count_dict


def estimate_eff_length(original_length, count_dict):
    """
    Calculate effective length of each transcript

    Parameters
    ----------
    original_length: array of int
        Original length of transcripts
    count_dict: dict
        key: length of transcripts
        value: count of transcripts with length = key

    Returns
    -------
    effective_length: array of int
    """
    # Calculate probility
    total_count = float(sum(count_dict.keys()))
    max_frag = max(count_dict.values())
    pdf = np.zeros(max_frag + 1)
    pdf[list(count_dict.keys())] = np.array(list(count_dict.values())) / max_frag

    # define the formula of effective length
    def eff_len(L, pmf):
        c_len = np.arange(0, min(L + 1, len(pmf)))
        effective_length = np.sum(pmf[c_len] * (L - c_len + 1))
        return effective_length

    # apply
    vfunc = np.vectorize(eff_len, excluded=['pmf'])
    return vfunc(L=original_length, pmf=pdf)


def reference_tpm(flux_simulator_pro, flux_simulator_lib, files_unpaired_reads, file_output):
    """
    Calculate the reference TMP by reading flux_simulator paramters

    The answer TPM is (simulation count - unpaired count) / effective length

    Parameters
    ----------
    flux_simulator_pro: str
        The path of flux_simulator.pro
    flux_simulator_lib: str
        The path of flux_simulator.pro
    files_unpaired_reads: list of str
        The path of read files
    file_output: str
        The path to save answer tsv
    """
    logger = logging.getLogger('QuantEval')
    logger.info(f"Calculate real TPM/count from {flux_simulator_pro} and {files_unpaired_reads}")
    # Read from pro
    flux_columns = ['locus', 'name', 'conding', 'length', 'molecular_fraction', 'molecular_count',
                    'fragment_fraction', 'fragment_count', 'read_fraction', 'read_count', 'covered',
                    'chi_square', 'variation']
    flux_simulator = pd.read_table(flux_simulator_pro, sep='\t', header=None, names=flux_columns)

    # Get sequences length distribution from lib
    count_dict = readLib(flux_simulator_lib)

    # Get read counts
    unpaired_read_count = readsCount(files_unpaired_reads)

    # Calculate effective length
    eff_length = estimate_eff_length(flux_simulator["length"].values, count_dict)
    flux_simulator["eff_length"] = eff_length

    # Merge
    flux_simulator = pd.merge(unpaired_read_count, flux_simulator, on='name', how='outer')
    flux_simulator = flux_simulator.fillna(value=0)

    # Calculate TPM, counts
    flux_simulator['answer_count'] = flux_simulator['read_count'] - flux_simulator['count']
    flux_simulator['read_per_nucleotide'] = flux_simulator['answer_count'] / flux_simulator['eff_length']
    flux_simulator['answer_tpm'] = 10 ** 6 * flux_simulator['read_per_nucleotide'] / flux_simulator['read_per_nucleotide'].sum()

    # output and save
    flux_simulator = flux_simulator.loc[:, ('name', 'answer_tpm', 'answer_count')]
    flux_simulator = flux_simulator.round({'answer_tpm': 3})
    flux_simulator.to_csv(file_output, sep='\t', index=False)


if __name__ == "__main__":
    data_folder = "data/yeast"
    """
    # test
    mRNAFilter(data_folder, refdir=f"{data_folder}/download", simdir=f"{data_folder}/simulation", shortest_length=500)
    splitInterleavedReads(f"{data_folder}/simulation/flux_simulator_yeast_low.fastq", f"{data_folder}/simlow")
    fasta_length_filter(f"{data_folder}/rnaspades/transcripts.fasta", "{data_folder}/simlow.rnaspades.fasta", 500)
    average_tpm({'kallisto': f"{data_folder}/explow.mRNA.kallisto.tsv",
                 'rsem': f"{data_folder}/explow.mRNA.rsem.tsv",
                 'salmon': f"{data_folder}/explow.mRNA.salmon.tsv"}, "tmp1.tsv")
    reference_tpm("data/yeast/simulation/flux_simulator_yeast_low.pro",
                  "data/yeast/simulation/flux_simulator_yeast_low.lib",
                  ["data/yeast/trimmed/simlow_r1.unpaired.fastq",
                   "data/yeast/trimmed/simlow_r2.unpaired.fastq"],
                  "tmp.tsv")
    """

    sys.exit(0)
