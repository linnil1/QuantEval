import pandas as pd
import numpy as np
import re
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from scripts.Match import Match


def add_prefix(df, prefix, include=[], exclude=["name"]):
    """
    Add prefix to column

    Parameters
    ----------
    df: pandas.DataFrame
    prefix: str
    include: list
        empty = all included
    exclude: list
        default: ["name"]
    """
    def addPrefix(col):
        """ A lambda-like function """
        if col in exclude:
            return col
        elif len(include) == 0:
            return prefix + col
        elif col in include:
            return prefix + col
        else:
            return col

    return df.rename(columns=addPrefix)


def read_expression(file_abundance):
    """Read expression data, included tpm and count"""
    table = pd.read_csv(file_abundance)
    table = add_prefix(table, 'xprs_')
    return table


def read_transrate(file_transrate):
    """Read transrate data"""
    # read
    transrate = pd.read_csv(file_transrate)

    # rename and select
    transrate = transrate.rename(columns={
        'contig_name': 'name',
        'length': 'length',
        'p_good': 'good',
        'p_bases_covered': 'bases_covered',
        'p_seq_true': 'seq_true',
        'score': 'score',
        'p_not_segmented': 'not_segmented'})
    transrate = transrate.loc[:, ('name', 'length', 'good', 'bases_covered',
                                  'seq_true', 'score', 'not_segmented')]

    transrate = add_prefix(transrate, "tr_", exclude=['name', 'length'])
    return transrate


def filter_blastn(file_blast_self, identity=70, evalue=1e-5, length=0):
    """
    Filter self-blast result with some critia

    Parameters
    ----------
    file_blast_self: str
        The path to tsv file generated from blastn
    identity: float
    evalue: float
    length: float

    Returns
    -------
    blast_df: pandas.Dataframe
        This is a filtered table.
        The columns are same as blast format.
    """

    # read
    columns = ['q_name', 'r_name', 'identity', 'm_length', 'mismatch', 'gap',
               'q_start', 'q_end', 'r_start', 'r_end', 'evalue', 'bitscore']
    blast_df = pd.read_table(file_blast_self, sep='\t', header=None, names=columns)

    # filter
    length_f    = blast_df.loc[:, 'm_length'] >= length
    identity_f  = blast_df.loc[:, 'identity'] >= identity
    evalue_f    = blast_df.loc[:, 'evalue']   <= evalue
    duplicate_f = blast_df.loc[:, 'q_name']   != blast_df.loc[:, 'r_name']

    # return
    return blast_df.loc[length_f & identity_f & evalue_f & duplicate_f, :]


def map_blastid(blast_df, map_q, map_r):
    """
    Map name to id

    Parameters
    ----------
    blast_df: pandas.Dataframe
    map_q: dict
        The map query name to id
    map_r: dict
        The map reference name to id

    Returns
    -------
    blast_df: pandas.Dataframe
    """
    blast_df['q_name'] = blast_df['q_name'].map(map_q)
    blast_df['r_name'] = blast_df['r_name'].map(map_r)
    return blast_df


def intersect_match(blast_df, q_df, r_df):
    """
    Connect two sequences together.
    The connection is based on blast result.

    Parameters:
    ----------
    blast_df: pandas.Dataframe
    q_df: pandas.DataFrame
        The query sequence dataframe
    r_df: pandas.DataFrame
        The reference sequence dataframe

    Returns
    -------
    matches: dict of Match
        key: combination name of two sequences
        value: Match
    """
    # Save each record into Match
    matches = {}
    for data in blast_df.itertuples():
        q_seq = q_df.loc[data.q_name]
        r_seq = r_df.loc[data.r_name]
        name = Match.getName(data)

        if name not in matches:
            matches[name] = Match(data, q_seq, r_seq)
        matches[name].extend(data)

    # Calculate the identity
    for m_name in matches.keys():
        matches[m_name].calculateIdentity()
    return matches


def reshape_component(component_label, component_size):
    """
    Reshape from component id pre sequences to
    sequences list per components

    Parameters
    ----------
    component_label: array of int
    component_size: int

    Returns
    -------
    components: array of list
    """
    components = [[] for _ in range(component_size)]
    for id, label in enumerate(component_label):
        components[label].append(id)
    return components


def construct_graph(matches, size, threshold=90):
    """
    Run connected component for sequences

    Parameters
    ----------
    matches: dict of Match
    threshold: float
        The threshold of global alignment for both accuracy and recovery

    Returns
    -------
    components: array of list
        The list is the sequences id in this component
    component_label: list of int
        The label of each transcript
    component_name: list of str
        The name of component
    """
    # link similar sequences
    links = []
    for match in matches.values():
        if match.q_global_identity > threshold or \
           match.r_global_identity > threshold:
            links.append((match.q_name, match.r_name))

    # run connected_components
    mat = coo_matrix((np.ones(len(links), dtype=np.bool), list(zip(*links))), shape=(size, size))
    component_size, component_label = connected_components(mat, directed=False)

    # reshape
    components = reshape_component(component_label, component_size)
    return components, component_label, range(component_size)


def read_gene(map_name_id, file_gtf):
    """
    Read gtf to preform gene-isoform matching

    Parameters
    ----------
    map_name_id: dict
        key: str
            the name of transcript
        value: int
            ID
    file_gtf: str
        The path to gtf data

    Returns
    -------
    components: array of list
        The list is the sequences id in this component
    """
    gene_name_id = {}  # Tmporary dict store name as key as value as index
    component_label = [-1] * len(map_name_id)

    # read data
    columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'header']
    ref_gtf = pd.read_table(file_gtf, sep='\t', header=None, names=columns)
    for i, data in ref_gtf.iterrows():
        regex = re.match(r'gene_id "(\S+)"; transcript_id "(\S+)"', data['header'])
        if regex:
            gene_name = regex.group(1)
            transcript_name = regex.group(2)

            # save to genes dictionary
            if gene_name not in gene_name_id:
                gene_name_id[gene_name] = len(gene_name_id)

            # Assign gene label to transcript
            # this also remove the duplicated transcript name
            component_label[map_name_id[gene_name]] = gene_name_id[gene_name]

    # reshape
    components = reshape_component(component_label, len(gene_name_id))
    return components, component_label, gene_name_id.keys()


def aggregate_abundance(sequence_df, components, column_abundance, prefix=""):
    """
    Aggregate the sequences TPM/counts in the same component

    Sum, Mean, Max, Percetage, Ratio are calcuated

    Parameters
    ----------
    sequence_df: pandas.DataFrame
        The table of sequences
    components: array of list
        The list is the sequences id in this component
    column_abundance: array of str
        The list of abundance columns name

    Returns
    -------
    sequence_df: pandas.DataFrame
        The table of sequences
    component_df: pandas.DataFrame
        The table of components
    """
    # set column name
    column_com = ["tot_" + name for name in column_abundance] + \
                 ["avg_" + name for name in column_abundance] + \
                 ["max_" + name for name in column_abundance]
    column_seq = ["contribute_" + name for name in column_abundance] + \
                 ["relative_" + name for name in column_abundance]

    # save to numpy is quicker
    n = len(column_abundance)
    ori_data = np.array(sequence_df[column_abundance])
    data_seq = np.zeros((len(ori_data),     len(column_seq)))
    data_com = np.zeros((len(components), len(column_com)))

    # Calculate sequences tpm and its components
    for i, comp in enumerate(components):
        data = ori_data[comp]
        data_com[i, 0:n] = data.sum(axis=0)
        data_com[i, n:n*2] = data.mean(axis=0)
        data_com[i, n*2:n*3] = data.max(axis=0)
        data_seq[comp, 0:n] = data / data.sum(axis=0)
        data_seq[comp, n:n*2] = data / data.max(axis=0)

    # Save to pandas
    sequence_df[column_seq] = data_seq
    component_df = pd.DataFrame(data_com, columns=column_com)
    component_df['size'] = [len(comp) for comp in components]

    # remove divided 0
    sequence_df.fillna(0)
    component_df.fillna(0)

    sequence_df = add_prefix(sequence_df, prefix, include=column_seq)
    component_df = add_prefix(component_df, prefix)
    return sequence_df, component_df


def matches_to_table(matches):
    """Turn Match object to dataframe"""
    matches_array = []
    columns = ['match_name', 'contig_index', 'ref_index', 'accuracy', 'recovery']
    for match in matches.values():
        matches_array.append((match.name, match.q_name, match.r_name,
                              match.q_global_identity, match.r_global_identity))
    return pd.DataFrame(matches_array, columns=columns)


def diff_ref_contig(merged_df):
    """Calculate the tpm/count/length difference"""
    # calculate difference in length
    merged_df['length_difference'] = (merged_df['contig_length'] - merged_df['ref_length']) / \
                                     (merged_df['contig_length'] + merged_df['ref_length']) * 100

    # calculate difference in tpm/count
    for matric in ["tpm", "count"]:
        column = list(filter(lambda a: str(a).startswith(f"contig_xprs_{matric}_"), merged_df.columns))
        column_abundance_error = [i.replace(f"_{matric}_", f"_{matric}_error_") for i in column]
        column_answer = f"ref_xprs_{matric}_answer"
        merged_df[column_abundance_error] = (merged_df[column].sub(merged_df[column_answer], axis=0)) / \
                                            (merged_df[column].add(merged_df[column_answer], axis=0)) * 100

        # remove divied 0
        merged_df[column_abundance_error].fillna(value=0)

    return merged_df


def round_expression(merged_df, decimals=3):
    """Round the tpm/count to given decimals"""
    columns = list(filter(lambda a: "xprs_" in a or "tr_" in a, merged_df.columns))
    merged_df[columns] = merged_df[columns].round(3)
    return merged_df
