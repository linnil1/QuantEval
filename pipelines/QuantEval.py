import argparse
import logging

import pandas as pd
from tqdm import tqdm

from quanteval_utils import *
from utils import setupLogger


def main(args):
    """
    Merge tpm/count by blast result
    """
    logger.debug("Read Expression")
    sequence_df = read_expression(args.abundance)
    column_abundance = sequence_df.columns.tolist()[1:]
    sequence_df.reset_index(level=0, inplace=True)

    # create a dict map name to id
    map_name_id = dict(zip(sequence_df['name'], sequence_df.index))
    n = len(map_name_id)
    if n != len(sequence_df):
        raise ValueError("The transcript name is duplicated")

    logger.debug("Read Transrate")
    transrate = read_transrate(args.transrate)
    sequence_df = pd.merge(sequence_df, transrate, on='name', how="outer")

    logger.debug("Read Blast")
    blast_df = filter_blastn(args.blast)
    blast_df = map_blastid(blast_df, map_name_id, map_name_id)

    # build connect component
    logger.debug("Run connect component for query sequences")
    name = "component"
    matches = intersect_match(blast_df, sequence_df, sequence_df)
    components, component_label, component_name = construct_graph(matches, size=n)

    logger.debug("Recalculate abundance")
    sequence_df, component_df = aggregate_abundance(sequence_df, components, column_abundance, name + "_")


    logger.debug("Merge the table of connected_components and sequences")
    sequence_df[name] = component_label
    component_df.reset_index(level=0, inplace=True)
    component_df[name] = component_name
    merged_df = pd.merge(sequence_df, component_df, how='left',
                         left_on=name, right_on='index', suffixes=["", '_' + name])

    # gene component (gene isoform)
    if args.gtf:
        name = "gene"
        logger.debug("Read GTF file")
        components, component_label, genes_name = read_gene(map_name_id, args.gtf)

        logger.debug("Recalculate abundance for gene")
        sequence_df, gene_df = aggregate_abundance(sequence_df, components, column_abundance, name + "_")

        logger.debug("Merge the table of connected_components and sequences and gene-isoform")
        sequence_df[name] = component_label
        merged_df[name] = component_label  # add gene label for sequences but also need to add in merged datafrmae
        gene_df.reset_index(level=0, inplace=True)
        gene_df[name] = genes_name
        merged_df = pd.merge(merged_df, gene_df, how='left',
                             left_on=name, right_on='index', suffixes=["", '_' + name])

    logger.debug(f"Save to {args.output}")
    merged_df.to_csv(args.output, sep='\t', index=False)


def match(args):
    """Match reference table and contig table"""
    logger.debug("Read reference and contig tsv")
    reference_df = pd.read_csv(args.reference, sep='\t')
    contig_df = pd.read_csv(args.contig, sep='\t')

    logger.debug("Read Blast")
    blast_df = filter_blastn(args.blast)
    map_qname_id = dict(zip(contig_df['name'], contig_df.index))
    map_rname_id = dict(zip(reference_df['name'], reference_df.index))
    blast_df = map_blastid(blast_df, map_qname_id, map_rname_id)

    logger.debug("Matches to table")
    matches = intersect_match(blast_df, contig_df, reference_df)
    matches_df = matches_to_table(matches)

    logger.debug("Merge tables")
    reference_df = add_prefix(reference_df, "ref_", exclude=[])
    contig_df = add_prefix(contig_df, "contig_", exclude=[])
    merged_df = pd.merge(matches_df, contig_df, on='contig_index', how='inner')
    merged_df = pd.merge(merged_df, reference_df, on='ref_index', how='inner')

    logger.debug("Calculate difference between contig and reference in each matches")
    merged_df = diff_ref_contig(merged_df)
    merged_df.to_csv(args.output, sep='\t', index=False)


if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description="QuantEval")
    subparsers = parser.add_subparsers(title="mode", dest="mode", required=True, help="Choose Mode")
    subparser = subparsers.add_parser("contig", aliases=['reference'])
    subparser.add_argument("--abundance", type=str, required=True, help="The file of abundance data from preprocessQuantEval in pipelines.ipy")
    subparser.add_argument("--transrate", type=str, required=True, help="The transrate/contig.tsv file")
    subparser.add_argument("--blast",     type=str, required=True, help="The self-blast result with format 6")
    subparser.add_argument("--gtf",       type=str, default="",    help="(Optional) The gtf file of mRNA")
    subparser.add_argument("--output",    type=str, required=True, help="The tsv of table")

    subparser = subparsers.add_parser("match")
    subparser.add_argument("--reference", type=str, required=True, help="The reference table")
    subparser.add_argument("--contig",    type=str, required=True, help="The contig table")
    subparser.add_argument("--blast",     type=str, required=True, help="The contig-mrna blast result with format 6")
    subparser.add_argument("--output",    type=str, required=True, help="The tsv of table")

    # logger
    logger = setupLogger()

    # main
    args = parser.parse_args()
    logger.info(f"Run QuanEval with mode {args.mode}")
    if args.mode == "match":
        match(args)
    else:
        main(args)
    logger.info(f"Done {args.output}")
