import click
import datetime
import logging
import sys
logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

from libs.fasta_file import FastaFile
from libs.fasta_analyzer import FastaAnalyzer
from libs.constants import *


@click.group()
def entry_point():
    pass


@entry_point.command("biggest")
@click.option("--name", "title", required=True)
def find_biggest_dist_cmd(title):
    logging.info("----- starting {} cmd-----".format(find_biggest_dist_cmd.name))
    logging.info("title: {}".format(title))
    input_filepath = "dna_files/{}.fasta".format(title)
    data = FastaFile.read(input_filepath)

    fasta_analyzer = FastaAnalyzer(data, 0)

    unique_seq_data, seq_map = fasta_analyzer.find_unique_seqs()
    fasta_analyzer.verify_seq_map(len(data), seq_map)
    sorted_seq_map = fasta_analyzer.sort_seq_map(seq_map)
    logging.info("unique seq: {}".format(len(unique_seq_data)))

    biggest_dist = fasta_analyzer.find_biggest_distance(sorted_seq_map)
    logging.info("biggest dist: {}".format(biggest_dist))


@entry_point.command("compute")
@click.option("--name", "title", required=True)
@click.option("--error-dist", "error_rate_dist", default=DEFAULT_LEV_ERROR_THRESH)
@click.option("--num-unique", "num_unique_thresh", default=DEFAULT_NUM_UNIQUE_COUNT_THRESH)
def compute_cmd(title, error_rate_dist, num_unique_thresh):
    logging.info("----- starting new job -----")
    logging.info("title: {}".format(title))
    logging.info("error_rate_dist: {}".format(error_rate_dist))

    input_filepath = "dna_files/{}.fasta".format(title)

    groups_file = "output/{}.groups.json".format(title)
    groups_fasta_file = "output/{}.groups.fasta".format(title)

    data = FastaFile.read(input_filepath)
    logging.info("total input seq: {}".format(len(data)))

    fasta_analyzer = FastaAnalyzer(data, error_rate_dist)

    unique_seq_data, seq_map = fasta_analyzer.find_unique_seqs()
    fasta_analyzer.verify_seq_map(len(data), seq_map)
    sorted_seq_map = fasta_analyzer.sort_seq_map(seq_map)
    logging.info("unique seq: {}".format(len(unique_seq_data)))

    logging.info("finding biggest dist")
    biggest_dist = fasta_analyzer.find_biggest_distance(sorted_seq_map)
    logging.info("biggest dist: {}".format(biggest_dist))

    top_count = sorted_seq_map[0]["count"]
    d0_filepath = "output/{}.d0_c{}_g{}_t{}.fasta".format(title, top_count, len(sorted_seq_map), len(data))
    FastaFile.write_seq_map_fasta(d0_filepath, title, sorted_seq_map)

    logging.info("starting to compute dist")
    group_data = fasta_analyzer.compute_groups(sorted_seq_map, error_rate_dist)
    logging.info("groups: {}".format(len(group_data)))

    errors = fasta_analyzer.verify_groups(group_data, error_rate_dist)
    logging.info("errors: {}".format(errors))

    FastaFile.write_family_fasta(title, "output", error_rate_dist, group_data, len(data))

    logging.info("done")


entry_point.add_command(compute_cmd)

if __name__ == "__main__":
    entry_point()

