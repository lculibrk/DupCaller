import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Trim fastq sequences and move barcodes to name"
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        type=str,
        help="folder where DupCallerCall results are stored",
    )
    parser.add_argument("-o", "--output", type=str, help="output filename")
    args = parser.parse_args()
    samples = args.input
    with open(args.output, "w") as output:
        output.write(
            f"sample\tpass_filter_reads\tunique_reads\tread_families\tduplication_rate\tread_family_efficiency\tsnv_number\tsnv_effective_coverage\tsnv_naive_burden\tsnv_least_square_burden\tsnv_least_square_burden_upper_uci\tsnv_least_square_burden_lci\tindel_number\tindel_effective_coverage\tindel_naive_burden\n"
        )
        for sample in samples:
            stats_file = f"{sample}/{sample}_stats.txt"
            snv_burden_file = f"{sample}/{sample}_snv_burden.txt"
            indel_burden_file = f"{sample}/{sample}_indel_burden.txt"
            with open(stats_file) as stats:
                lines = stats.readlines()
                uniq_reads = int(lines[0].strip("\n").split("\t")[1])
                pf_reads = int(lines[1].strip("\n").split("\t")[1])
                pf_read_family = int(lines[2].strip("\n").split("\t")[1])
                eff_cov = int(lines[3].strip("\n").split("\t")[1])
                dup_rate = float(lines[5].strip("\n").split("\t")[1])
                efficiency = float(lines[6].strip("\n").split("\t")[1])
            with open(snv_burden_file) as stats:
                lines = stats.readlines()
                snv_num = int(lines[0].strip("\n").split("\t")[1])
                naive_burden = float(lines[1].strip("\n").split("\t")[1])
                lsq_burden = float(lines[2].strip("\n").split("\t")[1])
                uci = float(lines[3].strip("\n").split("\t")[1])
                lci = float(lines[4].strip("\n").split("\t")[1])
            with open(indel_burden_file) as stats:
                lines = stats.readlines()
                indel_num = int(lines[0].strip("\n").split("\t")[1])
                indel_cov = int(lines[1].strip("\n").split("\t")[1])
                indel_naive_burden = float(lines[2].strip("\n").split("\t")[1])
            output.write(
                f"{sample}\t{pf_reads}\t{uniq_reads}\t{pf_read_family}\t{dup_rate}\t{efficiency}\t{snv_num}\t{eff_cov}\t{naive_burden}\t{lsq_burden}\t{uci}\t{lci}\t{indel_num}\t{indel_cov}\t{indel_naive_burden}\n"
            )
