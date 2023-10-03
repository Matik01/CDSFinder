import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Prokaryotic gene prediction tool")
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input FASTA file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output BED file")
    return parser.parse_args()


def fasta_file_parse(input_file_path):
    header_seq_dict = dict()
    with open(input_file_path, "r") as fasta:
        header = ''
        sequence = []
        for line in fasta:
            if line.startswith(">"):
                if header:
                    header_seq_dict[header] = "".join(sequence)
                    sequence = []
                header = line.strip()
            elif not line.startswith(">"):
                sequence.append(line.strip())

    return header_seq_dict


def seq_to_bed_convert(output_file_path, features_indexes):
    with open(output_file_path, "a") as f_out:
        for index in features_indexes:
            f_out.write(f"{index[3::]}\t{features_indexes[index][1]}\tStart codon\tTYPE=RBS_SITE\n")
        f_out.close()

def find_features(fasta, header, start_codon="ATG",
                  shine_dalgarno="AGGAGGT",
                  # pribnow_box="TATAAT"
                  ):
    header_feature_dict = dict()

    seq = fasta[header]
    header_feature_dict["SD_" + header] = greedy_search(seq, shine_dalgarno)
    # header_feature_dict["PB_" + header] = greedy_search(seq, pribnow_box) adding in next versions!
    header_feature_dict["SC_" + header] = greedy_search(seq, start_codon)

    return header_feature_dict


def greedy_search(seq, motif):
    kmer = len(motif)
    motif_score_matrix = dict()
    for i in range(0, len(seq) + 1, 1):
        target_motif = seq[i:i + kmer]
        motif_score = score(target_motif=target_motif, motif=motif)

        if kmer != 3 and motif_score > 4:
            motif_score_matrix[i] = motif_score
        if kmer == 3 and motif_score == 3:
            motif_score_matrix[i] = motif_score

    return motif_score_matrix


def score(target_motif, motif):
    if target_motif == motif:
        return len(motif)
    count = 0
    for i, j in zip(target_motif, motif):
        if i == j:
            count += 1
    return count


def gene_location_mapper(features):
    gene_location = dict()
    index_pb = int()
    for header in features:
        if header.startswith("SD_"):
            index_map = features[header]
            for index in index_map:
                index_pb = index

        elif header.startswith("SC_"):
            index_map = features[header]
            for index in index_map:
                if index_pb + 17 >= index > index_pb:
                    gene_location[header] = [index_pb, index]
                    index_pb = int()

    return gene_location


def main():
    args = parse_args()
    input_file_path = args.input
    output_file_path = args.output

    with open(output_file_path, "w") as f_out:
        f_out.write("Header\tPosition\tFeature\n")
        f_out.close()

    parsed_fasta = fasta_file_parse(input_file_path=input_file_path)

    for key in parsed_fasta:
        features = find_features(parsed_fasta, header=key)
        indexes = gene_location_mapper(features=features)
        seq_to_bed_convert(output_file_path=output_file_path, features_indexes=indexes)


if __name__ == "__main__":
    main()
