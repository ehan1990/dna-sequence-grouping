import json
from operator import itemgetter


class FastaFile:

    def __init__(self): pass

    @staticmethod
    def read(filepath):
        with open(filepath) as f:
            lines = f.readlines()
        parsed_data = FastaFile._parse(lines)
        return parsed_data

    @staticmethod
    def write_fasta(filepath, data):
        with open(filepath, "w") as f:
            for d in data:
                title = d["title"]
                seq = d["sequence"]
                f.write(title + "\n")
                f.write(seq + "\n")

    @staticmethod
    def write_json(filepath, data):
        with open(filepath, "w") as f:
            s = json.dumps(data, indent=2)
            f.write(s)

    @staticmethod
    def write_seq_map_fasta(filepath, title, data):
        with open(filepath, "w") as f:
            c = 1
            for item in data:
                f.write("> {}-{} seq appeared {} times\n".format(title, c, item["count"]))
                f.write(item["sequence"] + "\n")
                c += 1

    @staticmethod
    def write_group_fasta(filepath, group_data):
        with open(filepath, "w") as f:
            c = 1
            for item in group_data:
                count = item["source"]["count"]
                total_count = item["total_count"]
                sequence = item["source"]["sequence"]
                f.write("> seq appeared {} times, group appeared {} times\n".format(count, total_count))
                f.write(sequence + "\n")
                c += 1

    @staticmethod
    def write_family_fasta(title, root_dir, dist, group_data, input_total):
        top_source_count = group_data[0]["source"]["count"]
        filepath = "{}/{}_d{}_c{}_g{}_t{}.fasta".format(root_dir, title, dist, top_source_count, len(group_data), input_total)
        legend_filepath = "{}/{}_legend.txt".format(root_dir, title)
        legend_data = []

        with open(filepath, "w") as f:
            group_index = 1
            seq_index = 1
            for group in group_data:
                seqs = [group["source"]]
                total_count = group["total_count"]
                seqs.extend(group["sequences"])
                start_index = seq_index

                for item in seqs:
                    seq_count = item["count"]
                    f.write(">{}-{} seq appeared {} times, part of group {}, group appeared {}\n"
                            .format(title, seq_index, seq_count, group_index, total_count))
                    f.write(item["sequence"] + "\n")
                    end_index = seq_index
                    seq_index += 1

                legend_data.append({"group": group_index, "start": start_index, "end": end_index, "group_count": total_count})
                group_index += 1

        with open(legend_filepath, "w") as f:
            for item in legend_data:
                f.write("family {}, seq {}-{}, group appeared {}\n".
                        format(item["group"], item["start"], item["end"], item["group_count"]))


    @staticmethod
    def _parse(lines):
        parsed = []
        title = ""
        sequence = ""
        for line in lines:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == ">":
                if len(sequence) == 0:
                    title = line
                    continue
                parsed.append({"title": title, "sequence": sequence.lower()})
                title = line
                sequence = ""
            else:
                sequence += line
        parsed.append({"title": title, "sequence": sequence.lower()})
        return parsed
