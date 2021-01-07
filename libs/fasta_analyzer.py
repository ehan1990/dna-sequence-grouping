from operator import itemgetter
from libs import util
import logging
logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")

class FastaAnalyzer:

    __data = []
    __error_rate_dist = None
    __unique_data = []

    def __init__(self, data, error_rate_dist):
        self.__data = data
        self.__unique_data = []
        self.__result_data = []
        self.__error_rate_dist = error_rate_dist

    def get_result_data(self): return self.__result_data

    def find_biggest_distance(self, sorted_output_seq_map):
        biggest_dist = 0
        source = sorted_output_seq_map[0]
        source_seq = source["sequence"]
        for i in range(len(sorted_output_seq_map)):
            target_seq = sorted_output_seq_map[i]["sequence"]
            dist, norm = util.get_lev_distances(source_seq, target_seq)
            if dist > biggest_dist:
                biggest_dist = dist
        return biggest_dist

    def find_unique_seqs(self):
        seq_map = {}
        for item in self.__data:
            seq = item["sequence"]
            if seq_map.get(seq) is None:
                seq_map[seq] = 1
                self.__unique_data.append(item)
            else:
                seq_map[seq] += 1

        return self.__unique_data, seq_map

    def verify_seq_map(self, total_count, seq_map):
        sum = 0
        for key, val in seq_map.items():
            sum += val
        if total_count == sum:
           logging.info("verified seq map")
        else:
           logging.info("seq map contains invalid count")

    def sort_seq_map(self, seq_map):
        output_seq_map = []
        for key, val in seq_map.items():
            output_seq_map.append({"sequence": key, "count": val})
        sorted_output_seq_map = sorted(output_seq_map, key=itemgetter("count"), reverse=True)
        return sorted_output_seq_map

    def filter_low_unique_count_seq_map(self, sorted_output_seq_map, thresh):
        output = []
        for item in sorted_output_seq_map:
            if item["count"] > thresh:
                output.append(item)
        return output

    def compute_dist(self, sorted_seq_map):
        for i in range(len(sorted_seq_map)):
            if i+1 > len(sorted_seq_map):
                break
            logging.info("computing {}/{}".format(i+1, len(sorted_seq_map)))

            source = sorted_seq_map[i]
            source_seq = source["sequence"]
            source_data = {
                "source": source,
                "target_data": []
            }
            for target in sorted_seq_map[i+1:]:
                target_seq = target["sequence"]
                dist, norm_dist = util.get_lev_distances(source_seq, target_seq)
                is_subset = util.check_subset(source_seq, target_seq)
                similarity = round(1 - norm_dist, 3)

                target_data = {
                    "source": source,
                    "target": target,
                    "pair": "{}:{}".format(source_seq, target_seq),
                    "damerau_lev": dist,
                    "is_subset": is_subset,
                    "normalized_damerau_lev": norm_dist,
                    "similarity": similarity
                }
                source_data["target_data"].append(target_data)
            self.__result_data.append(source_data)

        all_target_data = []
        for result in self.__result_data:
            target_data = result["target_data"]
            all_target_data.extend(target_data)
        sorted_target_data = sorted(all_target_data, key=itemgetter("similarity"), reverse=True)
        return sorted_target_data

    def compute_groups(self, sorted_seq_map, error_dist_thresh):
        """
        [
            {
                "total_count": "5",
                "source": {"count": 45, "sequence": seq1}
                "sequences": [
                    {"count": 5, "sequence": seq2}
                ]
            }
        ]

        """
        output = []
        removed_seqs = {}

        for i in range(len(sorted_seq_map)):
            merged_count = len(removed_seqs.keys())
            logging.info("computing {}/{}, merged {}".format(i+1, len(sorted_seq_map), merged_count))
            source = sorted_seq_map[i]
            source_seq = source["sequence"]
            total_count = source["count"]
            similar_seq_data = []
            if i+1 > len(sorted_seq_map):
                break

            if removed_seqs.get(source_seq):
                continue

            j = 1
            for target in sorted_seq_map[i+1:]:
                target_seq = target["sequence"]

                if removed_seqs.get(target_seq):
                    continue

                is_subset = util.check_subset(source_seq, target_seq)
                dist, norm = util.get_lev_distances(source_seq, target_seq)

                if is_subset or dist <= error_dist_thresh:
                    total_count += target["count"]
                    similar_seq_data.append(target)
                    removed_seqs[target_seq] = True
                j += 1

            group = {
                "total_count": total_count,
                "source": source,
                "sequences": similar_seq_data
            }
            output.append(group)
        sorted_output = sorted(output, key=itemgetter("total_count"), reverse=True)

        return sorted_output

    def verify_groups(self, group_data, error_dist_thresh):
        errors = 0
        for i in range(len(group_data)):
            print("verifying {}/{}".format(i + 1, len(group_data)))
            item = group_data[i]
            source_seq = item["source"]["sequence"]
            seq_data = item["sequences"]

            for j in range(len(seq_data)):
                target_seq = seq_data[j]["sequence"]
                dist, norm = util.get_lev_distances(source_seq, target_seq)
                is_subset = util.check_subset(source_seq, target_seq)
                possible_error = util.check_possible_error(dist, error_dist_thresh)
                if is_subset or possible_error:
                    continue
                else:
                    print("error in {}/{}, {}:{}".format(i + 1, len(group_data), source_seq, target_seq))
                    errors += 1
        return errors


    def verify_group_count(self, sorted_group_data, ): pass

