from difflib import SequenceMatcher
from pyxdameraulevenshtein import damerau_levenshtein_distance, normalized_damerau_levenshtein_distance

from libs.constants import *


def remove_duplicates(data):
    unique_dna_list = set(data)
    return list(unique_dna_list)


def check_subset(source_seq, target_seq):
    """
    if source contains target, or vice-versa, return true
    """
    if source_seq in target_seq:
        return True
    if target_seq in source_seq:
        return True
    return False


def check_possible_error(lev_distance, error_rate_dist):
    if lev_distance <= error_rate_dist:
        return True
    return False


def get_lev_distances(source, target, round_to=DEFAULT_ROUND):
    dist = damerau_levenshtein_distance(source, target)
    norm_dist = normalized_damerau_levenshtein_distance(source, target)
    return dist, round(norm_dist, round_to)


def get_seq_matcher_ratio(source, target, round_to=DEFAULT_ROUND):
    s = SequenceMatcher(None, source, target)
    ratio = s.ratio()
    return round(ratio, round_to)


def get_fuzzy_token_ratio(source, target):
    try:
        if len(source) > len(target):
            substring = target
            longstring = source
        else:
            substring = source
            longstring = target
        res = fuzzysearch.find_near_matches(substring, longstring, max_l_dist=6)
    except Exception as e:
        print("error for {}:{}".format(source, target))
        print(e)
        return

    smallest_dist = None
    for item in res:
        if smallest_dist is None:
            smallest_dist = item.dist
            continue
        if item.dist < smallest_dist:
            smallest_dist = item.dist
    return smallest_dist


def find_seq_in_group_list(seq, group_list):
    for group in group_list:
        common = intersection(seq, group)
        if len(common) > 0:
            return group
    return []


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
