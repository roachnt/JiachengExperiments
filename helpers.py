from random import random


def print_run_ratio(bad_dict, good_dict):
    good_runs = 0
    bad_runs = 0
    for index in bad_dict:
        if bad_dict[index] == good_dict[index]:
            good_runs += 1
        else:
            bad_runs += 1
    print("Bad Runs:", bad_runs)
    print("Good Runs:", good_runs)


def get_run_ratio(bad_dict, good_dict):
    good_runs = 0
    bad_runs = 0
    for index in bad_dict:
        if bad_dict[index] == good_dict[index]:
            good_runs += 1
        else:
            bad_runs += 1
    return bad_runs, good_runs


def fuzzy(good_expression, gen_bad):
    if gen_bad:
        return good_expression*2*random()
    else:
        return good_expression
