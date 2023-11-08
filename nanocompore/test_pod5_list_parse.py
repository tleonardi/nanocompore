#!/usr/bin/env python3

from collections import defaultdict
import sys

def build_eventalign_fn_dict(pod5_list1, bam_list1, pod5_list2, bam_list2, label1, label2):
    """
    Build the eventalign_fn_dict from file lists and labels
    """
    pod5_list1 = pod5_list1.split(",")
    bam_list1 = bam_list1.split(",")

    pod5_list2 = pod5_list2.split(",")
    bam_list2 = bam_list2.split(",")
    if len(pod5_list1) == len(bam_list1) and len(pod5_list2) == len(bam_list2):           
        d = defaultdict(dict)
        d[label1]=build_condition_dict(pod5_list1, bam_list1, label1)
        d[label2] = build_condition_dict(pod5_list2, bam_list2, label2)
        return d
    else:
        pod5_list1_len = len(pod5_list1)
        bam_list1_len = len(bam_list1)
        pod5_list2_len = len(pod5_list2)
        bam_list2_len = len(bam_list2)
        sys.stderr.write(f"Mismatch in file list lengths:\npod5_list1 {pod5_list1_len}; bam_list1 {bam_list1_len}\npod5_list2 {pod5_list2_len}; bam_list2 {bam_list2_len}\n")

def build_condition_dict(pod5_list, bam_list, label):
    condition_list = defaultdict()
    for i, (pod5, bam) in enumerate(zip(pod5_list, bam_list)):
        condition_list[f"{label}_{i}"] = {'pod5':pod5, 'bam':bam}
    return condition_list


def build_eventalign_fn_dict_from_tsv(infile):
    fn_dict = defaultdict(dict)
    num_samples = set()
    with open(infile, 'r') as tsv:
        for i, line in enumerate(tsv):
            if line:
                try:
                    cond, sample, pod5, bam = line.strip().split('\t')
                    num_samples.add(sample)
                    fn_dict[cond][sample] = {'pod5':pod5, 'bam':bam}
                except:
                    sys.stderr.write(f"Error with entry {i} in input samples tsv file\n")
                    sys.exit()

    if len(num_samples) != i+1:
        sys.stderr.write(f"Not all sample labels are unique\n")
        sys.exit()

    return fn_dict


pod5_list1 = 'f1.pod5,f2.pod5'
bam_list1 = 'f1.bam,f2.bam'
pod5_list2 = 'f3.pod5,f4.pod5'
bam_list2 = 'f3.bam,f4.bam'
label1 = 'WT'
label2 = 'KD'

fn_dict = build_eventalign_fn_dict(pod5_list1, bam_list1, pod5_list2, bam_list2, label1, label2)

for condition in fn_dict:
    print(condition)
    for sample in fn_dict[condition]:
        print(sample, fn_dict[condition][sample]['pod5'], fn_dict[condition][sample]['bam'])



fn_dict = build_eventalign_fn_dict_from_tsv('samples.tsv')

for condition in fn_dict:
    print(condition)
    for sample in fn_dict[condition]:
        print(sample, fn_dict[condition][sample]['pod5'], fn_dict[condition][sample]['bam'])