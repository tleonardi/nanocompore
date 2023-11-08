#!/usr/bin/env python3

# Local package
from nanocompore.common import *
import nanocompore.Whitelist as Whitelist

import nanocompore as pkg

data_base_dir = '/hps/nobackup/birney/users/logan/nanocompore_demo_data/data'
kd_1_dir = f'{data_base_dir}/KD_1'
kd_2_dir = f'{data_base_dir}/KD_2'
wt_1_dir = f'{data_base_dir}/WT_1'
wt_2_dir = f'{data_base_dir}/WT_2'


pod5_list1 = f"{kd_1_dir}/pod5/output.pod5"
pod5_list2 = f"{wt_1_dir}/pod5/output.pod5"

bam_list1 = f"{kd_1_dir}/basecalled_pod5.bam"
bam_list2 = f"{wt_1_dir}/basecalled_pod5.bam"

fasta_fn = '/nfs/research/birney/users/logan/refs/human/gencode.v36.transcripts.fa'
fn_dict = build_eventalign_fn_dict(pod5_list1, bam_list1, pod5_list2, bam_list2, 'KD', 'WT')

whitelist_Obj = Whitelist.Whitelist(fn_dict,
                          fasta_fn,
                          min_coverage = 30,
                          min_ref_length = 100,
                          select_ref_id = [],
                          exclude_ref_id = [],
                          n_samples=2)

