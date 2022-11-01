#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   calculate_similarity.py
#         Author:   yujie
#    Description:   calculate_similarity.py
#        Version:   1.0
#           Time:   2022/10/27 11:59:35
#  Last Modified:   2022/10/27 11:59:35
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
# from humre import *  # 正则
from icecream import ic  # 打印
import argparse  # 命令行
import linecache  # 大文件行读取
import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
import time
import copy  # 深度拷贝
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
# \n\
# \n\
#       Filename:   calculate_similarity.py\n\
#         Author:   yujie\n\
#    Description:   calculate_similarity.py\n\
#        Version:   1.0\n\
#           Time:   2022/10/27 12:04:12\n\
#  Last Modified:   2022/10/27 12:04:12\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   GNU General Public License v3.0\n\
# \n\
# \n\
\n\
\npython3   calculate_similarity.py\n\
Function:\n\
1.常规使用\n\
1.1 -i [ ] -o [ ] \n\
2.其他使用\n\
2.1 -i [ ] -o [ ] \n\
\n\
# \n\
Path: E:\OneDrive\jshy信息部\Script\Denovo\calculate_similarity.py\n\
Path: /share/nas1/yuj/script/Denovo/assembly/calculate_similarity.py\n\
Version: 1.0\n\
# \n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i1', '--invcf', metavar='[vcf]', help='vcf', type=str, default='F:/2596-1/archive/samples.pop.snp.recode.vcf', required=False)
optional.add_argument(
    '-i2', '--infasta', metavar='[fa]', help='fa', type=str, default='F:/2596-1/archive/tags.consensus.fa', required=False)
optional.add_argument(
    '-g', '--ingroup', metavar='[group]', help='group', type=str, default='F:/2596-1/archive/group.txt', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/2596-1/相似度结果.log', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[help_information]')
args = parser.parse_args()

begin_time = time.time()
in_vcf_path = args.invcf
in_fa_path = args.infasta
in_group_path = args.ingroup
out_log_path = args.outfile

# 读取分组,生成所需变量
createVar = locals()
with open(in_group_path, 'r') as in_group_handle:
    f1_list = []
    for line in in_group_handle:
        if not line.startswith('#'):
            f1 = line.rstrip().split('\t')[0]
            f1_list.append(f1)
            createVar['dict_'+f1] = {}
            parent1 = line.rstrip().split('\t')[1]
            parent2 = line.rstrip().split('\t')[2]
            createVar['dict_'+f1][parent1] = 0
            createVar['dict_'+f1][parent2] = 0

# 先读取fa文件 计算出所有长度
fa_contents = linecache.getlines(in_fa_path)
seq_id = ''
dict_seq = {}
for line in fa_contents:
    if line.startswith('>'):
        seq_id = line.strip('\n').lstrip('>')
        dict_seq[seq_id] = 0
        seq = ''
    else:
        if line.startswith('A') or line.startswith('T') or line.startswith('G') or line.startswith('C'):
            seq += line.strip('\n')
            dict_seq[seq_id] = len(seq)


# 读取vcf并计算不同的
fa_len_count = 0
vcf_contents = linecache.getlines(in_vcf_path)

for line in vcf_contents:
    if line.startswith('#CH'):
        line_content = line.rstrip().split('\t')
        for f1 in f1_list:
            createVar['index_'+f1] = line_content.index(f1)
            for parent in createVar['dict_'+f1].keys():
                createVar['index_'+parent] = line_content.index(parent)

    if not line.startswith('#'):
        line_content = line.rstrip().split('\t')
        fa_len_count += dict_seq[line_content[0]]

        for f1 in f1_list:
            f1_type = line_content[createVar['index_'+f1]].split(':')[0]
            for parent in createVar['dict_'+f1].keys():
                parent_type = line_content[createVar['index_'+parent]].split(':')[
                    0]
                if f1_type != parent_type:
                    createVar['dict_'+f1][parent] += 1


# 写入文件
with open(out_log_path, 'w') as out_log_handle:
    print('子代\t亲本\t相似度')
    out_log_handle.write('子代\t亲本\t相似度\n')
    for f1 in f1_list:
        for parent in createVar['dict_'+f1].keys():
            print('{}\t{}\t{:.5f}'.format(f1, parent,
                                          1-(createVar['dict_'+f1][parent]/fa_len_count)))
            out_log_handle.write('{}\t{}\t{:.5f}\n'.format(f1, parent,
                                                           1-(createVar['dict_'+f1][parent]/fa_len_count)))


# 运行时间
print('Already Run {}s'.format(time.time()-begin_time))
