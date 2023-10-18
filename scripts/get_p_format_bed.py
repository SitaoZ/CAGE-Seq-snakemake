#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: get_p_format_gff.py
date: 2022/3/28 8:16 PM
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import re,os
from collections import defaultdict
from collections import OrderedDict


# ID=AT1G01010.1:five_prime_UTR:1;Parent=AT1G01010.1;
g_p = re.compile(r'ID=(\S+)\:five_prime_UTR:(\d);')
n_p = re.compile(r'Parent=(\S+);')
c_p = re.compile(r'peakCounts:(\S+)')

l_p = re.compile(r'Parent=(\S+);')

def tag_count(info):
    c_p_m = c_p.search(info)
    if c_p_m:
        c = c_p_m.group(1)
        return float(c)
    else:
        return 0

def each_tag(chrom, start, end, strand):
    """p@Chr1:222...444
       BED 0-based ---> GFF 1-based
    """
    
    tag = 'p@'+chrom + ":" + str(int(start) + 1) + ".." + str(int(end)) + "," + strand
    return tag

def each_position(chrom, start, end, strand):
    """Chr1:222...444
       BED 0-based ---> GFF 1-based
    """
    tag = chrom + ":" + str(int(start) + 1) + ".." + str(int(end)) + "," + strand
    return tag

def utr_tag(gene_name, i):
    return f'p{i}@{gene_name}'




def p_format_non_utr5(filePath, outPath='robust_mRNA.bed'):
    """ generate robust mRNA bed file"""
    outFile = open(outPath, 'w')
    pDict = defaultdict(list)
    gene_strand = dict()
    with open(filePath, 'r') as f:
        for line in f:
            line = line.strip()
            array = line.split('\t')
            l_p_m = l_p.search(array[8]) # transcript id
            if l_p_m :
                locus = l_p_m.group(1)
                # locus = '.'.join(locus.split('.')[:-1])
                if array[9:] in pDict[locus]:
                    continue
                gene_strand[locus] = array[6]
                pDict[locus].append(array[9:])

    for locus in pDict.keys():
        for peak in pDict[locus]:
            if gene_strand[locus] != peak[5]:
                # gene strand != peak strand
                continue
            c = tag_count(peak[3])
            # jbrowser
            # chrom start end starnd
            peak_id = each_tag(peak[0], peak[1], peak[2], peak[5])
            # for peak transcript density compare
            peak[3] = peak_id
            # peak start + 1 ---> GFF 1-based
            # peak[1] = int(peak[1]) + 1 
            print(*peak[:9], file=outFile, sep='\t')

def p_format_utr5(filePath, outPath='robust_utr5.bed'):
    outFile = open(outPath, 'w')
    pDict = defaultdict(list)
    gene_strand = dict()
    with open(filePath, 'r') as f:
        for line in f:
            line = line.strip()
            array = line.split('\t')
            # if 'Name' in array[8]:
            if array[8]:
                g_p_m = g_p.search(array[8])
                n_p_m = n_p.search(array[8])
                if g_p_m and n_p_m:
                    locus = g_p_m.group(1)
                    gene_name = n_p_m.group(1)
                    # tranript id to gene id 
                    # locus = '.'.join(locus.split('.')[:-1])
                    # gene_name = '.'.join(gene_name.split('.')[:-1])
                    key = locus + '+' + gene_name
                    if array[9:] in pDict[key]:
                        continue
                    # strand
                    gene_strand[key] = array[6]
                    pDict[key].append(array[9:])

    for key in pDict.keys():
        locus, gene_name = key.split('+')
        d = dict()
        # if len(pDict[key]) == 1:
        #     print(pDict[key])
        for index, peak in enumerate(sorted(pDict[key]), start=1):
            if gene_strand[key] != peak[5]:
                # gene strand != peak strand
                continue
            c = tag_count(peak[3])
            if c == 0:
                if gene_strand[key] == '+': 
                    c = len(pDict[key]) - index + 1 # forward strand c越大，排序的值越小, 正链的最开始应该是最大的索引
                else:
                    c = index  # backward strand (rigth to left)
            d[c] = peak
            if locus == 'AT1G02520':
                print(pDict[key])
                print(gene_strand[key])
        if locus == 'AT1G02520':
            print(d, '和输出robust_utr5.bed 文件的顺序相反')
        for i, c in enumerate(sorted(d, reverse=True), start=1):
            # jbrowser
            peak_id = utr_tag(key, i)
            # for peak transcript density compare
            d[c][3] = peak_id
            # peak start + 1 ---> GFF 1-based
            # d[c][1] = int(d[c][1]) + 1
            print(*d[c][:9], file=outFile, sep='\t')


def remove_mRNA_utr5(mRNA, utr5, outPath='robust.bed'):
    outFile = open(outPath,'w')
    utr5D = OrderedDict()
    with open(utr5,'r') as f:
        for line in f:
            line = line.strip()
            array = line.split('\t')
            key = '_'.join(array[:3])
            utr5D[key] = line
    with open(mRNA, 'r') as f:
        for line in f:
            line = line.strip()
            array = line.split('\t')
            key = '_'.join(array[:3])
            if utr5D.get(key):
                utr_array = utr5D[key].split('\t')
                print(*utr_array, file=outFile, sep='\t')
            else:
                print(*array, file=outFile, sep='\t')




if __name__ == '__main__':
    # Robust
    p_format_non_utr5('p_format_in_mRNA_fantom_prepare.csv', outPath='robust_mRNA.bed')
    p_format_utr5('p_format_in_fantom_prepare.csv', outPath='robust_utr5.bed')
    remove_mRNA_utr5('robust_mRNA.bed','robust_utr5.bed', outPath='robust.bed')
    os.system('cat robust.bed | sort | uniq | sort -k 1,1 -k2,2n > robust.sort.bed')
    # Permissive
    p_format_non_utr5('p_format_in_mRNA_fantom_prepare_permissive.csv', outPath='permissive_mRNA.bed')
    p_format_utr5('p_format_in_fantom_prepare_permissive.csv', outPath='permissive_utr5.bed')
    remove_mRNA_utr5('permissive_mRNA.bed', 'permissive_utr5.bed', outPath='permissive.bed')
    os.system('cat permissive.bed | sort | uniq | sort -k 1,1 -k2,2n > permissive.sort.bed')
