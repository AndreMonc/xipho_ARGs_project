# !/usr/bin/env python
# encoding: utf-8

"""
Goal: Produce peak-specific manhattan plots.
________________________________________________________________________________
Copyright 2025 Andre E. Moncrieff. All rights reserved.

"""

import argparse
import pandas
from collections import Counter
import numpy
import statistics
import argparse
import pandas
from collections import Counter
import statistics
import matplotlib.pyplot as plt
from qmplot import manhattanplot
import numpy as np
from scipy.stats import spearmanr
import math
import seaborn as sns
from sklearn.preprocessing import normalize


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--masterFile", required=True,
                        help="Enter the file name (+ .txt)",
                        type=str)
    #parser.add_argument("--fakerowRTHfile", required=True,
    #                    help="Enter the file name (+ .txt)",
    #                    type=str)
    args = parser.parse_args()
    return args


def read_in_file(txt_file_dataframe):
    raw_dataframe = pandas.read_csv(txt_file_dataframe,
                                    sep='\t',
                                    encoding = "ISO-8859-1")
    return raw_dataframe


def add_combined_stats(clean_df):
    #print(clean_df)
    # Add the JCR stat to my df
    jc_xin_bel_list = clean_df['TMRCAH_Xin_Bel'].tolist()
    jc_tap_xin_list = clean_df['TMRCAH_Tap_Xin'].tolist()
    
    jcr_stat_list = []
    for i in range(len(jc_xin_bel_list)):
        # print(tapajos_enrich_list[i])
        # print(type(tapajos_enrich_list[i]))
        jc_xin_bel = jc_xin_bel_list[i]
        jc_tap_xin = jc_tap_xin_list[i]
        if jc_tap_xin == 0:
            jcr_stat_list.append("NA")
        else:
            jcr_stat = jc_xin_bel/jc_tap_xin
            jcr_stat_list.append(jcr_stat)

    clean_df['JCR'] = jcr_stat_list
    clean_master = clean_df[clean_df["JCR"].astype(str).str.contains("NA") == False]
    clean_master["JCR"] = pandas.to_numeric(clean_master["JCR"])
    return clean_master


def fst_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
  
    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)  

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    
    # Add a chrom number column--Chromosome 4
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    five_SD = 0.14933548121106754 #(for all autosomes)
    print(window_df_v2)
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='Fst_belem_xingu',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='Fst_belem_xingu',
    xlabel='Chromosome 4',
    ylabel='Fst',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=five_SD,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("1_Fst_belem_xingu_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def xin_enr_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    print(window_df_v2)
    #window_df_v2.to_csv('xin_enr_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 7.074 #(for all autosomes)
    #window_df_v2.to_csv('test.txt', sep='\t', index=False) 
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='xipho_xingu_enrich',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='xin_enr',
    xlabel='Chromosome 4',
    ylabel='xingu_enr',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("2_xin_enr_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def bel_enr_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    print(window_df_v2)

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 8.389 #(for all autosomes)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='xipho_belem_enrich',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='bel_enr',
    xlabel='Chromosome 4',
    ylabel='bel_enr',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("3_bel_enr_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def xin_rth_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    print(window_df_v2)
    #window_df_v2.to_csv('xin_rth_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 3.816793893 #(for all autosomes)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='xin_RTH_inverse',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='xin_RTH',
    xlabel='Chromosome 4',
    ylabel='xin_RTH',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("4_xin_RTH_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def bel_rth_manhattan_plot(master_df):
    # Runs with an error first dataset
    # But then once I add a fake point at end, it works
    # something to do with the axes
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  
    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    print(window_df_v2)
    window_df_v2.to_csv('bel_rth_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 12.12165736 #(for all autosomes)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='bel_RTH_inverse',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='bel_rth',
    xlabel='Chromosome 4',
    ylabel='bel_rth',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("5_bel_rth_manhattan_Chr4_fake_final_point.pdf", format="pdf")
    #return clean_master2


def xin_bel_rth_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)
    #window_df_v2.to_csv('xin_bel_rth_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 2.1265 #(for all autosomes)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='TMRCAH_Xin_Bel_corr2_divbyTMRCAHall',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='xin_bel_rth',
    xlabel='Chromosome 4',
    ylabel='xin_bel_rth',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("6_xin_bel_rth_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def tap_xin_rth_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)
    #window_df_v2.to_csv('tap_xin_rth_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 1.526717557 #(for all autosomes)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='Tap_Xin_RTH_inverse',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='tap_xin_rth',
    xlabel='Chromosome 4',
    ylabel='tap_xin_rth',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("7_tap_xin_rth_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def jcr_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)
    #window_df_v2.to_csv('jcr_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 2.159139382 #(for all autosomes)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='JCR',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='jcr',
    xlabel='Chromosome 4',
    ylabel='jcr',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("8_jcr_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def Q95_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)
    #window_df_v2.to_csv('Q95_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 0.5 #(for all autosomes)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='Q95_bel_popA',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='Q95',
    xlabel='Chromosome 4',
    ylabel='Q95',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("9_Q95_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def Dxy_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)
    #window_df_v2.to_csv('Dxy_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 0.12081057375200367 #(mean for control autosome windows, gap version)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='dxy_Bel_Xin',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='Dxy',
    xlabel='Chromosome 4',
    ylabel='Dxy',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("10_Dxy_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def Pi_Xingu_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)
    #window_df_v2.to_csv('Pi_Xingu_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 0.122543512 #(mean for control autosome windows, gap version)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='pi_xingu',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='Pi_Xingu',
    xlabel='Chromosome 4',
    ylabel='Pi_Xingu',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("11_Pi_Xingu_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def Pi_Belem_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)
    #window_df_v2.to_csv('Pi_Xingu_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 0.106954281 #(mean for control autosome windows, gap version)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='pi_belem',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='Pi_Belem',
    xlabel='Chromosome 4',
    ylabel='Pi_Belem',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("12_Pi_Belem_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def recomb50site_manhattan_plot(master_df):
    clean_master = master_df.dropna(subset=['Fst_belem_xingu']).copy()
    clean_master.loc[clean_master['Fst_belem_xingu']<0, 'Fst_belem_xingu']=0
    clean_master2 = clean_master[clean_master['sites_summ_stats'] >= 50].copy() # filtering out windows with fewer than 50 sites
    #clean_master3 = clean_master2.dropna(subset=['xipho_xingu_enrich']).copy 
  

    chrom_specific_df  = clean_master2[clean_master2['chrom'].isin(['Chromosome_4_RagTag'])]
    #chrom_specific_df = clean_master3[clean_master3['chrom'].isin(['Chromosome_4_RagTag'])]
    window_df_v1 = chrom_specific_df[chrom_specific_df['start'] >= 69750000].copy()
    window_df_v2 = window_df_v1[window_df_v1['end'] <= 71750000].copy()

    #print(window_df_v2)
    #window_df_v2.to_csv('Pi_Xingu_workaround_df.txt', sep='\t', index=False) 

    start_list = window_df_v2['start'].tolist()
    end_list = window_df_v2['end'].tolist()

    position_list = []
    for i in range(len(start_list)):
        pos_val = (start_list[i] + end_list[i])/2
        position_list.append(pos_val)

    workaround_position_list = []
    for num in position_list:
        new_val = num - start_list[0]
        workaround_position_list.append(new_val)
    #print(workaround_position_list)

    window_df_v2['workaround_position_list'] = workaround_position_list

    # Add a chrom number column
    chrom_num = []
    for i in range(len(start_list)):
        chrom_val = int(4)
        chrom_num.append(chrom_val)
    window_df_v2['chrom_num'] = chrom_num

    pval = 2.88E-09 #(mean for control autosome windows, gap version, 92975 windows)
    
    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(
    window_df_v2,
    #chrom='final_chr_name',
    chrom='chrom_num',
    pos='workaround_position_list',
    #pos='position_list',
    pv='recombination_rate_50site_minimum',
    snp='ID',
    logp=False,
    ax=ax,
    marker='.',
    #color='#3B5488,#53BBD5',
    color='#000000,#000000',
    alpha=0.9,
    title='recomb50site',
    xlabel='Chromosome 4',
    ylabel='recomb50site',
    xtick_label_set=None,
    CHR="4",
    xticklabel_kws=None,
    suggestiveline=None,
    genomewideline=pval,
    sign_line_cols='#000000,#000000',
    hline_kws={"linestyle": "--", "lw": 1.3},
    #hline_kws=None,
    sign_marker_p=None,
    sign_marker_color='r',
    is_annotate_topsnp=False,
    highlight_other_SNPs_indcs=None,
    highlight_other_SNPs_color='r',
    highlight_other_SNPs_kwargs=None,
    text_kws=None,
    ld_block_size=50000)
    #plt.savefig("Fst_belem_xingu_manhattan_Chr4.png", dpi=600)
    plt.savefig("13_recomb50site_manhattan_Chr4.pdf", format="pdf")
    #return clean_master2


def main():
    #create args object
    args = parser()
    proto_master_df = read_in_file(args.masterFile)
    master_df = add_combined_stats(proto_master_df)
    print(master_df)
    
    #plots
    #fst_manhattan_plot(master_df)
    #xin_enr_manhattan_plot(master_df)
    #bel_enr_manhattan_plot(master_df)
    #xin_rth_manhattan_plot(master_df)
    #bel_rth_manhattan_plot(master_df)
    #xin_bel_rth_manhattan_plot(master_df)
    #tap_xin_rth_manhattan_plot(master_df)
    #jcr_manhattan_plot(master_df)
    #Q95_manhattan_plot(master_df)
    #Dxy_manhattan_plot(master_df)
    #Pi_Xingu_manhattan_plot(master_df)
    #Pi_Belem_manhattan_plot(master_df)
    #recomb50site_manhattan_plot(master_df)

    
if __name__ == '__main__':
    main()
