#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import math
import copy
import re

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scipy.stats import variation


# In[2]:


def split(obj, n):
    res = []
    if n == 0:
        res.append(obj)
        return res
    pos = 0.0
    length = float(len(obj)) / float(n)
    while pos < len(obj):
        res.append(obj[int(pos):int(pos + length)])
        pos += length
    return res


# In[3]:


def max_prefix(strs):
    if len(strs) == 0: return ""
    res = strs[0]
    for i in range(1, len(strs)):
        while(strs[i].startswith(res) == False):
            res = res[:-1]
            if len(res) == 0: return ""
    return res

def max_suffix(strs):
    rev = [s[::-1] for s in strs]
    return max_prefix(rev)

def remove_common(strs):
    res = copy.deepcopy(strs)
    prefix = len(max_prefix(res))
    suffix = len(max_suffix(res))
    for i in range(len(res)): res[i] = res[i][prefix:-suffix]
    return res


# In[4]:


def add_labels(bars, al = 1):
    max_height = 0
    if al == 1: als = 'center'
    else:
        if al == 0: als = 'left'
        else: als = 'right'
    for b in bars: 
        if b.get_height() > max_height: 
            max_height = b.get_height()
    for b in bars:
        height = b.get_height()
        if height > 0.3 * max_height: plt.text(b.get_x() + b.get_width()/2.0, height - max_height * 0.01, ('%.2f' % height).rstrip('0').rstrip('.'), ha=als, va='top', rotation='vertical')
        else: plt.text(b.get_x() + b.get_width()/2.0, height + max_height * 0.01, ('%.2f' % height).rstrip('0').rstrip('.'), ha=als, va='bottom', rotation='vertical') 
        
def bar_plot(title, x, y, axis = True, lab = False):
    N = math.floor(len(x) / 40)
    if N == 0: N = 1
    xl = split(x, N)
    yl = split(y, N)
    if (len(xl) >= 2): axis = True
    number = len(xl)
    index = 0
    if len(x) > 20: f = plt.figure(figsize = (20,5 * number))
    else: f = plt.figure()
    for i in range(0, len(xl)):
        index = index + 1
        s = f.add_subplot(number, 1, index)
        p = plt.bar(xl[i], yl[i], color = 'grey', alpha=0.5, edgecolor = 'black')
        if lab: add_labels(p)
        plt.title(title)
        if (axis == False): s.set_xticklabels(['']*len(xl[i]))
        else: s.set_xticklabels(xl[i], rotation = 45, ha = 'right')
        plt.tight_layout()
        
def double_bar_plot(title, x, y, z, legend_y, legend_z, axis = True, lab = False):
    N = math.floor(len(x) / 40)
    if N == 0: N = 1
    xl = split(x, N)
    yl = split(y, N)
    zl = split(z, N)
    if (len(xl) >= 2): axis = True
    number = len(xl)
    if len(x) > 20: f = plt.figure(figsize = (20,5 * number))
    else: f = plt.figure()
    for i in range(0, len(xl)):
        s = f.add_subplot(number, 1, i + 1)
        p1 = plt.bar(xl[i], yl[i], color = 'grey', alpha=0.5, edgecolor = 'black')
        p2 = plt.bar(xl[i], zl[i], color = 'firebrick', alpha=0.5, edgecolor = 'black')
        if lab: 
            add_labels(p1,0)
            add_labels(p2,2)
        plt.legend((p1[0],p2[0]), (legend_y, legend_z))
        plt.title(title)
        if (axis == False): s.set_xticklabels(['']*len(xl[i]))
        else: s.set_xticklabels(xl[i], rotation = 45, ha = 'right')
        plt.tight_layout()
        
def triple_bar_plot(title, x, y, z, u, legend_y, legend_z, legend_u, axis = True, lab = False):
    N = math.floor(len(x) / 40)
    if N == 0: N = 1
    xl = split(x, N)
    yl = split(y, N)
    zl = split(z, N)
    ul = split(u, N)
    if (len(xl) >= 2): axis = True
    number = len(xl)
    if len(x) > 20: f = plt.figure(figsize = (20,10 * number))
    else: f = plt.figure()
    for i in range(0, len(xl)):
        s = f.add_subplot(number, 1, i + 1)
        p1 = plt.bar(xl[i], yl[i], color = 'grey', alpha=0.5, edgecolor = 'black')
        p2 = plt.bar(xl[i], zl[i], color = 'firebrick', alpha=0.5, edgecolor = 'black')
        p3 = plt.bar(xl[i], ul[i], color = 'gold', alpha=0.5, edgecolor = 'black')
        if lab: 
            add_labels(p1,0)
            add_labels(p2,2)
            add_labels(p3,0)
        plt.legend((p1[0],p2[0],p3[0]), (legend_y, legend_z, legend_u))
        plt.title(title)
        if (axis == False): s.set_xticklabels(['']*len(xl[i]))
        else: s.set_xticklabels(xl[i], rotation = 45, ha = 'right')
        plt.tight_layout()

def corr_plot(x):
    lg = x.applymap(np.log2)
    f = plt.figure(figsize = (20.0, 20.0))
    plt.tight_layout()
    plt.matshow(lg.corr(), cmap=plt.cm.Reds)
    fsize = min(10, 150 / len(lg.columns))
    plt.gca().tick_params(width=min(1, 10.0 / len(lg.columns)))
    if len(lg.columns) <= 50: plt.xticks(range(lg.shape[1]), lg.columns, rotation=45, fontsize = fsize)
    else: plt.xticks(range(lg.shape[1]), lg.columns, rotation=90, fontsize = fsize)
    plt.yticks(range(lg.shape[1]), lg.columns, fontsize = fsize)
    cb = plt.colorbar()


# In[5]:


def report(stats, main, out):
    print("Generating report. Stats, main report, output file: " + stats + ", " + main + ", " + out)
    df = pd.read_csv(stats,sep='\t')
    df.loc[:,'File.Name'] = remove_common(df['File.Name'])
    
    quant = pd.read_csv(main,sep='\t')
    quant = quant[quant['Q.Value'] <= 0.01].reset_index(drop = True)
    quant.loc[:,'File.Name'] = remove_common(quant['File.Name'])
    quant_pg = quant[['File.Name','Protein.Group','PG.Normalised']].drop_duplicates().reset_index(drop = True)
    quant_gene = quant[quant['Proteotypic'] != 0][['File.Name','Genes','Gene.Normalised.Unique']].drop_duplicates().reset_index(drop = True)
    
    pr = quant.pivot_table(values = ['Precursor.Normalised'], columns = ['File.Name'], index = ['Precursor.Id'], aggfunc=sum)
    pg = quant_pg.pivot_table(values = ['PG.Normalised'], columns = ['File.Name'], index = ['Protein.Group'], aggfunc=sum)
    genes = quant_gene.pivot_table(values = ['Gene.Normalised.Unique'], columns = ['File.Name'], index = ['Genes'], aggfunc=sum)
    pr.columns, pg.columns, genes.columns = pr.columns.droplevel(), pg.columns.droplevel(), genes.columns.droplevel()
    totals = pr.sum(axis = 0, skipna = True)
    
    fnames = len(pg.columns)
    pr_ids = fnames - pr.count(axis = 1)
    pg_ids = fnames - pg.count(axis = 1)
    gene_ids = fnames - genes.count(axis = 1)
    
    skip_conditions = False
    try:
        df['Replicate'] = [re.findall(r'\d+', file)[-1] for file in df['File.Name']]
        df['Condition'] = [file[:file.rfind(re.findall(r'\d+', file)[-1])] + file[file.rfind(re.findall(r'\d+', file)[-1]) + len(re.findall(r'\d+', file)[-1]):] for file in df['File.Name']]
        conditions = df['Condition'].unique()
        df['Precursor.CV'] = 0
        df['Precursor.CV.20'] = 0
        df['Precursor.CV.10'] = 0
        df['PG.CV'] = 0
        df['PG.CV.20'] = 0
        df['PG.CV.10'] = 0
        df['Gene.CV'] = 0
        df['Gene.CV.20'] = 0
        df['Gene.CV.10'] = 0
        df['Precursor.N'] = 0
        df['PG.N'] = 0
        df['Gene.N'] = 0
        for condition in conditions:
            files = df['File.Name'][df['Condition'] == condition]
            if len([c for c in pr.columns if c in list(files)]) > 0:
                pr_cvs = np.ma.filled(variation(pr[files],axis = 1,nan_policy='omit'), float('nan'))
                pr_cvs[pr_cvs == 0] = float('nan')
                df.loc[df['Condition'] == condition,'Precursor.CV'] = np.nanmedian(pr_cvs)
                df.loc[df['Condition'] == condition,'Precursor.CV.20'] = len(pr_cvs[pr_cvs <= 0.2])
                df.loc[df['Condition'] == condition,'Precursor.CV.10'] = len(pr_cvs[pr_cvs <= 0.1])
                df.loc[df['Condition'] == condition,'Precursor.N'] = np.mean(pr[files].count()).astype(int)
            if len([c for c in pg.columns if c in list(files)]) > 0:
                pg_cvs = np.ma.filled(variation(pg[files],axis = 1,nan_policy='omit'), float('nan'))
                pg_cvs[pg_cvs == 0] = float('nan')
                df.loc[df['Condition'] == condition,'PG.CV'] = np.nanmedian(pg_cvs)
                df.loc[df['Condition'] == condition,'PG.CV.20'] = len(pg_cvs[pg_cvs <= 0.2])
                df.loc[df['Condition'] == condition,'PG.CV.10'] = len(pg_cvs[pg_cvs <= 0.1])
                df.loc[df['Condition'] == condition,'PG.N'] = np.mean(pg[files].count()).astype(int)
            if len([c for c in genes.columns if c in list(files)]) > 0:
                gene_cvs = np.ma.filled(variation(genes[files],axis = 1,nan_policy='omit'), float('nan'))
                gene_cvs[gene_cvs == 0] = float('nan')
                df.loc[df['Condition'] == condition,'Gene.CV'] = np.nanmedian(gene_cvs)
                df.loc[df['Condition'] == condition,'Gene.CV.20'] = len(gene_cvs[gene_cvs <= 0.2])
                df.loc[df['Condition'] == condition,'Gene.CV.10'] = len(gene_cvs[gene_cvs <= 0.1])
                df.loc[df['Condition'] == condition,'Gene.N'] = np.mean(genes[files].count()).astype(int)
    except: 
        print("Cannot infer conditions/replicates")
        skip_conditions = True

    
    with PdfPages(out) as pdf:        
        try:
            f = plt.figure(figsize = (20, 5))
            f.add_subplot(131)
            p1 = plt.hist(pr_ids, bins = fnames, cumulative = True, histtype = 'stepfilled', color = 'grey', alpha = 0.5, edgecolor = 'black')
            plt.title("Identification consistency: precursors, CDF")
            plt.xlabel("Missing values")
            plt.ylabel("IDs")
            f.add_subplot(132)
            p1 = plt.hist(pg_ids, bins = fnames, cumulative = True, histtype = 'stepfilled', color = 'grey', alpha = 0.5, edgecolor = 'black')
            plt.title("Identification consistency: protein groups, CDF")
            plt.xlabel("Missing values")
            plt.ylabel("IDs")
            f.add_subplot(133)
            if len(genes) > 0:
                p1 = plt.hist(gene_ids, bins = fnames, cumulative = True, histtype = 'stepfilled', color = 'grey', alpha = 0.5, edgecolor = 'black')
                plt.title("Identification consistency: unique genes, CDF")
                plt.xlabel("Missing values")
                plt.ylabel("IDs")
                pdf.savefig()
        except: pass

        try:
            f = plt.figure(figsize = (20, 20.0/3.0))
            f.add_subplot(131)
            hmp, ex, ey = np.histogram2d(quant['iRT'], quant['RT'], bins=250)
            plt.imshow(hmp.T, origin='lower', cmap = 'binary', extent = [ex[0], ex[-1], ey[0], ey[-1]], aspect = 'auto',interpolation='none')
            plt.title("Retention times heatmap, all runs")
            plt.xlabel("Library iRT")
            plt.ylabel("RT")
            f.add_subplot(132)
            hmp, ex, ey = np.histogram2d(quant['Predicted.RT'], quant['RT'], bins=250)
            plt.imshow(hmp.T, origin='lower', cmap = 'binary', extent = [ex[0], ex[-1], ey[0], ey[-1]], aspect = 'auto',interpolation='none')
            plt.title("Retention time accuracy heatmap, all runs")
            plt.xlabel("Predicted RT")
            plt.ylabel("RT")
            f.add_subplot(133)
            ratio = quant['Precursor.Normalised'] / quant['Precursor.Quantity']
            hmp, ex, ey = np.histogram2d(quant['RT'][ratio > 0], ratio[ratio > 0], bins=250)
            plt.imshow(hmp.T, origin='lower', cmap = 'binary', extent = [ex[0], ex[-1], ey[0], ey[-1]], aspect = 'auto',interpolation='none')
            plt.title("Normalisation factor heatmap, all runs")
            plt.xlabel("RT")
            plt.ylabel("Normalisation factor")
            pdf.savefig()
        except: pass
        
        try:
            corr_plot(pg)
            pdf.savefig()
        except: pass
        
        try:
            bar_plot("Total quantity, 1% FDR", df['File.Name'], df['Total.Quantity'])
            pdf.savefig()
            bar_plot("MS1 signal", df['File.Name'], df['MS1.Signal'])
            pdf.savefig()
            bar_plot("MS2 signal", df['File.Name'], df['MS2.Signal'])
            pdf.savefig()
            r = [x/y if y > 0 else 0 for x,y in zip(df['Total.Quantity'], df['MS2.Signal'])]
            bar_plot("Total quantity/MS2 signal ratio", df['File.Name'], r, lab = True)
            pdf.savefig()
            r = [x/y if y > 0 else 0 for x,y in zip(df['MS1.Signal'], df['MS2.Signal'])]
            bar_plot("MS1/MS2 signal ratio", df['File.Name'], r, lab = True)
            pdf.savefig()
            bar_plot("Precursors, 1% FDR", df['File.Name'], df['Precursors.Identified'], lab = True)
            pdf.savefig()
            if max(df['Proteins.Identified']) > 0:
                bar_plot("Unique proteins, 1% protein-level FDR", df['File.Name'], df['Proteins.Identified'], lab = True)
                pdf.savefig()
            bar_plot("Mean peak FWHM, in minutes", df['File.Name'], df['FWHM.RT'], lab = True)
            pdf.savefig()
            bar_plot("Mean peak FWHM, in MS2 scans", df['File.Name'], df['FWHM.Scans'], lab = True)
            pdf.savefig()
            bar_plot("Median RT prediction accuracy, minutes", df['File.Name'], df['Median.RT.Prediction.Acc'], lab = True)
            pdf.savefig()
            double_bar_plot("Median mass accuracy, MS2, ppm", df['File.Name'], df['Median.Mass.Acc.MS2'], df['Median.Mass.Acc.MS2.Corrected'], "Without correction","Corrected")
            pdf.savefig()
            double_bar_plot("Median mass accuracy, MS1, ppm", df['File.Name'], df['Median.Mass.Acc.MS1'], df['Median.Mass.Acc.MS1.Corrected'], "Without correction","Corrected")
            pdf.savefig()
            double_bar_plot("Peptide characteristics", df['File.Name'], df['Average.Peptide.Length'], df['Average.Peptide.Charge'], "Average length","Average charge")
            pdf.savefig()
            bar_plot("Average missed tryptic cleavages", df['File.Name'], df['Average.Missed.Tryptic.Cleavages'], lab = True)
            pdf.savefig()
        except: pass
        
        if skip_conditions == False:
            try:
                cvs = copy.deepcopy(df)[['Condition','Precursor.N','Precursor.CV','Precursor.CV.20','Precursor.CV.10','PG.N','PG.CV','PG.CV.20','PG.CV.10','Gene.N','Gene.CV','Gene.CV.20','Gene.CV.10']].drop_duplicates()
                cvs = cvs[cvs['Precursor.N'] > 0]
                cvs = cvs[cvs['Precursor.CV'] > 0.0]
                if len(cvs) > 0:
                    triple_bar_plot("Precursors, 1% FDR", cvs['Condition'], cvs['Precursor.N'], cvs['Precursor.CV.20'], cvs['Precursor.CV.10'], "Average", "CV < 20%", "CV < 10%", lab = True)
                    pdf.savefig()
                    bar_plot("Median precursor CV, 1% FDR", cvs['Condition'], cvs['Precursor.CV'], lab = True)
                    pdf.savefig()
                    cvs = cvs[cvs['PG.N'] > 0]
                    cvs = cvs[cvs['PG.CV'] > 0.0]
                    if len(cvs) > 0:
                        triple_bar_plot("Protein groups, 1% precursor-level FDR", cvs['Condition'], cvs['PG.N'], cvs['PG.CV.20'], cvs['PG.CV.10'], "Average", "CV < 20%", "CV < 10%", lab = True)
                        pdf.savefig()
                        bar_plot("Median protein group CV, 1% precursor-level FDR", cvs['Condition'], cvs['PG.CV'], lab = True)
                        pdf.savefig()
                        cvs = cvs[cvs['Gene.N'] > 0]
                        cvs = cvs[cvs['Gene.CV'] > 0.0]
                        if len(cvs) > 0:
                            triple_bar_plot("Unqiue genes, 1% precursor-level FDR", cvs['Condition'], cvs['Gene.N'], cvs['Gene.CV.20'], cvs['Gene.CV.10'], "Average", "CV < 20%", "CV < 10%", lab = True)
                            pdf.savefig()
                            bar_plot("Median gene CV, 1% precursor-level FDR", cvs['Condition'], cvs['Gene.CV'], lab = True)
                            pdf.savefig()
            except: pass


# In[ ]:


report(sys.argv[1], sys.argv[2], sys.argv[3])

