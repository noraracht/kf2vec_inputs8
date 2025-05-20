# Usage: python closeast2taxid.py /Users/nora/Documents/ml_metagenomics/cami_long_reads/tax_assignment /Users/nora/Documents/ml_metagenomics/65K_TOL/assembly_summary_refseq-after20190108-sampled_500_lineage.txt


import argparse
import os
import pandas as pd
import numpy as np
from pathlib import Path
import fnmatch
from collections import Counter
from scipy.stats import poisson

# Parse inputs
parser = argparse.ArgumentParser()
parser.add_argument('input_dir')
parser.add_argument('name_to_id')
parser.add_argument('length_file')

args = parser.parse_args()

input_dir  = args.input_dir
name_to_id = args.name_to_id
length_file = args.length_file

#print( f'Hello, {dist_mtrx } { name_to_id }!' )

formats = ['.csv']
files_names = [f for f in os.listdir(input_dir)  if True in (fnmatch.fnmatch(f, '*' + form) for form in formats)]
samples_names = [f.rsplit('.f', 1)[0] for f in files_names]

# Create species to taxid map
df_taxid = pd.read_csv(name_to_id, index_col=None, header=0, sep='\t')
df_length = pd.read_csv(length_file, index_col=None, header=None, sep='\t')

df_length[1] =df_length[1].astype(int)
my_map_dict_lens = pd.Series(df_length[1].values, index=df_length[0]).to_dict()

my_map_dict_sp = pd.Series(df_taxid['species'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_gen = pd.Series(df_taxid['genus'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_fam = pd.Series(df_taxid['family'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_ordr = pd.Series(df_taxid['order'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_cls = pd.Series(df_taxid['class'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_phyl = pd.Series(df_taxid['phylum'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_spkindm = pd.Series(df_taxid['superkingdom'].values, index=df_taxid["tree_label"]).to_dict()


# Get folder name
folder_name = os.path.basename(os.path.dirname(input_dir))
#print (folder_name)

# Create output files
rank_fname1 = "{}_v10_taxid.species".format("summary_subset")
rank_fname1_log = "{}_v10.log".format("summary_subset")

#rank_fname2 = "{}_taxid.genus".format("summary")
#rank_fname3 = "{}_taxid.family".format("summary")
#rank_fname4 = "{}_taxid.order".format("summary")
#rank_fname5 = "{}_taxid.class".format("summary")
#rank_fname6 = "{}_taxid.phylum".format("summary")

if os.path.isfile(os.path.join(input_dir, rank_fname1)):
    os.remove(os.path.join(input_dir, rank_fname1))
if os.path.isfile(os.path.join(input_dir, rank_fname1_log)):
    os.remove(os.path.join(input_dir, rank_fname1_log))

#if os.path.isfile(os.path.join(input_dir, rank_fname2)):
#    os.remove(os.path.join(input_dir, rank_fname2))
#if os.path.isfile(os.path.join(input_dir, rank_fname3)):
#    os.remove(os.path.join(input_dir, rank_fname3))
#if os.path.isfile(os.path.join(input_dir, rank_fname4)):
#    os.remove(os.path.join(input_dir, rank_fname4))
#if os.path.isfile(os.path.join(input_dir, rank_fname5)):
#    os.remove(os.path.join(input_dir, rank_fname5))
#if os.path.isfile(os.path.join(input_dir, rank_fname6)):
#    os.remove(os.path.join(input_dir, rank_fname6))


fo_sp = open(os.path.join(input_dir, rank_fname1), "a")
fo_sp_log = open(os.path.join(input_dir, rank_fname1_log), "a")

#fo_gen = open(os.path.join(input_dir, rank_fname2), "a")
#fo_fam = open(os.path.join(input_dir, rank_fname3), "a")
#fo_ordr = open(os.path.join(input_dir, rank_fname4), "a")
#fo_cls = open(os.path.join(input_dir, rank_fname5), "a")
#fo_phyl = open(os.path.join(input_dir, rank_fname6), "a")

fo_sp.write("@Version:0.9.0\n") 
fo_sp.write("@SampleID:marmgCAMI2_short_read_pooled_gold_standard_assembly\n") 
fo_sp.write("@@SEQUENCEID\tTAXID\n")

"""
fo_sp.write("@Version:0.9.0\n")
fo_sp.write("@SampleID:marmgCAMI2_short_read_pooled_gold_standard_assembly\n")
fo_sp.write("@@SEQUENCEID\tTAXID\n")

fo_sp.write("@Version:0.9.0\n")
fo_sp.write("@SampleID:marmgCAMI2_short_read_pooled_gold_standard_assembly\n")
fo_sp.write("@@SEQUENCEID\tTAXID\n")

fo_sp.write("@Version:0.9.0\n")
fo_sp.write("@SampleID:marmgCAMI2_short_read_pooled_gold_standard_assembly\n")
fo_sp.write("@@SEQUENCEID\tTAXID\n")

fo_sp.write("@Version:0.9.0\n")
fo_sp.write("@SampleID:marmgCAMI2_short_read_pooled_gold_standard_assembly\n")
fo_sp.write("@@SEQUENCEID\tTAXID\n")

fo_sp.write("@Version:0.9.0\n")
fo_sp.write("@SampleID:marmgCAMI2_short_read_pooled_gold_standard_assembly\n")
fo_sp.write("@@SEQUENCEID\tTAXID\n")
"""


#lens = eval(open('len.txt', 'r').read())
#print(lens)

for dist_mtrx in files_names:

    # Get the base name (myfile.txt)
    basename = os.path.basename(dist_mtrx)	
    #print (basename)

    # Get the directory path (/home/user/documents)
    dirname = os.path.dirname(dist_mtrx)

    #head_tail = os.path.split(dist_mtrx)
    #print(head_tail)


    # Read distance matrix
    df = pd.read_csv(os.path.join(input_dir, dist_mtrx), index_col=0, header=0, sep='\t')
    #print (df)

    out = df[df.iloc[0].sort_values(ascending=True).index]
    #print(out.head(n=2))
    """
    my_row = out.iloc[0].to_numpy()
    print (my_row[0])

    poiss_cut = np.array([1-poisson.cdf(q, my_row[0]) for q in my_row])
    print (poiss_cut)
    next_element = my_row[1:]
    print (next_element)
    curr_element = my_row[:-1]
    print (curr_element)
    diff = next_element-curr_element
    print (diff)
    with np.errstate(divide='ignore'):
    	grad = diff[1:]/diff[:-1]

    condition = poiss_cut >= 0.1

    index = np.argmax(condition)
    #if index > 0:
    #    sliced_array =  poiss_cut[index:]
    #print (sliced_array)
    """
    
    column_subset = []

    maxdelt = -1
    p = 0
    i = 0
    first = -1
    selected=5000
    query_label = basename.replace("apples_input_di_mtrx_query_", "").replace(".csv", "")
    #clen = lens[query_label]
    clen = my_map_dict_lens[query_label]
    #print("d","maxdelta","delta","delta/maxdelta","pvalue")
    for l in range (0, len(out.columns)):
        a=float(out.iat[0, l])
        if first == -1:
            first = a
        i = i + 1
        print(a,clen,maxdelt,a-p,(a-p)/maxdelt,a/(p+0.0000000001),1-poisson.cdf(a/100*clen,first/100*clen)) #,1-norm.cdf(a/100,first/100,1/sqrt(l)))
        if (maxdelt != -1 and (a-p)/maxdelt>3) or (maxdelt == -1 and first != a and a/p > 1.05) or 1-poisson.cdf(a*clen/100,first*clen/100)<0.1:
            if selected == 5000:
                selected = i
            #print("****")
            break
        if first !=a and a-p>maxdelt:
            maxdelt=a-p
        p=a

    column_subset = out.columns[:i-1] 
    #print (out)
 

    # Get closeast species
    #n1 = df.idxmin(axis=1).to_list()[0]
    #print (n1)
    #n1_root = "_".join(n1.split("_", 2)[:2])
    #print (n1_root)

    column_subset_sp = Counter([my_map_dict_sp[n1] for n1 in column_subset if my_map_dict_sp[n1] != -1])
    column_subset_gen = Counter([my_map_dict_gen[n1] for n1 in column_subset if my_map_dict_gen[n1] != -1])
    column_subset_fam = Counter([my_map_dict_fam[n1] for n1 in column_subset if my_map_dict_fam[n1] != -1])
    column_subset_ordr = Counter([my_map_dict_ordr[n1] for n1 in column_subset if my_map_dict_ordr[n1] != -1])
    column_subset_cls = Counter([my_map_dict_cls[n1] for n1 in column_subset if my_map_dict_cls[n1] != -1])
    column_subset_phyl = Counter([my_map_dict_phyl[n1] for n1 in column_subset if my_map_dict_phyl[n1] != -1])
    column_subset_spkindm = Counter([my_map_dict_spkindm[n1] for n1 in column_subset if my_map_dict_spkindm[n1] != -1])


    # change threshold below to 0.66 or 0.75
    # change 0.2 above to 0.1 or 0.05
    maj_thresh = 0.75
    print(column_subset_sp,column_subset_gen,column_subset_fam)
    if len(column_subset_sp)>0 and column_subset_sp.most_common(1)[0][1]/sum(column_subset_sp.values()) > maj_thresh:
        s=column_subset_sp.most_common(1)[0][0]
        print ("sp: ", s)
        fo_sp.write("{}\t{}\n".format(query_label, s))
        fo_sp_log.write("{}\t{}\t{}\t{}\n".format(query_label, i-1, "sp", s))
    elif len(column_subset_gen)>0 and column_subset_gen.most_common(1)[0][1]/sum(column_subset_gen.values()) > maj_thresh:
        s=column_subset_gen.most_common(1)[0][0]
        print ("gen: ", s)
        fo_sp.write("{}\t{}\n".format(query_label, s))
        fo_sp_log.write("{}\t{}\t{}\t{}\n".format(query_label, i-1, "gen", s))
    elif len(column_subset_fam)>0 and column_subset_fam.most_common(1)[0][1]/sum(column_subset_fam.values()) > maj_thresh:
        s=column_subset_fam.most_common(1)[0][0]
        print ("fam: ", s)
        fo_sp.write("{}\t{}\n".format(query_label, s))
        fo_sp_log.write("{}\t{}\t{}\t{}\n".format(query_label, i-1, "fam", s))
    elif len(column_subset_ordr)>0 and column_subset_ordr.most_common(1)[0][1]/sum(column_subset_ordr.values()) > maj_thresh:
        s=column_subset_ordr.most_common(1)[0][0]
        print ("ordr: ", s)
        fo_sp.write("{}\t{}\n".format(query_label, s))
        fo_sp_log.write("{}\t{}\t{}\t{}\n".format(query_label, i-1, "ordr", s))
    elif len(column_subset_cls)>0 and column_subset_cls.most_common(1)[0][1]/sum(column_subset_cls.values()) > maj_thresh:
        s=column_subset_cls.most_common(1)[0][0]
        print ("cls: ", s)
        fo_sp.write("{}\t{}\n".format(query_label, s))
        fo_sp_log.write("{}\t{}\t{}\t{}\n".format(query_label, i-1, "cls", s))
    elif len(column_subset_phyl)>0 and column_subset_phyl.most_common(1)[0][1]/sum(column_subset_phyl.values()) > maj_thresh:
        s=column_subset_phyl.most_common(1)[0][0]
        print ("phyl: ", s)
        fo_sp.write("{}\t{}\n".format(query_label, s))
        fo_sp_log.write("{}\t{}\t{}\t{}\n".format(query_label, i-1, "phyl", s))
    elif len(column_subset_spkindm)>0 and column_subset_spkindm.most_common(1)[0][1]/sum(column_subset_spkindm.values()) > maj_thresh:
        s=column_subset_spkindm.most_common(1)[0][0]
        print ("spkindm: ", s)
        fo_sp.write("{}\t{}\n".format(query_label, s))
        fo_sp_log.write("{}\t{}\t{}\t{}\n".format(query_label, i-1, "spkindm", s))
    else:
        print ("root: 1")
        fo_sp.write("{}\t1\n".format(query_label))
        fo_sp_log.write("{}\t{}\troot\t1\n".format(query_label, i-1))


    # Get taxid
    #print (my_map_dict[n1_root])

    #query_label = basename.replace("apples_input_di_mtrx_query_", "").replace(".csv", "")
    #print (query_label)
 
    """
    if my_map_dict_sp[n1] != -1:
        fo_sp.write("{}\t{}\n".format(query_label, my_map_dict_sp[n1]))
    if my_map_dict_gen[n1] != -1:
        fo_sp.write("{}\t{}\n".format(query_label, my_map_dict_gen[n1]))
    if my_map_dict_fam[n1] != -1:    
        fo_sp.write("{}\t{}\n".format(query_label, my_map_dict_fam[n1]))
    if my_map_dict_ordr[n1] != -1:    
        fo_sp.write("{}\t{}\n".format(query_label, my_map_dict_ordr[n1]))
    if my_map_dict_cls[n1] != -1:    
        fo_sp.write("{}\t{}\n".format(query_label, my_map_dict_cls[n1]))
    if my_map_dict_phyl[n1] != -1:    
        fo_sp.write("{}\t{}\n".format(query_label, my_map_dict_phyl[n1]))
    """

    #print(query_label, my_map_dict_sp[n1])


fo_sp.close()
fo_sp_log.close()

#fo_gen.close()
#fo_fam.close()
#fo_ordr.close()
#fo_cls.close()
#fo_phyl.close()
