#!/usr/bin/env python2
## Python script to assign NCBI taxonomy to Viral Clusters (VCs)
##(c) Martin Jahn
# dependencies: ete3, warnings, pandas, os

#################
# FUNCTION1: get taxid from NCBI accession. Note: Online queries run faster when America sleeps
#################
import os
import warnings
def accession_to_taxid(ACC):  
  warnings.filterwarnings("ignore")
  global taxid
  cmd = "curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="+ACC+"&rettype=fasta&retmode=xml' | grep TSeq_taxid | cut -d '>' -f 2 | cut -d '<' -f 1"
  accession_to_taxid.taxid = os.popen(cmd).read().strip()

#################
# FUNCTION2: get lineage from taxid
#################
import warnings
from ete3 import NCBITaxa
ncbi = NCBITaxa() ## load taxonomy  
def get_lineage(taxid , desired_ranks): 
  warnings.filterwarnings("ignore")
  global des_lin
  if taxid != "":   
    try:
      lineage = ncbi.get_lineage(int(taxid)) # get lineage tax of all tax. levels
      names = ncbi.get_taxid_translator(lineage)
      des_taxid = [k for k,v in ncbi.get_rank(names).items() if str(v) == desired_ranks] # get taxid of desired rank
      des_lin = ncbi.get_taxid_translator(des_taxid) # get lineage of desired rank
    except ValueError:
     pass

#################
# FUNCION3: calculate the proportion of the abundant cluster member 
#################
from collections import Counter
def prop_of_mayority(lista):
  global major_perc
  global major_id
  c= Counter(lista)
  major = [(i, float(c[i]) / float(len(lista)) * 100.0) for i, count in c.most_common(1)] # get percentage of most abundant list element
  major_perc = major[0][1]
  major_id = major[0][0]

  
#################
# LOOPY core: screen dataframe line by line and assemble taxonomic information
#################

import pandas as pd
import warnings

fii = open("clusters.genome.mci.I14","r") # read cluster file
"""Every enty needs to be prefixed with a database ID. Here samples of interest start with BCvir_ and ViralRefseq startswith RefSeq_  """
"""see file 'clusters.genome.mci.I14' as reference for how to format, one line per cluster """

req_rank = 'genus' # specify taxonomic rank to use
df = pd.DataFrame() # init cluster stat dataframe
df_memb = pd.DataFrame() # init member stat dataframe
cluster_counter = 0
member_counter = 0
unknown=''
major_name=''
yyy=''

for row in fii:
  warnings.filterwarnings("ignore")
  cluster_counter += 1
  print (("Begin to process cluster   %s") %cluster_counter)
  pool = []
  values = row.strip().split('\t')  # load cluster  
  nsample = len([i for i in values if i.startswith('BCvir_')]) # get sample number in cluster
  psample = float(nsample)/float(len(values))# Get proportion of sample members in cluster
  nrefseq = len([i for i in values if i.startswith('RefSeq_')]) # get refseq count in cluster  
  for memb in values:  # loop cluster members  
    member_counter +=1 # counter cluster members       
    df_memb = df_memb.append({'cluster_id': cluster_counter , 'member_id': memb},ignore_index=True)
    if nrefseq > 0: # any member in cluster is NCBI Refseq 
      unknown = False
    else:
      unknown = True # set cluster unknown    	
    if memb.startswith('RefSeq_'): # specific member is Refseq
      tmp = memb.replace("RefSeq_", "")
      accession_to_taxid(tmp) # retrieve taxid of entry
      get_lineage( (accession_to_taxid.taxid), req_rank) # get lineage of given rank
      pool.extend(des_lin)# get all cluster NCBI refseq of given rank to a list 
      tmp = ''	  
  lpool = len(pool)
  if unknown == True:
    print("cluster %s is new with no reference entry" % (cluster_counter))
    df = df.append({'cluster_id': cluster_counter , 'cluster_size': len(values), 'cluster_tax': "no_refseq",'majority_%' : "" , 'sample_members' : nsample , 'sample_proportion': psample}, ignore_index=True)    
  elif (unknown == False and lpool == 0):
    df = df.append({'cluster_id': cluster_counter , 'cluster_size': len(values), 'cluster_tax': "no_refseq_with_rank",'majority_%' : "", 'sample_members' : nsample , 'sample_proportion': psample}, ignore_index=True)
  elif (unknown == False and lpool > 0): 
    prop_of_mayority(pool)         
    if major_perc >= 75: #  specify here majority threshold 
      get_lineage( (major_id), req_rank) # get tax name of major hit
      df = df.append({'cluster_id': cluster_counter , 'cluster_size': len(values), 'cluster_tax': list(des_lin.values())[0].encode("ascii"), 'majority_%' : major_perc, 'sample_members' : nsample , 'sample_proportion': psample}, ignore_index=True)
      yyy = list(des_lin.values())[0].encode("ascii")
    else:
      df = df.append({'cluster_id': cluster_counter , 'cluster_size': len(values), 'cluster_tax': "no_refseq_with_majority", 'majority_%' : major_perc, 'sample_members' : nsample , 'sample_proportion': psample}, ignore_index=True)
  else:
    break # this should not happen
  print("cluster %s: members = %s  sample number = %s  sample proportion = %s tax_assignment %s "%(cluster_counter , len(values) , nsample , psample, yyy))
  yyy = ''
  lpool=''
  print(pool)
  pool = []
  des_lin = ''
  major_perc = ''
  member_counter = 0
  unknown=''
  
fii.close()

# Output files to specified path
df.to_csv('LM_final_cluster_summary.tsv', index=False, sep="\t")
df_memb.to_csv('LM_final_cluster_to_memb.map', index=False, sep="\t")

