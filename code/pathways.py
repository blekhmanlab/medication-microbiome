print("pathways.py: Importing libraries...", flush=True)
import json
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import pdb as pdb
import datetime
import matplotlib.pyplot as plt
import re
from sklearn.decomposition import PCA
import time
import pickle
import warnings
import skbio.stats.composition
import seaborn as sns
import scipy.stats
import statsmodels.stats.multitest
from sklearn.decomposition import SparsePCA
import statsmodels.api as sm
import statsmodels.formula.api as smf
from dateutil.relativedelta import relativedelta
import matplotlib
import icd10
import glob
import os
import networkx as nx
from sklearn.cross_decomposition import CCA
import sys
from sklearn.preprocessing import MultiLabelBinarizer


# normalize metacyc pathways hierarchy json file.
def to_lists(data, prefix=[]):
  if len(data) == 0:
    return "|".join(prefix)
  if type(data) == list:
    return [to_lists(child, prefix) for child in data]
  if type(data) == dict:
    return to_lists(data["children"], prefix + [data["id"]+"%"+re.sub("\([0-9]+\)","",re.sub("<.?a>","",data["text"])).strip()])

def flatten(items, seqtypes=(list, tuple, np.ndarray)):
  for i, x in enumerate(items):
    while i < len(items) and isinstance(items[i], seqtypes):
      items[i:i+1] = items[i]
  return items

class ProgressBar:
  def __init__(self, length):
    self.text =''
    self.length = length
    self.start_time = time.time()
  def update(self, progress, caption=''):
    self.progress = progress
    print('\b'*len(self.text), end='', flush=True)
    elapsed = time.time() - self.start_time
    estimated = elapsed * (self.length / max(1, self.progress))
    self.text = '%.2fs/%.2fs: (%.1f/%.1f) %s'%(elapsed, estimated, self.progress, self.length, caption)
    print(self.text, end='', flush=True)

# Metacyc Pathway Hierarchy
def load_pathway_hierarchy(table_pathways, table_results):
  # global pathways_flat, pathways_hierarchy, pathways_index, pathway_to_hierarchy, pathway_map, pathway_to_pathway_general, metabolite_processes, metabolites_lu, metabolite_to_pathway, table_metab_pathways, table_metab_pathways, metabs_positive, metabs_negative, include_pathways, pathways_positive, pathways_negative, results_enrichment_medication_metabolites
  global pathway_to_pathway_general, metabolite_processes, metabolites_lu, metabolite_to_pathway, table_metab_pathways, table_metab_pathways, metabs_positive, metabs_negative, include_pathways, pathways_positive, pathways_negative, results_enrichment_medication_metabolites
  try:
    pathway_to_pathway_general = pd.read_parquet("data/pathway_to_pathway_general.parquet")
    return
  except:
    pass
  print("load pathway hierarchy")
  with open('data/metacyc_pathway_hierarchy.json') as f:
    data = json.load(f)
  lists = to_lists(data)
  print("mapping pathways")
  pathways_flat = pd.DataFrame(flatten(lists))
  pathways_hierarchy = pd.DataFrame([x[0].split("|") for x in pathways_flat.values]).ffill(axis=1)
  pathways_index = pd.DataFrame({"pathway_name":table_pathways.index,"database_name":table_pathways.index})
  pathways_index["database_name"] = pathways_index["database_name"].str.replace("&|;","",regex=True)
  pathway_to_hierarchy = ({name:(pathways_flat[0][(pathways_flat[0].str.lower().str.contains((pathways_index.set_index("pathway_name").loc[name].item().lower()+":").lower().split(":")[0],regex=False)) & (pathways_flat[0].str.lower().str.contains((re.sub("&|;","",pathways_index.set_index("pathway_name").loc[name].item().lower())+":").lower().split(":")[1].strip(),regex=False))].values) for name in pathways_index["pathway_name"]})
  unmapped_pathways = pd.DataFrame([(x,len(pathway_to_hierarchy[x])) for x in pathway_to_hierarchy])
  print("mapping pathways 2")
  unmapped_pathways = unmapped_pathways[unmapped_pathways[1]==0]
  pathway_map = {source:source.split(":")[0]+": "+(list(pathways_flat[pathways_flat[0].str.contains(source.split(":")[0])].values.reshape(-1)) + [""])[0].split("%")[-1] for source in unmapped_pathways[0]}
  pathway_map["DTDPRHAMSYN-PWY: dTDP-beta-L-rhamnose biosynthesis"] = "DTDPRHAMSYN-PWY: dTDP-L-rhamnose biosynthesis"
  pathway_map["PWY-5431: aromatic compounds degradation via beta-ketoadipate"] = "PWY-5431: aromatic compounds degradation via 3-oxoadipate"
  pathways_index["database_name"] = pathways_index["database_name"].replace(pathway_map)
  pathway_to_hierarchy = ({name:(pathways_flat[0][(pathways_flat[0].str.lower().str.contains((pathways_index.set_index("pathway_name").loc[name].item().lower()+":").lower().split(":")[0],regex=False)) & (pathways_flat[0].str.lower().str.contains((re.sub("&|;","",pathways_index.set_index("pathway_name").loc[name].item().lower())+":").lower().split(":")[1].strip(),regex=False))].values) for name in pathways_index["pathway_name"]})
  unmapped_pathways = pd.DataFrame([(x,len(pathway_to_hierarchy[x])) for x in pathway_to_hierarchy])
  unmapped_pathways = unmapped_pathways[unmapped_pathways[1]==0]
  print("pathways to pathway general")
  pathway_to_pathway_general = pd.DataFrame([(x, pathway_to_hierarchy[x]) for x in pathway_to_hierarchy])
  pathway_to_pathway_general[1] = pathway_to_pathway_general[1].apply(lambda hierarchies: np.concatenate([hierarchies,["empty"]])).apply(lambda hierarchies: hierarchies[np.argmax([len(h) for h in hierarchies])])
  pathway_to_pathway_general = pathway_to_pathway_general.set_index(0)
  pathway_to_pathway_general["all"]= pathway_to_pathway_general[1]
  for level in range(5):
    pathway_to_pathway_general[level] = pathway_to_pathway_general["all"].str.split("|").apply(lambda items: items[-(level+1)] if len(items)>=(level+1) else items[0]).apply(lambda name: name.split("%")[-1])
  del pathways.pathways_flat, pathways.pathways_hierarchy, pathways.pathways_index, pathways.pathway_to_hierarchy, pathways.pathway_map
  pathway_to_pathway_general.columns = [str(x) for x in pathway_to_pathway_general.columns]
  pathway_to_pathway_general.index.name="pathway"
  pathway_to_pathway_general.to_parquet("data/pathway_to_pathway_general.parquet")
  #
  # Metabolite Pathway Enrichment Analysis:
  #

# pathways.pathway_to_pathway_general
# pathways.metabolite_to_pathway
# pathways.table_metab_pathways
# pathways.metabs_positive
# pathways.metabs_negative
# pathways.include_pathways
# pathways.pathways_positive
# pathways.pathways_negative

# pathways.metabolite_processes
# pathways.metabolites_lu

# pathways.table_metab_pathways


# df = pathways.metabolite_to_pathway.pathways.map(lambda xs: "".join(xs))
# cs = np.concatenate(pathways.metabolite_to_pathway.pathways.map(lambda xs: "".join(xs)).values).unique()

# pathways.pathway_to_pathway_general.to_parquet("code/data/pathways/pathway_to_pathway_general.parquet")
# df = pathways.metabolite_to_pathway.copy()
# df.pathways = df.pathways.map(lambda xs: ";".join(xs))
# df.to_parquet("code/data/pathways/metabolite_to_pathway.parquet")
# pd.DataFrame(pathways.table_metab_pathways.map(lambda xs: ";".join(xs))).to_parquet("code/data/pathways/table_metab_pathways.parquet")
# pd.DataFrame(pathways.metabs_positive,columns=["metabs"]).to_parquet("code/data/pathways/metabs_positive.parquet")
# pd.DataFrame(pathways.metabs_negative,columns=["metabs"]).to_parquet("code/data/pathways/metabs_negative.parquet")
# pd.DataFrame(pathways.include_pathways,columns=["pathway"]).to_parquet("code/data/pathways/include_pathways.parquet")
# pd.DataFrame(pathways.pathways_positive,columns=["pathways"]).to_parquet("code/data/pathways/pathways_positive.parquet")
# pd.DataFrame(pathways.pathways_negative,columns=["pathways"]).to_parquet("code/data/pathways/pathways_negative.parquet")
# pathways.results_enrichment_medication_metabolites.to_parquet("code/data/pathways/results_enrichment_medication_metabolites.parquet")

def load_metabolites(table_pathways, table_results, table_metab):
  print("loading metabolites tree")
  tree = ET.parse('data/hmdb_metabolites.xml')
  # Pathways
  print("Metabolite Pathways")
  metabolites = tree.findall(r"{http://www.hmdb.ca}metabolite")
  del tree
  progressbar = ProgressBar(len(metabolites))
  metabolites_lu = {}
  for mi,metabolite in enumerate(metabolites):
    progressbar.update(mi)
    synonyms = [x.text for x in metabolite.findall(r"{http://www.hmdb.ca}synonyms")[0].iter() if x.tag=='{http://www.hmdb.ca}synonym']
    pathways = np.array([[[x2.find(r"{http://www.hmdb.ca}name").text for x2 in x1.findall(r"{http://www.hmdb.ca}pathway")] for x1 in x.findall(r"{http://www.hmdb.ca}pathways")] for x in metabolite.findall(r"{http://www.hmdb.ca}biological_properties")]).reshape(-1)
    name =  metabolite.find(r"{http://www.hmdb.ca}name").text
    metabolites_lu[name] = {"synonyms":synonyms, "pathways":pathways}
  #
  # Ontology/Processes
  print("Ontology/Processes")
  metabolite_processes = {}
  progressbar = ProgressBar(len(metabolites))
  for mi,metabolite in enumerate(metabolites):
    progressbar.update(mi)
    name = metabolite.find(r"{http://www.hmdb.ca}name").text
    metabolite_processes[name] = []
    # print(name+":")
    for elem in metabolite.iter():
      if elem.text=="Naturally occurring process":
        # print(elem.text)
        for child in prev.iter():
          # print(child.tag, child.text)
          if child.tag==r"{http://www.hmdb.ca}term":
            if child.text not in ["Naturally occurring process","Biological process","Biochemical pathway","Biochemical process"]:
              metabolite_processes[name] += [child.text]
              # print(child.text)
      prev = elem
  del metabolites
  metabolite_processes = pd.DataFrame({key:[key,metabolite_processes[key]] for key in metabolite_processes.keys()}).T
  metabolite_processes = metabolite_processes.drop(0,axis=1)
  #
  metabolites_lu = pd.DataFrame(metabolites_lu).T
  metabolites_lu["processes"] = metabolite_processes[1]
  if "processes" in metabolites_lu.columns:
    metabolites_lu["all_pathways"] = metabolites_lu.apply(lambda row: np.unique(np.concatenate([row["pathways"],row["processes"]])), axis=1)
  else:
    metabolites_lu["all_pathways"] = metabolites_lu.apply(lambda row: np.unique(np.concatenate([row["pathways"]])), axis=1)
  print("Metabolite to Pathway")
  del metabolite_processes
  metabolite_to_pathway = metabolites_lu[["synonyms","pathways"]].copy()
  del metabolites_lu
  metabolite_to_pathway["synonyms"] = metabolite_to_pathway.apply(lambda row: [row.name] + row["synonyms"],axis=1)
  metabolite_to_pathway = metabolite_to_pathway.explode("synonyms").set_index("synonyms",drop=True)
  metabolite_to_pathway.index = metabolite_to_pathway.index.str.replace("[^a-zA-Z0-9]","",regex=True).str.lower()
  metabolite_to_pathway.columns = ["pathways"]
  synonyms = {
   '12dioxolithocholicacid':"",
   '3aminoisobutyrate2':"",
   '3oxodeoxycholicacidor3oxochenodeoxycholicacid':"",
   '3oxodeoxycholicor3oxochenodeoxycholicacid':"",
   '3oxolithocholicacid':"",
   '3oxoor3oxochenodeoxycholicacid':"",
   '5ohtryptophan':"5hydroxyltryptophan",
   '712dioxolithocholicacid':"",
   '7oxoor6oxolithocholicacid':"",
   'alloisolithocholicacid':"",
   'glycodehydrocholicacid':"",
   'omegamuricholicor3epicholicacid':"",
   'panthothenic':"panthothenic acid",
   'tauroalphaortaurobetamuricholicacid':"",
   'taurohyodeoxycholicacid':"",
   'transindole3acrylate':""
  }
  table_metab_pathways = pd.Series(table_metab.index).replace(synonyms).apply(lambda row: np.concatenate([np.array([],dtype=str), metabolite_to_pathway["pathways"].loc[row]]).tolist() if row in metabolite_to_pathway.index else "").apply(flatten)
  table_metab_pathways.index = table_metab.index
  table_metab_pathways = table_metab_pathways.loc[~table_metab_pathways.index.duplicated(keep="first")]
  pathway_to_pathway_general.to_parquet("code/data/pathways/pathway_to_pathway_general.parquet")
  df = metabolite_to_pathway.copy()
  df.pathways = df.pathways.map(lambda xs: ";".join(xs))
  df.to_parquet("code/data/pathways/metabolite_to_pathway.parquet")
  pd.DataFrame(table_metab_pathways.map(lambda xs: ";".join(xs))).to_parquet("code/data/pathways/table_metab_pathways.parquet")
  pd.DataFrame(metabs_positive,columns=["metabs"]).to_parquet("code/data/pathways/metabs_positive.parquet")
  pd.DataFrame(metabs_negative,columns=["metabs"]).to_parquet("code/data/pathways/metabs_negative.parquet")
  pd.DataFrame(include_pathways,columns=["pathway"]).to_parquet("code/data/pathways/include_pathways.parquet")
  pd.DataFrame(pathways_positive,columns=["pathways"]).to_parquet("code/data/pathways/pathways_positive.parquet")
  pd.DataFrame(pathways_negative,columns=["pathways"]).to_parquet("code/data/pathways/pathways_negative.parquet")
  results_enrichment_medication_metabolites.to_parquet("code/data/pathways/results_enrichment_medication_metabolites.parquet")

def load_saved():
  global pathway_to_pathway_general, metabolite_processes, metabolites_lu, metabolite_to_pathway, table_metab_pathways, table_metab_pathways, metabs_positive, metabs_negative, include_pathways, pathways_positive, pathways_negative, results_enrichment_medication_metabolites
  pathway_to_pathway_general = pd.read_parquet("data/pathways/pathway_to_pathway_general.parquet")
  metabolite_to_pathway = pd.read_parquet("data/pathways/metabolite_to_pathway.parquet")
  metabolite_to_pathway.pathways = metabolite_to_pathway.apply(lambda xs: xs.str.split(";"))
  table_metab_pathways = pd.read_parquet("data/pathways/table_metab_pathways.parquet")["compound"].str.split(";")
  metabs_positive = pd.read_parquet("data/pathways/metabs_positive.parquet")["metabs"]
  metabs_negative = pd.read_parquet("data/pathways/metabs_negative.parquet")["metabs"]
  include_pathways = pd.read_parquet("data/pathways/include_pathways.parquet")["pathway"]
  pathways_positive = pd.read_parquet("data/pathways/pathways_positive.parquet")["pathways"]
  pathways_negative = pd.read_parquet("data/pathways/pathways_negative.parquet")["pathways"]
  results_enrichment_medication_metabolites = pd.read_parquet("data/pathways/results_enrichment_medication_metabolites.parquet")


def pathway_enrichment(table_pathways, table_results, table_metab):
  global pathway_to_pathway_general, metabolite_processes, metabolites_lu, metabolite_to_pathway, table_metab_pathways, table_metab_pathways, metabs_positive, metabs_negative, include_pathways, pathways_positive, pathways_negative, results_enrichment_medication_metabolites
  table_metab_pathways = table_metab_pathways.loc[~table_metab_pathways.index.duplicated(keep="first")]
  print("Pathway Enrichment Analysis")
  binarizer = MultiLabelBinarizer(sparse_output=True)
  table_metab_pathways_binary = pd.DataFrame.sparse.from_spmatrix(binarizer.fit_transform(table_metab_pathways), index=table_metab_pathways.index, columns=binarizer.classes_)
  table_metab_pathways_binary = table_metab_pathways_binary.sparse.to_dense()
  metabs_positive = table_results[(table_results.sampletype=="metab") & (table_results.significant) & (table_results.coef>0)].drop_duplicates(["medication","microbe"]).groupby("medication").apply(lambda df: df["microbe"].values)
  metabs_negative = table_results[(table_results.sampletype=="metab") & (table_results.significant) & (table_results.coef<0)].drop_duplicates(["medication","microbe"]).groupby("medication").apply(lambda df: df["microbe"].values)
  metabs_agnostic = table_results[(table_results.sampletype=="metab") & (table_results.significant)].drop_duplicates(["medication","microbe"]).groupby("medication").apply(lambda df: df["microbe"].values)
  metabs_positive = metabs_positive.apply(lambda xs: [x.split("#")[-1] for x in xs])
  metabs_negative = metabs_negative.apply(lambda xs: [x.split("#")[-1] for x in xs])
  metabs_agnostic = metabs_agnostic.apply(lambda xs: [x.split("#")[-1] for x in xs])
  include_pathways = table_metab_pathways_binary.sum(axis=0).sort_values(ascending=False).head(100).index
  pathways_positive = metabs_positive.apply(lambda metabs: list(set(flatten([table_metab_pathways.loc[m] for m in metabs])).intersection(include_pathways)))
  pathways_negative = metabs_negative.apply(lambda metabs: list(set(flatten([table_metab_pathways.loc[m] for m in metabs])).intersection(include_pathways)))
  pathways_agnostic = metabs_agnostic.apply(lambda metabs: list(set(flatten([table_metab_pathways.loc[m] for m in metabs])).intersection(include_pathways)))
  results_enrichment_medication_metabolites = pd.DataFrame(np.nan, index=np.arange(len(include_pathways)*(len(metabs_positive)+len(metabs_negative)+len(metabs_agnostic))), columns=["medication","pathway","direction","pvalue"])
  ri = 0
  for direction in [1,0,-1]:
    df = {1:metabs_positive, -1:metabs_negative, 0:metabs_agnostic}[direction]
    for pathway in include_pathways:
      for med,metabs in df.items():
        total_metabolites_background = len(table_metab)   # total metabolites in background distribution.
        annotated_metabolites_to_pathway = table_metab_pathways_binary.loc[:,pathway].sum()
        number_metabolites_in_set = len(metabs)
        overlap_metabs_with_pathway = table_metab_pathways_binary.loc[metabs,pathway].sum()
        results_enrichment_medication_metabolites.iloc[ri,0] = med
        results_enrichment_medication_metabolites.iloc[ri,1] = pathway
        results_enrichment_medication_metabolites.iloc[ri,2] = direction
        results_enrichment_medication_metabolites.iloc[ri,3] = 1.0 - scipy.stats.hypergeom(total_metabolites_background, annotated_metabolites_to_pathway, number_metabolites_in_set).cdf(overlap_metabs_with_pathway)
        ri += 1
  results_enrichment_medication_metabolites["significant"] = statsmodels.stats.multitest.multipletests(results_enrichment_medication_metabolites["pvalue"], method='fdr_bh', alpha=0.1)[0]
  results_enrichment_medication_metabolites.loc[results_enrichment_medication_metabolites.direction==0,"significant"] = statsmodels.stats.multitest.multipletests(results_enrichment_medication_metabolites[results_enrichment_medication_metabolites.direction==0]["pvalue"], method='fdr_bh', alpha=0.1)[0]
  results_enrichment_medication_metabolites = results_enrichment_medication_metabolites.sort_values("pvalue")
  results_enrichment_medication_metabolites.head(50)
