print("Importing libraries...", flush=True)
import numpy as np
import pandas as pd
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
from textwrap import wrap
import statsmodels.api as sm
import statsmodels.formula.api as smf
import subprocess
import matplotlib.gridspec as gridspec
from dateutil.relativedelta import relativedelta
import matplotlib
import icd10
import glob
import os
import networkx as nx
from sklearn.cross_decomposition import CCA
import sys
import paperfigures
import itertools
import pathways
import importlib

start_time = time.time()
def get_program_time():
  return (time.time()-start_time)

filtered_taxa = pd.read_parquet("data/filtered_taxa.parquet")["name"]
metab_to_panel_quant = pd.read_parquet("data/metab_to_panel_quant.parquet")

class Results:
  def load_simple_associations(self, folder_name="data/all_associations", filter_taxa=True):
    # global all_results_startmed, lm_convergence, Xs, Ys, Ys_before, Ys_after, collection_lengths, microbe_clusters, correlation_threshold, tables_medication_started_decor, all_results_medications, table_pvalues, table_significant, table_coefs, table_results, table_results_synonyms
    print("load associations from %s..."%folder_name)
    all_results_startmed = pd.read_parquet("%s/all_results_startmed.parquet"%(folder_name))
    lm_convergence = pd.read_parquet("%s/lm_convergence.parquet"%(folder_name)).values.reshape(-1)
    table_stderror = pd.read_parquet("%s/stderror.parquet"%(folder_name))
    Xs = pd.read_parquet("%s/Xs.parquet"%(folder_name))
    Ys = pd.read_parquet("%s/Ys.parquet"%(folder_name))
    keys = Xs.index.to_frame()[["duration","sampletype","medicationtype"]].drop_duplicates().values
    Ys_before = {}
    Ys_after = {}
    tables_medication_started_decor={}
    all_results_startmed.columns = all_results_startmed.columns.str.replace("_indole_indole_","_@TEMP1_").str.replace("_indole_","_indole#").str.replace("_bile_","_bile#").str.replace("_scfa_","_scfa#").str.replace("@TEMP1","indole#indole")
    table_stderror.columns = table_stderror.columns.str.replace(r"^indole_indole_","@TEMP1_",regex=True).str.replace(r"^indole_","indole#",regex=True).str.replace(r"^bile_","bile#",regex=True).str.replace(r"^scfa_","scfa#",regex=True).str.replace("@TEMP1","indole#indole")
    for key in keys:
      Ys_before[tuple(key)] = pd.read_pickle("%s/Ys_before_%d_%s.pickle"%(folder_name,key[0],key[1],))
    for key in keys:
      Ys_after[tuple(key)] = pd.read_pickle("%s/Ys_after_%d_%s.pickle"%(folder_name,key[0],key[1],))
    with open("%s/collection_lengths.pickle"%(folder_name),"rb") as handle:
      collection_lengths = pickle.load(handle)
    microbe_clusters = None
    correlation_threshold = None
    has_microbe_clusters = os.path.exists("%s/microbe_clusters.parquet"%(folder_name))
    has_correlation_threshold = os.path.exists("%s/correlation_threshold.parquet"%(folder_name))
    if has_microbe_clusters:
      microbe_clusters = pd.read_parquet("%s/microbe_clusters.parquet"%(folder_name))
    if has_correlation_threshold:
      correlation_threshold = pd.read_parquet("%s/correlation_threshold.parquet"%(folder_name)).values.item()
    for key in range(len(np.unique(keys[:,0]))):
      print("loading tables_medication_started_decor: ", "%s/tables_medication_started_decor_%d.parquet"%(folder_name,key))
      tables_medication_started_decor[key] = pd.read_parquet("%s/tables_medication_started_decor_%d.parquet"%(folder_name,key))
    Xs = pd.read_parquet("%s/Xs.parquet"%(folder_name))
    keys = Xs.index.to_frame()[["duration","sampletype","medicationtype"]].drop_duplicates().values
    keys = [tuple(x) for x in keys]
    self.Xs = {key:Xs.reset_index().set_index(["duration","sampletype","medicationtype"]).loc[(key)] for key in keys}
    Ys = pd.read_parquet("%s/Ys.parquet"%(folder_name)).set_index(["duration","sampletype","medicationtype"])
    keys = Ys.index.to_frame()[["duration","sampletype","medicationtype"]].drop_duplicates().values
    keys = [tuple(x) for x in keys]
    self.Ys = {key:Ys.reset_index().set_index(["duration","sampletype","medicationtype"]).loc[(key)] for key in keys}
    # Xs = {tuple(key):Xs[(Xs.index.to_frame()[["duration","sampletype","medicationtype"]]==key).all(axis=1)].dropna(axis=1,how="all") for key in keys}
    # Ys = {tuple(key):Ys[(Ys.index.to_frame()[["duration","sampletype","medicationtype"]]==key).all(axis=1)].dropna(axis=1,how="all") for key in keys}
    if has_microbe_clusters:
      microbe_clusters = {tuple(key):microbe_clusters[(microbe_clusters.index.to_frame()[["duration","sampletype"]]==key).all(axis=1)].dropna(axis=1,how="all") for key in keys}
    print("process associations into results table...")
    all_results_medications = all_results_startmed[all_results_startmed.index.isin(np.concatenate([tables_medication_started_decor[0].columns, tables_medication_started_decor[0].columns+"/Dose", tables_medication_started_decor[0].columns+"/Stepwise1", tables_medication_started_decor[0].columns+"/Stepwise2", tables_medication_started_decor[0].columns+"/Stepwise3", tables_medication_started_decor[0].columns+"/Stepwise4"]))]
    table_pvalues = all_results_medications.loc[:,all_results_medications.columns.str.startswith("pvalue")]
    # print(table_pvalues)
    # print(lm_convergence)
    if lm_convergence.shape[0] != table_pvalues.shape[1]:
      print("Convergence array is of incorrect size: ",lm_convergence.shape[0], table_pvalues.shape[1])
      lm_convergence = np.ones(table_pvalues.shape[1]).astype(bool)
    table_pvalues = table_pvalues.loc[:, lm_convergence]
    table_pvalues.columns = table_pvalues.columns.str.replace("alpha_diversity","alphadiversity").str.replace("microbe_genus","genus").str.replace("microbe_species","species")
    table_pvalues.columns = table_pvalues.columns.str.split("_").str[1:].str.join("_")
    pvalues_flat = table_pvalues.values.reshape(-1)
    significant_notnull = statsmodels.stats.multitest.multipletests(pvalues_flat[~np.isnan(pvalues_flat)], method='fdr_bh', alpha=0.1)[0]
    significant = np.zeros(len(pvalues_flat)).astype("bool")
    significant[~np.isnan(pvalues_flat)] = significant_notnull
    print("N significant:",significant.sum())
    table_significant = pd.DataFrame(significant.reshape(table_pvalues.shape), index=table_pvalues.index, columns=table_pvalues.columns)
    table_coefs = all_results_medications.loc[:,all_results_medications.columns.str.startswith("coef")].fillna(0)
    table_coefs = table_coefs.loc[:, lm_convergence]
    table_coefs.columns = table_coefs.columns.str.split("_").str[1:].str.join("_")
    table_stderror = table_stderror.loc[:,lm_convergence]
    table_results = pd.DataFrame({"pvalue":table_pvalues.stack(),"coef":table_coefs.stack(),"significant":table_significant.stack(),"stderror":table_stderror.stack()})
    table_results = table_results.sort_values("pvalue").reset_index().rename({"level_0":"medication","level_1":"microbe"},axis=1)
    print("Determining medication classes...",end="")
    medications = table_results.medication.str.split("|").str[0].str.replace("/Dose","").unique()
    meds_ref = pd.read_parquet("data/medication_classes.parquet")
    medications = meds_ref.index.intersection(medications)
    medication_classes = meds_ref.loc[medications]
    print("Done")
    table_results.loc[table_results.medication.str.split("\\|").str[0].str.replace("(/Dose|/Stepwise1|/Stepwise2|/Stepwise3|/Stepwise4)","",regex=True).isin(medications),"medication_pharm_class"] = medication_classes.loc[table_results[table_results.medication.str.split("\\|").str[0].str.replace("(/Dose|/Stepwise1|/Stepwise2|/Stepwise3|/Stepwise4)","",regex=True).isin(medications)].medication.str.split("|").str[0].str.replace("(/Dose|/Stepwise1|/Stepwise2|/Stepwise3|/Stepwise4)","",regex=True)]["med_pharm_class"].reset_index(drop=True)
    table_results["microbe"] = table_results["microbe"].str.replace("alpha_diversity","alphadiversity").str.replace("microbe_genus","genus").str.replace("microbe_species","species")
    table_results["duration"] = table_results["microbe"].str.split("_").str[-3].astype("int")
    table_results["sampletype"] = table_results["microbe"].str.split("_").str[-2]
    table_results["modeltype"] = table_results["microbe"].str.split("_").str[-1]
    table_results["microbe"] = table_results["microbe"].str.split("_").str[0]
    table_results["medication_dose"] = table_results["medication"].str.endswith("/Dose")
    table_results["medication_stepwise"] = 0
    table_results.loc[table_results["medication"].str.endswith("/Stepwise1"),"medication_stepwise"] = 1
    table_results.loc[table_results["medication"].str.endswith("/Stepwise2"),"medication_stepwise"] = 2
    table_results.loc[table_results["medication"].str.endswith("/Stepwise3"),"medication_stepwise"] = 3
    table_results.loc[table_results["medication"].str.endswith("/Stepwise4"),"medication_stepwise"] = 4
    table_results["medication"] = table_results["medication"].str.split("/").str[0]
    if filter_taxa:
      table_results_unfiltered = table_results
      table_results = table_results[(table_results.sampletype.isin(["metab","pathway"])) | (table_results.microbe.isin(list(filtered_taxa)+["alphadiversity"]))].sort_values("pvalue")
      table_results = table_results[(~(table_results.sampletype=="metab")) | ((table_results.sampletype=="metab") & (table_results.microbe.str.startswith("quant#")))]
      table_results["microbe"].replace(metab_to_panel_quant[["compound_quant","compound_pfx"]].set_index("compound_quant",drop=True).to_dict()["compound_pfx"], inplace=True)
    for sampletype in table_results["sampletype"].unique():
      for duration in table_results["duration"].unique():
        for modeltype in table_results["modeltype"].unique():
          if len(table_results.loc[(~table_results["pvalue"].isna()) & (table_results["modeltype"]==modeltype) & (table_results["sampletype"]==sampletype) & (table_results["duration"]==duration),"significant"]) > 0:
            table_results.loc[(~table_results["pvalue"].isna()) & (table_results["modeltype"]==modeltype) & (table_results["sampletype"]==sampletype) & (table_results["duration"]==duration),"significant"] = statsmodels.stats.multitest.multipletests(table_results[(~table_results["pvalue"].isna()) & (table_results["modeltype"]==modeltype) & (table_results["sampletype"]==sampletype) & (table_results["duration"]==duration)]["pvalue"], method="fdr_bh", alpha=0.1)[0]
          if len(table_results_unfiltered.loc[(~table_results_unfiltered["pvalue"].isna()) & (table_results_unfiltered["modeltype"]==modeltype) & (table_results_unfiltered["sampletype"]==sampletype) & (table_results_unfiltered["duration"]==duration),"significant"]) > 0:
            table_results_unfiltered.loc[(~table_results_unfiltered["pvalue"].isna()) & (table_results_unfiltered["modeltype"]==modeltype) & (table_results_unfiltered["sampletype"]==sampletype) & (table_results_unfiltered["duration"]==duration),"significant"] = statsmodels.stats.multitest.multipletests(table_results_unfiltered[(~table_results_unfiltered["pvalue"].isna()) & (table_results_unfiltered["modeltype"]==modeltype) & (table_results_unfiltered["sampletype"]==sampletype) & (table_results_unfiltered["duration"]==duration)]["pvalue"], method="fdr_bh", alpha=0.1)[0]
    # table_results["coef"] /= np.log(2.0) # no need to do this anymore, since I've made the model in terms of log2. 
    table_results = table_results[~table_results["pvalue"].isna()]
    table_results_unfiltered = table_results_unfiltered[~table_results_unfiltered["pvalue"].isna()]
    self.all_results_startmed = all_results_startmed
    self.lm_convergence = lm_convergence
    # self.Xs = Xs
    # self.Ys = Ys
    self.Ys_before = Ys_before
    self.Ys_after = Ys_after
    self.collection_lengths = collection_lengths
    self.microbe_clusters = microbe_clusters
    self.correlation_threshold = correlation_threshold
    self.tables_medication_started_decor = tables_medication_started_decor
    self.all_results_medications = all_results_medications
    self.table_pvalues = table_pvalues
    self.table_significant = table_significant
    self.table_coefs = table_coefs
    self.table = table_results
    self.table_unfiltered = table_results_unfiltered

# load results table
print("(%.2f) Loading results..."%(get_program_time()), flush=True)
results = Results()
results.load_simple_associations("data/metab_associations_10_20_30_quant")

# load intermediate data files
if True:
  paperfigures.table_microbes = pd.read_parquet("data/table_microbes.parquet")
  paperfigures.table_microbes_clr = pd.read_parquet("data/table_microbes_clr.parquet")
  paperfigures.samples = pd.read_parquet("data/samples.parquet")
  paperfigures.table_demographics = pd.read_parquet("data/table_demographics.parquet")
  paperfigures.hospital_visits = pd.read_parquet("data/hospital_visits.parquet")
  paperfigures.study_intervals = pd.read_parquet("data/study_intervals.parquet")
  paperfigures.min_date = datetime.datetime.strptime(pd.read_parquet("data/min_date.parquet").values.item(), "%Y-%m-%d")
  paperfigures.medications = pd.read_parquet("data/medications.parquet")["medication"].values
  paperfigures.durations = pd.read_parquet("data/durations.parquet")
  paperfigures.medications_started_dose = pd.read_parquet("data/medications_started_dose.parquet")
  paperfigures.medications_started_starttime = pd.read_parquet("data/medications_started_starttime.parquet")
  paperfigures.table_microbes_relabund = pd.read_parquet("data/table_microbes_relabund.parquet")
  paperfigures.table_pathways = pd.read_parquet("data/table_pathways.parquet")
  table_pathways = paperfigures.table_pathways
  paperfigures.table_genus_metab_correl = pd.read_parquet("data/table_genus_metab_correl.parquet")
  paperfigures.matrix_genus = pd.read_parquet("data/matrix_genus.parquet")
  paperfigures.matrix_metab = pd.read_parquet("data/matrix_metab.parquet")
  paperfigures.table_invitro_reference = pd.read_parquet("data/table_invitro_reference.parquet")
  paperfigures.correspondences = pd.read_parquet("data/correspondences.parquet")
  paperfigures.table_pathways_general = pd.read_parquet("data/table_pathways_general.parquet")
  paperfigures.medication_classes = pd.read_parquet("data/medication_classes.parquet")
  paperfigures.medication_counts = pd.read_parquet("data/medication_counts.parquet")
  paperfigures.medication_class_counts = pd.read_parquet("data/medication_class_counts.parquet")
  paperfigures.results = results

pdb.set_trace()

# load pathways
if True:
  table_metab = pd.read_parquet("data/table_metab.parquet")
  importlib.reload(pathways)
  table_metab_reindexed = table_metab.copy()
  table_metab_reindexed.index = table_metab_reindexed.index.str.split("#").str[1]
  pathways.load_saved()
  pathways.pathway_enrichment(table_pathways, results.table, table_metab_reindexed)
  paperfigures.pathways = pathways

# build figures
def update_figures():
  importlib.reload(paperfigures)
  paperfigures.make_figures()

if __name__ == "__main__":
  update_figures()
