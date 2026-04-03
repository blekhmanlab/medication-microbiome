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

print("Importing rpy2...", flush=True)
import rpy2
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import pandas2ri
from sklearn.preprocessing import MultiLabelBinarizer
from rpy2.robjects.conversion import localconverter

start_time = time.time()
def get_program_time():
  return (time.time()-start_time)

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

def abbr(text, n=25):
  if len(text)>n-2:
    return text[:n-2]+'...'
  else:
    return text

def harmonize_medication(name):
  name = re.sub(r'\(IRB [0-9]+\)','', name)                   # remove IRB label.
  name = re.sub(r' ER ', ' ', name)                            # remove extended-release ER designation.
  name = re.sub(r'EXTENDED RELEASE', ' ', name)                            # remove extended-release ER designation.
  name = re.sub(r', HUMAN', ' ', name)                            # remove extended-release ER designation.
  name = re.sub(r'[><] ?[0-9]+ [A-Z][A-Z]', ' ', name)         # remove specifier (eg. > 50KG)
  name = re.sub(r' HFA ', ' ', name)                          # remove HFA designation (indicates propellant type in sprays)
  name = re.sub(r' FOR [A-Z0-9-]+( [A-Z0-9-]+)?', '', name)     # brand name
  name = re.sub(r' (IVPB|IV|IVBP)( |$)', ' ', name)     # IVPB -> IV
  name = re.sub(r' IN (NS|SWFI|SWFI/NS)?', ' ', name)     # IVPB -> IV
  name = re.sub(r' HCL ', ' ', name)     # IVPB -> IV
  name = re.sub(r' XL ', ' ', name)     # IVPB -> IV
  name = re.sub(r' XR ', ' ', name)     # IVPB -> IV
  name = re.sub(r' SOD SUCC ', ' ', name)     # IVPB -> IV
  name = re.sub(r' SODIUM SUCC ', ' ', name)     # IVPB -> IV
  name = re.sub(r'\([A-Za-z0-9-, ]*\)', ' ', name)     # IVPB -> IV
  name = re.sub(r'\(PF\)', ' ', name)     # IVPB -> IV
  name = re.sub(r'\(\s*\)', ' ', name)     # IVPB -> IV
  name = re.sub(r'\s+', ' ', name)                            # remove excess whitespace.
  regex = re.search(r' [0-9]', name)
  if regex is not None:
   name =  name[:regex.span()[0]]
  synonyms = {'AMINO ACID':'AMINO ACID', 'AMPHOTERICIN B LIPOSOME':'AMPHOTERICIN B LIPOSOMAL', 'BEVACIZUMAB-AWWB':'BEVACIZUMAB'}
  if name in synonyms.keys():
    name = synonyms[name]
  return name.strip()

print("Loading data tables...", flush=True)
directory = 'data/Pamer_DM/'
if os.path.isfile(directory+'PAMER_DM_ENC_RX.parquet'):
  meds = pd.read_parquet(directory+'PAMER_DM_ENC_RX.parquet')
else:
  meds = pd.read_csv(directory+'PAMER_DM_ENC_RX.txt', sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3)
  meds.to_parquet(directory+'PAMER_DM_ENC_RX.parquet')

if os.path.isfile(directory+'PAMER_DM_ENC_RX_MAR.parquet'):
  meds2 = pd.read_parquet(directory+'PAMER_DM_ENC_RX_MAR.parquet')
else:
  meds2 = pd.read_csv(directory+'PAMER_DM_ENC_RX_MAR_OLD_PreUpdate20230224.txt', sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3)
  meds2.to_parquet(directory+'PAMER_DM_ENC_RX_MAR.parquet')

meds_lookup = pd.read_csv(directory+"LU_CLARITY_RX.txt", sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3)
meds2 = pd.merge(meds2, meds_lookup, how="left", on="medication_id", suffixes=[None,"_clarity"])
meds2.rename({"ROUTE":"route"}, axis=1, inplace=True)
meds2 = meds2[meds2["mar_action"].isin(["Given","Rate Change","New Bag","Patch Applied","Restarted","Bolus","Full Dose Administered","Partial Dose"])]

print("(%.2f) Loading sequencing data..."%(get_program_time()), flush=True)
def load_sequencing_data(taxonomies=["G"], keep_taxa_list=[]):
  table_microbes_all = None
  table_microbes_relabund_all = None
  table_microbes_clr_all = None
  for taxonomy in taxonomies:
    table_microbes = pd.read_csv('data/bracken.csv')
    table_microbes = table_microbes[table_microbes['taxonomy_lvl']==taxonomy].pivot(columns='filename', values='new_est_reads', index='name').fillna(0)
    table_microbes.index = table_microbes.index.str.replace(r"\[|\]","", regex=True)
    table_microbes_relabund = table_microbes.copy()
    table_microbes_relabund.iloc[:,:] = (table_microbes_relabund / table_microbes_relabund.sum())
    table_microbes_clr = table_microbes.copy()
    table_microbes_clr.iloc[:,:] = skbio.stats.composition.clr(table_microbes_clr.T+1.0).T
    microbes_filter = ((table_microbes!=0).mean(axis=1)>0.3) | (table_microbes.index.isin(keep_taxa_list))
    table_microbes_clr = table_microbes_clr[microbes_filter]
    table_microbes_relabund = table_microbes_relabund[microbes_filter]
    table_microbes.loc["alphadiversity"] = table_microbes.apply(lambda row: scipy.stats.entropy(row), axis=0)
    table_microbes_all = pd.concat([table_microbes_all, table_microbes], axis=1).fillna(0)
    table_microbes_relabund_all = pd.concat([table_microbes_relabund_all, table_microbes_relabund], axis=1).fillna(0)
    table_microbes_clr_all = pd.concat([table_microbes_clr_all, table_microbes_clr], axis=1).fillna(0)
  return (table_microbes_all, table_microbes_relabund_all, table_microbes_clr_all)

table_maier2018_s3a = pd.read_csv("data/maier_2018_S3a.csv")
subset_species_maier2018 = table_maier2018_s3a.columns[4:].str.split("(").str[0].str.strip()
table_forslund_2021_s6c = pd.read_csv("data/forslund_2021_S6c.csv")
table_forslund_2021_s6c = table_forslund_2021_s6c[table_forslund_2021_s6c["Feature space"].str.startswith("mOTU, ")]
table_forslund_2021_s6c["Effector"] = table_forslund_2021_s6c["Effector"].str.split(" ").str[2].str.replace(r"\(|\)","", regex=True)
subset_species_forslund2021  = table_forslund_2021_s6c["Feature display name"].unique()

table_microbes, table_microbes_relabund, table_microbes_clr = load_sequencing_data(["G","S"], keep_taxa_list=list(subset_species_forslund2021) + list(subset_species_maier2018))
filtered_taxa = load_sequencing_data(["G","S"])[1].index

# microbial pathway abundances.
print("(%.2f) Loading pathways..."%(get_program_time()), flush=True)
if os.path.isfile("data/humann_pathway_abundances.parquet"):
  table_pathways = pd.read_parquet("data/humann_pathway_abundances.parquet")
else:
  files_pathways = glob.glob("data/humann_dfi_clin/*pathabundance.tsv")
  humann_abundances = [pd.read_csv(file,sep="\t").set_index("# Pathway") for file in files_pathways]
  columns = [file.columns[0] for file in humann_abundances]
  rows = [set(file.index) for file in humann_abundances]
  rows = list(set.union(*rows))
  table_pathways = pd.DataFrame(0, index=list(rows), columns=columns)
  for file in humann_abundances:
    table_pathways[file.columns[0]]=file
  table_pathways = table_pathways.fillna(0)
  table_pathways.to_parquet("data/humann_pathway_abundances.parquet")


print("(%.2f) Processing pathways..."%(get_program_time()), flush=True)
table_pathways_general = table_pathways.groupby(table_pathways.index.str.split("|").str[0],axis=0).sum()
table_pathways_columns_short = pd.concat([pd.Series(table_pathways.columns,name="pathway_filename"), pd.Series(table_pathways.columns.str.split("_mpa").str[0],name="prefix")],axis=1)
table_pathways_columns_short["microbe_filename"] = [list(table_microbes.columns[table_microbes.columns.str.contains(prefix) & table_microbes.columns.str.contains("genus")]) for prefix in table_pathways_columns_short.prefix]
table_pathways_columns_short["microbe_filename"] = table_pathways_columns_short["microbe_filename"].apply(lambda x: x[0] if len(x)>0 else None)
table_pathways_columns_short = pd.concat([table_pathways_columns_short, pd.DataFrame(data={"pathway_filename":None,"prefix":None,"microbe_filename":list(set(table_microbes.columns).difference(set(table_pathways_columns_short["microbe_filename"])))})])

table_pathways[None] = np.nan
table_pathways_general[None] = np.nan

pdb.set_trace()

# subset to species that were considered by maier et al. 2018.
subset_maier=True
subset_forslund=False
if subset_maier or subset_forslund:
  print("(%.2f) Creating subsets for maier/forslund studies..."%(get_program_time()), flush=True)
  if subset_maier:
    table_maier2018_s3a = pd.read_csv("data/maier_2018_S3a.csv")
    subset_species = table_maier2018_s3a.columns[4:].str.split("(").str[0].str.strip()
    table_invitro_reference = table_maier2018_s3a.drop(["prestwick_ID","drug_class","n_hit"],axis=1)
  if subset_forslund:
    table_forslund_2021_s6c = pd.read_csv("data/forslund_2021_S6c.csv")
    table_forslund_2021_s6c = table_forslund_2021_s6c[table_forslund_2021_s6c["Feature space"].str.startswith("mOTU, ")]
    table_forslund_2021_s6c["Effector"] = table_forslund_2021_s6c["Effector"].str.split(" ").str[2].str.replace(r"\(|\)","")
    subset_species  = table_forslund_2021_s6c["Feature display name"].unique()
    table_invitro_reference = table_forslund_2021_s6c[["Effector", "Feature display name", "FDR"]].rename({"Effector":"chemical_name"},axis=1).pivot_table(columns="Feature display name", values="FDR", index="chemical_name").fillna(1).reset_index()

table_diagnoses = pd.read_csv(directory+"/PAMER_DM_ENC_DX.txt", sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3)
table_drg = pd.read_csv(directory+"/PAMER_DM_ENC_DRG.txt", sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3)
table_procedures = pd.read_csv(directory+"/PAMER_DM_ENC_PROC_ICD.txt", sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3)

table_indole = pd.read_csv('data/indole_v3.csv')
table_bile = pd.read_csv('data/bile_v3.csv')
table_scfa = pd.read_csv('data/scfa_v3.csv')
table_indole["compound"] = "indole#" + table_indole["compound"]
table_bile["compound"] = "bile#" + table_bile["compound"]
table_scfa["compound"] = "scfa#" + table_scfa["compound"]

if True:
  table_metab = pd.concat([table_indole, table_bile, table_scfa])
  table_metab = table_metab[table_metab.type=='normalized'].dropna()[['compound','metabolomicsID','value','batch']].drop_duplicates(['compound','metabolomicsID']).pivot(columns='metabolomicsID', values='value', index='compound')
  table_metab = table_metab.groupby(table_metab.index.str.replace(r'(-| )','').str.lower().str.replace('indole3acetic$','indole3aceticacid').str.replace('indole3propionic$','indole3propionicacid')).apply(lambda x: x.ffill().bfill().median())
  table_metab.columns = table_metab.columns.str.replace(".","_",regex=False).str.replace("-","_",regex=False).str.replace("U_C_","UC_",regex=False)
  table_metab_annot = pd.DataFrame(table_metab.index.str.split("_").str[0], index=table_metab.index.str.split("_").str[1:].str.join(""))
  table_metab_annot.columns=["table"]
  samples = pd.read_csv('data/redcap_tbl_v3.csv')
  samples["metabolomicsID"] = samples["metabolomicsID"].str.replace(".","_",regex=False).str.replace("_$","")
  samples[samples["metabolomicsID"].str.contains("DFI").fillna(False)]["metabolomicsID"]
  table_metab.columns = table_metab.columns.str.strip()
  table_metab_quant = pd.read_csv("data/Quant_Metab_Data.csv")
  table_metab_quant = table_metab_quant.pivot(columns="#metabolomics_id",index="compound",values="value_mm")
  table_metab_quant.index = "quant#"+table_metab_quant.index.str.replace("-","").str.replace(" ","")
  table_metab_quant_mapping = pd.read_csv("data/Clin_Projects_Sample_IDs_mapping.csv")
  print("samples in sample list but not in normalized metabolites:", pd.Index(samples["metabolomicsID"]).difference(table_metab.columns))
  print("samples in normalized metabolites but not in sample list:", table_metab.columns.difference(pd.Index(samples["metabolomicsID"])))
  print("samples in sample list but not in quant metabolites:", pd.Index(samples["metabolomicsID"]).difference(table_metab_quant.columns))
  print("samples in quant metabolites but not in sample list:", table_metab_quant.columns.difference(pd.Index(samples["metabolomicsID"])))
  table_metab = table_metab.groupby(table_metab.columns,axis=1).mean()
  table_metab_quant = table_metab_quant.groupby(table_metab_quant.columns,axis=1).mean()
  table_metab = pd.concat([table_metab, table_metab_quant],axis=0)
  # healthy donors all have S in the metabolomics ID.

metab_to_panel_quant = pd.read_csv("data/Quant_Metab_Data.csv").drop_duplicates(["hmmf_panel","compound"])
metab_to_panel_quant["compound"] = metab_to_panel_quant["compound"].str.replace("-","").str.replace(" ","")
metab_to_panel_quant = metab_to_panel_quant.set_index("compound")[["hmmf_panel","pubchem_id"]]
metab_to_panel_quant["hmmf_panel"] = metab_to_panel_quant["hmmf_panel"].replace({"BileAcid":"bile", "PFBBr":"scfa", "Tryptophan":"indole"})
metab_to_panel_quant["pubchem_id"] = metab_to_panel_quant["pubchem_id"].fillna(-1).astype(int)
metab_to_panel_quant["compound_quant"] = "quant#"+metab_to_panel_quant.index
metab_to_panel_quant["compound_pfx"] = metab_to_panel_quant.hmmf_panel+"#"+metab_to_panel_quant.index

refdate = datetime.datetime(2000, 1, 1)
date_int = [(lambda date: datetime.datetime.strptime(str(date), '%Y-%m-%d') if str(date)!='nan' else None)(date) for date in samples.date_collected]
date_int = np.array([(date - refdate).days if date is not None else np.nan for date in date_int ])
samples['date_int'] = date_int

print("(%.2f) Processing encounters, medications, and demographics..."%(get_program_time()), flush=True)
table_demographics = pd.read_csv(directory+'PAMER_DM_PATIENT_DEMO.txt', sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3).set_index("mrn")
encs = pd.read_csv(directory+'PAMER_DM_ENC.txt', sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3)
encs = encs.drop(['min_service_date', 'max_service_date'], axis=1)
encs = encs.drop_duplicates()
encs = encs[~encs[['adm_date', 'disc_date']].isna().any(axis=1)]
encs["age_years"] = [relativedelta(datetime.datetime.strptime(row["adm_date"].split(" ")[0], "%Y-%m-%d"), datetime.datetime.strptime(table_demographics.loc[row["MRN"]]["birth_date"].split(" ")[0], "%Y-%m-%d")).years if row["MRN"] in table_demographics.index else np.nan for ri,row in encs.iterrows()]
encs['adm_date'] = encs['adm_date'].str.split('-| ')
encs['disc_date'] = encs['disc_date'].str.split('-| ')
encs['adm_date_int'] = [(datetime.datetime(int(x[0]), int(x[1]), int(x[2])) - datetime.datetime(2000, 1, 1)).days for x in encs['adm_date']] 
encs['disc_date_int'] = [(datetime.datetime(int(x[0]), int(x[1]), int(x[2])) - datetime.datetime(2000, 1, 1)).days for x in encs['disc_date']] 
min_date = encs.iloc[encs['adm_date_int'].argmin()]['adm_date']
min_date = datetime.datetime(int(min_date[0]), int(min_date[1]), int(min_date[2]))
min_date_int = encs['adm_date_int'].min()
encs['adm_date_int'] = encs['adm_date_int'] - min_date_int
encs['disc_date_int'] = encs['disc_date_int'] - min_date_int
samples['collect_date_int'] = samples['date_collected'].str.split('-')
samples['collect_date_int'] = [(datetime.datetime(int(row[0]), int(row[1]), int(row[2])) - datetime.datetime(2000, 1, 1)).days if not isinstance(row, float) else None for row in samples['collect_date_int']]
samples['collect_date_int'] = samples['collect_date_int'] - min_date_int

print("(%.2f) Processing medications..."%(get_program_time()), flush=True)
meds2['harmonized_generic'] = [harmonize_medication(row['med_name_generic'].upper()) if row['med_name_generic']==row['med_name_generic'] else harmonize_medication(row['medication_name']) for ri,row in meds2.iterrows()]
meds2['harmonized_generic_route'] = meds2['harmonized_generic'] + '_' + meds2['route']
meds2 = meds2[~meds2['take_med_dttm'].isna()]
meds2['year'] = meds2['take_med_dttm'].str.split('-').str[0].astype('int')
meds2['month'] = meds2['take_med_dttm'].str.split('-').str[1].astype('int')
meds2['day'] = meds2['take_med_dttm'].str.split('-| ').str[2].astype('int')
meds2['take_med_int'] = [(datetime.datetime(x['year'], x['month'], x['day']) - datetime.datetime(2000, 1, 1)).days for i,x in meds2.iterrows()]
meds2['take_med_int'] = meds2['take_med_int'] - min_date_int

print("(%.2f) Loading labs..."%(get_program_time()), flush=True)
if(os.path.exists("data/table_labs_agg.parquet")):
  table_labs = pd.read_parquet("data/table_labs_agg.parquet")
else:
  if os.path.isfile("data/labs.parquet"):
    table_labs = pd.read_parquet("data/labs.parquet")
  else:
    table_labs = pd.read_csv('labs_subset.txt', sep='|', encoding='windows-1252', engine='python', quotechar=None, quoting=3)
    table_labs.to_parquet("data/labs.parquet")
  # preprocess lab tests:
  print("(%.2f) Processing lab tests..."%(get_program_time()), flush=True)
  table_labs = table_labs[~table_labs["spec_take_time"].isna()]
  table_labs['spec_take_time'] = table_labs['spec_take_time'].str.split('-| ')
  table_labs["spec_take_time_int"] = [(datetime.datetime(int(x[0]), int(x[1]), int(x[2])) - datetime.datetime(2000, 1, 1)).days for x in table_labs['spec_take_time']]
  table_labs["spec_take_time_int"] = table_labs["spec_take_time_int"] - min_date_int
  # remove non-numeric lab entries and aggregate all labs taken the same day by using the median
  print("(%.2f) Aggregating lab tests taken on the same day..."%(get_program_time()), flush=True)
  table_labs = table_labs[pd.to_numeric(table_labs["ord_value"], errors='coerce').notnull()]
  table_labs["ord_value"] = pd.to_numeric(table_labs["ord_value"])
  table_labs = table_labs[["mrn","proc_name","component_name","spec_take_time_int","ord_value","reference_unit"]].groupby(["mrn","proc_name","component_name","spec_take_time_int"]).agg({"ord_value":"median", "reference_unit":pd.Series.mode})
  table_labs = table_labs.sort_values(["mrn","proc_name","component_name","spec_take_time_int"])
  table_labs["reference_unit"] = table_labs["reference_unit"].astype("str")
  table_labs.to_parquet("data/table_labs_agg.parquet")

print("(%.2f) Identifying hospital visits..."%(get_program_time()), flush=True)

del meds

try:
  hospital_visits = pd.read_pickle("data/hospital_visits.pickle")
except:
  print("(%.2f) Identifying hospital visits..."%(get_program_time()), flush=True)
  hospital_visits_list = []
  progressbar = ProgressBar(len(encs))
  for ei,encounter in encs[encs["encounter_EIO"]=="Inpatient"].reset_index().iterrows():
    progressbar.update(ei)
    collections = samples[(samples['mrn']==encounter['MRN']) & (samples['collect_date_int']>=encounter['adm_date_int']) & (samples['collect_date_int']<=encounter['disc_date_int'])]
    collections = collections.sort_values('collect_date_int')
    procedures = table_procedures[table_procedures["har"] == encounter["har"]]["icd_name"]
    diagnoses = table_diagnoses[table_diagnoses["har"] == encounter["har"]]["icd10_code"]
    visit = {'admitted':encounter['adm_date_int'], 'discharged':encounter['disc_date_int'],
      'mrn':encounter.MRN,
      'age_years':encounter["age_years"],
      'collection_date_int':list(collections['collect_date_int']),
      'sampleids':list(collections['sampleid']),
      'shotgun_collections':list(collections['shotgunSeq_id']),
      'ncollections':len(collections),
      'metabolomics_collections':list(collections['metabolomicsID']),
      'db':list(collections['db']),
      'diagnoses':list(diagnoses),
      'procedures':list(procedures)
    }
    hospital_visits_list += [visit]
  hospital_visits = pd.DataFrame(hospital_visits_list)
  last_length = np.nan
  nitrs = 0
  # combine overlapping hospital visits. there's probably a better way to do this
  # than using a loop, but it works for now.
  while len(hospital_visits) != last_length:
    last_length = len(hospital_visits)
    nitrs = nitrs+1
    hospital_visits = hospital_visits.sort_values(['mrn','admitted'])
    hospital_visits['group'] = ((hospital_visits['admitted']>hospital_visits['discharged'].shift()) | (hospital_visits['mrn'] != hospital_visits['mrn'].shift())).cumsum()
    hospital_visits = hospital_visits.reset_index(drop=True).groupby('group').agg(
      {'admitted':"min",
       'discharged':"max",
       'mrn':"min",
       'age_years':"min",
       'collection_date_int':"sum",
       'sampleids':"sum",
       'shotgun_collections':"sum",
       'metabolomics_collections':"sum",
       'db':"sum",
       'diagnoses':"sum",
       'procedures':"sum"})
  # a list of all of all unique inpatient hospital visits, after merging contiguous visits together.
  hospital_visits['ncollections'] = [len(np.unique(x)) for x in hospital_visits['sampleids']]
  hospital_visits["labs_ix"] = hospital_visits[["admitted","discharged","mrn"]].apply(lambda row: list(table_labs[(table_labs.index.get_level_values("mrn")==row["mrn"]) & (table_labs.index.get_level_values("spec_take_time_int")>=row["admitted"]) & (table_labs.index.get_level_values("spec_take_time_int")<=row["discharged"])].reset_index().values), axis=1)
  hospital_visits["labs_name"] = hospital_visits["labs_ix"].apply(lambda seq: [x[1]+" "+x[2] for x in seq])
  hospital_visits["labs_time_int"] = hospital_visits["labs_ix"].apply(lambda seq: [x[3] for x in seq])
  hospital_visits["labs_value"] = hospital_visits["labs_ix"].apply(lambda seq: [x[4] for x in seq])
  hospital_visits["labs_units"] = hospital_visits["labs_ix"].apply(lambda seq: [x[5] for x in seq])
  hospital_visits["meds_ix"] = hospital_visits[["admitted","discharged","mrn"]].apply(lambda row: list(meds2[(meds2["mrn"]==row["mrn"]) & (meds2["take_med_int"]>=row["admitted"]) & (meds2["take_med_int"]<=row["discharged"])].index), axis=1)
  hospital_visits.to_pickle("data/hospital_visits.pickle")

# all hospital visits, unique by (patient, admitted_date, discharged_date)
individuals = hospital_visits[(hospital_visits['discharged'] - hospital_visits['admitted']>0) & (hospital_visits['ncollections']>1)]  # inpatient visit length >= 1 day and >1 stool collection.
visits_filtered = pd.concat([row for row in individuals[['shotgun_collections', 'metabolomics_collections', 'collection_date_int', 'mrn', 'admitted', 'discharged', 'age_years','db', 'diagnoses', 'procedures']].apply(lambda visit: pd.DataFrame(
    {'shotgun_collection':visit['shotgun_collections'],
    'metabolomics_collection':visit['metabolomics_collections'],
    'collection_date_int':visit['collection_date_int'],
    'db':visit['db'],
    'diagnoses':" ".join(visit.diagnoses),
    'procedures':"|".join(visit.procedures),
    'visit':visit.name,
    'admitted':visit.admitted,
    'discharged':visit.discharged,
    'age_years':visit.age_years,
    'mrn':visit.mrn}).dropna(subset=["shotgun_collection", "metabolomics_collection"], how="all"), axis=1)], ignore_index=True)

visits_filtered['filename'] = visits_filtered['shotgun_collection'].apply(lambda filename: pd.NA if pd.isna(filename) else np.where(table_microbes.columns.str.contains(filename))[0])
visits_filtered.loc[~visits_filtered['filename'].isna(),'filename'] = visits_filtered['filename'][~visits_filtered['filename'].isna()].apply(lambda file_index: table_microbes.columns[file_index] if len(file_index)>0 else np.nan)
visits_filtered['filename_species'] = np.nan
visits_filtered.loc[~visits_filtered['filename'].isna(),'filename_species'] = visits_filtered['filename'][~visits_filtered['filename'].isna()].apply(lambda files: files[files.str.contains("species")]).apply(lambda files: files[0] if len(files)>0 else np.nan)
visits_filtered['filename_genus'] = np.nan
visits_filtered.loc[~visits_filtered['filename'].isna(),'filename_genus'] = visits_filtered['filename'][~visits_filtered['filename'].isna()].apply(lambda files: files[files.str.contains("genus")]).apply(lambda files: files[0] if len(files)>0 else np.nan)
visits_filtered['filename_metab'] = visits_filtered['metabolomics_collection'].apply(lambda filename: pd.NA if pd.isna(filename) else np.where(table_metab.columns.str.contains(filename))[0])
visits_filtered.loc[~visits_filtered['filename_metab'].isna(),'filename_metab'] = visits_filtered['filename_metab'][~visits_filtered['filename_metab'].isna()].apply(lambda file_index: table_metab.columns[file_index[0]] if len(file_index)>0 else np.nan)
visits_filtered.drop(["filename"], axis=1, inplace=True)

def prettify_medication_name(name, length=50):
  tokens = name.split("_")
  abbr_route = {
    "inhalation":"(Inhal)",
    "oral":"(Oral)",
    "intravenous":"(IV)",
    "nasal":"(Nas)",
    "injection":"(Inj)",
    "topical":"(Top)",
    "ophthalmic":"(Oph)",
    "subcutaneous":"(Subc)",
    "rectal":"(Rec)",
    "misc.(non-drug; combo route)":""
  }
  return abbr(tokens[0], length).title() + ((" "+abbr_route[tokens[1].lower()]) if tokens[1]!="unknown" else "")

def prettify_compound_medication_name(name, length=50, mode="all"):
  if mode=="trunc":
    len_first_name = len(prettify_medication_name(name.split("|")[0], length))
    return abbr("; ".join([prettify_medication_name(x, length) for x in name.split("|")]), max(length,len_first_name+2))
  if mode=="share":
    names = name.split("|")
    lengths = [int(np.ceil(length * (len(n) / len(name)))) for n in names]
    return "; ".join([prettify_medication_name(n, length) for n,length in zip(names,lengths)])
  if mode=="all":
    return "; ".join([prettify_medication_name(x, length) for x in name.split("|")])


all_results_startmed=None
lm_convergence=None
Xs=None
Ys=None
Ys_before=None
Ys_after=None
collection_lengths=None
microbe_clusters=None
correlation_threshold=None
tables_medication_started_decor=None
all_results_medications=None
table_pvalues=None
table_significant=None
table_coefs=None
table_results=None
table_results_synonyms=None

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
    meds_ref = meds2[["harmonized_generic_route","med_pharm_class","med_pharm_sub_class"]].drop_duplicates("harmonized_generic_route").set_index("harmonized_generic_route")
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

print("(%.2f) Loading results..."%(get_program_time()), flush=True)
results = Results()
results.load_simple_associations("data/metab_associations_10_20_30_quant")

# determine correlations between metabolites and microbial genera
if True:
  matrix_genus = results.table[(results.table.sampletype=="genus") & (results.table.modeltype=="binary") & (results.table.duration==0)].pivot(index="medication",columns="microbe",values="coef")
  matrix_metab = results.table[(results.table.sampletype=="metab") & (results.table.modeltype=="binary") & (results.table.duration==0)].pivot(index="medication",columns="microbe",values="coef")
  common_meds = matrix_genus.index.intersection(matrix_metab.index)
  matrix_genus = matrix_genus.loc[common_meds]
  matrix_metab = matrix_metab.loc[common_meds]
  all_pairs = list(itertools.product(matrix_genus.columns, matrix_metab.columns))
  table_genus_metab_correl = pd.DataFrame(index=np.arange(len((all_pairs))), columns=["genus","metab","statistic","pvalue"])
  progressbar = ProgressBar(len(table_genus_metab_correl))
  ti = 0
  for genus in matrix_genus.columns:
    for metab in matrix_metab.columns:
      progressbar.update(ti)
      test = scipy.stats.spearmanr(matrix_genus.loc[:,genus], matrix_metab.loc[:,metab])
      table_genus_metab_correl.iloc[ti]["genus"] = genus
      table_genus_metab_correl.iloc[ti]["metab"] = metab
      table_genus_metab_correl.iloc[ti]["statistic"] = test.statistic
      table_genus_metab_correl.iloc[ti]["pvalue"] = test.pvalue
      ti=ti+1

  table_genus_metab_correl = table_genus_metab_correl[~table_genus_metab_correl["statistic"].isna()].sort_values("pvalue")
  table_genus_metab_correl["significant"] = statsmodels.stats.multitest.multipletests(table_genus_metab_correl["pvalue"], method="fdr_bh", alpha=0.1)[0]

# load microbiome and metabolome pathways
if True:
  importlib.reload(pathways)
  table_metab_reindexed = table_metab.copy()
  table_metab_reindexed.index = table_metab_reindexed.index.str.split("#").str[1]
  print("(%.2f) Load pathway hierarchy..."%(get_program_time()), flush=True)
  pathways.load_pathway_hierarchy(table_pathways, results.table)
  print("(%.2f) Load metabolites..."%(get_program_time()), flush=True)
  pathways.load_metabolites(table_pathways, results.table, table_metab_reindexed)
  print("(%.2f) Pathway enrichment..."%(get_program_time()), flush=True)
  pathways.pathway_enrichment(table_pathways, results.table, table_metab_reindexed)

print("Identifying study intervals...", flush=True)
study_intervals = []
groups = visits_filtered.groupby('visit').filter(lambda x: len(x)>1).groupby('visit').apply(lambda x: (x.sort_values('collection_date_int')))
for mrn, visit_df in groups.groupby(level=0):
  visit_df = visit_df.drop_duplicates('collection_date_int')
  intervals = pd.DataFrame(columns=['mrn','date_collection1','date_collection2'])
  ixs = np.triu_indices(len(visit_df),1)
  intervals['date_collection1'] = list(visit_df['collection_date_int'].iloc[ixs[0]])
  intervals['shotgunSeq_species_id1'] = list(visit_df['filename_species'].iloc[ixs[0]])
  intervals['shotgunSeq_genus_id1'] = list(visit_df['filename_genus'].iloc[ixs[0]])
  intervals['metabolomics_id1'] = list(visit_df['filename_metab'].iloc[ixs[0]])
  intervals['date_collection2'] = list(visit_df['collection_date_int'].iloc[ixs[1]])
  intervals['shotgunSeq_species_id2'] = list(visit_df['filename_species'].iloc[ixs[1]])
  intervals['shotgunSeq_genus_id2'] = list(visit_df['filename_genus'].iloc[ixs[1]])
  intervals['metabolomics_id2'] = list(visit_df['filename_metab'].iloc[ixs[1]])
  intervals['diagnoses'] = visit_df['diagnoses'].iloc[0]
  intervals['procedures'] = visit_df['procedures'].iloc[0]
  intervals['mrn'] = visit_df['mrn'].iloc[0]
  intervals['age_years'] = visit_df['age_years'].iloc[0]
  intervals['admitted'] = visit_df['admitted'].iloc[0]
  intervals['discharged'] = visit_df['discharged'].iloc[0]
  intervals['visit'] = visit_df['visit'].iloc[0]
  intervals['db'] = visit_df['db'].iloc[0]
  # intervals["has_shotgun"] = (not pd.isna(intervals['shotgunSeq_id1'])) and (not pd.isna(intervals['shotgunSeq_id2']))
  # intervals["has_metab"] = (not pd.isna(intervals['metabolomics_id1'])) and (not pd.isna(intervals['metabolomics_id2']))
  study_intervals += [intervals]

study_intervals = pd.concat(study_intervals).reset_index(drop=True)
study_intervals["has_metagenomics_species"] = ~study_intervals[["shotgunSeq_species_id1", "shotgunSeq_species_id2"]].isna().any(axis=1)
study_intervals["has_metagenomics_genus"] = ~study_intervals[["shotgunSeq_genus_id1", "shotgunSeq_genus_id2"]].isna().any(axis=1)
study_intervals["has_metabolomics"] = ~study_intervals[["metabolomics_id1", "metabolomics_id2"]].isna().any(axis=1)
study_intervals["is_liver_disease"] = study_intervals["db"].isin(["LiverDisease", "LiverTransplant"])
study_intervals["is_heart_disease"] = study_intervals["db"].isin(["HeartTransplant"])
study_intervals["is_micu"] = study_intervals["db"].isin(["MICU"])
study_intervals["pathways_id1"] = study_intervals["shotgunSeq_genus_id1"].apply(lambda filename: list(table_pathways_columns_short[table_pathways_columns_short["microbe_filename"]==filename]["pathway_filename"]))
study_intervals["pathways_id1"] = study_intervals["pathways_id1"].apply(lambda filenames: filenames[0] if len(filenames)>0 else np.nan)
study_intervals["pathways_id2"] = study_intervals["shotgunSeq_genus_id2"].apply(lambda filename: list(table_pathways_columns_short[table_pathways_columns_short["microbe_filename"]==filename]["pathway_filename"]))
study_intervals["pathways_id2"] = study_intervals["pathways_id2"].apply(lambda filenames: filenames[0] if len(filenames)>0 else np.nan)
study_intervals["has_pathways"] = (~study_intervals["pathways_id1"].isna()) & (~study_intervals["pathways_id2"].isna())
study_intervals["is_readmission"] = study_intervals.apply(lambda row: visits_filtered[(visits_filtered["mrn"]==row.mrn) & (visits_filtered["discharged"] < row.admitted)].shape[0] > 0, axis=1)
study_intervals["days_since_admitted"] = study_intervals["date_collection1"] - study_intervals["admitted"]

# Determine the most common diagnoses
icd_lut = table_diagnoses.set_index("icd10_code")["dx_name"].drop_duplicates()
icd_general_lut = icd10.table_icd10[~icd10.table_icd10.index.str.contains("\\.")]
top_diagnoses = pd.DataFrame(pd.DataFrame(' '.join(list(study_intervals["diagnoses"])).split(" ")).value_counts())
top_diagnoses_general = pd.DataFrame(pd.DataFrame(' '.join(list(study_intervals["diagnoses"])).split(" "))[0].str.split(".").str[0].value_counts()).drop("")
top_diagnoses_general["Description"] = top_diagnoses_general.index.get_level_values(0).map(lambda x : icd10.table_icd10.loc[x]["description"])
top_procedures = table_procedures["icd_name"].value_counts().head(50)

# Determine the medications that are started between two stool collections.
print("Choosing medications...", flush=True)
durations = [(2,10), (10,20), (20,30)]
tables_medication_started = [None,None,None]
medication_column = ['harmonized_generic_route', 'PHARM_CLASS'][0]
print('Analyzing medications categorized by %s...'%medication_column, flush=True)

# Choose which medications to look at.
medications_invitro = []
if subset_maier or subset_forslund:
  medications_invitro =  meds2[meds2["route"].isin(["Oral","Intravenous"])][medication_column].str.split("_").str[0].drop_duplicates().str.lower()
  table_invitro_reference["chemical_name"] = table_invitro_reference["chemical_name"]
  correspondences = []
  for medname in table_invitro_reference["chemical_name"]:
    medname_abbrev = re.sub("( acid| salt| a| sulfate| chloride| sodium| maleate| bromide| fumarate| citrate| hydrochloride| dinitrate| monohydrate| hemifumarate| nitrate| dihydrate| tartrate| monohydrochloride)","",medname)
    for substr in medname.lower().split(" "):
      if (" "+medications_invitro+" ").str.contains(" "+substr+" ",regex=False).any():
        correspondences += [(medications_invitro[(" "+medications_invitro+" ").str.contains(" "+substr+" ",regex=False).fillna(False)].iloc[0], medname, medname_abbrev)]
  correspondences = pd.DataFrame(correspondences, columns=["uchicago_medication","maier_medication","maier_medication_harmonized"])
  medications_invitro = correspondences["uchicago_medication"].drop_duplicates().values
  all_medications = meds2[["harmonized_generic","harmonized_generic_route"]].value_counts()
  all_medications = all_medications[(all_medications.index.get_level_values("harmonized_generic").str.lower().isin(medications_invitro)) & (all_medications>40)]
  medications_invitro = all_medications.index.get_level_values("harmonized_generic_route")

medications_curated = pd.read_csv("data/model_include_medications.csv")
medications = meds2[(meds2["route"].isin(["Oral","Intravenous","Injection","Subcutaneous","Transdermal","Rectal","Intramuscular","Sublingual"])) & (~meds2[medication_column].isin(medications_curated["Medication"][~medications_curated["Exclude"].isna()]))][medication_column].value_counts().index[:150]
medications = medications.append(pd.Index(["FLUDROCORTISONE_Oral","DEXAMETHASONE_Oral"]))
print("Medications: ", list(medications))
medication_classes = meds2[["harmonized_generic_route","med_pharm_class","med_pharm_sub_class"]].drop_duplicates("harmonized_generic_route").set_index("harmonized_generic_route").loc[medications]

print("Determining medication intervals...", flush=True)
if True:
  medications_started = pd.DataFrame(0, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']].drop_duplicates()), columns=medications)
  progressbar = ProgressBar(len(medications_started)*3)
  for di,duration in enumerate(durations):
    medications_started = pd.DataFrame(0, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']].drop_duplicates()), columns=medications)
    effect_duration = duration
    for i,((mrn,date_int1,date_int2),row) in enumerate(medications_started.iterrows()):
      progressbar.update(i + (di*len(medications_started)))
      # medications the patient has taken in the last 90 days.
      medications_exclude = meds2[(meds2['mrn']==mrn) & (meds2['take_med_int'] < date_int1-2) & (meds2['take_med_int'] > date_int1-90)][medication_column].unique()
      medications_started_now = meds2[(meds2['mrn']==mrn)            # same patient.
        & (date_int2 - date_int1 < effect_duration[1])               # second stool collection is less than 30 days after first stool collection.
        & (date_int2 - date_int1 >= effect_duration[0])              # second stool collection is more than 7 days after first stool collection.
        & (meds2['take_med_int'] - date_int1 >= -2)                  # medication was taken up to 2 days before the first stool collection.
        & (date_int2 - meds2['take_med_int'] >= 2)                   # medication was taken up to 2 days before the second stool collection.
        & (~meds2[medication_column].isin(medications_exclude))]     # patient hasn't taken the medication before.
      medications_started_now = list(set(medications_started_now[medication_column]).intersection(set(medications)))
      medications_started.loc[mrn,date_int1,date_int2][medications_started_now] = 1.0
    medications_started = medications_started.fillna(0)
    tables_medication_started[di] = medications_started

tables_medication_started[0].to_parquet("data/table_medications_started_0.parquet")
tables_medication_started[1].to_parquet("data/table_medications_started_1.parquet")
tables_medication_started[2].to_parquet("data/table_medications_started_2.parquet")

progressbar = ProgressBar(len(medications_started))
medications_started_starttime = pd.DataFrame(np.nan, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']].drop_duplicates()), columns=medications)
medications_started_dose      = pd.DataFrame(0, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']].drop_duplicates()), columns=medications)
for i,((mrn,date_int1,date_int2),row) in enumerate(medications_started_starttime.iterrows()):
  progressbar.update(i)
  # medications the patient has taken in the last 90 days.
  medications_exclude = meds2[(meds2['mrn']==mrn) & (meds2['take_med_int'] < date_int1-2) & (meds2['take_med_int'] > date_int1-90)][medication_column].unique()
  medications_started_now = meds2[(meds2['mrn']==mrn)            # same patient.
    & (meds2['take_med_int'] - date_int1 >= -2)                  # medication was taken up to 2 days before the first stool collection.
    & (date_int2 - meds2['take_med_int'] >= 2)                   # medication was taken up to 2 days before the second stool collection.
    & (~meds2[medication_column].isin(medications_exclude))]     # patient hasn't taken the medication before.
  medication_starttime = medications_started_now[medications_started_now[medication_column].isin(medications)].groupby(medication_column)["take_med_int"].min()
  medication_dose = medications_started_now[medications_started_now[medication_column].isin(medications)].groupby(medication_column)["take_med_int"].nunique()
  medications_started_starttime.loc[mrn,date_int1,date_int2].loc[medication_starttime.index] = medication_starttime.values - date_int1
  medications_started_dose.loc[mrn,date_int1,date_int2].loc[medication_dose.index] = medication_dose.values


medications_started_dose = medications_started_dose.fillna(0)
n_intervals = (medications_started_dose!=0).sum().sort_values(ascending=False)

progressbar = ProgressBar(len(medications_started))
all_medications_taken = pd.DataFrame(0, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']].drop_duplicates()), columns=meds2["harmonized_generic_route"].value_counts().index)
for i,((mrn,date_int1,date_int2),row) in enumerate(all_medications_taken.iterrows()):
  progressbar.update(i)
  # medications the patient has taken in the last 90 days.
  medications_started_now = meds2[(meds2['mrn']==mrn)            # same patient.
    & (meds2['take_med_int'] - date_int1 >= -2)                  # medication was taken up to 2 days before the first stool collection.
    & (date_int2 - meds2['take_med_int'] >= 2)                   # medication was taken up to 2 days before the second stool collection.
    ]     # patient hasn't taken the medication before.
  medication_dose = medications_started_now.groupby(medication_column)["take_med_int"].nunique()
  all_medications_taken.loc[mrn,date_int1,date_int2].loc[medication_dose.index] = medication_dose.values

all_medications_taken = all_medications_taken.fillna(0)

progressbar = ProgressBar(len(medications_started))
medications_stopped = pd.DataFrame(np.nan, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']].drop_duplicates()), columns=medications)
for i,((mrn,date_int1,date_int2),row) in enumerate(medications_stopped.iterrows()):
  progressbar.update(i)
  # medications the patient has taken in the following 90 days.
  medication_taken_last_within_interval = meds2[(meds2['mrn']==mrn) & (meds2['take_med_int'] >= date_int1) & (meds2['take_med_int'] <= date_int2)][["harmonized_generic_route","take_med_int"]].groupby(["harmonized_generic_route"]).agg(list)["take_med_int"]
  medication_taken_first_after_interval = meds2[(meds2['mrn']==mrn) & (meds2['take_med_int'] >= date_int2)][["harmonized_generic_route","take_med_int"]].groupby(["harmonized_generic_route"]).agg(list)["take_med_int"]
  if len(medication_taken_first_after_interval) > 0 and len(medication_taken_last_within_interval) > 0:
    medication_taken_first_after_interval = medication_taken_first_after_interval.str[0]
    medication_taken_last_within_interval = medication_taken_last_within_interval.str[-1]
    medications_stopped_interval = (medication_taken_first_after_interval.fillna(10000000) - medication_taken_last_within_interval).dropna()
    medications_stopped_interval = medications_stopped_interval[medications_stopped_interval>30]
    medications_stopped_interval = medications_stopped_interval.loc[medications_stopped.columns.intersection(medications_stopped_interval.index)]
    medications_stopped_interval = medication_taken_last_within_interval.loc[medications_stopped_interval.index]
    medications_stopped.loc[mrn,date_int1,date_int2].loc[medications_stopped_interval.index] = medications_stopped_interval.values

# Determine correlations between variables and remove highly correlated covariates.
tables_medication_started_decor = [tables_medication_started[0].copy(), tables_medication_started[1].copy(), tables_medication_started[2].copy()]
if True:
  while True:
    tables_medication_started_decor_any = (tables_medication_started_decor[0] + tables_medication_started_decor[1] + tables_medication_started_decor[2])!=0
    corr = (tables_medication_started_decor_any).corr().fillna(0) * np.tri((tables_medication_started_decor_any).shape[1], k=-1)
    max_pair = corr.stack().index[np.argmax(np.abs(corr.values))]
    max_corr = corr.loc[max_pair[0], max_pair[1]]
    if max_corr > 0.7:
      for i in range(3):
        print("combining %s and %s. correlation: %.2f"%(max_pair[0], max_pair[1], max_corr))
        tables_medication_started_decor[i][max_pair[0] + "|" + max_pair[1]] = tables_medication_started_decor[i][[max_pair[0], max_pair[1]]].max(axis=1)
        tables_medication_started_decor[i] = tables_medication_started_decor[i].drop([max_pair[0], max_pair[1]], axis=1)
    else:
      break

print("Including diagnoses, procedures, and demographics...", flush=True)

# Determine patient diagnoses for each inpatient visit
use_general_diagnoses = False
if use_general_diagnoses:
  diagnoses_curated = pd.read_csv("data/model_include_diagnoses.csv")
  table_diagnoses_curated = table_diagnoses[~table_diagnoses.dx_name.isin(diagnoses_curated[~diagnoses_curated["Exclude"].isna()].dx_name)]
  include_diagnoses = table_diagnoses_curated["icd10_code"].str.split(".").str[0].value_counts().head(100).index
  table_diagnoses_interval = pd.DataFrame(0, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']]), columns=include_diagnoses)
  for i,(ri, row) in enumerate(table_diagnoses_interval.iterrows()):
    table_diagnoses_interval.iloc[i].loc[list(set([d.split(".")[0] for d in study_intervals.iloc[i]["diagnoses"].split(" ")]).intersection(include_diagnoses))] = 1
  table_diagnoses_interval.columns = [icd_general_lut.loc[name]["description"] for name in table_diagnoses_interval.columns]
else:
  diagnoses_curated = pd.read_csv("data/model_include_diagnoses.csv")
  table_diagnoses_curated = table_diagnoses[~table_diagnoses.dx_name.isin(diagnoses_curated[~diagnoses_curated["Exclude"].isna()].dx_name)]
  include_diagnoses = table_diagnoses_curated["icd10_code"].value_counts().head(100).index
  table_diagnoses_interval = pd.DataFrame(0, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']]), columns=include_diagnoses)
  for i,(ri, row) in enumerate(table_diagnoses_interval.iterrows()):
    table_diagnoses_interval.iloc[i].loc[list(set([d for d in study_intervals.iloc[i]["diagnoses"].split(" ")]).intersection(include_diagnoses))] = 1
  table_diagnoses_interval.columns = table_diagnoses_curated.drop_duplicates("icd10_code").set_index("icd10_code").loc[table_diagnoses_interval.columns]["dx_name"]

# Identify procedures for each inpatient visit
procedures_curated = pd.read_csv("data/model_include_procedures.csv")
include_procedures = procedures_curated[(procedures_curated["Exclude"].isna()) & (~procedures_curated["Category"].isna())]["icd_name"]
table_procedures_interval = pd.DataFrame(0, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']]), columns=include_procedures)
for i,(ri, row) in enumerate(table_procedures_interval.iterrows()):
  table_procedures_interval.iloc[i].loc[list(set(study_intervals.iloc[i]["procedures"].split("|")).intersection(include_procedures))] = 1

table_procedures_interval.columns = procedures_curated.set_index("icd_name").loc[table_procedures_interval.columns]["Category"]
table_procedures_interval = table_procedures_interval.groupby(table_procedures_interval.columns,axis=1).any().astype(float)

# Identify demographics for each inpatient visit
table_demographics_interval = pd.DataFrame(0, index=pd.MultiIndex.from_frame(study_intervals[['mrn','date_collection1','date_collection2']]), columns=["age"])
table_demographics_interval["age"] = study_intervals["age_years"].values
table_demographics_interval["sex"] = pd.merge(table_demographics_interval, table_demographics[["sex"]], how="left", on="mrn")["sex"].values
table_demographics_interval["sex"] = (table_demographics_interval["sex"] == "Female").astype("int")

# Regression model: log-fold change in relative microbiome abundances, 
#   Predictors: Indicator matrix of whether a medication was started within this interval.
#   Groups: the patient indicator matrix.

def spearmanr_pval(x,y):
  return scipy.stats.spearmanr(x,y)[1]

def draw_grid(grid, title=''):
  ordering_medication = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid), no_plot=True, color_threshold=-np.inf)['leaves']
  ordering_microbe = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid.T), no_plot=True, color_threshold=-np.inf)['leaves']
  grid = grid.iloc[ordering_medication,ordering_microbe]
  grid.index = [abbr(x, 80) for x in grid.index]
  fig, ax = plt.subplots(1,1,tight_layout=True)
  sns.heatmap(grid, xticklabels=grid.columns, yticklabels=grid.index, cmap='PiYG', center=0, ax=ax, cbar=None)
  ax.set_title(title, font='serif', fontsize='small')
  ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', font='serif', fontsize='xx-small')
  ax.set_yticklabels(ax.get_yticklabels(), font='serif', fontsize='xx-small')

# save associations to a file.
def save_simple_associations(all_results_startmed, results_stderror, lm_convergence, Xs, Ys, Ys_before, Ys_after, collection_lengths, tables_medication_started_decor, folder_name="data/all_associations"):
  if not os.path.isdir("%s"%(folder_name)):
    os.mkdir("%s"%(folder_name))
  all_results_startmed.to_parquet("%s/all_results_startmed.parquet"%(folder_name))
  results_stderror.to_parquet("%s/stderror.parquet"%(folder_name))
  pd.DataFrame(lm_convergence, columns=["converged"]).to_parquet("%s/lm_convergence.parquet"%(folder_name))
  # pdb.set_trace()
  pd.concat([pd.concat([Xs[key]], keys=[key], names=["duration","sampletype","medicationtype"]) for key in Xs]).to_parquet("%s/Xs.parquet"%(folder_name))
  pd.concat([pd.concat([Ys[key]], keys=[key], names=["duration","sampletype","medicationtype"]).reset_index() for key in list(Ys.keys())],ignore_index=True).to_parquet("%s/Ys.parquet"%(folder_name))
  for key in Ys_before.keys():
    Ys_before[key].to_pickle("%s/Ys_before_%d_%s.pickle"%(folder_name,key[0],key[1]))
  for key in Ys_after.keys():
    Ys_after[key].to_pickle("%s/Ys_after_%d_%s.pickle"%(folder_name,key[0],key[1]))
  with open("%s/collection_lengths.pickle"%(folder_name),"wb") as handle:
    pickle.dump(collection_lengths, handle, protocol=pickle.HIGHEST_PROTOCOL)
  # pd.concat([pd.concat([microbe_clusters[key]], keys=[key], names=["duration","sampletype"]) for key in microbe_clusters]).to_parquet("%s/microbe_clusters.parquet"%(folder_name))
  # pd.DataFrame([correlation_threshold], columns=["correlation_threshold"]).to_parquet("%s/correlation_threshold.parquet"%(folder_name))
  for key in range(len(tables_medication_started_decor)):
    tables_medication_started_decor[key].to_parquet("%s/tables_medication_started_decor_%d.parquet"%(folder_name,key))

print("(%.2f) Running models..."%(get_program_time()), flush=True)

def run_models(modeldir="data/metab_associations_10_20_30_quant"):
  if not os.path.exists(modeldir):
    os.mkdir(modeldir)
  cca_results = [{},{},{},{}]
  all_results_startmed = None
  results_stderror = None
  lm_convergence = []
  all_models = []
  Xs = {}
  Ys = {}
  Ys_before={}
  Ys_after={}
  collection_lengths={}
  microbe_clusters = {}
  correlation_threshold = 0.9
  do_lm = True
  cov_matrices = {}
  list_sampletypes = ["metab","genus","species","pathway"]
  list_medicationtypes = ["binary","dose"]
  results_stderror = pd.read_csv(modeldir+"/partial_results_stderr.csv",index_col=0)
  all_results_startmed = pd.read_csv(modeldir+"/partial_all_results_startmed.csv",index_col=0)
  lm_convergence = list(pd.read_csv(modeldir+"/partial_lm_convergence.csv",index_col=0).values.reshape(-1))
  for di,duration in enumerate(durations):
    for sampletype in list_sampletypes:
      for medicationtype in list_medicationtypes:
        key = (di,sampletype,medicationtype)
        print("\n",key, flush=True)
        if medicationtype=="binary":
          X = tables_medication_started_decor[di].loc[[(row['mrn'], row['date_collection1'], row['date_collection2']) for ri,row in study_intervals.iterrows()]]
        else:
          if di != 0:
            print("skipping.")
            continue
          if medicationtype=="dose":
            X1 = medications_started_dose.loc[[(row['mrn'], row['date_collection1'], row['date_collection2']) for ri,row in study_intervals.iterrows()]].copy()
            X2 = medications_started_dose.loc[[(row['mrn'], row['date_collection1'], row['date_collection2']) for ri,row in study_intervals.iterrows()]].copy()
            X1.columns += ["/Dose"]
            X2 = (X2!=0).astype("float")
            X = pd.concat([X1, X2], axis=1)
          elif medicationtype=="null":
            X = pd.DataFrame()
        X = X.copy()
        n_covariates = len(X.columns)
        X = pd.concat([X,table_diagnoses_interval, table_procedures_interval, table_demographics_interval], axis=1)
        X["is_liver_disease"] = study_intervals["is_liver_disease"].values.astype("float")
        X["is_heart_disease"] = study_intervals["is_heart_disease"].values.astype("float")
        X["is_micu"] = study_intervals["is_micu"].values.astype("float")
        X["is_readmission"] = study_intervals["is_readmission"].astype("float")
        X["days_since_admitted"] = study_intervals["days_since_admitted"].astype("float")
        # log change in relative abundance.
        if sampletype == "genus":
          X = X.loc[study_intervals["has_metagenomics_genus"].values]
          Y = np.log2((0.01+table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_genus"]]['shotgunSeq_genus_id2']]) / (0.01+table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_genus"]]['shotgunSeq_genus_id1']].values)).T
          alpha1 = table_microbes[study_intervals[study_intervals["has_metagenomics_genus"]]['shotgunSeq_genus_id1']].apply(lambda row: scipy.stats.entropy(row))
          alpha2 = table_microbes[study_intervals[study_intervals["has_metagenomics_genus"]]['shotgunSeq_genus_id2']].apply(lambda row: scipy.stats.entropy(row))
          Y_mask = (((table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_genus"]]['shotgunSeq_genus_id1'].drop_duplicates(keep="first")]>0).sum(axis=1)>5) & ((table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_genus"]]['shotgunSeq_genus_id2'].drop_duplicates(keep="first")]>0).sum(axis=1)>5))
          Y = Y.loc[:,Y_mask]
          Y["alphadiversity"] = alpha2.values - alpha1.values
        elif sampletype == "species":
          X = X.loc[study_intervals["has_metagenomics_species"].values]
          Y = np.log2((0.01+table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_species"]]['shotgunSeq_species_id2']]) / (0.01+table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_species"]]['shotgunSeq_species_id1']].values)).T
          alpha1 = table_microbes[study_intervals[study_intervals["has_metagenomics_species"]]['shotgunSeq_species_id1']].apply(lambda row: scipy.stats.entropy(row))
          alpha2 = table_microbes[study_intervals[study_intervals["has_metagenomics_species"]]['shotgunSeq_species_id2']].apply(lambda row: scipy.stats.entropy(row))
          Y_mask = (((table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_species"]]['shotgunSeq_species_id1'].drop_duplicates(keep="first")]>0).sum(axis=1)>5) & ((table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_species"]]['shotgunSeq_species_id2'].drop_duplicates(keep="first")]>0).sum(axis=1)>5))
          Y = Y.loc[:,Y_mask]
          Y["alphadiversity"] = alpha2.values - alpha1.values
        elif sampletype=="metab":
          X = X.loc[study_intervals["has_metabolomics"].values]
          Y = np.log2((0.01+table_metab[study_intervals[study_intervals["has_metabolomics"]]['metabolomics_id2']]) / (0.01+table_metab[study_intervals[study_intervals["has_metabolomics"]]['metabolomics_id1']].values)).T
          Y_mask = (((table_metab[study_intervals[study_intervals["has_metabolomics"]]['metabolomics_id1'].drop_duplicates(keep="first")]>0).sum(axis=1)>5) & ((table_metab[study_intervals[study_intervals["has_metabolomics"]]['metabolomics_id2'].drop_duplicates(keep="first")]>0).sum(axis=1)>5))
          Y = Y.loc[:,Y_mask]
        elif sampletype=="pathway":
          X = X.loc[study_intervals["has_pathways"].values]      
          Y = np.log2((0.01+table_pathways_general[study_intervals[study_intervals["has_pathways"]]['pathways_id2']]) / (0.01+table_pathways_general[study_intervals[study_intervals["has_pathways"]]['pathways_id1']].values)).T
          Y_mask = (((table_pathways_general[study_intervals[study_intervals["has_pathways"]]['pathways_id1'].drop_duplicates(keep="first")]>0).sum(axis=1)>5) & ((table_pathways_general[study_intervals[study_intervals["has_pathways"]]['pathways_id2'].drop_duplicates(keep="first")]>0).sum(axis=1)>5))
          Y = Y.loc[:,Y_mask]
        mask = X.sum(axis=1)>0    # remove medications that are all zero.
        Y = Y.loc[mask.values]
        X = X.loc[mask.values]
        X = X.loc[:,(X.var(axis=0)>0).values] # remove medications that are all zero or all one. 
        Y = Y.loc[:,((Y.var(axis=0)>0)).values]
        Xs[key] = X
        Ys[key] = Y
        X = X.loc[:,X.apply(lambda col: col.value_counts().max() < len(col)-10)]
        X_corr = (X.corr()-np.eye(X.shape[1]))
        X_corr_abs = (X.corr()-np.eye(X.shape[1])).abs()
        X = X.loc[:,(X_corr_abs * np.tri(X.shape[1])).max()<correlation_threshold]  # remove highly correlated covariates.
        if sampletype=="genus":
          Ys_before[key] = table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_genus"]]['shotgunSeq_genus_id1']].loc[:,mask.values]
          Ys_after[key] = table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_genus"]]['shotgunSeq_genus_id2']].loc[:,mask.values]
        elif sampletype=="species":
          Ys_before[key] = table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_species"]]['shotgunSeq_species_id1']].loc[:,mask.values]
          Ys_after[key] = table_microbes_relabund[study_intervals[study_intervals["has_metagenomics_species"]]['shotgunSeq_species_id2']].loc[:,mask.values]
        elif sampletype=="metab":
          Ys_before[key] = table_metab[study_intervals[study_intervals["has_metabolomics"]]['metabolomics_id1']].loc[:,mask.values]
          Ys_after[key] = table_metab[study_intervals[study_intervals["has_metabolomics"]]['metabolomics_id2']].loc[:,mask.values]
        elif sampletype=="pathway":
          Ys_before[key] = table_pathways_general[study_intervals[study_intervals["has_pathways"]]['pathways_id1']].loc[:,mask.values]
          Ys_after[key] = table_pathways_general[study_intervals[study_intervals["has_pathways"]]['pathways_id2']].loc[:,mask.values]
        collection_lengths[key] = (X.index.get_level_values("date_collection2") - X.index.get_level_values("date_collection1")).values
        print("n models: ", len(Y.columns), "n tests: ", len(Y.columns)*n_covariates)
        if do_lm:
          groups = X.index.get_level_values(0)
          progressbar = ProgressBar(len(Y.columns))
          for mi,microbe in enumerate(Y.columns):
            colname_pvalue = 'pvalue_'+microbe+'_'+str(di)+'_'+sampletype+"_"+medicationtype
            colname_coef = 'coef_'+microbe+'_'+str(di)+'_'+sampletype+"_"+medicationtype
            if (all_results_startmed is not None) and (colname_pvalue in all_results_startmed):
              print("skipping duplicate.")
              continue
            # paused at mi = 183
            progressbar.update(mi)
            outcome = Y[microbe].reset_index(drop=True)
            X_sub = X.reset_index(drop=True)[~outcome.isna()]
            Y_sub = outcome[~outcome.isna()]
            group_sub = groups[~outcome.isna()]
            X_sub = X_sub.loc[:,X_sub.var()!=0]
            if np.linalg.eigvalsh(X_sub.cov()).min() <= 1e-10:
              X_subset = X_sub.iloc[:,0]
              for i in range(1, X_sub.shape[1]):
                X_newcol = pd.concat([X_subset, X_sub.iloc[:,i]],axis=1)
                # s = np.linalg.svd(X_newcol.cov(), compute_uv=False)
                # if s.min() > 1e-10:
                # print(i,s)
                if np.linalg.eigvalsh(X_newcol.cov()).min() > 1e-10:
                  X_subset = X_newcol
              X_sub = X_subset
            try:
              model = sm.MixedLM(Y_sub,  sm.add_constant(X_sub), group_sub).fit()
            except: 
              model = sm.MixedLM(Y_sub,  (X_sub), group_sub).fit()
            results = pd.DataFrame({colname_pvalue:model.pvalues, colname_coef:model.params})
            key_filename_str = "%d_%s_%s_%s"%(di,sampletype,medicationtype,microbe)
            key_filename_str = re.sub("[^a-zA-Z0-9_]","-",key_filename_str)
            model.cov_params().to_csv(modeldir+"/covariance_%s"%key_filename_str)
            results_stderror = pd.concat([results_stderror, model.bse], axis=1)
            all_results_startmed = pd.concat([all_results_startmed, results], axis=1)
            lm_convergence += [model.converged]
        print("(%.2f) Saving partial results. Done with %d,%s,%s..."%(get_program_time(), di, sampletype, medicationtype), flush=True)
        results_stderror.to_csv(modeldir+"/partial_results_stderr.csv")
        all_results_startmed.to_csv(modeldir+"/partial_all_results_startmed.csv")
        pd.DataFrame(lm_convergence).to_csv(modeldir+"/partial_lm_convergence.csv")
  results_stderror.columns = all_results_startmed.columns[all_results_startmed.columns.str.startswith("coef_")].str.split("coef_").str[1]
  print("(%.2f) Saving simple associations..."%(get_program_time()), flush=True)
  save_simple_associations(all_results_startmed, results_stderror, lm_convergence, Xs, Ys, Ys_before, Ys_after, collection_lengths, tables_medication_started_decor, modeldir)
  print("done.")

def update_figures():
  importlib.reload(paperfigures)
  paperfigures.table_microbes = table_microbes
  paperfigures.table_microbes_clr = table_microbes_clr
  paperfigures.samples = samples
  paperfigures.table_demographics = table_demographics
  paperfigures.hospital_visits = hospital_visits
  paperfigures.study_intervals = study_intervals
  paperfigures.min_date = min_date
  paperfigures.meds2 = meds2
  paperfigures.medications = medications
  paperfigures.abbr = abbr
  paperfigures.prettify_medication_name = prettify_medication_name
  paperfigures.prettify_compound_medication_name = prettify_compound_medication_name
  paperfigures.results = results
  paperfigures.durations = durations
  paperfigures.medications_started_dose = medications_started_dose
  paperfigures.medications_started_starttime = medications_started_starttime
  paperfigures.table_microbes_relabund = table_microbes_relabund
  paperfigures.table_pathways = table_pathways
  paperfigures.pathways = pathways
  paperfigures.table_genus_metab_correl = table_genus_metab_correl
  paperfigures.matrix_genus = matrix_genus
  paperfigures.matrix_metab = matrix_metab
  paperfigures.table_invitro_reference = table_invitro_reference
  paperfigures.correspondences = correspondences
  paperfigures.table_pathways_general = table_pathways_general
  paperfigures.make_figures()
  paperfigures.make_supplementary_figures()

def compute_estimates_return_to_baseline():
  # We use Fieller's theorem to estimate a 90% confidence interval for the ratio of coefficients"
  #   https://en.wikipedia.org/wiki/Fieller%27s_theorem
  #   We use the formulas from section 2.3 (page 6) from here: https://fbe.unimelb.edu.au/__data/assets/pdf_file/0005/2704649/2037_Joe-Hirschberg_Fieller-Examples-WP-version.pdf
  #   Confidence  Intervals for Ratios: Econometric Examples  in  Stata and R. Jenny Lye and Joe Hirschberg. Feb 2018. Revised July 2018 .
  time_return_baseline = pd.DataFrame(columns=["estimate","low","high","low_fieller","high_fieller","intercept","slope"])
  for medication in medications:
    di = 0
    sampletype="genus"
    modeltype="dose"
    microbe="alphadiversity"
    table = pd.read_csv("data/metab_associations_10_20_30_quant/covariance_%d_%s_%s_%s"%(di,sampletype,modeltype,microbe),index_col=0)
    if medication not in table.index or medication+"/Dose" not in table.index:
      continue
    covariance = table.loc[[medication,medication+"/Dose"],[medication,medication+"/Dose"]]
    v11 = covariance.iloc[0,0]
    v12 = covariance.iloc[0,1]
    v21 = covariance.iloc[1,0]
    v22 = covariance.iloc[1,1]
    # time to return back to baseline is intercept / slope = (b/a)
    a = results.table[(results.table.medication == medication) & (results.table.modeltype=="dose") & (results.table.sampletype==sampletype) & (results.table.microbe==microbe) & (results.table.medication_dose==False)]["coef"] # intercept
    b = results.table[(results.table.medication == medication) & (results.table.modeltype=="dose") & (results.table.sampletype==sampletype) & (results.table.microbe==microbe) & (results.table.medication_dose==True)]["coef"]  # slope
    if len(a)==0 or len(b) == 0:
      continue
    a = a.item()
    b = b.item()
    # interested in a / b.
    estimate = a / b
    estimate2 = estimate*estimate
    nsamples = results.Xs[(di,sampletype,modeltype)].shape[0]
    npredictors = results.Xs[(di,sampletype,modeltype)].shape[1]
    alpha = 0.10    # 90% confidence interval
    s = 0.1
    t = scipy.stats.t.ppf(q=1.0 - (alpha/2.0), df=nsamples-npredictors)
    # Delta method:
    radius = t * (1/b) * np.sqrt(v11 + (estimate2*v22) - (2*estimate*v12))
    interval = np.array([-estimate-np.abs(radius), -estimate+np.abs(radius)])
    # Fieller's method:
    c0 = (b*b) - (t*t*v22)
    c1 = 2*((t*v12) - (a*b))
    c2 = (a*a) - (t*t * v11)
    fieller = [(-c1/(2*c0)) - (np.sqrt((c1*c1) - (4*c0*c2))/(2*c0)), (-c1/(2*c0)) + (np.sqrt((c1*c1) - (4*c0*c2))/(2*c0))]
    fieller = [-fieller[1], -fieller[0]]
    fieller = [np.min(fieller), np.max(fieller)]
    # g = (t*t * s*s * v22) / (b*b)
    # mL = (1/(1-g)) * ((a/b) - (g*v12/v22) - (t*s/b) * np.sqrt(v11 - (2*(a/b)*v12) + ((a*a)/(b*b))*v22 - (g*(v11 - (v12*v12/v22)))) )
    # mU = (1/(1-g)) * ((a/b) - (g*v12/v22) + (t*s/b) * np.sqrt(v11 - (2*(a/b)*v12) + ((a*a)/(b*b))*v22 - (g*(v11 - (v12*v12/v22)))) )
    # interval = (mL, mU)
    time_return_baseline.loc[medication] = [-a/b, interval[0], interval[1], fieller[0], fieller[1], a, b]
  time_return_baseline["significant"] = ~((0 < time_return_baseline[["low","high"]].max(axis=1)) & (0 > time_return_baseline[["low","high"]].min(axis=1)))
  return time_return_baseline

if __name__ == "__main__":
  run_models()
