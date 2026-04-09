from sklearn.decomposition import PCA
import matplotlib.pyplot as pl
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
import plots
import itertools
import matplotlib.colors
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pypdf
import io
import reportlab
import reportlab.pdfgen
import reportlab.lib.pagesizes
import reportlab.pdfbase
import reportlab.pdfbase.pdfmetrics
import reportlab.pdfbase.ttfonts
import reportlab.pdfgen.canvas
from sklearn.preprocessing import MultiLabelBinarizer
import pdf2image
import pathways

from util import *


durations = [(2, 10), (10, 20), (20, 30)]
# initialize fonts
font_variants = ("Arial","Arial Italic","Arial Bold")
folder = '/System/Library/Fonts/Supplemental/'
for variant in font_variants:
  reportlab.pdfbase.pdfmetrics.registerFont(reportlab.pdfbase.ttfonts.TTFont(variant, os.path.join(folder, variant+'.ttf'))) 

def dot_heatmap(data, ax_heatmap, ax_legend, cmap1, cmap2, class_labels=None, size=10.0, min_size=2, exponent=0.15, legend_values=[-0.3, -0.05, -0.001, 0.001, 0.05, 0.3], legend_labels=["-0.3", "-0.05", "-0.001", "0.001", "0.05", "0.3"], cmap_skip=0.2, text_fontsize=8, classlabel_lw=1, classlabel_line_offset=0.25, classlabel_text_offset=0.6, legend_kws={}, marker="o", onlyplot=False):
  grid = data.copy()
  grid[grid==0]=np.nan
  # hm = sns.heatmap(grid[microbe_columns[ordering_microbe]].astype(float), ax=ax_heatmap, cbar=False,cmap=cmap, center=0, linewidths=0.4, linecolor='#d0d0d0')
  grid_sc = grid.reset_index(drop=True).T.reset_index(drop=True).T.astype(float).stack().reset_index(name="corr")
  top = cmap1
  bottom = cmap2
  newcolors = np.vstack((top(np.linspace(0, 1.0-cmap_skip, 128)), bottom(np.linspace(cmap_skip, 1, 128))))
  newcmp = ListedColormap(newcolors, name='OrangeBlue')
  vmin = grid.min().min()
  vmax = grid.max().max()
  vmax = max(-vmin,np.abs(vmax))
  vmin = min(vmin,-vmax)
  cmap_norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=False)
  cmap_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=newcmp.copy())
  # legend
  legend_corrs = legend_values
  if len(legend_values) > 0:
    legend_minval = np.min([grid_sc["corr"].abs().min(), np.min(legend_values)])
  else:
    legend_minval = grid_sc["corr"].abs().min()
  legend_sizes = pow(2 + (pow(np.abs(legend_corrs) - legend_minval, exponent) * size), 2.0)+min_size
  # print("legend sizes: ", legend_sizes)
  legend_colors = cmap_map.to_rgba(legend_corrs)
  cmap_map.to_rgba(0)
  # scatterplot
  sizes = pow(2 + (pow(grid_sc["corr"].abs() - legend_minval, exponent) * size), 2.0)+min_size
  size_to_value = lambda s: pow((pow(x-min_size,1/2.0)-2)/size,1.0/exponent)
  sc = ax_heatmap.scatter(grid_sc["level_1"], grid_sc["level_0"], marker=marker, c=cmap_map.to_rgba(grid_sc["corr"]), s=sizes, zorder=1, linewidths=0)
  if ax_legend is not None:
    leg = ax_legend.legend([matplotlib.lines.Line2D([0], [0], ls="", color=legend_colors[i], markerfacecolor=legend_colors[i], ms=np.sqrt(legend_sizes[i]), marker=sc.get_paths()[0]) for i in range(len(legend_corrs))], legend_labels, numpoints=1, markerscale=1, **legend_kws)
    ax_legend.add_artist(leg)
  if onlyplot:
    return
  ax_heatmap.set_axisbelow(True)
  ax_heatmap.xaxis.grid(zorder=0, alpha=0.2)
  ax_heatmap.yaxis.grid(zorder=0, alpha=0.2)
  ax_heatmap.set_title("Most Significant Genus-Medication Associations")
  ax_heatmap.set_xticks(np.arange(len(grid.columns)))
  ax_heatmap.set_xticklabels(grid.columns, fontsize=8, rotation=45, ha="right", rotation_mode="anchor")
  ax_heatmap.set_yticks(np.arange(len(grid.index)))
  ax_heatmap.set_yticklabels(grid.index, fontsize=8)
  ax_heatmap.set_xlabel("Genus", fontsize=11)
  ax_heatmap.set_ylabel("Medication", fontsize=11)
  # line1 = matplotlib.lines.Line2D([], [], color="white", marker='o', markerfacecolor="red")
  # line2 = matplotlib.lines.Line2D([], [], color="white", marker='o', markerfacecolor="green")
  # line3 = matplotlib.lines.Line2D([], [], color="white", marker='o', markersize=5,  markerfacecolor="slategray")
  # line4 = matplotlib.lines.Line2D([], [], color="white", marker='o', markersize=10, markerfacecolor="slategray")
  # handles, labels = sc.legend_elements(prop="sizes", alpha=0.6, func=lambda s: s/10.0)
  # legend_values = [float(re.split("{|}",x)[1]) for x in labels]
  # for handle,label,value in zip(handles,labels,legend_values):
  #   handle.set_color(newcmp(value))
  #   handle.set_linewidth(0)
  # ax_legend.legend(handles, labels, numpoints=1, markerscale=1, bbox_to_anchor=(1, 0.5), loc="center left", **legend_kws)
  # ax_legend.axis("off")
  # handles, labels = legend_elts = sc.legend_elements("sizes", num=3)
  coord_y = -0.5
  text_offset = 0
  if class_labels is not None:
    for class_label in class_labels:
      height = class_label[1]
      line = ax_heatmap.plot([grid.shape[1]+classlabel_line_offset+text_offset, grid.shape[1]+classlabel_line_offset+text_offset], [coord_y+0.1, coord_y+height-0.25], clip_on=False, color="black", lw=classlabel_lw)
      ax_heatmap.text(grid.shape[1]+classlabel_text_offset+text_offset, coord_y+height-0.6, class_label[0].title(), horizontalalignment='left', verticalalignment='center', weight="bold", color="#505050", fontsize=text_fontsize)
      print([-4, -4], [coord_y, coord_y+height])
      coord_y = coord_y+height
  ax_heatmap.set_xlim(0-1, grid.shape[1])
  ax_heatmap.set_ylim(0-1, grid.shape[0])


from matplotlib.colors import ListedColormap, LinearSegmentedColormap
def make_supplementary_figure_pcaplots():
  table_microbes_genus = table_microbes.loc[:,table_microbes_clr.columns.str.endswith("genus.txt")].copy().T
  table_pathways_temp = table_pathways_general.copy().T
  table_pathways_temp = table_pathways_temp.loc[:,(table_pathways_temp!=0).mean(axis=0)>0.05]
  table_microbes_genus.index = table_microbes_genus.index.str.split("_R1").str[0]
  table_pathways_temp.index = table_pathways_temp.index.str.split("_mpa").str[0]
  common_indices = table_microbes_genus.index.intersection(table_pathways_temp.index)
  table_pathways_temp = table_pathways_temp.loc[common_indices]
  diversity = table_microbes_genus.apply(lambda row: scipy.stats.entropy(row), axis=1)
  table_microbes_genus_clr = table_microbes_genus.copy()
  table_microbes_genus_clr.iloc[:,:] = skbio.stats.composition.clr(table_microbes_genus.T+1.0).T
  pca_microbes = PCA(n_components=4).fit(table_microbes_genus_clr)
  pca_microbes.explained_variance_ratio_
  pca_microbes_proj = pca_microbes.transform(table_microbes_genus_clr)
  table_pathways_clr = table_pathways_temp.copy().drop(index=[None],errors="ignore")
  table_pathways_clr.iloc[:,:] = skbio.stats.composition.clr(table_pathways_clr.T+1.0).T
  pca_pathways = PCA(n_components=4).fit(table_pathways_clr)
  pca_pathways.explained_variance_ratio_
  pca_pathways_proj = pca_pathways.transform(table_pathways_clr)
  pathways_diversity = diversity.loc[common_indices]
  # table_microbes_genus = table_microbes.loc[:,table_microbes_clr.columns.str.endswith("genus.txt")].copy().T
  shotgun_to_sample = table_microbes.columns.to_series().apply(lambda filename: samples["shotgunSeq_id"].apply(lambda id: str(id) in filename))
  shotgun_to_sample = shotgun_to_sample.apply(lambda row: samples.loc[np.where(row.values)].assign(filename=row.name), axis=1)
  shotgun_to_sample = pd.concat(shotgun_to_sample.values).set_index("filename")
  shotgun_to_sample = shotgun_to_sample[~shotgun_to_sample["mrn"].isna()]
  shotgun_to_sample["mrn"] = shotgun_to_sample["mrn"].astype("int")
  shotgun_to_sample.index = shotgun_to_sample.index.str.split("_R1").str[0]
  # include_taxa = table_microbes.sum(axis=0).sort_values(ascending=False).head(8).index
  # table_microbes_barplot = table_microbes.loc[:,include_taxa]
  # table_microbes_barplot["Other"] = table_microbes.loc[:,~table_microbes.columns.isin(include_taxa)].sum(axis=1).values
  # table_microbes_barplot = table_microbes_barplot.div(table_microbes_barplot.sum(axis=1),axis=0)
  diversity = table_microbes_genus.apply(lambda row: scipy.stats.entropy(row), axis=1)
  # table_microbes_genus.index.map(lambda filename: samples[])
  demo_samples = table_microbes_genus.index.to_series().apply(lambda filename: samples[samples["shotgunSeq_id"].fillna("NAN").apply(lambda id: id in filename)])
  # demo_samples = table_microbes_genus.index.to_series().apply(lambda filename: samples[samples["shotgunSeq_id"].str.contains(filename, na=False)])
  demo_mrn = demo_samples.apply(lambda df: ((df.iloc[0].mrn) if (len(df)>0) else 0))
  demo_mrn[demo_mrn.isna()] = 0
  demo_mrn = demo_mrn.astype("int")
  demo_sex = demo_mrn.apply(lambda mrn: table_demographics.loc[mrn]["sex"] if mrn in table_demographics.index else None)
  demo_birthdate = demo_mrn.apply(lambda mrn: table_demographics.loc[mrn]["birth_date"].split(" ")[0] if mrn in table_demographics.index else None)
  demo_birthyear = demo_birthdate.str.split("-").str[0]
  demographics = pd.concat([demo_mrn.rename("mrn"), demo_sex.rename("sex"), demo_samples.rename("samples"), demo_birthyear.rename("birthyear")], axis=1)
  demographics["age"] = demographics.apply(lambda row: float(str(row.samples.iloc[0].date_collected).split("-")[0]) - float(np.nan if row.birthyear is None else row.birthyear) if len(row.samples)>0 else None, axis=1)
  demographics["cohorts"] = demographics["samples"].apply(lambda df: " ".join([str(x) for x in np.sort(np.unique(df["db"]))]))
  # days hospitalized
  collection_to_visit = hospital_visits.copy()
  collection_to_visit = collection_to_visit[["mrn","admitted","discharged","collection_date_int","shotgun_collections"]]
  collection_to_visit["keep_indices"] = collection_to_visit.apply(lambda row: ~((pd.Series(row["shotgun_collections"]).duplicated(keep="first")) | pd.Series(row["shotgun_collections"]).isna()).values, axis=1)
  collection_to_visit["collection_date_int"] = collection_to_visit.apply(lambda row: np.array(row["collection_date_int"])[row["keep_indices"]], axis=1)
  collection_to_visit["shotgun_collections"] = collection_to_visit.apply(lambda row: np.array(row["shotgun_collections"])[row["keep_indices"]], axis=1)
  collection_to_visit = collection_to_visit[collection_to_visit["shotgun_collections"].str.len()>0]
  collection_to_visit["shotgun_collections_dates"] = collection_to_visit.apply(lambda row: list(zip(row["collection_date_int"].astype(int), row["shotgun_collections"])), axis=1)
  collection_to_visit = collection_to_visit.explode("shotgun_collections_dates")
  collection_to_visit["collection_date_int"] = collection_to_visit["shotgun_collections_dates"].str[0]
  collection_to_visit["shotgun_collections"] = collection_to_visit["shotgun_collections_dates"].str[1]
  collection_to_visit["days_since_admitted"] = collection_to_visit["collection_date_int"] - collection_to_visit["admitted"]
  collection_to_visit = collection_to_visit.drop_duplicates("shotgun_collections").set_index("shotgun_collections",drop=False)
  demographics["days_since_admitted"] = demographics.samples.apply(lambda row: collection_to_visit.loc[collection_to_visit.index.intersection(row["shotgunSeq_id"].values)].days_since_admitted.values)
  demographics["days_since_admitted"] = demographics["days_since_admitted"].apply(lambda row: row[0] if len(row)>0 else np.nan)
  year = demographics.apply(lambda row: row.samples.iloc[0]["date_collected"] if len(row.samples)>0 else np.nan, axis=1)
  year[~year.isna()] = year[~year.isna()].apply(lambda x: x.split("-")[0])
  demographics["year"] = year.astype(float)
  sys.setrecursionlimit(100000)
  table_microbes_genus_only = table_microbes_genus.drop("alphadiversity",axis=1)
  include_taxa = table_microbes_genus_only.sum(axis=0).sort_values(ascending=False).head(10).index
  fig = plt.figure(figsize=(7.5,7.5))
  spec = gridspec.GridSpec(ncols=3, nrows=5, figure=fig, height_ratios=[0.45,0.45,0.45,0.45,0.45], width_ratios=[1,1,0.05], left=0.075,top=0.9,bottom=0.02,right=0.8,wspace=0.2,hspace=0.2)
  axes = np.array([[fig.add_subplot(spec[i,j]) for i in range(5)] for j in range(3)])
  for ax in axes.reshape(-1):
    ax.tick_params(
    axis='both',which='both',bottom=False,top=False,labelbottom=False,left=False,labelleft=False)
  # Alpha Diversity
  sc = axes[0,0].scatter(pca_microbes_proj[:,0], pca_microbes_proj[:,1], c=diversity, s=1)
  # axes[0,0].set_title("Alpha Diversity", font="serif")
  fig.suptitle("Sample Genus and Pathway Abundances (n=3,124)",fontsize=8)
  axes[0,0].set_title("PC 1 (%.2f%%)"%(pca_microbes.explained_variance_ratio_[0]*100.0), fontsize=8)
  axes[0,0].set_ylabel("PC 2 (%.2f%%)"%(pca_microbes.explained_variance_ratio_[1]*100.0), fontsize=8)
  axes[1,0].scatter(pca_pathways_proj[:,0], pca_pathways_proj[:,1], c=pathways_diversity, s=1)
  axes[1,0].set_title("PC 1 (%.2f%%)"%(pca_pathways.explained_variance_ratio_[0]*100.0), fontsize=8)
  axes[1,0].set_ylabel("PC 2 (%.2f%%)"%(pca_pathways.explained_variance_ratio_[1]*100.0), fontsize=8)
  axes[0,0].text(0.5, 1.25, "Genus Abundance", transform = axes[0,0].transAxes, horizontalalignment='center', verticalalignment='center', rotation=0, fontweight="bold", fontsize=8)
  axes[1,0].text(0.5, 1.25, "Pathway Abundance", transform = axes[1,0].transAxes, horizontalalignment='center', verticalalignment='center', rotation=0, fontweight="bold", fontsize=8)
  cbar = fig.colorbar(sc, cax=axes[2,0], orientation='vertical', fraction=0.02, aspect=40, label="Alpha Diversity")
  cbar.ax.tick_params(labelsize=6) 
  cbar.set_label(label="Alpha Diversity", fontsize=8)
  # Age
  sc = axes[0,1].scatter(pca_microbes_proj[:,0], pca_microbes_proj[:,1], c=demographics.age, s=1, cmap="gnuplot")
  # axes[0,1].set_title("Age")
  axes[1,1].scatter(pca_pathways_proj[:,0], pca_pathways_proj[:,1], c=demographics.loc[common_indices].age, s=1, cmap="gnuplot")
  cbar = fig.colorbar(sc, cax=axes[2,1], orientation='vertical', fraction=0.02, aspect=40, label="Age")
  cbar.ax.tick_params(labelsize=6) 
  cbar.set_label(label="Age", fontsize=8)
  # Cohort
  cohort_labeled = demographics["cohorts"].copy()
  top_cohorts = demographics["cohorts"].value_counts().head(5).index
  top_cohorts = ["HeartTransplant","LiverDisease","LiverTransplant","MICU"]
  cohort_labeled[~cohort_labeled.isin(top_cohorts)] = "Other"
  labels, colormap, colors = plots.make_labels(cohort_labeled)
  colormap["HeartTransplant"] = "orange"
  colormap["LiverDisease"] = "dodgerblue"
  colormap["LiverDisease LiverTransplant"] = "teal"
  colormap["LiverTransplant"] = "limegreen"
  colormap["MICU"] = "purple"
  colormap["Other"] = "#ffffff00"
  colors = np.array([colormap[x] for x in cohort_labeled])
  # axes[0,2].set_title("Disease Cohort")
  sc = axes[0,2].scatter(pca_microbes_proj[:,0], pca_microbes_proj[:,1], c=colors, s=1)
  sc = axes[1,2].scatter(pca_pathways_proj[:,0], pca_pathways_proj[:,1], c=np.array(colors)[demographics.index.isin(common_indices)], s=1)
  patches = [mpatches.Patch(color=colormap[label], linewidth=2, edgecolor="black", label=label) for label in top_cohorts]
  axes[2,2].legend(handles=patches,title_fontsize=8,fontsize=8,loc=(0,0.2),title="Cohort")
  axes[2,2].axis("off")
  # Days Hospitalized
  demographics.days_since_admitted[demographics.days_since_admitted>500] = np.nan
  cmap = matplotlib.cm.get_cmap('Spectral', 256)
  newcmp = ListedColormap(cmap(np.power(np.linspace(0, 1, 256), 0.6)))
  sc = axes[0,3].scatter(pca_microbes_proj[:,0], pca_microbes_proj[:,1], c=(demographics.days_since_admitted), s=1, cmap=newcmp)
  axes[1,3].scatter(pca_pathways_proj[:,0], pca_pathways_proj[:,1], c=(demographics.loc[common_indices].days_since_admitted), s=1, cmap=newcmp)
  cbar = fig.colorbar(sc, cax=axes[2,3], orientation='vertical', fraction=0.02, aspect=40, label="Days After Admission")
  cbar.ax.tick_params(labelsize=6) 
  cbar.set_label(label="Days After Admission", fontsize=8)
  sc = axes[0,4].scatter(pca_microbes_proj[:,0][~demographics.year.isna()], pca_microbes_proj[:,1][~demographics.year.isna()], c=(demographics[~demographics.year.isna()].year.values), s=1, cmap=cmap)
  axes[1,4].scatter(pca_pathways_proj[:,0][~demographics.loc[common_indices].year.isna()], pca_pathways_proj[:,1][~demographics.loc[common_indices].year.isna()], c=(demographics.loc[common_indices][~demographics.loc[common_indices].year.isna()].year.values.astype(int)), s=1, cmap=cmap)
  cbar = fig.colorbar(sc, cax=axes[2,4], orientation='vertical', fraction=0.02, aspect=40, label="Admission Year")
  cbar.ax.tick_params(labelsize=6) 
  cbar.ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([2020,2021,2022,2023]))
  cbar.set_label(label="Admission Year", fontsize=8)
  fig.savefig("out/figure_supp_pcaplots.pdf")
  fig.savefig("out/figure_supp_pcaplots.png",dpi=400)

def make_figure_1b():
  table_microbes_genus = table_microbes.loc[:,table_microbes_clr.columns.str.endswith("genus.txt")].copy().T
  shotgun_to_sample = table_microbes.columns.to_series().apply(lambda filename: samples["shotgunSeq_id"].apply(lambda id: str(id) in filename))
  shotgun_to_sample = shotgun_to_sample.apply(lambda row: samples.loc[np.where(row.values)].assign(filename=row.name), axis=1)
  shotgun_to_sample = pd.concat(shotgun_to_sample.values).set_index("filename")
  shotgun_to_sample = shotgun_to_sample[~shotgun_to_sample["mrn"].isna()]
  shotgun_to_sample["mrn"] = shotgun_to_sample["mrn"].astype("int")
  shotgun_to_sample.index = shotgun_to_sample.index.str.split("_R1").str[0]
  # include_taxa = table_microbes.sum(axis=0).sort_values(ascending=False).head(8).index
  # table_microbes_barplot = table_microbes.loc[:,include_taxa]
  # table_microbes_barplot["Other"] = table_microbes.loc[:,~table_microbes.columns.isin(include_taxa)].sum(axis=1).values
  # table_microbes_barplot = table_microbes_barplot.div(table_microbes_barplot.sum(axis=1),axis=0)
  diversity = table_microbes_genus.apply(lambda row: scipy.stats.entropy(row), axis=1)
  # table_microbes_genus.index.map(lambda filename: samples[])
  demo_samples = table_microbes_genus.index.to_series().apply(lambda filename: samples[samples["shotgunSeq_id"].fillna("NAN").apply(lambda id: id in filename)])
  # demo_samples = table_microbes_genus.index.to_series().apply(lambda filename: samples[samples["shotgunSeq_id"].str.contains(filename, na=False)])
  demo_mrn = demo_samples.apply(lambda df: ((df.iloc[0].mrn) if (len(df)>0) else 0))
  demo_mrn[demo_mrn.isna()] = 0
  demo_mrn = demo_mrn.astype("int")
  demo_sex = demo_mrn.apply(lambda mrn: table_demographics.loc[mrn]["sex"] if mrn in table_demographics.index else None)
  demo_birthdate = demo_mrn.apply(lambda mrn: table_demographics.loc[mrn]["birth_date"].split(" ")[0] if mrn in table_demographics.index else None)
  demo_birthyear = demo_birthdate.str.split("-").str[0]
  demographics = pd.concat([demo_mrn.rename("mrn"), demo_sex.rename("sex"), demo_samples.rename("samples"), demo_birthyear.rename("birthyear")], axis=1)
  age2 = demographics.apply(lambda row: (str(row.samples.iloc[0].date_collected).split("-")[0]) if len(row.samples)>0 else None, axis=1).apply(lambda x: np.nan if x=="None" else x).astype(float)
  age1 = demographics.apply(lambda row: float(np.nan if row.birthyear is None else row.birthyear) if len(row.samples)>0 else None, axis=1)
  demographics["age"] = age2-age1
  # demographics["age"] = demographics.apply(lambda row: float(str(row.samples.iloc[0].date_collected).split("-")[0]) - float(np.nan if row.birthyear is None else row.birthyear) if len(row.samples)>0 else None, axis=1)
  demographics["cohorts"] = demographics["samples"].apply(lambda df: " ".join([str(x) for x in np.sort(np.unique(df["db"]))]))
  # days hospitalized
  collection_to_visit = hospital_visits.copy()
  collection_to_visit = collection_to_visit[["mrn","admitted","discharged","collection_date_int","shotgun_collections"]]
  collection_to_visit["keep_indices"] = collection_to_visit.apply(lambda row: ~((pd.Series(row["shotgun_collections"]).duplicated(keep="first")) | pd.Series(row["shotgun_collections"]).isna()).values, axis=1)
  collection_to_visit["collection_date_int"] = collection_to_visit.apply(lambda row: np.array(row["collection_date_int"])[row["keep_indices"]], axis=1)
  collection_to_visit["shotgun_collections"] = collection_to_visit.apply(lambda row: np.array(row["shotgun_collections"])[row["keep_indices"]], axis=1)
  collection_to_visit = collection_to_visit[collection_to_visit["shotgun_collections"].str.len()>0]
  collection_to_visit["shotgun_collections_dates"] = collection_to_visit.apply(lambda row: list(zip(row["collection_date_int"].astype(int), row["shotgun_collections"])), axis=1)
  collection_to_visit = collection_to_visit.explode("shotgun_collections_dates")
  collection_to_visit["collection_date_int"] = collection_to_visit["shotgun_collections_dates"].str[0]
  collection_to_visit["shotgun_collections"] = collection_to_visit["shotgun_collections_dates"].str[1]
  collection_to_visit["days_since_admitted"] = collection_to_visit["collection_date_int"] - collection_to_visit["admitted"]
  collection_to_visit = collection_to_visit.drop_duplicates("shotgun_collections").set_index("shotgun_collections",drop=False)
  demographics["days_since_admitted"] = demographics.samples.apply(lambda row: collection_to_visit.loc[collection_to_visit.index.intersection(row["shotgunSeq_id"].values)].days_since_admitted.values)
  demographics["days_since_admitted"] = demographics["days_since_admitted"].apply(lambda row: row[0] if len(row)>0 else np.nan)
  sys.setrecursionlimit(100000)
  table_microbes_genus_only = table_microbes_genus.drop("alphadiversity",axis=1)
  include_taxa = table_microbes_genus_only.sum(axis=0).sort_values(ascending=False).head(10).index
  # alphadiversity = table_microbes_genus_only["alphadiversity"]
  table_microbes_barplot = table_microbes_genus_only.loc[:,include_taxa]
  table_microbes_barplot["Other"] = table_microbes_genus_only.loc[:,~table_microbes_genus_only.columns.isin(include_taxa)].sum(axis=1).values
  table_microbes_barplot = table_microbes_barplot.div(table_microbes_barplot.sum(axis=1),axis=0)
  if False:
    ordering = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(table_microbes_barplot, method="complete", metric="cityblock"), no_plot=True, color_threshold=-np.inf)['leaves']
  elif False:
    top_microbe = table_microbes_barplot.sum(axis=0).drop("Other").sort_values(ascending=False).index[0]
    ordering = table_microbes_barplot[top_microbe].argsort()
  else:
    # top_microbe = "alphadiversity"
    # ordering = table_microbes_barplot[top_microbe].argsort()
    ordering = diversity.argsort()
  table_microbes_barplot = table_microbes_barplot.iloc[ordering,:]
  diversity_ord = diversity[ordering]
  demo_sex_ord = demo_sex[ordering]
  demographics_ord = demographics.iloc[ordering]
  yresolution = 1000
  image = np.zeros((yresolution, table_microbes_barplot.shape[0], 3))
  labels, colormap, colors = plots.make_labels(table_microbes_barplot.columns)
  colors = [matplotlib.colors.to_rgb(color) for color in colors]
  colormap = {key:matplotlib.colors.to_rgb(colormap[key]) for key in colormap}
  bottom = np.zeros(table_microbes_barplot.shape[0])
  for taxon, abundances in table_microbes_barplot.items():
    color = colormap[taxon]
    for si,sample in enumerate(abundances):
      image[int(np.floor(bottom[si]*yresolution*0.9999)):int(np.floor((bottom[si]+sample)*yresolution*0.9999)),si] = color
    bottom += abundances.values
  patches = [mpatches.Patch(color=colormap[label], label=label) for label in labels]
  patches_2 = [mpatches.Patch(color=(0.8,0.2,0.2), label="Male"), mpatches.Patch(color=(0.2,0.2,0.8), label="Female")]
  fig = plt.figure(constrained_layout=True, figsize=(7.5,2.9))
  axh1 = 0.05
  axhpad = 0.02
  axwpad = 0.10
  axright = 0.02
  axes = [None, None, None, None]
  axes[3] = fig.add_axes([axwpad, 1.0-(4*axhpad)-(4*axh1)-0.36, (1.0 - (axwpad + axright)), 0.36])
  axes[2] = fig.add_axes([axwpad, 1.0-(3*axhpad)-(3*axh1), (1.0 - (axwpad + axright)), axh1])
  axes[1] = fig.add_axes([axwpad, 1.0-(2*axhpad)-(2*axh1), (1.0 - (axwpad + axright)), axh1])
  axes[0] = fig.add_axes([axwpad, 1.0-axhpad-axh1, (1.0 - (axwpad + axright)), axh1])
  # axend = 1.0-(4*axhpad)-(4*axh1)-0.40
  # spec = gridspec.GridSpec(ncols=1, nrows=4, figure=fig, height_ratios=[0.2,0.2,0.2,2], hspace=0.001)
  # axes = [None, None, None, None]
  # axes[0] = fig.add_subplot(spec[0])
  # axes[1] = fig.add_subplot(spec[1])
  # axes[2] = fig.add_subplot(spec[2])
  # axes[3] = fig.add_subplot(spec[3], sharex=axes[1])
  width_bars = (1.0 - (axwpad + axright)) * 7.5
  height_bars = axh1 * 2.9
  height_graph = 0.36 * 2.9
  aspect_bars = len(diversity_ord) * (height_bars / width_bars)
  aspect_graph = aspect_bars / yresolution * (height_graph/height_bars)     # 0.8
  axes[0].set_title("Sample Microbiome Composition",fontsize=8)
  axes[0].imshow(np.array(diversity_ord.values).reshape(1,-1), aspect=aspect_bars, interpolation="none")
  axes[0].set_xticklabels([])
  axes[0].set_yticks([])
  axes[0].set_yticklabels([])
  axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  axes[0].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
  axes[0].set_ylabel("α-diversity", fontsize=8, rotation=0, horizontalalignment='right', verticalalignment='center')
  image_sex = np.array(demo_sex_ord.values)
  image_sex = np.zeros((1, demo_sex_ord.shape[0], 3))
  image_sex[0,demo_sex_ord.values=="Male"]=np.array([0.8,0.2,0.2])
  image_sex[0,demo_sex_ord.values=="Female"]=np.array([0.2,0.2,0.8])
  image_sex[0,demo_sex_ord.values==None]=np.array([1,1,1])
  axes[1].imshow(image_sex, aspect=aspect_bars, interpolation="none")
  axes[1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  axes[1].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
  axes[1].set_ylabel("sex", fontsize=8, rotation=0, horizontalalignment='right', verticalalignment='center')
  min_age = demographics_ord["age"].min()
  max_age = demographics_ord["age"].max()
  cmap = matplotlib.cm.get_cmap("cool")
  cmap.set_bad("white")
  age_normalized = (demographics_ord["age"] - min_age) / max_age
  image_age = np.array([cmap(x) for x in age_normalized]).reshape(1,-1,4)
  # image_age = np.array((demographics_ord["age"]).values).reshape(1,-1)
  axes[2].imshow(image_age, cmap="cool", aspect=aspect_bars, interpolation="none")
  axes[2].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  axes[2].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
  axes[2].set_ylabel("age", fontsize=8, rotation=0, horizontalalignment='right', verticalalignment='center')
  axes[3].imshow(image, aspect=aspect_graph, interpolation_stage="rgba", interpolation="none")
  # axes[3].set_yticklabels([])
  axes[3].set_ylabel("Genus Relative\nAbundance", fontsize=8)
  axes[3].set_yticks([0,200,400,600,800,1000])
  axes[3].set_xticklabels(axes[3].get_xticklabels(), fontsize=6)
  axes[3].set_yticklabels(["0.0","0.2","0.4","0.6","0.8","1.0"][::-1], fontsize=6)
  axes[3].set_xlabel("Sample", fontsize=8)
  fig.legend(handles=patches, loc='upper left', bbox_to_anchor=[axwpad, 0, 0.40, axes[3].get_position().y0-0.14], fontsize=6.5, ncol=3)
  fig.legend(handles=patches_2, loc='upper left', bbox_to_anchor=(axwpad + 0.55, 0, 0.10, axes[3].get_position().y0-0.14), fontsize=6.5, ncol=1, title="sex", title_fontsize=8)
  cax = fig.add_axes([axwpad + 0.67, 0.08, 0.015, 0.14])
  cbar_age = fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=demographics_ord["age"].min(), vmax=demographics_ord["age"].max(), clip=False), cmap=matplotlib.cm.get_cmap("cool").copy()), cax = cax, orientation='vertical', label="Age")
  cbar_age.ax.tick_params(labelsize=6) 
  cbar_age.set_label('age',size=8)
  cax = fig.add_axes([axwpad + 0.75, 0.08, 0.015, 0.14])
  cbar_diversity = fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=diversity_ord.min(), vmax=diversity_ord.max(), clip=False), cmap=matplotlib.cm.get_cmap("viridis").copy()), cax = cax, orientation='vertical', label="α-diversity")
  cbar_diversity.ax.tick_params(labelsize=6) 
  cbar_diversity.set_label('α-diversity',size=8)
  fig.savefig("out/figure_1b.pdf")
  # plt.show()
  plt.close()

def make_figure_1c():
  # fig,axes = plt.subplots(1,1, figsize=(2.35,3.25), tight_layout=True)
  # axes = np.array([axes]).reshape(-1)
  fig = plt.figure(figsize=(2.35,3.10))
  spec = gridspec.GridSpec(figure=fig, ncols=1, nrows=1, left=0.15, bottom=0.2, right=0.95, top=0.92)
  axes = [fig.add_subplot(spec[0])]
  ypos = pd.DataFrame(np.arange(study_intervals["mrn"].nunique()), index=study_intervals.sort_values("date_collection1").drop_duplicates("mrn")["mrn"].values, columns=["ypos"])
  points_x = []
  points_y = []
  lines_x = []
  lines_y = []
  for i,(ri,row) in enumerate(study_intervals.sort_values("date_collection1").iterrows()):
    if row['discharged']-row['admitted'] < 200:
      lines_x += [[row['admitted'],row['discharged']]]
      lines_y += [[ypos.loc[row["mrn"]]["ypos"],ypos.loc[row["mrn"]]["ypos"]]]
      points_x += [row['date_collection1'], row['date_collection2']]
      points_y += [ypos.loc[row["mrn"]]["ypos"],ypos.loc[row["mrn"]]["ypos"]]
  axes[0].scatter(points_x, points_y, c='black', s=1, zorder=2, edgecolors='none')
  axes[0].set_xlabel('Year', fontsize=8)
  axes[0].set_ylabel('Patient', fontsize=8)
  # axes[0].set_xticks(np.linspace(encs[encs['adm_date_int']!=encs['disc_date_int']]['adm_date_int'].min(), encs[encs['adm_date_int']!=encs['disc_date_int']]['disc_date_int'].max(), 5).astype('int'))
  ax_date1 = datetime.datetime((min_date + datetime.timedelta(days=study_intervals["date_collection1"].min()-1)).year, 1, 1)    ## ERROR HERE
  ax_date2 = datetime.datetime((min_date + datetime.timedelta(days=int(study_intervals['discharged'].max())-1)).year+1, 1, 1)
  tick_dates = [ax_date1.replace(year=ax_date1.year+yi) for yi in range((ax_date2.year+1) - ax_date1.year)]
  axes[0].set_xlim((tick_dates[0] - min_date).days, (tick_dates[1] - min_date).days)
  axes[0].set_xticks([(date-min_date).days for date in tick_dates])
  axes[0].set_xticklabels([date.strftime('%Y') for date in tick_dates], fontsize=6)
  axes[0].xaxis.grid(linestyle=(0, (1, 2)))
  # axes[0].set_xticks(np.linspace(study_intervals['admitted'].min(), study_intervals['discharged'].max(), 5).astype('int'))
  # axes[0].set_xticklabels([(datetime.timedelta(days=int(x))+min_date).strftime('%b-%d-%Y') for x in axes[0].get_xticks()], fontsize=8)
  fig.suptitle('Inpatient Visits With 2+ Stool Collections', fontsize=8)
  axes[0].plot(np.array(lines_x).T, np.array(lines_y).T, c='gray', lw=0.5, zorder=1)
  axes[0].set_yticks([])
  # plt.show()
  # axes[0].set_yticklabels(axes[0].get_yticklabels(), font='serif', fontsize=12)
  fig.savefig("out/figure_1c.pdf")

def make_figure_1d():
  top_medications = medication_counts.loc[medications]["count"].sort_values(ascending=False).head(20)
  fig = plt.figure(figsize=(2.5,3.10))
  spec = gridspec.GridSpec(figure=fig, ncols=1, nrows=1, left=0.47, bottom=0.2, right=0.90, top=0.92)
  ax = fig.add_subplot(spec[0])
  ax.barh(np.arange(len(top_medications))[::-1], top_medications.values, 0.8, color="gray")
  ax.set_yticks(np.arange(len(top_medications))[::-1])
  ax.set_yticklabels(top_medications.index.map(lambda name: abbr(";".join([prettify_medication_name(subname.lower(), length=50) for subname in name.split("|")]), 60)), fontsize=6)
  # ax.set_xticklabels(ax.get_xticklabels(), fontsize=6)
  ax.set_xticks([0,600,1200,1800])
  ax.set_xticklabels(ax.get_xticks(), fontsize=6,rotation=90)
  fig.suptitle("Most Commonly Given Medications", fontsize=8)
  ax.set_ylabel("Medication", fontsize=8)
  ax.set_xlabel("Count (Patient-Days)", fontsize=8, loc="center")
  ax.xaxis.set_label_coords(0, -0.15)
  fig.savefig("out/figure_1d.pdf")

def make_figure_1e():
  variance_explained = pd.DataFrame(0, index=["r2"], columns=[])
  keys = results.Xs.keys()
  for key in keys:
    print(key)
    Y = results.Ys[tuple(key)]
    X = results.Xs[tuple(key)]
    X["const"] = 1.0
    coefs = results.all_results_startmed.loc[:,results.all_results_startmed.columns.str.endswith("%s_%s_%s"%(key[0],key[1],key[2]))]
    covariates = [x for x in X.columns if x in coefs.index]
    # coefs_aligned = coefs.loc[covariates,"coef_%s_%d_%s_%s"%(Y.columns[1], key[0], key[1], key[2])].fillna(0)
    X_aligned = X.loc[:,covariates].fillna(0)
    # X_aligned @ coefs_aligned
    # coefs_aligned = coefs.loc[covariates].fillna(0)
    coefs_aligned = coefs.loc[covariates,coefs.columns.intersection(["coef_%s_%d_%s_%s"%(col, key[0], key[1], key[2]) for col in Y.columns[1:]])].fillna(0)
    Y_covariates = Y.columns.intersection(coefs_aligned.columns.str.split("_").str[1])
    coefs_aligned = coefs_aligned.loc[:,"coef_"+Y_covariates+"_%s_%s_%s"%(key[0],key[1],key[2])]
    # coefs_aligned = coefs.loc[covariates,:].fillna(0)
    key_str = "_".join([str(k) for k in key])
    Y_subset = Y[Y_covariates]
    variance_explained[key_str] = 1.0 - (((X_aligned @ coefs_aligned) - Y_subset.values).var().sum() / Y_subset.var().sum())
    print("Total Variance Explained: ", variance_explained)
  variance_explained = variance_explained.T
  variance_explained["duration"] = variance_explained.index.str.split("_").str[0]
  variance_explained["Phenotype"] = variance_explained.index.str.split("_").str[1].str.title().str.replace("Metab","Metabolite")
  variance_explained["modeltype"] = variance_explained.index.str.split("_").str[2]
  variance_explained["model_descr"] = (variance_explained["duration"]+"_"+variance_explained["modeltype"]).apply(lambda x: {"0_binary":"2-10\nDays", "1_binary":"10-20\nDays", "2_binary":"20-30\nDays", "0_dose":"Duration", "0_null":"Null"}[x])
  variance_explained = variance_explained.sort_values(["model_descr","Phenotype"], key=lambda x: x.map({"Null":-1, "2-10\nDays":0, "10-20\nDays":1, "20-30\nDays":2, "Duration":3, "Species":0, "Genus":1, "Pathway":2,"Metabolite":3}))
  variance_explained = variance_explained[variance_explained.model_descr!="Null"]
  # columns = list(variance_explained.columns[variance_explained.columns.str.contains("0_.*_binary",regex=True)]) + list(variance_explained.columns[variance_explained.columns.str.contains("1_.*_binary",regex=True)]) + list(variance_explained.columns[variance_explained.columns.str.contains("2_.*_binary",regex=True)]) + list(variance_explained.columns[variance_explained.columns.str.contains("dose",regex=True)])
  # variance_explained = variance_explained[columns]
  for include_genus in [True,False]:
    fig = plt.figure(figsize=(2.35, 3.10))
    spec = gridspec.GridSpec(figure=fig, ncols=1, nrows=1, left=0.2, bottom=0.2, right=0.98, top=0.92)
    axes = fig.add_subplot(spec[0])
    # fig,axes = plt.subplots(1,1,constrained_layout=True, figsize=(2.4,2.3))
    sns.barplot(variance_explained[variance_explained["Phenotype"].isin(["Species","Pathway","Metabolite"]+(["Genus"] if include_genus else []))], x="model_descr", y="r2", hue="Phenotype", ax=axes)
    fig.suptitle("Model Variance Explained By Fixed Effects", fontsize=8)
    axes.set_xlabel("Model", fontsize=8)
    axes.set_xticklabels(axes.get_xticklabels(),fontsize=6)
    axes.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    axes.set_yticklabels(axes.get_yticklabels(),fontsize=6)
    axes.set_ylabel("Explained Variance ($R^{2}$)", fontsize=8)
    axes.legend(prop=dict(size=6),ncol=1)
    # axes.get_legend().set(size=6)
    fig.savefig("out/figure_supplementary_1e.pdf" if include_genus else "out/figure_1e.pdf")
  # plt.show()
  plt.close()


# results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["species","genus"]))].groupby("medication_pharm_class")["coef"].apply(lambda xs: np.abs(xs).sum()).sort_values(ascending=False).head(50)
# aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["species","genus"]))][["medication_pharm_class","coef"]].groupby(["medication_pharm_class"]).apply(lambda xs: xs.abs().sum()).sort_values("coef",ascending=False)

def make_figure_2a():
  figure_normalize = "effectsize"
  df_medication_counts = medication_counts["count"].loc[medications].sort_values(ascending=False)
  # medication_counts = meds2.drop_duplicates(["take_med_int","harmonized_generic_route"]).groupby("harmonized_generic_route").size().loc[medications].sort_values(ascending=False)
  # medication_counts = medications_started_dose.sum()
  for key in set(results.table.medication.unique()) - set(df_medication_counts.index):
    keymeds = key.split("|")
    df_medication_counts.loc[key] = df_medication_counts.loc[keymeds].sum()
  # medication_counts = meds2.drop_duplicates(["take_med_int","harmonized_generic_route"]).groupby("harmonized_generic_route").size().loc[medications].sort_values(ascending=False)
  figure_dims = (7.5,3.25)
  fig = plt.figure(figsize=figure_dims)
  # fig, axes = plt.subplots(2,1, tight_layout=True, figsize=figure_dims)
  # spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, left=0.3, right=0.95, top=0.925)
  axes = [fig.add_axes([0.27,0.92 - 0.8,0.18,0.78]), fig.add_axes([0.8,0.92 - (0.8*(2/3)),0.18,0.78*(2/3)])]
  # axes = [fig.add_subplot(spec[0]), fig.add_subplot(spec[1])]
  axes = np.array([axes]).reshape(-1)
  for li,level in enumerate(["medication","medication_pharm_class"]):
    num_display_rows = {"medication":30, "medication_pharm_class":20}[level]  
    if level=="medication":
      if figure_normalize == "count":
        aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary")][["medication","sampletype","duration"]].value_counts()
        # medication_counts = medications_started_dose.sum()
        # medication_counts = meds2.drop_duplicates(["take_med_int","harmonized_generic_route"]).groupby("harmonized_generic_route").size().loc[medications].sort_values(ascending=False)
        aggregated = aggregated.div(df_medication_counts.loc[aggregated.index.get_level_values("medication")].values)
      elif figure_normalize == "effectsize":
        aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary")][["medication","sampletype","duration","coef"]].groupby(["medication","sampletype","duration"]).apply(lambda xs: xs.abs().mean())
      elif figure_normalize == "pvalue":
        aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary")][["medication","sampletype","duration","pvalue"]].groupby(["medication","sampletype","duration"]).apply(lambda xs: -np.log(xs).sum())
        # results.table[(results.table.modeltype=="binary")][["medication","sampletype","duration","coef"]].groupby(["medication","sampletype","duration"]).apply(lambda xs: xs.abs().mean())
        # results.table[(results.table.modeltype=="binary") & (results.table.significant)][["medication","sampletype","duration","coef"]].groupby(["medication","sampletype","duration"]).apply(lambda xs: xs.sum())
      else:
        aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary")][["medication","sampletype","duration"]].value_counts()
      figure_number_of_associations = aggregated.reset_index().pivot(index="medication",columns=["sampletype","duration"]).fillna(0)
    elif level=="medication_pharm_class":
      if figure_normalize == "count":
        aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary")].explode("medication_pharm_class")[["medication_pharm_class","sampletype","duration"]].value_counts()
        # medication_class_counts = meds2.drop_duplicates(["take_med_int","pharm_class"]).groupby("pharm_class").size().loc[aggregated.index.get_level_values("medication_pharm_class")]
        medication_class_counts = medication_class_counts["count"].loc[aggregated.index.get_level_values("medication_pharm_class")]
        aggregated = aggregated.div(medication_class_counts.values)
      elif figure_normalize == "effectsize":
        aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary")].explode("medication_pharm_class")[["medication_pharm_class","sampletype","duration","coef"]].groupby(["medication_pharm_class","sampletype","duration"]).apply(lambda xs: xs.abs().mean())
      elif figure_normalize == "pvalue":
        aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary")].explode("medication_pharm_class")[["medication_pharm_class","sampletype","duration","pvalue"]].groupby(["medication_pharm_class","sampletype","duration"]).apply(lambda xs: -np.log(xs).sum())
      else:
        aggregated = results.table[(results.table.significant) & (results.table.modeltype=="binary")].explode("medication_pharm_class")[["medication_pharm_class","sampletype","duration"]].value_counts()
      figure_number_of_associations = aggregated.reset_index().pivot(index="medication_pharm_class",columns=["sampletype","duration"]).fillna(0)
      figure_number_of_associations.index = figure_number_of_associations.index.str.title()
      num_display_rows = 20
      figure_dims = (6,2.333)
    figure_number_of_associations.columns = figure_number_of_associations.columns.droplevel(0)
    # figure_number_of_associations = figure_number_of_associations.reindex([("genus",0),("genus",1),("genus",2),("genus",3), ("metab",0),("metab",1),("metab",2),("metab",3), ("pathway",0),("pathway",1),("pathway",2),("pathway",3)],axis=1)
    figure_number_of_associations = figure_number_of_associations.reindex([("genus",0),("genus",1),("genus",2), ("metab",0),("metab",1),("metab",2), ("pathway",0),("pathway",1),("pathway",2)],axis=1)
    figure_number_of_associations["total"] = figure_number_of_associations.sum(axis=1)
    figure_number_of_associations = figure_number_of_associations.sort_values("total",ascending=True).fillna(0)
    figure_number_of_associations.to_csv("out/supp_table_top_%s.csv"%(level))
    figure_number_of_associations = figure_number_of_associations.tail(num_display_rows)
    ax = axes[li]
    bottom = figure_number_of_associations.iloc[:,0]*0
    # colors = ["#b6d7a8","#93c47d","#6aa84f","green", "#ffe599","#ffd966","#f1c232","orange", "#9fc5e8","#6fa8dc","#3d85c6","blue"]
    colors = ["#b6d7a8","#93c47d","#6aa84f", "#ffe599","#ffd966","#f1c232", "#9fc5e8","#6fa8dc","#3d85c6"]
    for ci,column in enumerate(figure_number_of_associations.columns[:-1]):
      if level=="medication":
        _ = ax.barh(figure_number_of_associations.index.map(lambda x: prettify_compound_medication_name(x,30,mode="trunc")), figure_number_of_associations[column], 0.9, left=bottom, color=colors[ci], label="%s (%d-%d days)"%({"genus":"genus","metab":"metabolite","pathway":"pathway"}[column[0]], durations[column[1]][0], durations[column[1]][1] ))
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
      else:
        _ = ax.barh(figure_number_of_associations.index, figure_number_of_associations[column], 0.9, left=bottom, color=colors[ci], label="%s (%d-%d days)"%({"genus":"genus","metab":"metabolite","pathway":"pathway"}[column[0]], durations[column[1]][0], durations[column[1]][1] ))
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
      bottom += figure_number_of_associations[column]
    ax.set_ylim(0-0.5-0.5, num_display_rows)
    if level=="medication":
      if figure_normalize == "effectsize":
        ax.set_title("Medications with the Largest Effects", fontsize=8)
      else:
        ax.set_title("Medications with the Most Associations", fontsize=8)
      ax.set_ylabel("Medication", fontsize=8)
    else:
      if figure_normalize == "effectsize":
        ax.set_title("Medication Classes\nwith the Largest Effects", fontsize=8)
      else:
        ax.set_title("Medications with the Most Associations", fontsize=8)
      ax.set_ylabel("Medication Class", fontsize=8)
    if figure_normalize == "count":
      ax.set_xlabel("Proportion of Significant Associations", fontsize=8)
    elif figure_normalize == "effectsize":
      ax.set_xlabel("Mean Log2 Relative Abundance     ", fontsize=8)
    else:
      ax.set_xlabel("Number of Significant Associations", fontsize=8)    
    ax.set_xticklabels(ax.get_xticklabels(),fontsize=6)
  labels = np.concatenate([["%s (%d-%d days)"%(si, durations[di][0], durations[di][1]) for di in range(len(durations))] for si in ["genus","metabolite","pathway"]])
  patches = [mpatches.Patch(color=colors[i], linewidth=2, edgecolor="black", label=labels[i]) for i in range(len(labels))]
  fig.legend(handles=patches, fontsize=6,loc='upper center', bbox_to_anchor=(0.725, 0.25), ncol=3)
  fig.savefig("out/figure_2a.pdf")
  fig.savefig("out/figure_2a.png")
  # plt.show()
  plt.close()

def make_figure_2b():
  fig = plt.figure(figsize=(2,3.3))
  # spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[1,0.05], wspace=0.05, left=0.7, bottom=0.3, right=0.99, top=0.9)
  # axes = np.array([fig.add_subplot(spec[0]), fig.add_subplot(spec[1])])
  # axes = np.array([axes]).reshape(-1)
  axes = [fig.add_axes([0.75,0.4,0.2,0.53]), fig.add_axes([0.2,0.15,0.6,0.025])]
  grid_medications = results.table[(results.table["microbe"]=="alphadiversity") & (results.table["significant"]) & (results.table["modeltype"]=="binary")]["medication"].unique()
  grid = results.table[(results.table["microbe"]=="alphadiversity") & (results.table["significant"]) & (results.table["sampletype"]=="species") & (results.table["modeltype"]=="binary")].pivot(columns="medication", index=["duration"], values="coef").T
  grid_filtered = grid.fillna(0)
  grid_filtered.to_csv("out/table_alphadiversity.csv")
  # ordering_medication = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid_filtered,method="average"), no_plot=True, color_threshold=-np.inf)['leaves']
  ordering_medication = grid_filtered.mean(axis=1).argsort()
  grid = grid.iloc[ordering_medication,:]
  cmap = matplotlib.cm.get_cmap("PiYG").copy()
  cmap.set_bad(color="#d0d0d0")  
  hm = sns.heatmap(grid, ax=axes[0], cmap=cmap, cbar=False, center=0, linewidths=0.1, linecolor='#a0a0a0')
  fig.suptitle("\n".join(wrap("Associations with α-Diversity",45)), fontsize=8)
  axes[0].set_yticks(np.arange(grid.shape[0])+0.5, grid.index)
  axes[0].set_yticklabels(grid.index.str.lower().map(lambda name: ";".join([prettify_medication_name(subname, length=50) for subname in name.split("|")])), fontsize=6)
  axes[0].set_xticks([0.5,1.5,2.5])
  axes[0].set_xticklabels(["2-10 days", "10-20 days", "20-30 days"], rotation=45, fontsize=6, ha="right", rotation_mode="anchor")
  axes[0].set_xlabel("Time Window      ", fontsize=8)
  axes[0].set_ylabel("Medication", fontsize=8)
  cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=grid.min().min(), vmax=grid.max().max(), clip=False), cmap=cmap.copy()), cax=axes[1], orientation="horizontal", label="Change in alpha-diversity\n(Δ Shannon Index)")
  cbar.ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(6))
  cbar.ax.tick_params(labelsize=6) 
  cbar.set_label(label="Change in alpha-diversity\n(Δ Shannon Index)", fontsize=8)
  fig.savefig("out/figure_2b.pdf")
  plt.close()
  # cbar_kws={"label":"Change in alpha-diversity\n(Δ Shannon Index)"}
  # for di,duration in enumerate(["2-10 days", "10-30 days", "30-90 days"]):
  #   axes[0].text((di*2)+0.5, grid.shape[0]+6, duration, horizontalalignment='left', verticalalignment='center', weight="bold", color="#505050", fontsize=8)
  # plt.show()

def make_figure_2c():
  temporal_collections = pd.DataFrame()
  temporal_collections["admitted"] = study_intervals.groupby(["visit"]).apply(lambda df: df["admitted"].iloc[0])
  temporal_collections["date_collections"] = study_intervals.groupby(["visit"]).apply(lambda df: list(df.date_collection1) + list(df.date_collection2))
  temporal_collections["filename_collections"] = study_intervals.groupby("visit").apply(lambda df: list(df.shotgunSeq_genus_id1) + list(df.shotgunSeq_genus_id2))
  temporal_collections["unique_ixs"] = temporal_collections["date_collections"].apply(lambda x: np.unique(x, return_index=True)[1])
  temporal_collections["filename_collections"] = temporal_collections.apply(lambda df: np.array(df["filename_collections"])[df["unique_ixs"]],axis=1)
  temporal_collections["date_collections"] = temporal_collections.apply(lambda df: np.array(df["date_collections"])[df["unique_ixs"]],axis=1)
  temporal_collections["notnull_ixs"] = temporal_collections["filename_collections"].apply(lambda xs: np.where([not (pd.isnull(x) or x=="nan") for x in xs])[0])
  temporal_collections["filename_collections"] = temporal_collections.apply(lambda df: np.array(df["filename_collections"])[df["notnull_ixs"]],axis=1)
  temporal_collections["date_collections"] = temporal_collections.apply(lambda df: np.array(df["date_collections"])[df["notnull_ixs"]],axis=1)
  temporal_collections["filename_collections"].apply(lambda xs: np.where([x for x in xs if (pd.isnull(x) or x=="nan")])[0]).str.len().value_counts()
  temporal_collections = temporal_collections.drop("unique_ixs",axis=1)
  temporal_collections = temporal_collections.drop("notnull_ixs",axis=1)
  colorlist = ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#ffe119', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#d0d0d0', '#000000']
  selected_medications = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype=="genus") & (results.table.microbe=="alphadiversity")].sort_values("coef",key=abs,ascending=False).head(8).medication.unique()
  df_values = temporal_collections["filename_collections"].apply(lambda xs: table_microbes.loc[:,xs]).apply(lambda df: df.loc["alphadiversity"].values)
  df_dates  = temporal_collections["date_collections"]
  df_first_collection_value = df_values.apply(lambda xs: xs[0] if len(xs)>0 else np.nan)
  df_first_collection_value.name = "first_collection_value"
  df_first_collection_date = df_dates.apply(lambda xs: xs[0] if len(xs)>0 else np.nan)
  df_first_collection_date.name = "first_collection_date"
  medications_taken = ((medications_started_dose!=0).groupby(study_intervals["visit"].values)).apply(lambda df: df.any())
  medications_starttime = medications_started_starttime.add(medications_started_starttime.index.get_level_values("date_collection1"),axis=0).groupby(study_intervals["visit"].values).min()
  medications_starttime.columns = "time_"+medications_starttime.columns
  df = pd.concat([df_dates,df_values,df_first_collection_date,df_first_collection_value,medications_taken,medications_starttime],axis=1)
  df = df[~df_first_collection_date.isna()]
  df = df.explode(["date_collections","filename_collections"])
  df["admitted"] = hospital_visits.loc[df.index]["admitted"]
  df["discharged"] = hospital_visits.loc[df.index]["discharged"]
  # fig,ax=plt.subplots(1,1,figsize=(3,3),constrained_layout=True)
  fig=plt.figure(figsize=(2.5,3))
  spec = gridspec.GridSpec(ncols=2, nrows=4, figure=fig, left=0.15, bottom=0.15, right=0.95, top=0.875, hspace=0.75, wspace=0.35)
  axes = [fig.add_subplot(spec[i,j]) for i in range(4) for j in range(2)]
  axes = np.array(axes).reshape(-1)
  # axes = np.array([fig.add_subplot(spec[0]), fig.add_subplot(spec[1])])
  linewidths = [1]
  linestyles = ["-",(0,(1,0.5))]
  markers = ['o','v','^','s']
  for mi,(medication,linestyle) in enumerate(zip(selected_medications,itertools.product(linewidths,linestyles,markers))):
    if False:
      xs = df["date_collections"] - df["first_collection_date"]
      ys = df["filename_collections"] - df["first_collection_value"]
      mask = df[medication]
      dfplot = pd.DataFrame({"value":ys[mask].groupby(xs[mask]).mean(), "count":ys[mask].groupby(xs[mask]).size()})
    else:
      medlist = [x for x in medication.split("|")]
      # xs = df["date_collections"] - df["time_"+medication]
      xs = df["date_collections"] - df.loc[:,["time_"+x for x in medication.split("|")]].min(axis=1)
      ys = df["filename_collections"] - df["first_collection_value"]
      mask = df[medlist].any(axis=1)
      dfplot = pd.DataFrame({"value":ys[mask].groupby(xs[mask]).mean(), "count":ys[mask].groupby(xs[mask]).size()})
    dfplot = dfplot.reindex(np.arange(0,40)).interpolate(method="linear",limit_area="inside")
    # dfplot["rolling"] = (dfplot["value"]-dfplot["value"].iloc[0]).rolling(5,min_periods=1).mean()
    dfplot["rolling"] = (dfplot["value"]).rolling(5,min_periods=1,center=False).mean()
    dfplot["rolling"] = dfplot["rolling"]- dfplot["rolling"].iloc[0]
    # dfplot = dfplot.reindex(np.arange(0,35,0.1)).interpolate(method="linear",limit_area="inside")
    # points = dfplot["rolling"].reset_index().values.reshape(-1,1,2)
    # segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # lc = matplotlib.collections.LineCollection(segments, linewidths=dfplot["count"].fillna(0),color=colorlist[mi])
    # dfcount = dfcount.reindex(np.arange(0,35)).interpolate(method="linear",limit_area="inside").rolling(5,min_periods=1).mean()
    # ax.add_collection(lc)
    n_individuals = (medications_started_dose[medication]!=0).reset_index().groupby("mrn")[medication].apply(any).sum()
    n_intervals = (medications_started_dose[medication]!=0).sum()
    axes[mi].set_title(prettify_compound_medication_name(medication),fontsize=6)
    axes[mi].plot(dfplot.index,dfplot["rolling"],zorder=1,label=prettify_compound_medication_name(medication),lw=1,c="black")
    # ax.scatter(dfplot.index,dfplot["rolling"],zorder=1,s=dfplot["count"])
  # xs = df["date_collections"] - df["first_collection_date"]
  xs = df["date_collections"] - df["admitted"]
  ys = df["filename_collections"] - df["first_collection_value"]
  dfcount = ys[mask].groupby(xs[mask]).size()
  dfplot = ys.groupby(xs).mean()
  dfplot = dfplot.reindex(np.arange(0,40)).interpolate(method="linear",limit_area="inside")
  for ai,ax in enumerate(axes):
    ax.plot(dfplot.index,dfplot.rolling(5,min_periods=1).mean(),color="#c0c0c0",lw=2,zorder=0)
    ax.set_yticks([-2,0,2])
    ax.set_xticks([0,10,20,30])
    if ai==len(axes)-1 or ai==len(axes)-2:
      ax.set_xticklabels(["0","10","20","30"],fontsize=6)
      # ax.set_xlabel("# Days Since Starting Medication",fontsize=6)
    else:
      ax.set_xticklabels([])
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=6)
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=6)
    ax.set_xlim(0,35)
    ax.set_ylim(-2.2,2.2)
  fig.text(0.02,0.5,"Change in Alpha Diversity", fontsize=8, rotation=90,ha="left",va="center")
  fig.text(0.5,0.02,"Days Since Starting Medication", fontsize=8, ha="center",va="bottom")
  fig.suptitle("Patient Alpha Diversity Over Time",fontsize=8)
  # ax.set_title("Patient Alpha Diversity Over Time",fontsize=8)
  # ax.set_ylabel("Change in Alpha Diversity",fontsize=8)
  # ax.set_xlabel("Days Since Starting Medication",fontsize=8)
  # ax.legend(fontsize=7,ncol=2,loc='upper center', bbox_to_anchor=(0.42, -0.3),handlelength=1,columnspacing=0.6)
  fig.savefig("out/figure_2c.pdf")
  plt.close()

def make_figure_2d():
  grid_duration = 0
  grid = (results.table[(results.table["duration"]==grid_duration) & (results.table["sampletype"]=="species") & (results.table["modeltype"]=="binary")].set_index(["medication","microbe","duration","sampletype"])["significant"] * results.table[(results.table["duration"]==grid_duration) & (results.table["sampletype"]=="species") & (results.table["modeltype"]=="binary")].set_index(["medication","microbe","duration","sampletype"])["coef"]).unstack(level=0).T
  grid_all_durations = (results.table[(results.table["sampletype"]=="species") & (results.table["modeltype"]=="binary")].set_index(["medication","microbe","duration","sampletype"])["significant"] * results.table[(results.table["duration"]==grid_duration) & (results.table["sampletype"]=="species") & (results.table["modeltype"]=="binary")].set_index(["medication","microbe","duration","sampletype"])["coef"]).unstack(level=0).T
  grid_all_durations.to_csv("out/table_results_species.csv")
  grid_significant = (results.table[(results.table["duration"]==grid_duration) & (results.table["sampletype"]=="species") & (results.table["modeltype"]=="binary")].set_index(["medication","microbe","duration","sampletype"])["significant"]).unstack(level=0).T
  selected_medications = results.table[(results.table.significant) & (results.table.sampletype=="species") & (results.table["modeltype"]=="binary")]["medication"].value_counts().head(50).index
  # selected_microbes = results.table[(results.table.significant) & (results.table.sampletype=="species")]["microbe"].value_counts().head(30).index
  # selected_microbes = table_microbes_relabund.index[table_microbes_relabund.mean(axis=1)>0.001].intersection(results.table)
  selected_microbes = table_microbes_relabund.loc[results.table[(results.table.sampletype=="species") & (results.table.microbe!="alphadiversity")].microbe.unique()]
  selected_microbes = selected_microbes.loc[~selected_microbes.index.str.contains(" sp. ")].mean(axis=1).sort_values(ascending=False).head(25).index
  # table_microbes_relabund.loc[results.table[results.table.sampletype=="species"].microbe.unique()]
  # grid = grid.loc[selected_medications,selected_microbes]
  grid = grid.loc[:,selected_microbes]
  grid = grid.loc[(grid!=0).sum(axis=1).sort_values(ascending=False).head(35).index,:]
  grid = grid.loc[(grid!=0).sum(axis=1)>0, (grid!=0).sum(axis=0)>1]
  grid.columns = grid.columns.get_level_values(0)
  # grid = grid.loc[((grid!=0).sum(axis=1)>5),((grid!=0).sum(axis=0)>7)]
  grid_classes = medication_classes.loc[grid.index.str.split("|").str[0]].reset_index()
  # grid_classes_ix_medname = grid_classes.set_index("index")
  # ordering_X = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid), no_plot=True, color_threshold=-np.inf)['leaves']
  ordering_meds = grid_classes.sort_values("med_pharm_class").index
  grid_classes_grouped = grid.fillna(0).groupby(grid_classes.med_pharm_class.values, axis=0)
  grid_classes_means = grid_classes_grouped.mean()
  ordering_med_classes = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid_classes_means), no_plot=True, color_threshold=-np.inf)['leaves']
  ordering_med_names = np.concatenate([grid_classes_grouped.groups[grid_classes_means.index[ci]].values for ci in ordering_med_classes])
  microbe_columns = (grid!=0).sum().sort_values(ascending=False).index
  ordering_microbe = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid[microbe_columns].fillna(0).T), no_plot=True, color_threshold=-np.inf)['leaves']
  grid = grid.loc[ordering_med_names,:]
  grid = grid[microbe_columns[ordering_microbe]]
  grid = grid.loc[:,grid.columns.sort_values()]
  grid.index = grid.index.map(lambda name: abbr(";".join([prettify_medication_name(subname.lower(), length=30) for subname in name.split("|")]), 60))
  grid = grid.T
  fig = plt.figure(figsize=(5.5,3.3))
  spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, left=0.275, bottom=0.41, right=(5.35 / 5.5), top=0.93)
  axes = [fig.add_subplot(spec[0])]
  axes = np.array([axes]).reshape(-1)
  class_labels = [(grid_classes_means.index[class_ix],len(grid_classes_grouped.groups[grid_classes_means.index[class_ix]].values)) for class_ix in ordering_med_classes]
  # grid.index = grid.index.map(lambda name: abbr(";".join([prettify_medication_name(subname.lower(), length=50) for subname in name.split("|")]), 60))
  dot_heatmap(grid, axes[0], axes[0], matplotlib.cm.get_cmap('Blues_r', 128), matplotlib.cm.get_cmap('Oranges', 128), None, legend_values=[-2.4,-1.6,-0.8,0.8,1.6,2.4], legend_labels=["-2.4","-1.6","-0.8","0.8","1.6","2.4"], min_size=0.2, size=0.2, exponent=2, text_fontsize=5, classlabel_lw=0.75, classlabel_line_offset=0.5, classlabel_text_offset=0.8, legend_kws={"title_fontsize":8,"fontsize":6,"title":"Log-fold\nChange\nin Relative\nAbundance","ncol":2,"bbox_to_anchor":(-0.125, -0.1, 0, 0.05), "loc":"upper right","columnspacing":0.9,"handlelength":1})
  axes[0].set_xlabel("Medication",fontsize=8)
  axes[0].set_ylabel("Species",fontsize=8)
  axes[0].set_xticklabels(grid.columns,fontsize=6,rotation=40,ha="right",rotation_mode="anchor")
  axes[0].set_yticklabels(axes[0].get_yticklabels(),fontsize=6,rotation=0,ha="right",rotation_mode="anchor")
  axes[0].set_title("")
  fig.suptitle("Significant 2-10 Day Species-Level Associations",fontsize=8)
  # axes[1].set_title("Log-fold Change in Relative Abundance", fontsize=8)
  # fig.savefig("out/results_heatmap_microbiome_subset.pdf")
  fig.savefig("out/figure_2d.pdf")
  # plt.show()
  plt.close()

def make_figure_2e():
  microbe = "Bacteroides"
  selected_medications = results.table[(results.table.microbe=="Bacteroides") & (results.table.modeltype=="binary")].sort_values(["medication","duration"]).pivot(index="medication",columns="duration",values="significant").fillna(False).any(axis=1)
  selected_medications = selected_medications[selected_medications].index
  results_pvalue = results.table[(results.table.microbe=="Bacteroides") & (results.table.modeltype=="binary") & (results.table.medication.isin(selected_medications))].sort_values(["medication","duration"]).pivot(index="medication",columns="duration",values="pvalue")
  results_sig = results.table[(results.table.microbe=="Bacteroides") & (results.table.modeltype=="binary") & (results.table.medication.isin(selected_medications))].sort_values(["medication","duration"]).pivot(index="medication",columns="duration",values="significant")
  results_subset = results.table[(results.table.microbe=="Bacteroides") & (results.table.modeltype=="binary") & (results.table.medication.isin(selected_medications))].sort_values(["medication","duration"]).pivot(index="medication",columns="duration",values="coef")
  results_duration = results.table[(results.table.microbe=="Bacteroides") & (results.table.modeltype=="dose") & (results.table.medication.isin(selected_medications))].sort_values(["medication","duration"]).pivot(index="medication",columns="medication_dose",values="coef")
  results_duration_sig = results.table[(results.table.microbe=="Bacteroides") & (results.table.modeltype=="dose") & (results.table.medication.isin(selected_medications))].sort_values(["medication","duration"]).pivot(index="medication",columns="medication_dose",values="significant")
  results_duration_pvalue = results.table[(results.table.microbe=="Bacteroides") & (results.table.modeltype=="dose") & (results.table.medication.isin(selected_medications))].sort_values(["medication","duration"]).pivot(index="medication",columns="medication_dose",values="pvalue")
  ordering_meds = results_subset.mean(axis=1).sort_values().index
  results_subset = results_subset.loc[ordering_meds]
  results_duration = results_duration.loc[ordering_meds]
  results_pvalue = results_pvalue.loc[ordering_meds] 
  results_sig = results_sig.loc[ordering_meds] 
  results_subset = results_subset.loc[ordering_meds] 
  results_duration_sig = results_duration_sig.loc[ordering_meds] 
  results_duration_pvalue = results_duration_pvalue.loc[ordering_meds] 
  grid = pd.concat([results_subset, results_duration],axis=1)
  grid.columns = ["2-10 Days", "10-20 Days", "20-30 Days", "intercept","slope"]
  grid.to_csv("out/supp_table_bacteroides.csv")
  fig = plt.figure(figsize=(2.5,3))
  axes = [fig.add_axes([0.62, 0.3, 0.18, 0.575]), fig.add_axes(fig.add_axes([0.84, 0.3, 0.1, 0.575]))]
  ax_colorbar = fig.add_axes([0.1, 0.15, 0.45, 0.025])
  cmap = matplotlib.cm.get_cmap("PiYG").copy()
  cmap.set_bad(color="#c0c0c0")
  vmax = results.table[(results.table.microbe=="Bacteroides") & (results.table.medication.isin(selected_medications))]["coef"].abs().max()
  sns.heatmap(results_subset, ax=axes[0], cmap=cmap,center=0,cbar_ax=ax_colorbar,cbar_kws={"orientation":"horizontal"},vmin=-vmax, vmax=vmax)
  cbar = axes[0].collections[0].colorbar
  cbar.ax.tick_params(labelsize=6)
  cbar.ax.set_xlabel("Effect Size (Log-Fold\nRelative Abundance)",fontsize=8)
  results_duration[True]*=5
  sns.heatmap(results_duration, ax=axes[1], cmap=cmap,center=0,cbar=False)
  axes[0].set_yticks(np.arange(results_subset.shape[0])+0.5,results_subset.index.map(prettify_compound_medication_name),fontsize=6)
  axes[1].set_yticks([])
  axes[0].set_xticks([0.5,1.5,2.5])
  axes[1].set_xticks([0.5,1.5])
  axes[1].set_xticklabels(["intercept","slope"],fontsize=6,rotation=90)
  axes[1].set_ylabel("")
  axes[1].set_xlabel("")
  axes[0].set_xlabel("")
  axes[0].set_xticklabels(["%d-%d Days"%(d[0],d[1]) for d in durations],fontsize=6,rotation=90)
  axes[0].set_ylabel("Medication",fontsize=8)
  for x in range(results_subset.shape[1]):
    for y in range(results_subset.shape[0]):
      if results_sig.fillna(False).iloc[y,x]:
        axes[0].text(x+0.5,y+1-.25,"**",fontsize=6,ha="center",va="center")
      elif results_pvalue.fillna(1).iloc[y,x] < 0.1:
        axes[0].text(x+0.5,y+1-.25,"*",fontsize=6,ha="center",va="center")
  for x in range(results_duration.shape[1]):
    for y in range(results_duration.shape[0]):
      if results_duration_sig.fillna(False).iloc[y,x]:
        axes[1].text(x+0.5,y+1-.25,"**",fontsize=6,ha="center",va="center")
      elif results_duration_pvalue.fillna(1).iloc[y,x] < 0.1:
        axes[1].text(x+0.5,y+1-.25,"*",fontsize=6,ha="center",va="center")
  fig.text(0.7,0.07,"Presence/\nAbsence",fontsize=7,fontweight="bold",ha="center",va="top")
  fig.text(0.89,0.11,"Duration",fontsize=7,fontweight="bold",ha="center",va="top")
  fig.suptitle("\n".join(wrap("Medications Associated with Bacteroides Abundance",35)),fontsize=8)
  fig.savefig("out/figure_2e.pdf")
  plt.close()
  # plt.show()

def make_figure_2f():
  microbe = "Bifidobacterium"
  medication = "LACTULOSE_Oral"
  results_subset = results.table_unfiltered[(results.table_unfiltered.sampletype.isin(["genus","species"])) & (results.table_unfiltered.medication==medication) & (results.table_unfiltered.microbe.str.contains("Bifidobacterium")) & (results.table_unfiltered.modeltype=="binary")].sort_values(["medication","duration"])
  results_subset = results_subset.groupby("microbe").filter(lambda df: df.significant.sum()>=0).sort_values(["microbe","pvalue"])
  mean_coefs = results_subset.groupby("microbe")["coef"].mean()
  results_subset["_mean_coef"] = results_subset["microbe"].map(mean_coefs)
  results_subset = results_subset.sort_values(by=["_mean_coef", "microbe"]).drop(columns="_mean_coef")
  fig = plt.figure(figsize=(2.5,3))
  ax_durations = np.sort(results_subset["duration"].unique())
  ax_heights = np.array([len(results_subset[results_subset.duration==di]) / len(results_subset) for di in ax_durations])
  # ax_heights = np.array([0.25,ax_heights[0],0.01,ax_heights[1],0.01,ax_heights[2],0.01,ax_heights[3],0.15])
  ax_bottom = 0.15
  ax_top = 0.875
  ax_heights = np.concatenate([np.array([[ax_heights[i],0.01] for i in range(len(ax_heights))]).reshape(-1)])[:-1]
  ax_heights /= ax_heights.sum()
  ax_heights *= (ax_top - ax_bottom)
  ax_heights = np.concatenate([[ax_bottom], ax_heights, [1-ax_top]])
  ax_ypos = np.cumsum(ax_heights)
  ax_bottom = ax_ypos
  ax_left = 0.5
  ax_width = 0.925 - ax_left
  axes = [fig.add_axes([ax_left,ax_ypos[0],ax_width,ax_heights[1]]), fig.add_axes([ax_left,ax_ypos[2],ax_width,ax_heights[3]]), fig.add_axes([ax_left,ax_ypos[4],ax_width,ax_heights[5]])][::-1]
  label_colors = ["#EFFFD6","#82B536","#37471F","orange"]
  colors = results_subset.apply(lambda row: label_colors[row.duration], axis=1)
  for di in ax_durations:
    ax = axes[di]
    ax.text(1.05,0.5,"%d-%d Days"%(durations[di][0], durations[di][1]), fontsize=7, fontweight="bold", color="gray", rotation=90,ha="left",va="center", transform=ax.transAxes)
    results_subset_dur = results_subset[results_subset.duration==di]
    for ri,row in results_subset_dur.reset_index(drop=True).iterrows():
      ax.scatter(row["coef"],ri,c="black",edgecolors='black',zorder=1, s=2, lw=0.25)
      region = scipy.stats.norm(loc=row["coef"],scale=row["stderror"]).interval(0.9)
      ax.plot([region[0],region[1]],[ri,ri],lw=0.5,c="black",zorder=0)
      ax.plot([region[0],region[0]],[ri-0.3,ri+0.3],lw=0.5,c="black",zorder=0)
      ax.plot([region[1],region[1]],[ri-0.3,ri+0.3],lw=0.5,c="black",zorder=0)
      ax.set_yticks(np.arange(len(results_subset_dur)))
      ax.set_yticklabels(results_subset_dur["microbe"].str.replace("Bifidobacterium ","B. ") + results_subset_dur["significant"].apply(lambda sig: "*" if sig else ""),fontsize=6)
      # if di == ax_durations[2]:
      #   pass
      #   # ax.set_xticklabels(ax.get_xticklabels(), fontsize=6)
      # else:
      #   ax.set_xticklabels([])
  xlim = (np.min([axes[0].get_xlim()[0], axes[1].get_xlim()[0], axes[2].get_xlim()[0]]), np.max([axes[0].get_xlim()[1], axes[1].get_xlim()[1], axes[2].get_xlim()[1]]))
  # axes[-1].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4,steps=[1, 2, 4, 5, 10]))
  axes[-1].set_xticklabels(axes[-1].get_xticklabels(), fontsize=6)
  for ax in axes[:-1]:
    ax.set_xticklabels([])
  fig.suptitle("\n".join(wrap("Associations Between Oral Lactulose and %s"%microbe,30)),fontsize=8)
  axes[-1].set_xlabel("Log-Fold Relative\nAbundance",fontsize=8)
  axes[1].set_ylabel("Taxon",fontsize=8)
  for ax in axes:
    ax.set_xlim(xlim)
    ax.axvline(0, color="gray", zorder=0, linestyle="dotted", lw=0.5)
  fig.savefig("out/figure_2f.pdf")
  plt.close()

def make_figure_3a():
  grid = (results.table[(results.table["duration"]==1) & (results.table["sampletype"]=="pathway") & (results.table.modeltype=="binary")].set_index(["medication","microbe","duration","sampletype"])["significant"] * results.table[(results.table["duration"]==1) & (results.table["sampletype"]=="pathway") & (results.table.modeltype=="binary")].set_index(["medication","microbe","duration","sampletype"])["coef"]).unstack(level=0).T
  grid_all = (results.table[(results.table["sampletype"]=="pathway") & (results.table.modeltype=="binary")].set_index(["medication","microbe","duration","sampletype"])["significant"] * results.table[(results.table["duration"]==1) & (results.table["sampletype"]=="pathway") & (results.table.modeltype=="binary")].set_index(["medication","microbe","duration","sampletype"])["coef"]).unstack(level=0).T.fillna(0)
  grid_all.to_csv("out/supp_table_pathways.csv")
  grid_significant = (results.table[(results.table["duration"]==1) & (results.table["sampletype"]=="pathway") & (results.table["modeltype"]=="binary")].set_index(["medication","microbe","duration","sampletype"])["significant"]).unstack(level=0).T
  grid = grid.loc[grid_significant.mean(axis=1).sort_values(ascending=False).index.get_level_values(0)[:45],grid_significant.mean(axis=0).sort_values(ascending=False).index.get_level_values(0)[:30]]
  # grid = grid.loc[grid.sum(axis=1)>5,grid.sum(axis=0)>1]
  grid.columns = grid.columns.get_level_values(0)
  ordering_ax0 = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid.fillna(0), method="centroid"), no_plot=True, color_threshold=-np.inf)['leaves']
  ordering_ax1 = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid.fillna(0).T, method="centroid"), no_plot=True, color_threshold=-np.inf)['leaves']
  grid = grid.iloc[ordering_ax0,ordering_ax1]
  grid = grid.T
  grid_pathway_taxa = grid.index.map(lambda x: table_pathways[table_pathways.index.str.contains(x,regex=False)].sum(axis=1)).map(lambda df: pd.Series(df.values,index=df.index.str.split("(.g__|.s__)",regex=True).str[2]))
  grid_pathway_taxa = grid_pathway_taxa.map(lambda df: df.groupby(df.index).sum()).map(lambda df: pd.concat([df,pd.Series([0],index=["unclassified"])]))
  grid_pathway_taxa = grid_pathway_taxa.map(lambda df: df[~df.index.isna()].idxmax())
  grid = grid.iloc[grid_pathway_taxa.argsort()[::-1]]
  grid_pathway_taxa = grid_pathway_taxa[grid_pathway_taxa.argsort()[::-1]]
  # grid_pathway_taxa = grid_pathway_taxa.map(lambda x: "" if x=="unclassified" else "  ("+x+")")
  grid.columns = grid.columns.map(lambda name: abbr(";".join([prettify_medication_name(subname.lower(), length=50).title() for subname in name.split("|")]), 60))
  grid.index = grid.index.str.split(":").str[1].str.replace("&alpha;", "α").str.replace("&beta;", "β").map(lambda x: abbr(x,90))
  # grid.index = grid.index + grid_pathway_taxa
  fig = plt.figure(figsize=(7.5,4.5))
  spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, wspace=0.05, left=0.42, bottom=0.42, right=0.89, top=0.925)
  axes = [fig.add_subplot(spec[0])]
  axes = np.array([axes]).reshape(-1)
  dot_heatmap(grid.astype(float), axes[0], axes[0], matplotlib.cm.get_cmap('Blues_r', 128), matplotlib.cm.get_cmap('Oranges', 128), size=0.5, exponent=0.4, legend_values=[-8,-4,-2,2,4,8], legend_labels=["-8.0","-4.0","-2.0","2.0","4.0","8.0"], legend_kws={"title_fontsize":8,"fontsize":8,"title":"Effect Size","ncol":2,"bbox_to_anchor":(-0.75, -0.5, 0.5, 0.5), "loc":"upper left"})
  for ti,tx in enumerate(grid_pathway_taxa):
    axes[0].text(grid.shape[1]+0.2,ti,tx,fontsize=6,ha="left",va="center",c="#505050")
  axes[0].set_yticklabels(axes[0].get_yticklabels(),fontsize=6)
  axes[0].set_xticklabels(axes[0].get_xticklabels(),fontsize=6,rotation=90,va="center")
  axes[0].set_xlabel("Medication", fontsize=8)
  axes[0].set_ylabel("Pathway", fontsize=8)
  axes[0].set_title("")
  axes[0].text(-1,-3,"Medication",rotation=90,ha="right",va="top",fontweight="bold",fontsize=7)
  axes[0].text(grid.shape[1]+0.25,-2,"Genus Driving\nPathway\nAbundance",ha="left",va="top",fontweight="bold",fontsize=7)
  fig.suptitle("Most Significant 2-10 Day Pathway Associations",fontsize=8)
  fig.savefig("out/results_heatmap_microbiome_pathway_subset.pdf")
  fig.savefig("out/figure_3a.pdf")
  plt.close()

def make_figure_3b():
  npathwaycols = 50         # 50
  nupsetcols = 15            # 15
  height_ratios = [1,0.4]   # [1,0.4]
  gridspec_left = 0.45      # 0.28
  gridspec_right = 0.97     # 0.97
  gridspec_top = 0.90       # 0.87
  gridspec_bottom = 0.4    # 0.47
  figsize = (4.5,4)      # (7.5,3.25)
  results_steroids = results.table[(results.table.significant) & (results.table.sampletype=="pathway") & results.table.medication_pharm_class.str.contains("ANTIB")].copy()
  # results_steroids = results.table[(results.table.modeltype=="binary") & (results.table.significant) & (results.table.sampletype=="pathway") & (results.table.medication_pharm_class==df.index[0])].copy()
  results_steroids["pathway_general"] = pathways.pathway_to_pathway_general["1"].loc[results_steroids.microbe].values
  # include_pathways = results_steroids.groupby("medication").head(100)["pathway_general"].drop_duplicates()
  # include_pathways = results_steroids["pathway_general"].value_counts().head(30).index
  include_pathways = results_steroids.head(npathwaycols)["pathway_general"].drop_duplicates()
  df = pd.DataFrame(results.table.groupby("medication_pharm_class")["medication"].unique())
  df["count"] = df["medication"].str.len()
  df = df.sort_values("count",ascending=False)
  df_upset = results_steroids[["medication","microbe"]].groupby("microbe").agg(set).reset_index()
  mlb = MultiLabelBinarizer()
  df_upset = df_upset.join(pd.DataFrame(mlb.fit_transform(df_upset.pop('medication')),columns=mlb.classes_,index=df_upset.index))
  df_upset = df_upset.set_index(pd.MultiIndex.from_frame(df_upset.iloc[:,1:].astype(bool))).drop(df_upset.columns[1:],axis=1)
  df_upset_ix = df_upset.index.names
  df_upset = df_upset.groupby(df_upset.index).size()
  df_upset.index = pd.MultiIndex.from_frame(pd.DataFrame(np.stack(df_upset.index), columns=df_upset_ix))
  pairs = list(itertools.product(df_upset.index.names,df_upset.index.names))
  counts = [df_upset[(df_upset.index.get_level_values(pair[0]).values) & (df_upset.index.get_level_values(pair[1]).values)].sum() for pair in pairs]
  pair_counts = pd.DataFrame(pairs)
  pair_counts["counts"] = counts
  pair_counts = pair_counts.sort_values("counts",ascending=False)
  df_upset = df_upset.sort_values(ascending=False).head(nupsetcols)
  # grid = results_steroids.pivot_table(columns=["pathway_general"], values="significant", index=["medication"], aggfunc=any).T.astype(float)
  grid = results_steroids.pivot_table(columns=["pathway_general"], values="significant", index=["medication"], aggfunc=any) * results_steroids.pivot_table(columns=["pathway_general"], values="coef", index=["medication"], aggfunc=np.mean)
  grid = np.sign(grid.fillna(0))
  # grid = grid.T
  grid_pval = -np.log10(results_steroids.pivot_table(columns=["pathway_general"], values="pvalue", index=["medication"], aggfunc=min).T.astype(float)).T
  grid.to_csv("out/supp_table_antibiotics_pathways.csv")
  grid = grid.loc[:,include_pathways]
  ordering_X = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid.fillna(0)), no_plot=True, color_threshold=-np.inf)['leaves']
  ordering_X = grid.sum(axis=1).argsort()
  grid = grid.iloc[ordering_X]
  grid = grid.loc[df_upset.index.names]
  grid = grid.loc[:,grid.sum().sort_values(ascending=False).index]
  fig,ax = plt.subplots(1,1,figsize=figsize,constrained_layout=True)
  fig = plt.figure(figsize=figsize)
  spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, left=0.64, bottom=0.26, right=0.98, top=0.93)
  ax = fig.add_subplot(spec[0,0])
  # spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, width_ratios=[1,0.75], height_ratios=height_ratios, wspace=0.015, hspace=0.10, left=gridspec_left, bottom=gridspec_bottom, right=gridspec_right, top=gridspec_top)
  # axes = [fig.add_subplot(spec[0:2,0]), fig.add_subplot(spec[1,1]), fig.add_subplot(spec[0,1])]
  # axes = np.array([axes]).reshape(-1)
  colors = [(0, 0, 0), (0, 0, 0), (0, 0, 0)]  # R -> G -> B
  n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
  cmap_name = 'my_list'
  cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=10)
  cmap = matplotlib.cm.get_cmap("coolwarm").copy()
  base_cmap = matplotlib.cm.get_cmap('RdYlGn', 512)
  cmap = ListedColormap(base_cmap(np.linspace(0.1, 0.8, 256)))
  dot_heatmap(grid.loc[grid.index,grid.columns].astype(float).T, ax, ax, cmap, cmap, min_size=9, size=0.2, exponent=1, legend_values=[1,-1], legend_labels=["positive","negative"], legend_kws={"title_fontsize":8,"fontsize":8,"title":"Direction of Effect","ncol":1,"bbox_to_anchor":(-1.1, -0.3, 0.5, 0.5), "loc":"lower left"}, cmap_skip=0.4, marker="o")
  # sns.heatmap(grid, axes[0]=axes[0], cmap=cmap, vmin=0, linewidths=0.25, linecolor="black", cbar_kws={"label":"Log-fold change in Relative Abundance"})
  fig.suptitle("Associations Between Antibiotics and Microbial Pathway Categories",fontsize=8)
  ax.set_title("",fontsize=7)
  ax.set_xlabel("Medication",fontsize=8)
  ax.set_ylabel("Microbial Pathway",fontsize=8)
  # ax.set_yticks(np.arange(len(grid.index))+0.5)
  # ax.set_yticklabels(grid.index.str.split(":").str[1].str.replace("&alpha;","α").str.replace("&beta;","β"), fontsize=5)
  # ax.set_xticklabels(grid.columns.to_frame().apply(lambda x: "\n".join(wrap(prettify_compound_medication_name(x[0]),35))), fontsize=6, rotation=40,ha="right",rotation_mode="anchor")
  # ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
  ax.set_xticks(np.arange(len(grid.index)))
  ax.set_yticks(np.arange(len(grid.columns)))
  ax.set_xticklabels(grid.index.map(lambda x: "\n".join(wrap(prettify_compound_medication_name(x),35))), fontsize=6, rotation=45,ha="right",rotation_mode="anchor")
  ax.set_yticklabels(grid.columns.str.replace("&alpha;","α").str.replace("&beta;","β"), fontsize=6, rotation=0,ha="right",rotation_mode="anchor")
  # ax.set_xticks(np.arange(len(grid.columns))+0.5)
  # ax.set_yticklabels(grid.index.str.replace("&alpha;","α").str.replace("&beta;","β"), axis=1), fontsize=6.5)
  # axes[1].set_yticks([])
  # axes[1].set_xticks([])
  # ax.set_ylim(axes[0].get_ylim()[0]+0.5, axes[0].get_ylim()[1]-0.5)
  # axes[1].set_ylim(axes[0].get_ylim())
  np.where(np.ones([4,6]))
  df_upset = df_upset.T
  df_upset_grid = np.stack(df_upset.index.values)
  # axes[1].set_title("Shared Pathway\nAssociations",fontsize=7)
  # axes[1].scatter(np.where(np.ones([len(df_upset.index[0]),len(df_upset.index)]))[1], np.where(np.ones([len(df_upset.index[0]),len(df_upset.index)]))[0],c="#e0e0e0",s=9)
  # axes[1].scatter(np.where(df_upset_grid)[0], np.where(df_upset_grid)[1],c="black",s=9)
  lines_y = np.stack([df_upset_grid.argmax(axis=1), df_upset_grid.shape[1]-1-df_upset_grid[:,::-1].argmax(axis=1)])
  lines_x = np.tile(np.arange(df_upset_grid.shape[0]),[2,1])
  # axes[1].plot(lines_x,lines_y,c="black",lw=1)
  # axes[1].set_xlim(axes[1].get_xlim()[0]-0.5,axes[1].get_xlim()[1]+0.5)
  # axes[1].set_xticks(np.arange(df_upset_grid.shape[0]), df_upset,fontsize=6)
  # axes[2].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
  # axes[1].set_xticks(np.arange(len(grid.index)))
  # axes[1].set_xticklabels(grid.index.map(lambda x: "\n".join(wrap(prettify_compound_medication_name(x),35))), fontsize=6, rotation=45,ha="right",rotation_mode="anchor")
  # axes[2].bar(np.arange(df_upset_grid.shape[0]), df_upset.values, color="gray")
  # axes[2].set_xlim(axes[1].get_xlim())
  # axes[2].set_xticklabels([])
  # axes[2].set_xlabel("# Pathway\nAssociations",fontsize=7)
  # axes[2].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5,steps=[1, 2, 4, 5, 10]))
  # axes[2].set_yticklabels(axes[2].get_yticklabels(),fontsize=6)
  # axes[2].yaxis.set_label_position("right")
  # axes[2].yaxis.tick_right()
  fig.savefig("out/figure_3b.pdf")
  plt.close()

# set(results_steroids[(results_steroids.medication=="VANCOMYCIN_Intravenous")].microbe).intersection(results_steroids[(results_steroids.medication=="CEFAZOLIN_Intravenous")].microbe)

def make_figure_3c():
  npathwaycols = 50         # 50
  nupsetcols = 15            # 15
  height_ratios = [1,0.4]   # [1,0.4]
  gridspec_left = 0.37      # 0.28
  gridspec_right = 0.97     # 0.97
  gridspec_top = 0.90       # 0.87
  gridspec_bottom = 0.4    # 0.47
  results_steroids = results.table[(results.table.significant) & (results.table.sampletype=="pathway") & results.table.medication_pharm_class.str.contains("ANTIB")].copy()
  results_steroids["pathway_general"] = pathways.pathway_to_pathway_general["1"].loc[results_steroids.microbe].values
  include_pathways = results_steroids.head(npathwaycols)["pathway_general"].drop_duplicates()
  df_upset = results_steroids[["medication","microbe"]].groupby("microbe").agg(set).reset_index()
  mlb = MultiLabelBinarizer()
  df_upset = df_upset.join(pd.DataFrame(mlb.fit_transform(df_upset.pop('medication')),columns=mlb.classes_,index=df_upset.index))
  df_upset = df_upset.set_index(pd.MultiIndex.from_frame(df_upset.iloc[:,1:].astype(bool))).drop(df_upset.columns[1:],axis=1)
  df_upset_ix = df_upset.index.names
  df_upset = df_upset.groupby(df_upset.index).size()
  df_upset.index = pd.MultiIndex.from_frame(pd.DataFrame(np.stack(df_upset.index), columns=df_upset_ix))
  pairs = list(itertools.product(df_upset.index.names,df_upset.index.names))
  counts = [df_upset[(df_upset.index.get_level_values(pair[0]).values) & (df_upset.index.get_level_values(pair[1]).values)].sum() for pair in pairs]
  pair_counts = pd.DataFrame(pairs)
  pair_counts["counts"] = counts
  pair_counts = pair_counts.sort_values("counts",ascending=False)
  df_upset = df_upset.sort_values(ascending=False).head(nupsetcols)
  grid = results_steroids.pivot_table(columns=["pathway_general"], values="significant", index=["medication"], aggfunc=any) * results_steroids.pivot_table(columns=["pathway_general"], values="coef", index=["medication"], aggfunc=np.mean)
  grid = np.sign(grid.fillna(0))
  # grid = grid.T
  grid_pval = -np.log10(results_steroids.pivot_table(columns=["pathway_general"], values="pvalue", index=["medication"], aggfunc=min).T.astype(float)).T
  grid = grid.loc[:,include_pathways]
  ordering_X = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid.fillna(0)), no_plot=True, color_threshold=-np.inf)['leaves']
  ordering_X = grid.sum(axis=1).argsort()
  grid = grid.iloc[ordering_X]
  grid = grid.loc[df_upset.index.names]
  df_upset = results_steroids[["medication","microbe"]].groupby("microbe").agg(set).reset_index()
  mlb = MultiLabelBinarizer()
  df_upset = df_upset.join(pd.DataFrame(mlb.fit_transform(df_upset.pop('medication')),columns=mlb.classes_,index=df_upset.index))
  df_upset = df_upset.set_index(pd.MultiIndex.from_frame(df_upset.iloc[:,1:].astype(bool))).drop(df_upset.columns[1:],axis=1)
  df_upset_ix = df_upset.index.names
  df_upset = df_upset.groupby(df_upset.index).size()
  df_upset.index = pd.MultiIndex.from_frame(pd.DataFrame(np.stack(df_upset.index), columns=df_upset_ix))
  pairs = list(itertools.product(df_upset.index.names,df_upset.index.names))
  counts = [df_upset[(df_upset.index.get_level_values(pair[0]).values) & (df_upset.index.get_level_values(pair[1]).values)].sum() for pair in pairs]
  totals = [df_upset[(df_upset.index.get_level_values(medication).values)].sum() for medication in df_upset.index.names]
  pair_counts = pd.DataFrame(pairs)
  pair_counts["counts"] = counts
  pair_counts = pair_counts.sort_values("counts",ascending=False)
  df_upset = df_upset.sort_values(ascending=False).head(nupsetcols)
  fig = plt.figure(figsize=(3,4))
  spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, width_ratios=[1,0.67], height_ratios=[1,1], wspace=0.05, hspace=0.10, left=0.15, bottom=0.27, right=gridspec_right, top=0.93)
  axes = [fig.add_subplot(spec[1,0]), fig.add_subplot(spec[1,1]), fig.add_subplot(spec[0,0])]
  # axes = axes[::-1]
  # axes[1].set_ylim(axes[1].get_ylim())
  np.where(np.ones([4,6]))
  df_upset = df_upset.T
  df_upset_grid = np.stack(df_upset.index.values)
  fig.suptitle("Shared Pathway Associations",fontsize=8)
  axes[0].scatter(np.where(np.ones([len(df_upset.index[0]),len(df_upset.index)]))[0], np.where(np.ones([len(df_upset.index[0]),len(df_upset.index)]))[1],c="#e0e0e0",s=9)
  axes[0].scatter(np.where(df_upset_grid)[1], np.where(df_upset_grid)[0],c="black",s=9)
  lines_x = np.stack([df_upset_grid.argmax(axis=1), df_upset_grid.shape[1]-1-df_upset_grid[:,::-1].argmax(axis=1)])
  lines_y = np.tile(np.arange(df_upset_grid.shape[0]),[2,1])
  axes[0].plot(lines_x,lines_y,c="black",lw=1)
  axes[0].set_xticks(np.arange(df_upset_grid.shape[0]), df_upset,fontsize=6)
  # axes[1].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
  axes[0].set_xticks(np.arange(len(grid.index)))
  axes[0].set_xticklabels(grid.index.map(lambda x: "\n".join(wrap(prettify_compound_medication_name(x),35))), fontsize=6, rotation=45,ha="right",rotation_mode="anchor")
  axes[0].set_xlim(-0.5, len(df_upset.index.names)-0.5)
  axes[0].set_yticklabels([])
  axes[1].barh(np.arange(df_upset_grid.shape[0]), df_upset.values, color="gray")
  axes[1].set_ylim(axes[0].get_ylim())
  axes[1].set_yticks([])
  axes[1].set_xlabel("# Pathway\nAssociations",fontsize=7)
  axes[1].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5,steps=[1, 2, 4, 5, 10]))
  axes[1].set_xticklabels(axes[1].get_xticklabels(),fontsize=6)
  grid = pair_counts.pivot(columns=0,index=1,values="counts")
  grid.to_csv("out/supp_table_num_shared_antibiotic_pathways.csv")
  # grid_totals = 
  vals = grid.values.astype(float)
  vals[np.triu_indices(vals.shape[0])] = np.nan
  grid.iloc[:,:] = vals.T
  # sns.heatmap(grid,ax=axes[2])
  cmap = matplotlib.cm.get_cmap("gnuplot2_r").copy()
  dot_heatmap(grid, axes[2], axes[0], cmap, cmap, min_size=0.02, size=0.0004, exponent=2, legend_values=[40,80,120,160], legend_labels=["40","80","120","160"], legend_kws={"title_fontsize":8,"fontsize":8,"title":"# Shared Associations","ncol":2,"bbox_to_anchor":(0.8, -0.87, 0.37, 0.5), "loc":"lower left","handletextpad":0.25,"handlelength":2.0,"columnspacing":2.0}, cmap_skip=0.0, marker="o")
  axes[2].set_xlim(axes[0].get_xlim())
  axes[2].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
  axes[2].yaxis.tick_right()
  axes[2].set_yticklabels(grid.index.map(lambda x: "\n".join(wrap(prettify_compound_medication_name(x),35))), fontsize=6)
  axes[2].set_title("")
  axes[2].set_ylabel("")
  axes[2].set_xlabel("")
  axes[2].text(-0.02,0.5,"Medication Pairs",fontsize=7,ha="right",va="center",rotation=90,transform=axes[2].transAxes)
  axes[0].text(-0.02,0.5,"Medication Subsets",fontsize=7,ha="right",va="center",rotation=90,transform=axes[0].transAxes)
  axes[0].tick_params(axis='y',which='both',left=False,labelleft=False)
  axes[0].set_xlabel("Medication",fontsize=8)
  # axes[0].set_title("Medication Subsets",fontsize=7)
  # axes[1].yaxis.set_label_position("right")
  # axes[1].yaxis.tick_right()
  fig.savefig("out/figure_3c.pdf")
  plt.close()

# if True:
#   grid = pd.pivot(results.table[(results.table.sampletype=="metab") & (results.table.significant) & (results.table.modeltype=="binary") & (results.table.duration==0)], columns="medication",index="microbe",values="coef")
#   grid_small = grid
#   # if False:
#   # grid_small = grid.loc[grid.abs().sum(axis=1).sort_values(ascending=False).head(50).index, grid.abs().sum(axis=0).sort_values(ascending=False).head(60).index]
#   # # grid_small = grid.loc[(~grid.isna()).sum(axis=1).sort_values(ascending=False,key=abs).head(50).index,(~grid.isna()).sum(axis=0).sort_values(ascending=False,key=abs).head(30).index].replace([np.inf, -np.inf], np.nan)
#   # grid_small = grid_small.loc[(~grid_small.isna()).sum(axis=1)>0,(~grid_small.isna()).sum(axis=0)>0]
#   # ranking0 = grid.abs().sum(axis=0).sort_values(ascending=False)
#   # ranking1 = grid.abs().sum(axis=1).sort_values(ascending=False)
#   # threshold = 1
#   ranking0 = ((0.5*grid.abs().max(axis=0).abs().rank()) + (0.5*(~grid.isna()).sum(axis=0).abs().rank())).sort_values(ascending=False)
#   ranking1 = ((0.5*grid.abs().max(axis=1).abs().rank()) + (0.5*(~grid.isna()).sum(axis=1).abs().rank())).sort_values(ascending=False)
#   threshold = 1
#   grid_small = grid.loc[ranking1.head(50).index,ranking0.head(50).index]
#   grid_small = grid_small.loc[(~grid_small.isna()).sum(axis=1)>threshold,(~grid_small.isna()).sum(axis=0)>threshold]

def make_figure_4a():
  grid = pd.pivot(results.table[(results.table.sampletype=="metab") & (results.table.significant) & (results.table.modeltype=="binary") & (results.table.duration==0)], columns="medication",index="microbe",values="coef")
  grid_pathways = pathways.results_enrichment_medication_metabolites[(pathways.results_enrichment_medication_metabolites.significant) & (pathways.results_enrichment_medication_metabolites.direction==0)].pivot(columns="medication",index="pathway",values="pvalue")
  pd.concat([grid,grid_pathways]).fillna(0).to_csv("out/supp_table_metabolites.csv")
  grid_small = grid
  # if False:
  # grid_small = grid.loc[grid.abs().sum(axis=1).sort_values(ascending=False).head(50).index, grid.abs().sum(axis=0).sort_values(ascending=False).head(60).index]
  # # grid_small = grid.loc[(~grid.isna()).sum(axis=1).sort_values(ascending=False,key=abs).head(50).index,(~grid.isna()).sum(axis=0).sort_values(ascending=False,key=abs).head(30).index].replace([np.inf, -np.inf], np.nan)
  # grid_small = grid_small.loc[(~grid_small.isna()).sum(axis=1)>0,(~grid_small.isna()).sum(axis=0)>0]
  # ranking0 = grid.abs().sum(axis=0).sort_values(ascending=False)
  # ranking1 = grid.abs().sum(axis=1).sort_values(ascending=False)
  # threshold = 1
  ranking0 = ((0.5*grid.abs().max(axis=0).abs().rank()) + (0.5*(~grid.isna()).sum(axis=0).abs().rank())).sort_values(ascending=False)
  ranking1 = ((0.5*grid.abs().max(axis=1).abs().rank()) + (0.5*(~grid.isna()).sum(axis=1).abs().rank())).sort_values(ascending=False)
  threshold = 0
  grid_small = grid.loc[ranking1.head(50).index,ranking0.head(50).index]
  grid_small = grid_small.loc[(~grid_small.isna()).sum(axis=1)>0,(~grid_small.isna()).sum(axis=0)>1]
  grid_bile = grid_small.loc[grid_small.index.str.split("#").str[0]=="bile"]
  grid_scfa = grid_small.loc[grid_small.index.str.split("#").str[0]=="scfa"]
  grid_indole = grid_small.loc[grid_small.index.str.split("#").str[0]=="indole"]
  grid_pathways = pathways.results_enrichment_medication_metabolites[(pathways.results_enrichment_medication_metabolites.significant) & (pathways.results_enrichment_medication_metabolites.direction==0)].pivot(columns="medication",index="pathway",values="pvalue")
  grid_pathways = -np.log10(grid_pathways.loc[:,grid_small.columns.intersection(grid_pathways.columns)])
  grid_pathways = grid_pathways.sort_index(key=lambda ix: ix.str.len()).drop_duplicates().sort_index()
  grid_pathways[grid_pathways==np.inf] = 4
  # grid_pathways = (~grid_pathways.loc[:,grid_indole.columns.intersection(grid_pathways.columns)].drop_duplicates().isna()).astype(float) * 1.0
  grids = [grid_bile, grid_scfa, grid_indole,grid_pathways]
  for g in grids:
    ordering_y = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(g.fillna(0)), no_plot=True, color_threshold=-np.inf)['leaves']
    ordering_x = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(g.fillna(0).T), no_plot=True, color_threshold=-np.inf)['leaves']
    g = g.iloc[ordering_y, ordering_x]
  grid_small = pd.concat(grids)
  fig,ax = plt.subplots(1,1,figsize=(7.5,5.5))
  grid_small_mask1 = grid_small.copy()
  grid_small_mask2 = grid_small.copy()
  grid_small_mask1.loc[grid_pathways.index] = pd.NA
  grid_small_mask2.loc[~grid_small_mask1.index.isin(grid_pathways.index)] = pd.NA
  # grid_small.T.reset_index(drop=True).T.reset_index(drop=True).T.astype(float).stack().reset_index(name="corr")
  colors = ["black", "black"] 
  single_color_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("single_color_map", colors)
  dot_heatmap(grid_small_mask1.T, ax, ax, matplotlib.cm.get_cmap('Blues_r', 128), matplotlib.cm.get_cmap('Oranges', 128), None, cmap_skip=0.4, size=0.03, min_size=2, exponent=1.5,legend_values=[-8,-4,-2,2,8,16], legend_labels=["-8.0","-4.0","-2.0","2.0","8.0","16.0"], legend_kws={"title_fontsize":8,"fontsize":8,"title":"Metabolite Effect Size","ncol":2,"bbox_to_anchor":(0.15, -0.68, 0.5, 0.5), "loc":"lower left", "columnspacing":1.0, "handlelength":1.0})
  dot_heatmap(grid_small_mask2.T, ax, ax, matplotlib.cm.get_cmap('plasma_r', 128), matplotlib.cm.get_cmap('plasma_r', 128), None, cmap_skip=0.2, size=0.7, min_size=3, exponent=1, onlyplot=True, legend_values=[2,4,6,8], legend_labels=["2.0","4.0","6.0",">8.0"], legend_kws={"title_fontsize":8,"fontsize":8,"title":"Pathway -$\mathregular{log_{10}}$ p","ncol":2,"bbox_to_anchor":(0.5,-0.64,0.5,0.5), "loc":"lower right", "columnspacing":1.0, "handlelength":1.0})
  fig.suptitle("Significant Medication-Metabolite Associations After 2-10 Days",fontsize=8)
  fig.subplots_adjust(top=0.93,left=0.21,right=0.99,bottom=0.5)
  ax.set_title(None)
  ax.set_xticklabels(grid_small.index.str.split("#").str[-1], rotation=45, ha="right", va="center",fontsize=6,rotation_mode="anchor")
  ax.set_yticks(np.arange(len(grid_small.columns)))
  ax.set_yticklabels(grid_small.columns.map(lambda x: prettify_compound_medication_name(x,mode="trunc",length=30)),fontsize=6)
  ax.set_ylabel("Medication",fontsize=8)
  ax.set_xlabel("Metabolite",fontsize=8,ha="center")
  grid_sizes = [g.shape[0] for g in grids]
  grid_sizes_cum = np.cumsum(grid_sizes)
  ax.plot([-0.3,grid_sizes_cum[0]-0.7], [grid_small.shape[1] + 0.5,grid_small.shape[1] + 0.5], color="black", lw=1, clip_on=False)
  ax.plot([grid_sizes_cum[0]-0.3,grid_sizes_cum[1]-0.7], [grid_small.shape[1] + 0.5,grid_small.shape[1] + 0.5], color="black", lw=1, clip_on=False)
  ax.plot([grid_sizes_cum[1]-0.3,grid_sizes_cum[2]-0.7], [grid_small.shape[1] + 0.5,grid_small.shape[1] + 0.5], color="black", lw=1, clip_on=False)
  ax.plot([grid_sizes_cum[2]-0.3,grid_sizes_cum[3]-0.7], [grid_small.shape[1] + 0.5,grid_small.shape[1] + 0.5], color="black", lw=1, clip_on=False)
  ax.text(grid_sizes_cum[0]/2.0, grid_small.shape[1] + 0.5, "Bile Acid", va="bottom", ha="center", fontsize=6)
  ax.text((grid_sizes_cum[0]+grid_sizes_cum[1])/2.0, grid_small.shape[1] + 0.5, "SCFA", va="bottom", ha="center", fontsize=6)
  ax.text((grid_sizes_cum[1]+grid_sizes_cum[2])/2.0, grid_small.shape[1] + 0.5, "Indole", va="bottom", ha="center", fontsize=6)
  ax.text((grid_sizes_cum[2]+grid_sizes_cum[3])/2.0, grid_small.shape[1] + 0.5, "Metabolomic Pathway", va="bottom", ha="center", fontsize=6)
  fig.savefig("out/figure_4a.pdf")
  plt.close()

def make_supplementary_figure_metabolites():
  grid = pd.pivot(results.table[(results.table.sampletype=="metab") & (results.table.significant) & (results.table.modeltype=="binary") & (results.table.duration==0)], columns="medication",index="microbe",values="coef")
  grid_small =grid
  # grid_small = grid.loc[(~grid.isna()).sum(axis=1).sort_values(ascending=False,key=abs).head(50).index,(~grid.isna()).sum(axis=0).sort_values(ascending=False,key=abs).head(30).index].replace([np.inf, -np.inf], np.nan)
  grid_small = grid_small.loc[(~grid_small.isna()).sum(axis=1)>0,(~grid_small.isna()).sum(axis=0)>0]
  grid_bile = grid_small.loc[grid_small.index.str.split("#").str[0]=="bile"]
  grid_scfa = grid_small.loc[grid_small.index.str.split("#").str[0]=="scfa"]
  grid_indole = grid_small.loc[grid_small.index.str.split("#").str[0]=="indole"]
  grid_pathways = pathways.results_enrichment_medication_metabolites[(pathways.results_enrichment_medication_metabolites.significant) & (pathways.results_enrichment_medication_metabolites.direction==0)].pivot(columns="medication",index="pathway",values="pvalue")
  grid_pathways = -np.log10(grid_pathways.loc[:,grid_indole.columns.intersection(grid_pathways.columns)].drop_duplicates())
  grid_pathways[grid_pathways==np.inf] = 4
  # grid_pathways = (~grid_pathways.loc[:,grid_indole.columns.intersection(grid_pathways.columns)].drop_duplicates().isna()).astype(float) * 1.0
  grids = [grid_bile, grid_scfa, grid_indole,grid_pathways]
  for g in grids:
    ordering_y = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(g.fillna(0)), no_plot=True, color_threshold=-np.inf)['leaves']
    ordering_x = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(g.fillna(0).T), no_plot=True, color_threshold=-np.inf)['leaves']
    g = g.iloc[ordering_y, ordering_x]
  grid_small = pd.concat(grids)
  fig,ax = plt.subplots(1,1,figsize=(7.5,10))
  grid_small_mask1 = grid_small.copy()
  grid_small_mask2 = grid_small.copy()
  grid_small_mask1.loc[grid_pathways.index] = pd.NA
  grid_small_mask2.loc[~grid_small_mask1.index.isin(grid_pathways.index)] = pd.NA
  # grid_small.T.reset_index(drop=True).T.reset_index(drop=True).T.astype(float).stack().reset_index(name="corr")
  colors = ["black", "black"] 
  single_color_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("single_color_map", colors)
  dot_heatmap(grid_small_mask1.T, ax, ax, matplotlib.cm.get_cmap('Blues_r', 128), matplotlib.cm.get_cmap('Oranges', 128), None, size=0.1, min_size=1, exponent=1.5,legend_values=[-6,-4,-2,2,4,6], legend_labels=["-6.0","-4.0","-2.0","2.0","4.0","6.0"], legend_kws={"title_fontsize":8,"fontsize":8,"title":"Metabolite Effect Size","ncol":2,"bbox_to_anchor":(-0.4, -0.5, 0.5, 0.5), "loc":"upper left"})
  dot_heatmap(grid_small_mask2.T, ax, ax, matplotlib.cm.get_cmap('plasma_r', 128), matplotlib.cm.get_cmap('plasma_r', 128), None, cmap_skip=0.2, size=1, min_size=3, exponent=1.2, onlyplot=True, legend_values=[1, 2, 3, 4], legend_labels=["1.0","2.0","3.0",">4.0"], legend_kws={"title_fontsize":8,"fontsize":8,"title":"Pathway -$\mathregular{log_{10}}$ p","ncol":2,"bbox_to_anchor":(-0.4, -0.75, 0.5, 0.25), "loc":"lower left"})
  fig.suptitle("Significant Medication-Metabolite Associations After 2-10 Days",fontsize=8)
  fig.subplots_adjust(top=0.96,left=0.2,right=0.99,bottom=0.3)
  ax.set_title(None)
  ax.set_xticklabels(grid_small.index.str.split("#").str[-1], rotation=90,fontsize=6)
  ax.set_yticks(np.arange(len(grid_small.columns)))
  ax.set_yticklabels(grid_small.columns.map(lambda x: prettify_compound_medication_name(x,mode="trunc",length=15)),fontsize=6)
  ax.set_ylabel("Medication",fontsize=8)
  ax.set_xlabel("Metabolite",fontsize=8,ha="center")
  grid_sizes = [g.shape[0] for g in grids]
  grid_sizes_cum = np.cumsum(grid_sizes)
  ax.plot([-0.3,grid_sizes_cum[0]-0.7], [grid_small.shape[1] + 0.5,grid_small.shape[1] + 0.5], color="black", lw=1, clip_on=False)
  ax.plot([grid_sizes_cum[0]-0.3,grid_sizes_cum[1]-0.7], [grid_small.shape[1] + 0.5,grid_small.shape[1] + 0.5], color="black", lw=1, clip_on=False)
  ax.plot([grid_sizes_cum[1]-0.3,grid_sizes_cum[2]-0.7], [grid_small.shape[1] + 0.5,grid_small.shape[1] + 0.5], color="black", lw=1, clip_on=False)
  ax.plot([grid_sizes_cum[2]-0.3,grid_sizes_cum[3]-0.7], [grid_small.shape[1] + 0.5,grid_small.shape[1] + 0.5], color="black", lw=1, clip_on=False)
  ax.text(grid_sizes_cum[0]/2.0, grid_small.shape[1] + 0.5, "Bile Acid", va="bottom", ha="center", fontsize=6)
  ax.text((grid_sizes_cum[0]+grid_sizes_cum[1])/2.0, grid_small.shape[1] + 0.5, "SCFA", va="bottom", ha="center", fontsize=6)
  ax.text((grid_sizes_cum[1]+grid_sizes_cum[2])/2.0, grid_small.shape[1] + 0.5, "Indole", va="bottom", ha="center", fontsize=6)
  ax.text((grid_sizes_cum[2]+grid_sizes_cum[3])/2.0, grid_small.shape[1] + 0.5, "Metabolomic Pathway", va="bottom", ha="center", fontsize=6)
  fig.savefig("out/figure_supplementary_metabolites_all.pdf")
  plt.close()

def figure_4b_binomial_test(df):
  grid = df["medication"].apply(len).unstack().fillna(0)
  mu = grid.mean().mean()
  print("mean", grid.mean().mean())
  print(grid)
  for size in df["size"].unique():
    binomial_prob = scipy.stats.binomtest(size,len(medications),grid.mean().mean()/len(medications),alternative="greater")
    print("Probability of %d medications shared between genus and metabolite is: %4f. p=%4f"%(size, binomial_prob.pvalue,grid.mean().mean()/len(medications)))
  print(df.to_string())
  print(df["size"].value_counts())

def make_figure_4b():
  pair_to_max_coef = None
  pair_to_min_pval = None
  def max_abs(xs):
    return xs[np.abs(np.array(xs)).argmax()]
  nrows=10
  alpha=0
  pair_to_max_coef = None
  pair_to_min_pval = None
  column="genus"
  if alpha > 0:
    results_binary_metab = results.table[(results.table.modeltype.isin(["binary"])) & (results.table.sampletype=="metab")].copy()
    results_binary_column = results.table[(results.table.modeltype.isin(["binary"])) & (results.table.sampletype==column)].copy()
    results_binary_metab["sig100"] = statsmodels.stats.multitest.multipletests(results_binary_metab["pvalue"], method="fdr_bh", alpha=alpha)[0]
    results_binary_column["sig100"] = statsmodels.stats.multitest.multipletests(results_binary_column["pvalue"], method="fdr_bh", alpha=alpha)[0]
    df = pd.concat([results_binary_column[results_binary_column.sig100][["medication","microbe","pvalue"]].groupby("medication").agg(set), 
    results_binary_metab[results_binary_metab.sig100][["medication","microbe","pvalue"]].groupby("medication").agg(set)], axis=1)
  else:
    df = pd.concat([results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype==column)][["medication","microbe","pvalue"]].groupby("medication").agg(set), 
    results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype=="metab")][["medication","microbe","pvalue"]].groupby("medication").agg(set)], axis=1)
  df.columns = ["genus","genus_pvalue","metab","metab_pvalue"]
  df["metab_pvalue"][df["metab_pvalue"].isna()]=np.array(set())
  df["genus_pvalue"][df["genus_pvalue"].isna()]=np.array(set())
  df["minp"] = df[["genus_pvalue","metab_pvalue"]].apply(lambda row: (row["genus_pvalue"] | row["metab_pvalue"]),axis=1).apply(lambda row: min(row))
  df.explode("genus").explode("metab")
  minp = df.explode("genus").explode("metab").reset_index().groupby(["genus","metab"])["minp"].agg(min)
  df = df.explode("genus").explode("metab").reset_index().groupby(["genus","metab"])["medication"].agg(set)
  df = pd.DataFrame(df)
  df["minp"] = minp
  df["neg_minp"] = -minp
  df["size"] = df["medication"].str.len()
  df = df.sort_values("size",ascending=False)
  df = df.sort_values("minp",ascending=True)
  df = df.sort_values(["size","neg_minp"],ascending=False)
  df[["medication","minp","size"]].to_csv("out/supp_table_metabolite_genus_pairs.csv")
  # figure_4b_binomial_test(df)
  # nrows = 15
  selected_genus = pd.Series(df.index.get_level_values(0)[:nrows]).drop_duplicates().values
  selected_metab = pd.Series(df.index.get_level_values(1)[:nrows]).drop_duplicates().values
  selected_meds = pd.Series(np.concatenate(df["medication"].head(nrows).apply(list).values)).drop_duplicates().values
  adj = pd.DataFrame(False, index = np.concatenate([selected_genus, selected_metab, selected_meds]), columns = np.concatenate([selected_genus, selected_metab, selected_meds]))
  # edges1 = df.head(nrows).explode("medication").reset_index()[["genus","metab"]].drop_duplicates()
  edges2 = df.head(nrows).explode("medication").reset_index()[["genus","medication"]].drop_duplicates()
  edges3 = df.head(nrows).explode("medication").reset_index()[["medication","metab"]].drop_duplicates()
  # edges1.columns = [0,1]
  edges2.columns = [0,1]
  edges3.columns = [0,1]
  edges = pd.concat([edges2,edges3])
  for ei,edge in edges.iterrows():
    adj.loc[edge[0],edge[1]] = True
    adj.loc[edge[1],edge[0]] = True
  graph = nx.Graph(adj.astype(bool))
  labels = {i:i for i in adj.index}
  colors = {x:"#a0a0a0" for x in selected_meds} | {x:"#a0ffa0" for x in selected_genus} | {x:"#a0a0ff" for x in selected_metab}
  df1 = df.head(nrows).explode("medication").reset_index()[["genus","medication"]]
  df2 = df.head(nrows).explode("medication").reset_index()[["metab","medication"]]
  df1.columns = [0,1]
  df2.columns = [0,1]
  df3 = pd.concat([df1,df2])
  results.table[results.table.modeltype=="binary"]
  # if pair_to_max_coef is None:
  pair_to_max_coef = results.table[results.table.modeltype=="binary"][["medication","microbe","coef"]].groupby(["medication","microbe"])["coef"].apply(lambda df: max_abs(df.values))
  # if pair_to_min_pval is None:
  pair_to_min_pval = results.table[results.table.modeltype=="binary"][["medication","microbe","pvalue"]].groupby(["medication","microbe"])["pvalue"].apply(lambda df: np.min(df.values))
  df3["coef"] = df3.apply(lambda row: pair_to_max_coef[(row[1],row[0])],axis=1)
  df3["pval"] = df3.apply(lambda row: pair_to_min_pval[(row[1],row[0])],axis=1)
  df3 = df3.set_index([0,1])
  max_col = np.abs(max_abs(df3["coef"].values))
  pd.concat([df.head(nrows).explode("medication").reset_index()[["genus","medication"]], df.head(nrows).explode("medication").reset_index()[["metab","medication"]]])
  fig, ax = plt.subplots(figsize=(3.7,3.3))
  fig.subplots_adjust(top=0.93,left=0.005,right=0.995,bottom=0.005)
  cmap_skip = 0.2
  top = matplotlib.cm.get_cmap('Blues_r', 128)
  bottom = matplotlib.cm.get_cmap('Oranges', 128)
  newcolors = np.vstack((top(np.linspace(0, 1.0-cmap_skip, 128)), bottom(np.linspace(cmap_skip, 1, 128))))
  cmap = ListedColormap(newcolors, name='OrangeBlue')
  # cmap = matplotlib.cm.get_cmap("PiYG").copy()
  # cmap.set_bad(color="#b0b0b0")
  selected_genus = pd.Series(df.index.get_level_values(0)[:nrows]).values
  selected_metab = pd.Series(df.index.get_level_values(1)[:nrows]).values
  selected_meds = pd.Series(np.concatenate(df["medication"].head(nrows).apply(list).values)).drop_duplicates().values
  components = list(nx.connected_components(graph))
  yoff = 0
  for comp in components:
    selected_genus_subset = list(comp & set(selected_genus))
    selected_metab_subset = list(comp & set(selected_metab))
    selected_meds_subset = list(comp & set(selected_meds))
    height = np.max([len(selected_genus_subset), len(selected_metab_subset), len(selected_meds_subset)])
    pos = {}
    for ti,med in enumerate(selected_meds_subset):
      pos[med] = (1,height - ((1.0+ti)*height/(1.0+len(selected_meds_subset))) + yoff)
      ax.text(pos[med][0], pos[med][1],(prettify_compound_medication_name(med,26,mode="trunc")),ha="center",va="center",fontsize=5.5, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round"),zorder=1)
      ax.text(pos[med][0], pos[med][1],(prettify_compound_medication_name(med,26,mode="trunc")),ha="center",va="center",fontsize=5.5, zorder=2)
    for ti,metab in enumerate(selected_metab_subset):
      if metab not in pos:
        pos[metab] = (0.1,height - ((1.0+ti)*height/(1.0+len(selected_metab_subset))) + yoff)
      ax.text(0.08-0.02, pos[metab][1],metab.split("#")[-1],ha="right",va="center",fontsize=5.5, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round"))
      # for med in selected_meds_subset:
        # if df3.loc[genus,med].ax.plot([pos[genus][0], pos[med][0]],[pos[genus][1], pos[med][1]], c="black")
    for ti,genus in enumerate(selected_genus_subset):
      if genus not in pos:
        pos[genus] = (1.9,height - ((1.0+ti)*height/(1.0+len(selected_genus_subset))) + yoff)
      ax.text(1.92+0.02, pos[genus][1],genus,ha="left",va="center",fontsize=5.5, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round"))
    df3_component = df3[(df3.index.get_level_values(0).isin(selected_metab_subset+selected_genus_subset)) & (df3.index.get_level_values(1).isin(selected_meds_subset))]
    for ri,row in df3_component.reset_index().iterrows():
      xbegin = pos[row[0]][0]
      xend = pos[row[1]][0]
      ybegin = pos[row[0]][1]
      yend = pos[row[1]][1]
      linex = np.array([xbegin, (0.8*xbegin)+(0.2*xend), (0.75*xbegin)+(0.25*xend), (0.2*xbegin)+(0.8*xend), xend])
      liney = np.array([ybegin, (1.0*ybegin)+(0.0*yend), (0.0*ybegin)+(1.0*yend), (0.0*ybegin)+(1.0*yend), yend])
      liney = liney[np.argsort(linex)]
      linex = linex[np.argsort(linex)]
      linex_interp = np.linspace(linex.min(), linex.max(), 100)
      spl = scipy.interpolate.make_interp_spline(linex, liney, k=2)
      liney_interp = spl(linex_interp)
      f_cubic = scipy.interpolate.interp1d(linex, liney, kind='slinear')
      liney_interp = f_cubic(linex_interp)
      liney_interp = scipy.ndimage.gaussian_filter1d(liney_interp, 8.0,mode="nearest")
      # ax.plot(linex, liney,color=cmap(matplotlib.colors.Normalize(vmin=-max_col, vmax=max_col)(df3.loc[(row[0],row[1])].coef.iloc[0])), lw=0.1 + 0.3*(-np.log10((df3.loc[(row[0],row[1])].pval.iloc[0]))), zorder=0)
      ax.plot(linex_interp, liney_interp,color=cmap(matplotlib.colors.Normalize(vmin=-max_col, vmax=max_col)(df3.loc[(row[0],row[1])].coef.iloc[0])), lw=(0.1 + 0.3*(-np.log10((df3.loc[(row[0],row[1])].pval.iloc[0])))).clip(max=4.0), zorder=0, alpha=0.5)
    #   df3.loc[(row[0],row[1])].coef.iloc[0]
    yoff += height
  ax.text(1.00, yoff + 0.5,"Medication",ha="center",va="center",weight="bold",fontsize=6, zorder=2)
  ax.text(-0.25, yoff + 0.5,"Metabolite",ha="center",va="center",weight="bold",fontsize=6, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round"))
  ax.text(2.15, yoff + 0.5,"Genus",ha="center",va="center",weight="bold",fontsize=6, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round"))
  ax.set_title("Medications Associated With Genus-Metabolite Pairs",fontsize=8)
  ax.set_xlim(-1.1,2.8)
  ax.set_ylim(0, yoff + 1.5)
  ax.set_xticks([])
  ax.set_yticks([])
  ax.axis("off")      ## turn off axis borders
  fig.savefig("out/figure_4b.pdf")
  # plt.show()

def make_figure_4c():
  table_genus_metab_correl.to_csv("out/supp_table_genus_metab_correlations.csv")
  table_genus_metab_correl_coefs = table_genus_metab_correl.pivot(index="genus",columns="metab",values="statistic").astype(float)
  table_genus_metab_correl_pvalues = table_genus_metab_correl.pivot(index="genus",columns="metab",values="pvalue").astype(float)
  table_genus_metab_correl_sig = table_genus_metab_correl.pivot(index="genus",columns="metab",values="significant").astype(float)
  cmap = matplotlib.cm.get_cmap("BrBG").copy()
  cmap.set_bad(color="#b0b0b0")
  select_genus = table_genus_metab_correl[table_genus_metab_correl.significant]["genus"].value_counts().head(24)
  select_metab = table_genus_metab_correl[table_genus_metab_correl.significant]["metab"].value_counts().head(35)
  print(select_metab)
  print(select_genus)
  grid = table_genus_metab_correl_coefs.loc[select_genus.index, select_metab.index].astype(float) * table_genus_metab_correl_sig.loc[select_genus.index, select_metab.index].astype(float).replace(0,np.nan)
  ordering_0 = np.array(scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid.fillna(0),method="ward"), no_plot=True, color_threshold=-np.inf)['leaves'])
  ordering_1 = np.array(scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(grid.fillna(0).T,method="ward"), no_plot=True, color_threshold=-np.inf)['leaves'])
  grid = grid.iloc[ordering_0,ordering_1]
  fig = plt.figure(figsize=(3.7,3.3))
  spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, left=0.35,top=0.9,bottom=0.3,right=0.99)
  ax = fig.add_subplot(spec[0,0])
  fig.suptitle("\n".join(wrap("Correlations in Effect Sizes Between Metabolites and Microbial Genera",60)),fontsize=8)
  ax_colorbar = fig.add_axes([0.02, 0.25, 0.30, 0.025])
  # vmax = results.table[(results.table.microbe=="Bacteroides") & (results.table.medication.isin(selected_medications))]["coef"].abs().max()
  # sns.heatmap(results_subset, ax=axes[0], cmap=cmap,center=0,cbar_ax=ax_colorbar,cbar_kws={"orientation":"horizontal"},vmin=-vmax, vmax=vmax)
  hm = sns.heatmap(grid.T,ax=ax,cmap=cmap,center=0,cbar_ax=ax_colorbar,cbar_kws={"orientation":"horizontal"}, linewidth=0.1, linecolor="gray")
  cbar = ax.collections[0].colorbar
  cbar.ax.tick_params(labelsize=6)
  cbar.ax.set_xlabel("Spearman\nCorrelation",fontsize=8)
  # cbar = ax.collections[0].colorbar
  # cbar.ax.tick_params(labelsize=6)
  # cbar.ax.set_ylabel("Spearman Correlation",fontsize=6)
  ax.set_yticks(np.arange(len(grid.columns))+0.5)
  ax.set_xticks(np.arange(len(grid.index))+0.5)
  ax.set_yticklabels(grid.columns.str.split("#").str[-1],fontsize=6)
  ax.set_xticklabels(grid.index,fontsize=6)
  fig.savefig("out/figure_4c.pdf")
  # plt.show()

import sklearn.neighbors
def make_supplementary_figure_metabolite_genus_correlations():
  v1 = "alphadiversity"
  v2 = "scfa#butyrate"
  fig,ax=plt.subplots(1,1,figsize=(7.5,3))
  ax.scatter(matrix_genus.loc[:,v1], matrix_metab.loc[:,v2],s=3,c="gray")
  scatterpoints = pd.concat([matrix_genus.loc[:,v1], matrix_metab.loc[:,v2]],axis=1)
  mu = scatterpoints.mean(axis=0)
  Sigma = np.cov(scatterpoints, rowvar=False)
  Sigma_inv = np.linalg.inv(Sigma)
  df_text = pd.concat([matrix_genus.loc[:,v1], matrix_metab.loc[:,v2]],axis=1)
  df_text[["orig0","orig1"]] = df_text[[v1,v2]]
  df_text["text"] = df_text.index.map(lambda x: prettify_compound_medication_name(x,length=20,mode="trunc"))
  points = df_text[[v1,v2]]
  aniso_scale=np.array([1,1.5])
  nn = sklearn.neighbors.NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(points)
  dists, idxs = nn.kneighbors(points)
  adj_dists = np.linalg.norm((points.iloc[idxs[:,1]]*aniso_scale) - (points.values*aniso_scale),axis=1)
  df_text = df_text.iloc[::20]
  # df_text = df_text[adj_dists>np.quantile(adj_dists,.92)]
  for i in range(200):
    points = df_text[[v1,v2]]
    nn = sklearn.neighbors.NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(points)
    dists, idxs = nn.kneighbors(points)
    adj_dists = np.linalg.norm((points.iloc[idxs[:,1]]*aniso_scale) - (points.values*aniso_scale),axis=1, ord=1)
    speed = 0.6
    vector = (points.iloc[idxs[:,1]] - points.values) * (np.abs(adj_dists)<1).reshape(-1,1) * -speed
    diff = points - mu
    v = Sigma_inv @ diff.T
    vector.iloc[:,:] += ((1.0/(0.15+v.T.values))*0.01) * (np.linalg.norm(diff,axis=1)<20).reshape(-1,1)
    # vector += points.values * (np.linalg.norm(points,axis=1)<0.8).reshape(-1,1) * speed
    vector += (points.values[:,0]<-1.2).reshape(-1,1) @ np.array([speed,0.0]).reshape(1,-1)
    vector += (points.values[:,0]>1.2).reshape(-1,1) @ np.array([-speed,0.0]).reshape(1,-1)
    vector += (points.values[:,1]<-2).reshape(-1,1) @ np.array([0.0,speed]).reshape(1,-1)
    vector += (points.values[:,1]>1.5).reshape(-1,1) @ np.array([0.0,-speed]).reshape(1,-1)
    df_text[[v1,v2]] += vector.values
  ix = 0
  for ti,text in df_text.reset_index(drop=True).iterrows():
    if np.linalg.norm(np.array([text["orig0"],text["orig1"]]) - np.array([text[0], text[1]])) > 0.1:
      ax.annotate(text["text"].replace(", ","\n"), xy=(text["orig0"],text["orig1"]), xytext=(text[0], text[1]), xycoords="data", textcoords="data", arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.05" if ix%2 == 0 else "arc3,rad=-.05",color="black",alpha=1), fontsize=6, ha="center", va="center")
    else:
      ax.annotate(text["text"].replace(", ","\n"), xy=(text["orig0"],text["orig1"]), xytext=(text[0], text[1]), xycoords="data", textcoords="data", fontsize=6, ha="center", va="center")
    ix = ix+1
  ax.set_xlabel("Effect size with "+v1.replace("alphadiversity", "alpha diversity"),fontsize=7)
  ax.set_ylabel("Effect size with "+v2.split("#")[-1],fontsize=7)
  ax.set_xlim(-1.7, 1.7)
  ax.set_xticklabels(ax.get_xticklabels(),fontsize=6)
  ax.set_yticklabels(ax.get_yticklabels(),fontsize=6)
  ax.text(0.05,0.05,"Spearman correl: %.2f"%table_genus_metab_correl.set_index(["genus","metab"]).loc[(v1,v2)]["statistic"], transform = ax.transAxes, fontsize=6, fontweight="bold",ha="left",va="bottom")
  ax.set_title("Effect size correlation between butyrate and microbial diversity",fontsize=8)
  fig.savefig("out/figure_supp_effect_size_correlations_microbe_metabolite.pdf")
  fig.savefig("out/figure_supp_effect_size_correlations_microbe_metabolite.png", dpi=400)
  plt.close()
  # plt.show()

def make_supplementary_figure_invitro_overlap():
  table_invitro_subset = table_invitro_reference
  table_invitro_subset.columns = table_invitro_subset.columns.str.split("(").str[0].str.strip()
  table_invitro_subset = table_invitro_subset.set_index("chemical_name").loc[correspondences.maier_medication.values]
  table_invitro_subset = table_invitro_subset.loc[:, table_microbes_relabund.index.intersection(table_invitro_subset.columns)]
  table_coefs_0_microbe = results.table_coefs.loc[:,results.table_coefs.columns.str.contains("_species")]
  table_pvalues_0_microbe = results.table_pvalues.loc[:,results.table_pvalues.columns.str.contains("_species")]
  table_pvalues_0_microbe = table_pvalues_0_microbe.loc[~table_pvalues_0_microbe.index.str.contains("/Dose"),:]
  table_coefs_0_microbe = table_coefs_0_microbe.loc[~table_coefs_0_microbe.index.str.contains("/Dose"),:]
  table_pvalues_0_microbe = table_pvalues_0_microbe.mask((table_coefs_0_microbe.values>=0), 1.0)
  table_pvalues_0_microbe = table_pvalues_0_microbe.loc[:,~table_pvalues_0_microbe.columns.str.contains("dose")]
  table_pvalues_0_microbe.columns = table_pvalues_0_microbe.columns.str.replace("_0_species","").str.replace("_dose","").str.replace("_binary","")
  table_pvalues_0_microbe.index = table_pvalues_0_microbe.index.str.split("_").str[0].str.lower()
  maier_medications_intersection = table_pvalues_0_microbe.index.intersection([x.lower() for x in table_invitro_subset.index])
  maier_microbes_intersection = table_pvalues_0_microbe.columns.str.lower().intersection([x.lower() for x in table_invitro_subset.columns])
  table_pvalues_0_microbe.columns = table_pvalues_0_microbe.columns.str.lower()
  table_invitro_subset.columns = table_invitro_subset.columns.str.lower()
  table_invitro_subset.index = table_invitro_subset.index.str.lower()
  table_pvalues_0_microbe = table_pvalues_0_microbe.loc[maier_medications_intersection,maier_microbes_intersection]
  table_invitro_subset = table_invitro_subset.loc[maier_medications_intersection,maier_microbes_intersection]
  table_invitro_subset = table_invitro_subset.loc[:,~table_invitro_subset.columns.duplicated()]
  table_pvalues_0_microbe = table_pvalues_0_microbe.groupby(table_pvalues_0_microbe.columns,axis=1).min()
  table_pvalues_0_microbe = table_pvalues_0_microbe.groupby(table_pvalues_0_microbe.index).min().sort_index().reindex(sorted(table_pvalues_0_microbe.columns.unique()), axis=1)
  table_invitro_subset = table_invitro_subset.groupby(table_invitro_subset.index).min().sort_index().reindex(sorted(table_invitro_subset.columns.unique()), axis=1)
  threshs = 10.0**(np.arange(-3,-0.8,0.05))
  fdr_thresh = 0.1
  fisher_pvals = []
  for thresh in threshs:
    df_overlap = pd.DataFrame(index=["overall_sig_dfi","overall_nsig_dfi"], columns=["overall_sig_maier","overall_nsig_maier"])
    df_overlap.loc["overall_sig_dfi","overall_sig_maier"] = (((table_invitro_subset<thresh)) & (table_pvalues_0_microbe<thresh)).sum().sum()
    df_overlap.loc["overall_sig_dfi","overall_nsig_maier"] = (((table_invitro_subset>=thresh)) & (table_pvalues_0_microbe<thresh)).sum().sum()
    df_overlap.loc["overall_nsig_dfi","overall_sig_maier"] = (((table_invitro_subset<thresh)) & (table_pvalues_0_microbe>=thresh)).sum().sum()
    df_overlap.loc["overall_nsig_dfi","overall_nsig_maier"] = (((table_invitro_subset>=thresh)) & (table_pvalues_0_microbe>=thresh)).sum().sum()
    agreement = (((table_invitro_subset<thresh)) == (table_pvalues_0_microbe<thresh)).sum().sum()
    print(thresh, scipy.stats.fisher_exact(df_overlap), agreement)
    fisher_pvals += [scipy.stats.fisher_exact(df_overlap).pvalue]
  # fig,ax=plt.subplots(1,1,tight_layout=True)
  # ax.scatter(threshs,fisher_pvals, s=2, color="blue")
  # ax.plot(threshs,fisher_pvals, color="black")
  # plt.show()
  thresh = 0.1
  df_overlap = pd.DataFrame(index=["overall_sig_dfi","overall_nsig_dfi"], columns=["overall_sig_maier","overall_nsig_maier"])
  df_overlap.loc["overall_sig_dfi","overall_sig_maier"] = (((table_invitro_subset<thresh)) & (table_pvalues_0_microbe<thresh)).sum().sum()
  df_overlap.loc["overall_sig_dfi","overall_nsig_maier"] = (((table_invitro_subset>=thresh)) & (table_pvalues_0_microbe<thresh)).sum().sum()
  df_overlap.loc["overall_nsig_dfi","overall_sig_maier"] = (((table_invitro_subset<thresh)) & (table_pvalues_0_microbe>=thresh)).sum().sum()
  df_overlap.loc["overall_nsig_dfi","overall_nsig_maier"] = (((table_invitro_subset>=thresh)) & (table_pvalues_0_microbe>=thresh)).sum().sum()
  print(thresh, scipy.stats.fisher_exact(df_overlap))
  fisher_pvals += [scipy.stats.fisher_exact(df_overlap).pvalue]
  fig,ax = plt.subplots(1,1,tight_layout=True,figsize=(4,5))
  ax.set_title("Overlap With 2-10 Day Associations\n(Fisher's exact test p=%.5f)"%(fisher_pvals[np.sum(threshs<0.01)]),fontsize=8)
  grid_agree = ((table_invitro_subset<thresh) == (table_pvalues_0_microbe<thresh)).astype("float")
  grid_sig = (table_invitro_subset<thresh).astype("float")
  grid_nsig = (table_invitro_subset>=thresh).astype("float")
  grid = (grid_agree*grid_sig) - (grid_agree*grid_nsig)
  cmap = matplotlib.cm.get_cmap("tab20c").copy()
  cmap.set_bad(color="#909090")
  cmap.set_under(color="darkturquoise")
  cmap.set_over(color="orange")
  grid = grid.replace(0,np.nan)
  sns.heatmap(grid, ax=ax, cbar=False, cmap=cmap, linewidths=0.5, linecolor='#606060',vmin=-0.5,vmax=0.5)
  # sns.heatmap((table_invitro_subset<thresh).astype("float"), ax=axes[1], cbar=False, cmap="Blues")
  # sns.heatmap((table_pvalues_0_microbe<thresh).astype("float"), ax=axes[2], cbar=False, cmap="Blues")
  ax.set_xticks(np.arange(len(table_invitro_subset.columns))+0.5)
  ax.set_xticklabels(table_invitro_subset.columns.str.capitalize(), fontsize=6, rotation=45, ha="right", rotation_mode="anchor")
  # ax.set_title("\n".join(wrap("Green (both significant); Red (both not significant); White (disagreement)",40)), fontsize=8)
  ax.set_yticks(np.arange(len(table_invitro_subset.index))+0.5)
  ax.set_yticklabels(table_invitro_subset.index.str.capitalize(), fontsize=6)
  fig.savefig("out/supplementary_figure_invitro_overlap.pdf")
  fig.savefig("out/supplementary_figure_invitro_overlap.png",dpi=400)
  plt.close()


def flatten(xss):
  return [x for xs in xss for x in xs]


# df1.apply(lambda xs: np.concatenate([x for x in xs]))
def make_figure_4d():
  # for each medication, get the list of all genera, metabolites and pathways it is associated with, along with their coefficients.
  # df = pd.concat([
  #   results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype=="genus")][["medication","microbe","coef"]].groupby("medication").agg(list), 
  #   results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype=="metab")][["medication","microbe","coef"]].groupby("medication").agg(list), 
  #   results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype=="pathway")][["medication","microbe","coef"]].groupby("medication").agg(list)], axis=1)
  # df.columns = ["genus","genus_coef","metab","metab_coef","pathway","pathway_coef"]
  # # df["metab_coef"][df["metab_coef"].isna()]=[]
  # # df["genus_coef"][df["genus_coef"].isna()]=np.array([])
  # # df["pathway_coef"][df["pathway_coef"].isna()]=np.array([])
  # # df["minp"] = df[["genus_pvalue","metab_pvalue"]].apply(lambda row: (row["genus_pvalue"] | row["metab_pvalue"]),axis=1).apply(lambda row: min(row))
  # # get list of 
  # df2 = df.explode(["genus","genus_coef"]).explode(["metab","metab_coef"]).explode(["pathway","pathway_coef"])[["metab","metab_coef","pathway","pathway_coef"]].drop_duplicates()
  # df3 = df.explode(["genus","genus_coef"]).explode(["metab","metab_coef"]).explode(["pathway","pathway_coef"])[["metab","metab_coef","genus","genus_coef"]].drop_duplicates()
  # df = pd.concat([pd.crosstab(df2["metab"],df2["pathway"]),pd.crosstab(df3["metab"],df3["genus"])],axis=1)
  metab_connectivity_1 = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["metab"]))][["medication","microbe"]].drop_duplicates()["microbe"].value_counts().head(16)
  df1 = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["metab"]))][["medication","microbe"]].drop_duplicates().rename({"microbe":"metab"},axis=1)
  df2 = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["genus","pathway"]))][["medication","microbe"]].drop_duplicates().rename({"microbe":"microbe"},axis=1)
  df1 = df1.groupby("metab").apply(lambda df: list(df.medication.drop_duplicates()))
  df2 = df2.groupby("medication").apply(lambda x: list(x.microbe.drop_duplicates()))
  metab_connectivity_2 = df1.apply(lambda xs: np.unique(flatten([df2.loc[x] for x in xs]))).str.len().sort_values(ascending=False).head(16)
  metab_connectivity_all_1 = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["metab"]))][["medication","microbe"]].drop_duplicates()["microbe"].value_counts()
  metab_connectivity_all_2 = df1.apply(lambda xs: np.unique(flatten([df2.loc[x] for x in xs]))).str.len().sort_values(ascending=False)
  grid_csv = pd.concat([metab_connectivity_all_1, metab_connectivity_all_2],axis=1)
  grid_csv.columns = ["degree1", "degree2"]
  grid_csv.to_csv("out/supp_table_metab_connectivity.csv")
  fig,axes = plt.subplots(1,2,constrained_layout=True, figsize=(3.8,2))
  fig.suptitle("Number of Graph Connections for Metabolites",fontsize=8)
  axes[0].set_title("Degree 1",fontsize=6)
  axes[1].set_title("Degree 2",fontsize=6)
  axes[0].set_ylabel("Metabolite",fontsize=6)
  axes[0].set_xlabel("# Medications",fontsize=6)
  axes[1].set_xlabel("# Genus/Pathway  ",fontsize=6)
  axes[0].barh(np.arange(len(metab_connectivity_1))[::-1], metab_connectivity_1.values, 0.8, color="gray")
  axes[0].set_yticks(np.arange(len(metab_connectivity_1))[::-1])
  axes[0].set_yticklabels(metab_connectivity_1.index.str.split("#").str[-1],fontsize=6)
  axes[0].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4,steps=[1, 2, 4, 5, 10]))
  # axes[0].xaxis.set_minor_locator(matplotlib.ticker.MaxNLocator(20,steps=[1, 2, 4, 5, 10]))
  axes[0].set_xticklabels(axes[0].get_xticklabels(),fontsize=6)
  axes[1].barh(np.arange(len(metab_connectivity_2))[::-1], metab_connectivity_2.values, 0.8, color="gray")
  axes[1].set_yticks(np.arange(len(metab_connectivity_2))[::-1])
  axes[1].set_yticklabels(metab_connectivity_2.index.str.split("#").str[-1],fontsize=6)
  axes[1].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4,steps=[1, 2, 4, 6]))
  # axes[1].xaxis.set_minor_locator(matplotlib.ticker.MaxNLocator(20,steps=[1, 2, 4, 5, 10]))
  axes[1].set_xticklabels(axes[1].get_xticklabels(),fontsize=6)
  fig.savefig("out/figure_4d.pdf")
  plt.close()

def make_figure_4e():
  df = pd.concat([
    results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype=="genus")][["medication","microbe","coef"]].groupby("medication").agg(list), 
    results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype=="metab")][["medication","microbe","coef"]].groupby("medication").agg(list), 
    results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype=="pathway")][["medication","microbe","coef"]].groupby("medication").agg(list)], axis=1)
  df.columns = ["genus","genus_coef","metab","metab_coef","pathway","pathway_coef"]
  # df["metab_coef"][df["metab_coef"].isna()]=[]
  # df["genus_coef"][df["genus_coef"].isna()]=np.array([])
  # df["pathway_coef"][df["pathway_coef"].isna()]=np.array([])
  # df["minp"] = df[["genus_pvalue","metab_pvalue"]].apply(lambda row: (row["genus_pvalue"] | row["metab_pvalue"]),axis=1).apply(lambda row: min(row))
  df2 = df.explode(["genus","genus_coef"]).explode(["metab","metab_coef"]).explode(["pathway","pathway_coef"])[["metab","metab_coef","pathway","pathway_coef"]].drop_duplicates()
  df3 = df.explode(["genus","genus_coef"]).explode(["metab","metab_coef"]).explode(["pathway","pathway_coef"])[["metab","metab_coef","genus","genus_coef"]].drop_duplicates()
  df = pd.concat([pd.crosstab(df2["metab"],df2["pathway"]),pd.crosstab(df3["metab"],df3["genus"])],axis=1)
  metab_connectivity_1 = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["metab"]))][["medication","microbe"]].drop_duplicates()["microbe"].value_counts().head(20)
  metab_connectivity_2 = df.sum(axis=1).sort_values(ascending=False).head(19)
  adj = pd.concat([df,df.T]).fillna(0).astype(bool)
  adj_num = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["metab","genus","pathway"]))].pivot_table(index="medication",columns="microbe",values="coef", aggfunc = lambda xs: max(xs, key=abs))
  adj_num = pd.concat([adj_num,adj_num.T])
  adj_pval = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["metab","genus","pathway"]))].pivot_table(index="medication",columns="microbe",values="pvalue", aggfunc = lambda xs: min(xs, key=abs))
  adj_pval = pd.concat([adj_pval,adj_pval.T])
  selected_metabs = list(metab_connectivity_2.head(3).index)
  selected_meds = adj_num.loc[selected_metabs].abs().mean().sort_values(ascending=False).head(12).index
  all_metab = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["metab"]))]["microbe"].unique()
  all_genus = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["genus"]))]["microbe"].unique()
  all_pathway = results.table[(results.table.significant) & (results.table.modeltype=="binary") & (results.table.sampletype.isin(["pathway"]))]["microbe"].unique()
  selected_t2_1 = adj_num.loc[selected_meds,all_metab].abs().mean().sort_values(ascending=False).head(15).index
  selected_t2_2 = adj_num.loc[selected_meds,all_genus].abs().mean().sort_values(ascending=False).head(16).index
  selected_t2_3 = adj_num.loc[selected_meds,all_pathway].abs().mean().sort_values(ascending=False).head(0).index
  selected_t2 = list(selected_t2_2) + list(selected_t2_3)
  # selected_t2 = list(selected_t2_2)
  adj_sub = adj_num.loc[selected_metabs + list(selected_meds) + list(selected_t2),selected_metabs + list(selected_meds) + list(selected_t2)]
  adj_sub_pval = adj_pval.loc[selected_metabs + list(selected_meds) + list(selected_t2),selected_metabs + list(selected_meds) + list(selected_t2)]
  graph = nx.Graph(adj_sub)
  pos = {}
  fig,ax=plt.subplots(1,1,figsize=(3.45,2))
  fig.subplots_adjust(top=0.93,left=0.005,right=0.995,bottom=0.005)
  fig.suptitle("Interaction Graph for Subset of Metabolites",fontsize=8)
  for i in range(len(selected_metabs)):
    pos[selected_metabs[i]] = (-0.5,(i+1) / (len(selected_metabs)+1))
    ax.text(pos[selected_metabs[i]][0], pos[selected_metabs[i]][1],selected_metabs[i].split("#")[1],ha="right",va="center",fontsize=5.5, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round",alpha=0.7))
  for i in range(len(selected_meds)):
    pos[selected_meds[i]] = (0.6,(i+1) / (len(selected_meds)+1))
    ax.text(pos[selected_meds[i]][0], pos[selected_meds[i]][1],"\n".join(wrap(prettify_compound_medication_name(selected_meds[i]),32)),ha="center",va="center",fontsize=5.5, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round",alpha=0.7))
  for i in range(len(selected_t2)):
    pos[selected_t2[i]] = (2,(i+1) / (len(selected_t2)+1))
    label = selected_t2[i]
    if ":" in label:
      label = label.split(":")[1].strip().replace("&alpha;", "α").replace("&beta;", "β").replace("&delta;", "Δ")
      label = "\n".join(wrap(label,35))
    ax.text(pos[selected_t2[i]][0], pos[selected_t2[i]][1],label,ha="left",va="center",fontsize=5.5, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round",alpha=0.7))
  df = adj_sub.stack().reset_index()
  df["pvalue"] = adj_sub_pval.stack().reset_index()[0]
  df["linewidth"] = (-np.log10(df["pvalue"]))
  df["linewidth"] = (0.4 + ((df["linewidth"]-(df["linewidth"].min()))*1.0)).clip(upper=4.0)
  df = df[df["level_0"].isin(selected_meds)]
  max_col = df[0].abs().max()
  cmap_skip = 0.2
  top = matplotlib.cm.get_cmap('Blues_r', 128)
  bottom = matplotlib.cm.get_cmap('Oranges', 128)
  newcolors = np.vstack((top(np.linspace(0, 1.0-cmap_skip, 128)), bottom(np.linspace(cmap_skip, 1, 128))))
  cmap = ListedColormap(newcolors, name='OrangeBlue')
  # cmap = matplotlib.cm.get_cmap("PiYG").copy()
  cmap.set_bad(color="#b0b0b0")
  for ri,row in df.iterrows():
    textpos = [0,0]
    if row.iloc[1] in selected_metabs:
      xbegin = pos[row.iloc[0]][0]-0.3
      xend = pos[row.iloc[1]][0]
      ybegin = pos[row.iloc[0]][1]
      yend = pos[row.iloc[1]][1]
      # textpos = [pos[row.iloc[0]][0]-0.4, pos[row.iloc[1]][0]]
      # _ = ax.plot([pos[row.iloc[0]][0]-0.4, pos[row.iloc[1]][0]], [pos[row.iloc[0]][1], pos[row.iloc[1]][1]], color = cmap(matplotlib.colors.Normalize(vmin=-max_col, vmax=max_col)(row[0])))
    else:
      xbegin = pos[row.iloc[0]][0]+0.6
      xend = pos[row.iloc[1]][0]
      ybegin = pos[row.iloc[0]][1]
      yend = pos[row.iloc[1]][1]
      # textpos = [pos[row.iloc[0]][0]+0.4, pos[row.iloc[1]][0]]
      # _ = ax.plot([pos[row.iloc[0]][0]+0.4, pos[row.iloc[1]][0]], [pos[row.iloc[0]][1], pos[row.iloc[1]][1]], color = cmap(matplotlib.colors.Normalize(vmin=-max_col, vmax=max_col)(row[0])))
    # linex = np.array([xbegin, (0.9*xbegin)+(0.1*xend), (0.9*xbegin)+(0.1*xend), xend])
    # liney = np.array([ybegin, (0.0*ybegin)+(1.0*yend), (0.0*ybegin)+(1.0*yend), yend])
    linex = np.array([xbegin, (0.8*xbegin)+(0.2*xend), (0.2*xbegin)+(0.8*xend), xend])
    liney = np.array([ybegin, ybegin, yend, yend])
    liney = liney[np.argsort(linex)]
    linex = linex[np.argsort(linex)]
    linex_interp = np.linspace(linex.min(), linex.max(), 100)
    spl = scipy.interpolate.make_interp_spline(linex, liney, k=2)
    liney_interp = spl(linex_interp)
    f_cubic = scipy.interpolate.interp1d(linex, liney, kind='slinear')
    liney_interp = f_cubic(linex_interp)
    liney_interp = scipy.ndimage.gaussian_filter1d(liney_interp, 8.0,mode="nearest")
    # ax.plot(linex, liney,color=cmap(matplotlib.colors.Normalize(vmin=-max_col, vmax=max_col)(df3.loc[(row[0],row[1])].coef.iloc[0])), lw=0.1 + 0.3*(-np.log10((df3.loc[(row[0],row[1])].pval.iloc[0]))), zorder=0)
    ax.plot(linex_interp, liney_interp,color=cmap(matplotlib.colors.Normalize(vmin=-max_col, vmax=max_col)(row[0])), zorder=0, alpha=0.5, lw=row["linewidth"])
  # _ = nx.draw_networkx_nodes(graph,pos=pos,node_shape="s",node_color=[colors[x] for x in adj_sub.index], ax=ax, node_size=20)
  # _ = nx.draw_networkx_edges(graph,pos=pos, ax=ax, edgelist=graph.edges, edge_color=[adj.loc[e[0],e[1]] for e in graph.edges],edge_cmap=cmap, edge_vmin=-adj.abs().max().max(), edge_vmax=adj.abs().max().max())
  yoff = np.max([pos[p][1] for p in pos])
  ax.text(0.6, yoff + 0.06,"Medication",ha="center",va="center",weight="bold",fontsize=6, zorder=2)
  ax.text(2.3, yoff + 0.06,"Genus",ha="center",va="center",weight="bold",fontsize=6, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round"))
  ax.text(-1.0, yoff + 0.06,"Metabolite",ha="center",va="center",weight="bold",fontsize=6, bbox=dict(facecolor='white', edgecolor='none', boxstyle="round"))
  ax.set_ylim(ax.get_ylim()[0], yoff+0.15)
  ax.set_xlim(-1.8,3.35)
  ax.set_xticks([])
  ax.set_yticks([])
  ax.axis("off")      ## turn off axis borders
  # _ = nx.draw_networkx_labels(graph, pos=pos, ax=ax, font_size=6)
  fig.savefig("out/figure_4e.pdf")
  plt.close()

def make_spplementary_figure_effect_size_correlations_durations():
  sig01 = results.table[(results.table.duration.isin([0,1])) & (results.table.modeltype=="binary") & (results.table.significant) & (results.table.sampletype=="genus")].groupby(["medication","microbe"]).filter(lambda df: len(df)==2).groupby(["medication","microbe"]).apply(lambda df: pd.Series(df.sort_values("duration").coef.values))
  sig02 = results.table[(results.table.duration.isin([0,2])) & (results.table.modeltype=="binary") & (results.table.significant) & (results.table.sampletype=="genus")].groupby(["medication","microbe"]).filter(lambda df: len(df)==2).groupby(["medication","microbe"]).apply(lambda df: pd.Series(df.sort_values("duration").coef.values))
  coefsd01 = sig01
  coefsd02 = sig02
  coefsd01_annot = [abbr(prettify_compound_medication_name(x["medication"]),30)+", "+x["microbe"] for ix,x in sig01.index.to_frame().iterrows()]
  coefsd02_annot = [abbr(prettify_compound_medication_name(x["medication"]),30)+", "+x["microbe"] for ix,x in sig02.index.to_frame().iterrows()]
  # coefsd01 = [results.table[(results.table["medication"]==x.split("*")[0]) & (results.table["microbe"]==x.split("*")[1]) & (results.table["duration"].isin([0,1])) & (results.table["modeltype"]=="binary") & (~results.table["medication_dose"])].sort_values("duration") for x in list(set(sig0).intersection(sig1))]
  # coefsd02 = [results.table[(results.table["medication"]==x.split("*")[0]) & (results.table["microbe"]==x.split("*")[1]) & (results.table["duration"].isin([0,2])) & (results.table["modeltype"]=="binary") & (~results.table["medication_dose"])].sort_values("duration") for x in list(set(sig0).intersection(sig2))]
  # coefsd01_annot = [abbr(prettify_compound_medication_name(x.iloc[0]["medication"]),30)+", "+x.iloc[0]["microbe"] for x in coefsd01]
  # coefsd02_annot = [abbr(prettify_compound_medication_name(x.iloc[0]["medication"]),30)+", "+x.iloc[0]["microbe"] for x in coefsd02]
  # coefsd01 = [x["coef"].values for x in coefsd01]
  # coefsd02 = [x["coef"].values for x in coefsd02]
  data = pd.DataFrame({"x":[c[0] for c in coefsd01.values],"y":[c[1] for c in coefsd01.values]})
  model1 = sm.RLM(data["y"], sm.add_constant(data["x"]), M=sm.robust.norms.HuberT(1))
  results1 = model1.fit()
  data = pd.DataFrame({"x":[c[0] for c in coefsd02.values],"y":[c[1] for c in coefsd02.values]})
  model2 = sm.RLM(data["y"], sm.add_constant(data["x"]), M=sm.robust.norms.HuberT(1))
  results2 = model2.fit()
  coefsd01 = coefsd01.values
  coefsd02 = coefsd02.values
  df_text_0 = pd.DataFrame(coefsd01)[[0,1]]
  df_text_0["text"] = coefsd01_annot
  df_text_0[["orig0","orig1"]] = df_text_0[[0,1]]
  points = df_text_0[[0,1]]
  aniso_scale=np.array([1,9])
  nn = sklearn.neighbors.NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(points)
  dists, idxs = nn.kneighbors(points)
  adj_dists = np.linalg.norm((points.iloc[idxs[:,1]]*aniso_scale) - (points.values*aniso_scale),axis=1)
  df_text_0 = df_text_0[adj_dists>np.quantile(adj_dists,.92)]
  for i in range(200):
    points = df_text_0[[0,1]]
    nn = sklearn.neighbors.NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(points)
    dists, idxs = nn.kneighbors(points)
    adj_dists = np.linalg.norm((points.iloc[idxs[:,1]]*aniso_scale) - (points.values*aniso_scale),axis=1, ord=1)
    speed = 0.6
    vector = (points.iloc[idxs[:,1]] - points.values) * (np.abs(adj_dists)<10).reshape(-1,1) * -speed
    vector += points.values * (np.linalg.norm(points,axis=1)<0.8).reshape(-1,1) * speed
    vector += (points.values[:,0]<-2.5).reshape(-1,1) @ np.array([speed,0.0]).reshape(1,-1)
    vector += (points.values[:,0]>2.5).reshape(-1,1) @ np.array([-speed,0.0]).reshape(1,-1)
    vector += (points.values[:,1]<-2.5).reshape(-1,1) @ np.array([0.0,speed]).reshape(1,-1)
    vector += (points.values[:,1]>2.5).reshape(-1,1) @ np.array([0.0,-speed]).reshape(1,-1)
    df_text_0[[0,1]] += vector.values
  df_text_1 = pd.DataFrame(coefsd02)[[0,1]]
  df_text_1["text"] = coefsd02_annot
  df_text_1[["orig0","orig1"]] = df_text_1[[0,1]]
  points = df_text_1[[0,1]]
  aniso_scale=np.array([1,9])
  nn = sklearn.neighbors.NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(points)
  dists, idxs = nn.kneighbors(points)
  adj_dists = np.linalg.norm((points.iloc[idxs[:,1]]*aniso_scale) - (points.values*aniso_scale),axis=1)
  df_text_1 = df_text_1[adj_dists>np.quantile(adj_dists,.70)]
  for i in range(200):
    points = df_text_1[[0,1]]
    nn = sklearn.neighbors.NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(points)
    dists, idxs = nn.kneighbors(points)
    adj_dists = np.linalg.norm((points.iloc[idxs[:,1]]*aniso_scale) - (points.values*aniso_scale),axis=1, ord=1)
    speed = 0.6
    vector = (points.iloc[idxs[:,1]] - points.values) * (np.abs(adj_dists)<10).reshape(-1,1) * -speed
    vector += points.values * (np.linalg.norm(points,axis=1)<0.8).reshape(-1,1) * speed
    vector += (points.values[:,0]<-2.5).reshape(-1,1) @ np.array([speed,0.0]).reshape(1,-1)
    vector += (points.values[:,0]>2.5).reshape(-1,1) @ np.array([-speed,0.0]).reshape(1,-1)
    vector += (points.values[:,1]<-2.5).reshape(-1,1) @ np.array([0.0,speed]).reshape(1,-1)
    vector += (points.values[:,1]>2.5).reshape(-1,1) @ np.array([0.0,-speed]).reshape(1,-1)
    df_text_1[[0,1]] += vector.values
  fig, axes = plt.subplots(1,2,constrained_layout=True,figsize=(7.5,3.75))
  fig.suptitle("Correlation Between Significant Effect Sizes Across Time Intervals",fontsize=8)
  minx = np.concatenate([[c[0] for c in coefsd01], [c[0] for c in coefsd02]]).min()
  maxx = np.concatenate([[c[0] for c in coefsd01], [c[0] for c in coefsd02]]).max()
  miny = np.concatenate([[c[1] for c in coefsd01], [c[1] for c in coefsd02]]).min()
  maxy = np.concatenate([[c[1] for c in coefsd01], [c[1] for c in coefsd02]]).max()
  minv = min(minx,miny)-1.0
  maxv = max(maxx,maxy)+1.0
  minv = min(minv, -maxv)
  maxv = max(maxv, -minv)
  axes[0].set_aspect("equal")
  axes[1].set_aspect("equal")
  # axes[0].scatter(coefs_duration0[(coefs_duration0!=0) & (coefs_duration1!=0)].values, coefs_duration1[(coefs_duration0!=0) & (coefs_duration1!=0)].values, s=5, c="black")
  axes[0].scatter([c[0] for c in coefsd01], [c[1] for c in coefsd01], s=5, c="black")
  axes[0].set_xlabel("Effect Size at %d-%d Days"%(durations[0][0], durations[0][1]),fontsize=8)
  axes[0].set_ylabel("Effect Size at %d-%d Days"%(durations[1][0], durations[1][1]),fontsize=8)
  axes[1].set_aspect("equal")
  # axes[1].scatter(coefs_duration0[(coefs_duration0!=0) & (coefs_duration2!=0)].values, coefs_duration2[(coefs_duration0!=0) & (coefs_duration2!=0)].values, s=5, c="black")
  axes[1].scatter([c[0] for c in coefsd02], [c[1] for c in coefsd02], s=5, c="black")
  axes[1].set_xlabel("Effect Size at %d-%d Days"%(durations[0][0], durations[0][1]),fontsize=8)
  axes[1].set_ylabel("Effect Size at %d-%d Days"%(durations[2][0], durations[2][1]),fontsize=8)
  axes[0].set_xlim(minv, maxv)
  axes[1].set_xlim(minv, maxv)
  axes[0].set_ylim(minv, maxv)
  axes[1].set_ylim(minv, maxv)
  axes[0].axline((0,results1.params[0]), slope=results1.params[1], linestyle="--", color="gray", lw=1, zorder=0)
  axes[0].text(maxv*0.8, minv+1.1, "slope: %.2f"%(results1.params[1]),ha="right",va="bottom",fontweight="bold",fontsize=8)
  axes[1].text(maxv*0.8, minv+1.1, "slope: %.2f"%(results2.params[1]),ha="right",va="bottom",fontweight="bold",fontsize=8)
  axes[1].axline((0,results2.params[0]), slope=results2.params[1], linestyle="--", color="gray", lw=1, zorder=0)
  axes[0].grid(color="gray", alpha=0.2, linewidth=0.5, linestyle="--")
  axes[1].grid(color="gray", alpha=0.2, linewidth=0.5, linestyle="--")
  axes[0].tick_params(axis='both', which='major', labelsize=6)
  axes[0].tick_params(axis='both', which='minor', labelsize=6)
  axes[1].tick_params(axis='both', which='major', labelsize=6)
  axes[1].tick_params(axis='both', which='minor', labelsize=6)
  for ti,text in df_text_0.reset_index(drop=True).iterrows():
    if np.linalg.norm(np.array([text["orig0"],text["orig1"]]) - np.array([text[0], text[1]])) > 0.1:
      axes[0].annotate(text["text"].replace(", ","\n"), xy=(text["orig0"],text["orig1"]), xytext=(text[0], text[1]), xycoords="data", textcoords="data", arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2",color="gray",alpha=0.5), fontsize=6, ha="center", va="center")
    else:
      axes[0].annotate(text["text"].replace(", ","\n"), xy=(text["orig0"],text["orig1"]), xytext=(text[0], text[1]), xycoords="data", textcoords="data", fontsize=6, ha="center", va="center")
  for ti,text in df_text_1.reset_index(drop=True).iterrows():
    if np.linalg.norm(np.array([text["orig0"],text["orig1"]]) - np.array([text[0], text[1]])) > 0.1:
      axes[1].annotate(text["text"].replace(", ","\n"), xy=(text["orig0"],text["orig1"]), xytext=(text[0], text[1]), xycoords="data", textcoords="data", arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2",color="gray",alpha=0.5), fontsize=6, ha="center", va="center")
    else:
      axes[1].annotate(text["text"].replace(", ","\n"), xy=(text["orig0"],text["orig1"]), xytext=(text[0], text[1]), xycoords="data", textcoords="data", fontsize=6, ha="center", va="center")
  fig.savefig("out/figure_supp_effect_size_correlations_durations.pdf")
  fig.savefig("out/figure_supp_effect_size_correlations_durations.png",dpi=400)

def make_supplementary_figure_combined_pvalues():
  sampletypes = ["genus","species","pathway","metab"]
  fig,axes = plt.subplots(1,5,constrained_layout=True,figsize=(7.5,9.5))
  df = pd.DataFrame(results.table[(results.table.modeltype=="binary")].groupby("medication").apply(lambda df: scipy.stats.combine_pvalues(df["pvalue"])[1]).sort_values(),columns=["all"])
  for sampletype in sampletypes:
    df[sampletype] = pd.DataFrame(results.table[(results.table.sampletype==sampletype) & (results.table.modeltype=="binary")].groupby("medication").apply(lambda df: scipy.stats.combine_pvalues(df["pvalue"])[1]).sort_values(),columns=["all"])
  df = df.iloc[::-1]
  df = -np.log10(df)
  m = df[df!=np.inf].max()
  df = df.replace(np.inf,m)
  cmap = matplotlib.cm.get_cmap("plasma").copy()
  cmap.set_bad(color="#c0c0c0")
  for si,sampletype in enumerate(df.columns):
    colors = cmap((df[sampletype].rank() / df.shape[0]) * 0.75)
    axes[si].barh(np.arange(df.shape[0]), df.iloc[:,si].values, color=colors)
    axes[si].set_yticks([])
    axes[si].set_ylim(-1,df.shape[0]+1)
    axes[si].xaxis.set_tick_params(labelsize=6)
    axes[si].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4,steps=[1, 2, 4, 5, 10]))
    axes[si].set_xlabel("-$\mathregular{log_{10}}$ p", fontsize=8)
    for i in range(5,df.shape[0],5):
      axes[si].axhline(i-0.5,color="gray",alpha=0.5,linestyle="--",lw=0.5)
  axes[0].set_yticks(np.arange(df.shape[0]))
  axes[0].set_yticklabels(df.index.map(lambda x: prettify_compound_medication_name(x,length=100)),fontsize=6)
  fig.suptitle("Overall Impact of Medications on Microbiome Composition (Combined p, Fisher's method)",fontsize=8)
  axes[0].set_title("Total",fontsize=8)
  axes[1].set_title("Genus",fontsize=8)
  axes[2].set_title("Species",fontsize=8)
  axes[3].set_title("Pathway",fontsize=8)
  axes[4].set_title("Metabolite",fontsize=8)
  plt.savefig("out/figure_supplementary_overall_pvalues.pdf")
  plt.savefig("out/figure_supplementary_overall_pvalues.png",dpi=400)
  plt.close()

def make_figure1_panels():
  make_figure_1b()
  make_figure_1c()
  make_figure_1d()
  make_figure_1e()

def make_figure2_panels():
  make_figure_2a()
  make_figure_2b()
  make_figure_2c()
  make_figure_2d()
  make_figure_2e()
  make_figure_2f()

def make_figure3_panels():
  make_figure_3a()
  make_figure_3b()
  make_figure_3c()

def make_figure4_panels():
  make_figure_4a()
  make_figure_4b()
  make_figure_4c()
  make_figure_4d()
  make_figure_4e()

def make_supplementary_figures():
  make_supplementary_figure_metabolite_genus_correlations()
  make_spplementary_figure_effect_size_correlations_durations()
  make_supplementary_figure_combined_pvalues()
  make_supplementary_figure_invitro_overlap()

def arrange_figure_on_page(page, filename, tx, ty):
  reader = pypdf.PdfReader(filename)
  page_figure = reader.pages[0]
  transformation = pypdf.Transformation().translate(tx=72*tx, ty=72*ty)
  page_figure.add_transformation(transformation)
  mb = page.mediabox
  page_figure.mediabox = pypdf.generic.RectangleObject((mb.left, mb.bottom, mb.right, mb.top))
  page_figure.cropbox = pypdf.generic.RectangleObject((mb.left, mb.bottom, mb.right, mb.top))
  page_figure.trimbox = pypdf.generic.RectangleObject((mb.left, mb.bottom, mb.right, mb.top))
  page_figure.bleedbox = pypdf.generic.RectangleObject((mb.left, mb.bottom, mb.right, mb.top))
  page_figure.artbox = pypdf.generic.RectangleObject((mb.left, mb.bottom, mb.right, mb.top))
  page.merge_page(page_figure, expand=True)

def make_figure(filename, figures):
  writer = pypdf.PdfWriter()
  blank_page = writer.add_blank_page(width=7.5 * 72, height=10 * 72)
  for figure in figures:
    arrange_figure_on_page(blank_page, figure[0], figure[1], figure[2])
  with open(filename, "wb") as fp:
    writer.write(fp)
 
def add_labels(filename, labels):
  packet = io.BytesIO()
  can = reportlab.pdfgen.canvas.Canvas(packet, pagesize=reportlab.lib.pagesizes.letter)
  can.setFont("Arial Bold", 12)
  for label in labels:
    can.drawString(label[1], label[2], label[0])
  can.save()
  packet.seek(0)
  new_pdf = pypdf.PdfReader(packet)
  existing_pdf = pypdf.PdfReader(open(filename, "rb"))
  output = pypdf.PdfWriter()
  page = existing_pdf.pages[0]
  page.merge_page(new_pdf.pages[0])
  output.add_page(page)
  output_stream = open(filename, "wb")
  output.write(output_stream)
  output_stream.close()

def make_figure_panels():
  make_figure1_panels()
  make_figure2_panels()
  make_figure3_panels()
  make_figure4_panels()
  make_supplementary_figures()

def place_figures():
  make_figure("out/figure_1.pdf", [["out/figure_1a.pdf", 0, 6.5], ["out/figure_1b.pdf", 0, 3.35], ["out/figure_1c.pdf", 0, 0], ["out/figure_1d.pdf", 2.5, 0], ["out/figure_1e.pdf", 5, 0]])
  add_labels("out/figure_1.pdf", [["a.", 0.1 * 72,9.85 * 72],["b.", 0.1 * 72,6.35 * 72],["c.", 0.1 * 72,3.15 * 72],["d.", 2.5 * 72,3.15 * 72],["e.", 5 * 72,3.15 * 72]])
  make_figure("out/figure_2.pdf", [["out/figure_2a.pdf", 0, 6.65], ["out/figure_2b.pdf", 0, 3.25], ["out/figure_2d.pdf", 2, 3.25], ["out/figure_2c.pdf", 0, 0.10], ["out/figure_2e.pdf", 2.5, 0.10], ["out/figure_2f.pdf", 5, 0.10]])
  add_labels("out/figure_2.pdf", [["a.", 0.1 * 72,9.85 * 72],["b.", 0.1 * 72,6.55 * 72],["d.", 2.1 * 72,6.55 * 72],["c.", 0.1 * 72,3.1 * 72],["e.", 2.5 * 72,3.1 * 72],["f.", 5 * 72,3.1 * 72]])
  make_figure("out/figure_3.pdf", [["out/figure_3a.pdf", 0, 5.5], ["out/figure_3b.pdf", 0, 1.5], ["out/figure_3c.pdf", 4.5, 1.5]])
  add_labels("out/figure_3.pdf", [["a.", 0.1 * 72,9.85 * 72],["b.", 0.1 * 72,5.45 * 72],["c.", 4.6 * 72,5.45 * 72]])
  make_figure("out/figure_4.pdf", [["out/figure_4a.pdf", 0, 4.5], ["out/figure_4b.pdf", 0, 2], ["out/figure_4c.pdf", 3.75, 2], ["out/figure_4d.pdf", 0, 0], ["out/figure_4e.pdf", 4, 0]])
  add_labels("out/figure_4.pdf", [["a.", 0.1 * 72,9.85 * 72],["b.", 0.1 * 72,5.35 * 72],["c.", 3.85 * 72,5.35 * 72],["d.", 0.1 * 72,1.8 * 72],["e.", 3.85 * 72,1.8 * 72]])

poppler_path = r"/usr/bin"    ## UPDATE THIS TO YOUR POPPLER INSTALL DIRECTORY

def convert_to_png():
  for n in [1,2,3,4]:
    pages = pdf2image.convert_from_path("out/figure_%d.pdf"%n, 500, poppler_path=poppler_path)
    for count, page in enumerate(pages):
      page.save(f'out/figure_%d.png'%n, 'PNG')

def make_figures():
  make_figure_panels()
  place_figures()
  # convert_to_png()


