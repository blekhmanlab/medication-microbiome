import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def plot_stackedbar_p(df, labels, colors, title, subtitle, ax):
    fields = df.columns.tolist()
    
    # figure and axis
    # fig, ax = plt.subplots(1, figsize=(12, 10))
# plot bars
    left = len(df) * [0]
    # print(df.index)
    for idx, name in enumerate(fields):
        ax.barh(np.arange(len(df.index)), df[name], left = left, color=colors[idx], linewidth=0, tick_label=df.index)
        left = left + df[name]
# title and subtitle
    ax.set_title(title, loc='left')
    # ax.text(0, ax.get_yticks()[-1] + 0.75, subtitle)
# legend
    ax.legend(labels, ncol=3, frameon=False, fontsize='x-small', loc='best')
# remove spines
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
# format x ticks
    xticks = np.arange(0,1.1,0.1)
    xlabels = ['{}%'.format(i) for i in np.arange(0,101,10)]
    # ax.xticks(xticks, xlabels)
# adjust limits and draw grid lines
    ax.set_ylim(-0.5, ax.get_yticks()[-1] + 0.5)
    ax.xaxis.grid(color='gray', linestyle='dashed')


def distsq(x,y):
  return (x[0]-y[0])**2.0 + (x[1]-y[1])**2.0 + (x[2]-y[2])**2.0

def choose_colors(n, itr=1000):
  colorlist = ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#ffe119', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#d0d0d0', '#000000']
  if n < len(colorlist):
    return colorlist[:n]
  n -= len(colorlist)
  delta = 0.01
  e0 = np.array([1,0,0])
  e1 = np.array([0,1,0])
  e2 = np.array([0,0,1])
  colors = np.random.rand(n,3)
  x0 = colors

  # maximize the minimum difference between colors.
  def f(x):
    # print('f')
    x = x.reshape(-1,3)
    v = 1
    for i in range(x.shape[0]):
      for j in range(x.shape[0]):
        if i<=j:
          continue
        # print(x[0])
        # print('v',np.linalg.norm(x[i]-x[j]))
        v = np.min([v, distsq(x[i],x[j])])
        # print('vvv')
    # v is the minimum difference between colors.
    # print(v)
    if(np.any(x<0) or np.any(x>1)):
      return float('inf')
    return -v

  if itr == 0:
    return colorlist + list(x0.reshape(n,3))
  res = minimize(f, x0, method='Nelder-Mead', tol=0.001, options={'maxiter':itr})
  # print(res.x.reshape(n,3))
  return colorlist + list(res.x.reshape(n,3))

def make_labels(labels):
  colorlist = choose_colors(len(labels),0)
  unique = pd.Series(labels).value_counts().index
  nlabels = len(unique)
  colormap = {unique[i]:colorlist[i] for i in range(nlabels)}
  return (unique, colormap, [colormap[l] for l in labels])

