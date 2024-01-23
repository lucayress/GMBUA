# -*- coding: utf-8 -*-
"""
Synthetic data results - Raincloud plots

Reference:
https://github.com/RainCloudPlots/RainCloudPlots.git
Tutorial:
https://github.com/pog87/PtitPrince/blob/27260b4337af269f356218f8ba9ca15819082fd1/tutorial_python/raincloud_tutorial_python.ipynb

@author: luciano Ayres
Created on Thu Oct 26 18:06:36 2023
"""

import numpy as np
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import matplotlib.collections as clt
import h5py
import sys
import ptitprince as pt

# print(os.environ['PATH'])
sns.set(style="whitegrid",font_scale=1, font="Times New Roman", rc={"text.usetex": True})

sys.path.append('../ptitprince/')
sys.path.append('C:\\Users\\lucia\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64')
data_dir = 'test_data/synthetic4/'

#%% Functions

def create_dataframe(file_name):
  # Load file
  mat_contents = h5py.File(file_name,'r')
  # mat_contents.keys()
  
  # Read variables
  K = int(mat_contents.get('K')[:].item())
  SRE_FCLSU_k = mat_contents.get('SRE_FCLSU_k')[:]
  SRE_collaborative_k = mat_contents.get('SRE_collaborative_k')[:]
  SRE_group_k = mat_contents.get('SRE_group_k')[:]
  SRE_elitist_k = mat_contents.get('SRE_elitist_k')[:]
  SRE_fractional_k = mat_contents.get('SRE_fractional_k')[:]
  SRE_SUnCNN_k = mat_contents.get('SRE_SUnCNN_k')[:]
  SRE_MUA_k = mat_contents.get('SRE_MUA')[:]
  SRE_MUA_k10 = mat_contents.get('SRE_MUA10')[:]
  SRE_MUA_k20 = mat_contents.get('SRE_MUA20')[:]

  # Create dataframe
  names = ['FCLSU','Collaborative','Group','Elitist','Fractional','SUnCNN', r'GMBUA, $K=10$', r'GMBUA, $K=20$', r'GMBUA, $K=30$']
  scores = [SRE_FCLSU_k, SRE_collaborative_k, SRE_group_k, SRE_elitist_k, SRE_fractional_k, SRE_SUnCNN_k, SRE_MUA_k10, SRE_MUA_k20, SRE_MUA_k]
  df = pd.DataFrame()
  for i in range(len(names)):
    aux = []
    for j in range(K):
      aux.append(names[i])
      algo = pd.DataFrame(aux, columns=[''])
      score = pd.DataFrame(scores[i], columns=['SRE (dB)'])
      df_aux = algo.join(score)
    df = pd.concat([df, df_aux], ignore_index=True)
  return df

def create_dataframe2(file_name):
  # Load file
  mat_contents = h5py.File(file_name,'r')
  # mat_contents.keys()
  
  # Read variables
  K = int(mat_contents.get('K')[:].item())
  SRE_FCLSU_k = mat_contents.get('SRE_FCLSU_k')[:]
  SRE_collaborative_k = mat_contents.get('SRE_collaborative_k')[:]
  SRE_group_k = mat_contents.get('SRE_group_k')[:]
  SRE_elitist_k = mat_contents.get('SRE_elitist_k')[:]
  SRE_fractional_k = mat_contents.get('SRE_fractional_k')[:]
  SRE_SUnCNN_k = mat_contents.get('SRE_SUnCNN_k')[:]
  SRE_MUA_k = mat_contents.get('SRE_MUA')[:]

  # Create dataframe
  names = ['FCLSU','Collaborative','Group','Elitist','Fractional','SUnCNN','GMBUA']
  scores = [SRE_FCLSU_k, SRE_collaborative_k, SRE_group_k, SRE_elitist_k, SRE_fractional_k, SRE_SUnCNN_k, SRE_MUA_k]
  df = pd.DataFrame()
  for i in range(len(names)):
    aux = []
    for j in range(K):
      aux.append(names[i])
      algo = pd.DataFrame(aux, columns=[''])
      score = pd.DataFrame(scores[i], columns=['SRE (dB)'])
      df_aux = algo.join(score)
    df = pd.concat([df, df_aux], ignore_index=True)
  return df

#%% synthetic - 20dB SNR - R = 30, K = 30

image_name = 'synthetic4'
fname =  'synthetic4_compare_all_noise20_R30.mat'
df = create_dataframe2(data_dir + 'results/' + fname)

R = 30
SNR = 20

# Raincloud plot
dx = 'SRE (dB)'; dy = ''; ort = "h"; pal = "Set2"; sigma = .2
f, ax = plt.subplots(figsize=(5, 5))

ax=pt.RainCloud(x = dy, y = dx, data = df, palette = pal, bw = sigma,
                 width_viol = .6, ax = ax, orient = ort, move = .2, alpha = .65)
plt.xlim(2, 12.0)
plt.xlabel('SRE' + r'($\hat{\mathbf{Z}}$)' + ' (dB)')
plt.xticks(np.arange(2, 12+1, step=1))
# plt.title('K = ' + str(K))
plt.show()
# Save figure
fig_name = image_name + '_rainclouds' + '_noise' + str(SNR) + '_R' + str(R)
# f.savefig(data_dir + '/prints/' + fig_name + '.pdf', format="pdf", bbox_inches="tight")
plt.close()

#%% synthetic - 30dB SNR - R = 30, K = 30

image_name = 'synthetic4'
fname =  'synthetic4_compare_all_noise30_R30.mat'
df = create_dataframe2(data_dir + 'results/' + fname)

R = 30
SNR = 30

# Raincloud plot
dx = 'SRE (dB)'; dy = ''; ort = "h"; pal = "Set2"; sigma = .2
f, ax = plt.subplots(figsize=(5, 5))

ax=pt.RainCloud(x = dy, y = dx, data = df, palette = pal, bw = sigma,
                 width_viol = .6, ax = ax, orient = ort, move = .2, alpha = .65)
plt.xlim(3, 13.0)
plt.xlabel('SRE' + r'($\hat{\mathbf{Z}}$)' + ' (dB)')
plt.xticks(np.arange(3, 13+1, step=1))
# plt.title('K = ' + str(K))
plt.show()
# Save figure
fig_name = image_name + '_rainclouds' + '_noise' + str(SNR) + '_R' + str(R)
# f.savefig(data_dir + '/prints/' + fig_name + '.pdf', format="pdf", bbox_inches="tight")
plt.close()

#%% synthetic - 20dB SNR - R = 30, K = 10, 20, 30

image_name = 'synthetic4'
fname =  'synthetic4_compare_all_noise20_R30.mat'
df = create_dataframe(data_dir + 'results/' + fname)

R = 30
SNR = 20

# Raincloud plot
dx = 'SRE (dB)'; dy = ''; ort = "h"; pal = "Set2"; sigma = .2
f, ax = plt.subplots(figsize=(5, 5))

ax=pt.RainCloud(x = dy, y = dx, data = df, palette = pal, bw = sigma,
                 width_viol = .6, ax = ax, orient = ort, move = .2, alpha = .65)
plt.xlim(2, 12.0)
plt.xlabel('SRE' + r'($\hat{\mathbf{Z}}$)' + ' (dB)')
plt.xticks(np.arange(2, 12+1, step=1))
# plt.title('K = ' + str(K))
plt.show()
# Save figure
fig_name = image_name + '_rainclouds' + '_noise' + str(SNR) + '_R' + str(R) + '_K'
# f.savefig(data_dir + '/prints/' + fig_name + '.pdf', format="pdf", bbox_inches="tight")
plt.close()

#%% synthetic - 30dB SNR - R = 30, K = 10, 20, 30

image_name = 'synthetic4'
fname =  'synthetic4_compare_all_noise30_R30.mat'
df = create_dataframe(data_dir + 'results/' + fname)

R = 30
SNR = 30

# Raincloud plot
dx = 'SRE (dB)'; dy = ''; ort = "h"; pal = "Set2"; sigma = .2
f, ax = plt.subplots(figsize=(5, 5))

ax=pt.RainCloud(x = dy, y = dx, data = df, palette = pal, bw = sigma,
                 width_viol = .6, ax = ax, orient = ort, move = .2, alpha = .65)
plt.xlim(3, 13.0)
plt.xlabel('SRE' + r'($\hat{\mathbf{Z}}$)' + ' (dB)')
plt.xticks(np.arange(3, 13+1, step=1))
# plt.title('K = ' + str(K))
plt.show()
# Save figure
fig_name = image_name + '_rainclouds' + '_noise' + str(SNR) + '_R' + str(R) + '_K'
# f.savefig(data_dir + '/prints/' + fig_name + '.pdf', format="pdf", bbox_inches="tight")
plt.close()

# No raindrops

#%% synthetic - 20dB SNR - R = 30, K = 10, 20, 30

image_name = 'synthetic4'
fname =  'synthetic4_compare_all_noise20_R30.mat'
df = create_dataframe(data_dir + 'results/' + fname)

R = 30
SNR = 20

# Raincloud plot
dx = 'SRE (dB)'; dy = ''; ort = "h"; pal = "Set2"; sigma = .2
f, ax = plt.subplots(figsize=(5, 4))

ax=pt.RainCloud(x = dy, y = dx, data = df, palette = pal, bw = sigma,
                 linewidth = 0.6, width_viol = .6, width_box = .2, 
                 ax = ax, orient = ort, move = .0, alpha = .8, offset = None,
                 boxplot = True, raindrops = True)
plt.xlim(2, 12.0)
plt.xlabel('SRE' + r'($\hat{\mathbf{Z}}$)' + ' (dB)',fontsize = 10)
plt.xticks(np.arange(2, 12+1, step=1))
# plt.title('K = ' + str(K))
plt.show()
# Save figure
fig_name = image_name + '_rainclouds' + '_noise' + str(SNR) + '_R' + str(R) + '_K'
f.savefig(data_dir + '/prints/' + fig_name + '_' '.pdf', format="pdf", bbox_inches="tight")
plt.close()

#%% synthetic - 20dB SNR - R = 30, K = 10, 20, 30 - no raindrops

image_name = 'synthetic4'
fname =  'synthetic4_compare_all_noise20_R30.mat'
df = create_dataframe(data_dir + 'results/' + fname)

R = 30
SNR = 20

# Raincloud plot
dx = 'SRE (dB)'; dy = ''; ort = "h"; pal = "Set2"; sigma = .2
f, ax = plt.subplots(figsize=(5, 4))

ax=pt.RainCloud(x = dy, y = dx, data = df, palette = pal, bw = sigma,
                 linewidth = 0.6, width_viol = .6, width_box = .2, 
                 ax = ax, orient = ort, move = .0, alpha = .8, offset = None,
                 boxplot = True, raindrops = False)
plt.xlim(2, 12.0)
plt.xlabel('SRE' + r'($\hat{\mathbf{Z}}$)' + ' (dB)',fontsize = 10)
plt.xticks(np.arange(2, 12+1, step=1))
# plt.title('K = ' + str(K))
plt.show()
# Save figure
fig_name = image_name + '_rainclouds' + '_noise' + str(SNR) + '_R' + str(R) + '_K'
f.savefig(data_dir + '/prints/' + fig_name + '__' '.pdf', format="pdf", bbox_inches="tight")
plt.close()