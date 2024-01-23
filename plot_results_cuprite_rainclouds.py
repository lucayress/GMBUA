# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 18:06:36 2023

@author: lucia
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

sns.set(style="whitegrid",font_scale=1)
sys.path.append('../ptitprince/')
data_dir = 'real_data/cuprite/'

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
  SRE_MUA = mat_contents.get('SRE_MUA')[:]

  SRE_MUA_k = SRE_MUA

  # Create dataframe

  names = ['FCLSU','Collab.','Group','Elitist','Fractional','GMBUA']
  scores = [SRE_FCLSU_k, SRE_collaborative_k, SRE_group_k, SRE_elitist_k, SRE_fractional_k, SRE_MUA_k]
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

#%% synthetic - R = 30, K = 30

image_name = 'cuprite'
fname =  'cuprite_compare_noise0_R30_K30.mat'
df = create_dataframe(data_dir + 'results/' + fname)

R = 30
SNR = 0

# Raincloud plot
dx = 'SRE (dB)'; dy = ''; ort = "h"; pal = "Set2"; sigma = .2
f, ax = plt.subplots(figsize=(5, 5))

ax=pt.RainCloud(x = dy, y = dx, data = df, palette = pal, bw = sigma,
                 width_viol = .6, ax = ax, orient = ort, move = .2, alpha = .65)
plt.xlim(31, 41.0)
plt.xlabel('SRE' + r'($\hat{\mathbf{Y}}$)' + ' (dB)')
plt.xticks(np.arange(31, 41+1, step=1))
# plt.title('K = ' + str(K))
plt.show()
# Save figure
fig_name = image_name + '_rainclouds' + '_noise' + str(SNR) + '_R' + str(R)
f.savefig(data_dir + '/prints/' + fig_name + '.pdf', format="pdf", bbox_inches="tight")
plt.close()