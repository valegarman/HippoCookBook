#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
import numpy as np
from scipy.stats import zscore
from joblib import load
import sys, os
from pathlib import Path
import argparse


# --- Configure paths ---

parser = argparse.ArgumentParser()
parser.add_argument('--folder', required=True, help='Path to the folder containing .mat files')
parser.add_argument('--path', required=True, help='Path to Hipp-LFP-embedding folder')
args = parser.parse_args()

folder = args.folder
path = args.path
fnames = ['wave_ripple_mean', 'wave_theta_mean']  # .mat filenames (without extension)

# --- Load and preprocess LFP waveform data ---
data = {}
for fname in fnames:
    filepath = os.path.join(folder, f'{fname}.mat')
    data_ = loadmat(filepath)[fname]
    key = 'sw' if 'ripple' in fname else 'theta'
    data[key] = data_


# Pad ripple waveform with one sample at beginning and end
data['sw'] = np.hstack((
     data['sw'][:, 0][:, None],
     data['sw'],
     data['sw'][:, -1][:, None]
 ))
data['sw'] = data['sw'].T            # now shape (627, 64) - 64 features
data['theta'] = data['theta'].T      # now shape (187, 64)


# Z-score across samples for each channel
for key in data:
    data[key] = zscore(data[key], axis=1)

# --- Load embedding models and project data ---

sys.path.append(os.path.abspath(path))
from hipp_embedding import embedding, trajectory

pca_model_single = load(f'{path}/models/pca')
pca_model = {
    'sw': pca_model_single,
    'theta': pca_model_single
}
embedding_model = load(f'{path}/models/iso')

# Debug info: stampa le shape dei dati e delle PCA
print("Shape di data['sw']:", data['sw'].shape)
print("Shape di data['theta']:", data['theta'].shape)

print("Numero componenti PCA sw:", pca_model['sw'].n_components_)
print("Numero componenti PCA theta:", pca_model['theta'].n_components_)

projections = embedding.project(data, pca_model, embedding_model)

# --- Define the anatomical trajectory and identify labeled control points ---

train_data_path = f'{path}/data/trajectory_points.npz'
traj, ctrlpts_proj, ctrlpts_labels = trajectory.define_trajectory(pca_model, embedding_model, input_data=train_data_path)

layeris = [i for i, label in enumerate(ctrlpts_labels) if '_interp' not in label]
layer_labels = [ctrlpts_labels[i] for i in layeris]
meanproj_layer = np.nanmean(ctrlpts_proj, axis=0)

# RGB color mapping for anatomical layers
layer_colors = {
    'pyr': np.array([0.949, 0.439, 1.000]),
    'rad': np.array([0.212, 0.655, 0.086]),
    'lm' : np.array([0.867, 0.710, 0.239]),
    'hf' : np.array([0.000, 0.000, 0.000]),
    'om' : np.array([0.349, 1.000, 0.749]),
    'mm' : np.array([0.376, 0.067, 1.000]),
    'im' : np.array([0.596, 0.851, 0.184]),
    'gr' : np.array([0.557, 0.545, 0.161])
}

# --- Plot 2D embedding with anatomical layers and channel projections ---
plt.figure(figsize=(5, 5))
plt.plot(traj[:, 1], traj[:, 0], color='k', lw=3)

for i in layeris:
    label = ctrlpts_labels[i]
    color = layer_colors[label]
    plt.plot(meanproj_layer[i, 1], meanproj_layer[i, 0],
             marker='o', linestyle='none',
             markerfacecolor=color, markeredgecolor='black',
             markeredgewidth=1.5, markersize=12)

cmap = plt.cm.jet
colors = cmap(np.linspace(0, 1, len(projections)))
for chi, proj_ in enumerate(projections):
    plt.plot(proj_[1], proj_[0], 'v', mew=3, markersize=8,
             color=colors[chi], label=f'Ch {chi+1}')

for i, label in enumerate(layer_labels):
    plt.text(0.75, 0.9 - i * 0.04, label,
             transform=plt.gca().transAxes,
             fontsize=14, fontweight='bold',
             color=layer_colors[label],
             verticalalignment='top')

for chi in range(len(projections)):
    plt.text(1.05, 0.95 - chi * 0.035, f'Ch {chi+1}',
             transform=plt.gca().transAxes,
             fontsize=10, color=colors[chi],
             verticalalignment='top')

plt.grid(True)
plt.axis('equal')
plt.xlabel('Feature component 2', fontsize=13)
plt.ylabel('Feature component 1', fontsize=13)
plt.tight_layout()

outpath1 = os.path.join(folder, "lfp_embedding_2d_projection.png")
plt.savefig(outpath1, dpi=300)
print(f"Saved: {outpath1}")

plt.close()


# --- Compute linearized coordinate along anatomical axis ---
precision = 2
pyr_trajis = np.zeros((len(ctrlpts_proj)))
for proji, proj_ in enumerate(ctrlpts_proj[:, 0]):
    traji = np.argmin(np.sum((proj_ - traj)**2, axis=1))
    pyr_trajis[proji] = traji
offset = np.mean(pyr_trajis / precision)

trajis = np.zeros((len(projections)))
for proji, proj_ in enumerate(projections):
    traji = np.argmin(np.sum((proj_ - traj)**2, axis=1))
    trajis[proji] = traji
trajis = trajis / precision - offset

# --- Plot ripple waveforms with their anatomical coordinate ---
plt.figure(figsize=(5, 8))
_plotoffsets = []
for chi, lfp in enumerate(data['sw']):
    offset = 3
    plt.plot(lfp - chi * offset, lw=2)
    _plotoffsets.append(-chi * offset)
    plt.text(625, -chi * offset, f'{trajis[chi]:.2f}', fontweight='bold')

plt.yticks(_plotoffsets, 1 + np.arange(len(_plotoffsets)))
plt.xticks([len(lfp) / 2], [0])
plt.xlabel('Time (samples)')
plt.ylabel('Channels')
plt.grid(True)
plt.tight_layout()
plt.savefig("lfp_embedding_ripples.png", dpi=300)
plt.close()

# --- Save embedding results to .mat file ---
p = Path(folder)
folder_name = p.name
filename = f"{folder_name}.lfp_embedding.channelinfo.mat"

lfp_embeding = {
    'data_projection_1d': trajis.T,
    'data_projection_2d': projections,
    'model_trajectory': traj,
    'training_data_points': meanproj_layer[layeris],
    'training_data_labels': layer_labels,
    'all_training_data_points': meanproj_layer,
    'all_training_data_labels': ctrlpts_labels,
    'precision': precision,
}

savemat(filename, {'lfp_embeding': lfp_embeding})