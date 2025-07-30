import numpy as np
import os
import joblib
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap

def fit_embedding(input_data_path, output_folder,
                  pca_ncomp={'sw': 4, 'theta': 4},
                  n_neighbors=15, 
                  n_components=2):
    """
    Fits a low-dimensional embedding of LFP features using PCA followed by ISOMAP,
    and saves the resulting models into the output folder.

    Parameters
    ----------
    input_data_path : str
        Path to the compressed .npz file containing the input data.
    output_folder : str
        Folder where the fitted models will be saved.
    pca_ncomp : dict
        Dictionary defining the number of PCA components for each network pattern ('sw', and 'theta' in original data).
    n_neighbors : int
        Number of neighbors for ISOMAP.
    n_components : int
        Number of output dimensions for ISOMAP.
    """

    # --- Load input data ---
    # This file contains (z-scored) mean waveforms obtained for different ephys network patterns
    data = np.load(input_data_path)
    pattern_labels = data.files  # Expected: 'sw', 'theta'

    # --- Prepare and fit PCA models ---
    pca_models = {}
    for pattern_label in pattern_labels:
        if pattern_label in pca_ncomp:
            pca = PCA(n_components=pca_ncomp[pattern_label])
            pca.fit(data[pattern_label])
            pca_models[pattern_label] = pca
            pca_models[pattern_label].scores = pca.transform(data[pattern_label])  # keep scores attached

    # --- Concatenate PCA scores across patterns in specified order ---
    scores_concat = np.concatenate([pca_models[label].scores for label in ['sw', 'theta']],axis=1)

    # --- Fit ISOMAP on concatenated PCA-reduced features ---
    embedding_model = Isomap(n_neighbors=n_neighbors, n_components=n_components)
    embedding_model.fit(scores_concat)

    # --- Create output folder if it doesn't exist ---
    os.makedirs(output_folder, exist_ok=True)

    # --- Save models ---
    joblib.dump(pca_models, os.path.join(output_folder, 'pca'))
    joblib.dump(embedding_model, os.path.join(output_folder, 'iso'))

def project(data, pca_model, embedding_model):
    """
    Projects LFP features into a prebuilt low-dimensional embedding.

    Parameters
    ----------
    data : dict
        Dictionary containing LFP network patterns for multiple points 
        (in the original implementation, 'sw' and 'theta').
        Each key maps to an array of shape (n_points, n_time_samples).
    pca_model : dict
        Dictionary containing fitted PCA models for each LFP network pattern.
        Keys should match the network pattern labels.
    embedding_model : object
        Prebuilt embedding model used to transform the PCA-reduced data into low-dimensional coordinates.
        In the original work, this is a fitted sklearn.manifold.Isomap object.

    Returns
    -------
    embedding_proj : ndarray
        Array of shape (n_points, n_dim) containing the low-dimensional coordinates
        for each input point. Points with invalid data are filled with NaNs.
    """

    # --- Ensure correct feature order ---
    pattern_labels = ['sw', 'theta']  # Hardcoded to avoid order-dependent bugs

    # --- Identify valid points based on first pattern ---
    valid_points = ~np.isnan(data[pattern_labels[0]][:, 0])
    n_points = data[pattern_labels[0]].shape[0]

    # --- Initialize output ---
    embedding_proj = np.full((n_points, embedding_model.n_components), np.nan)

    # --- PCA projections ---
    pca_scores = {}
    for pattern_label in pattern_labels:
        pca_scores[pattern_label] = pca_model[pattern_label].transform(data[pattern_label][valid_points])

    # --- Concatenate PCA scores in correct order ---
    pca_scores_concat = np.concatenate([pca_scores[label] for label in pattern_labels],axis=1)

    # --- Embed the valid points ---
    embedding_proj[valid_points] = embedding_model.transform(pca_scores_concat)

    return embedding_proj
