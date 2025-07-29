import numpy as np
from scipy.interpolate import interp1d
from hipp_embedding.embedding import project

def interpolate_trajectory(ctrlpts_proj, precision=2):
    """
    Interpolates a smooth trajectory through the average of projected control points 
    from multiple sessions.

    Parameters
    ----------
    ctrlpts_proj : ndarray
        Array of shape (n_sessions, n_control_points, n_dim), containing the low-dimensional
        projection of each control point in each recording session.
    precision : float
        Spatial resolution in control point units. Determines how many interpolation
        points to insert between each pair of control points (higher = more dense).

    Returns
    -------
    traj : ndarray
        Smooth interpolated trajectory of shape (n_interp_points, n_dim).
    """

    # --- Compute the mean trajectory across sessions ---
    ctrlpts_proj_mean = np.nanmean(ctrlpts_proj, axis=0)  # shape: (n_control_points, n_dim)

    # --- Create interpolated spacing along the control point index axis ---
    x = np.arange(ctrlpts_proj_mean.shape[0])
    x_new = []

    for i in range(len(ctrlpts_proj_mean) - 1):
        dist = np.linalg.norm(ctrlpts_proj_mean[i + 1] - ctrlpts_proj_mean[i])
        n_interp = max(int(np.round(precision * dist)), 1)
        x_new_segment = np.linspace(i, i + 1, n_interp, endpoint=False)
        x_new.extend(x_new_segment)

    x_new.append(len(ctrlpts_proj_mean) - 1)
    x_new = np.array(x_new)

    # --- Interpolate each dimension independently ---
    n_dim = ctrlpts_proj_mean.shape[1]
    traj = np.zeros((len(x_new), n_dim))

    for dim in range(n_dim):
        f = interp1d(x, ctrlpts_proj_mean[:, dim], kind='quadratic')
        traj[:, dim] = f(x_new)

    return traj


def extrapolate_trajectory(traj, nextrappoints=[50, 20]):
    """
    Extrapolates a trajectory by linearly extending it before the first point 
    and after the last point.

    Parameters
    ----------
    traj : ndarray
        Trajectory of shape (n_points, n_dim).
    nextrappoints : list of two ints
        Number of points to extrapolate before and after the trajectory.

    Returns
    -------
    traj_extended : ndarray
        Extended trajectory of shape (n_points + sum(nextrappoints), n_dim).
    """

    # Before
    dx0 = traj[1] - traj[0]
    extrap_before = traj[0] + np.outer(np.arange(-nextrappoints[0], 0), dx0)

    # After
    dx1 = traj[-1] - traj[-2]
    extrap_after = traj[-1] + np.outer(np.arange(1, nextrappoints[1] + 1), dx1)

    traj_extended = np.vstack((extrap_before, traj, extrap_after))

    return traj_extended


def define_trajectory(pca_model, embedding_model, input_data='Hipp-LFP-embedding/data/trajectory_points.npz'):
    """
    Builds the full embedding trajectory by projecting waveform data from all sessions,
    then interpolating and extrapolating it.

    Parameters
    ----------
    pca_model : dict
        Dictionary of fitted PCA models, keyed by pattern label (in this implementation: 'theta' and 'sw').
    embedding_model : object
        Fitted low-dimensional embedding model (in this implementation: sklearn.manifold.Isomap).
    input_data : str
        Path to a .npz file containing 'theta', 'sw', and 'layers'.

    Returns
    -------
    traj : ndarray
        Interpolated and extrapolated trajectory through the low-dimensional space.
        Shape: (n_points_interp, n_dim).
    ctrlpts_proj : ndarray
        Raw projected control points per session before interpolation.
        Shape: (n_mice, n_control_points, n_dim).
    layer_labels : ndarray
        Array of anatomical or interpolated labels corresponding to the control points.
        Shape: (n_control_points,).
    """

    # --- Load waveform data ---
    wfs = np.load(input_data)  # see data/README.md for structure
    theta = wfs['theta']       # shape: (n_mice, n_points, n_time)
    sw = wfs['sw']             # shape: (n_mice, n_points, n_time)
    ctrlpts_labels = wfs['layers'] 

    n_mice, n_points = theta.shape[:2]
    n_dim = embedding_model.n_components

    ctrlpts_proj = np.zeros((n_mice, n_points, n_dim))

    # --- Project each session into embedding space ---
    for mousei in range(n_mice):
        wfs_ = {'theta': theta[mousei], 'sw': sw[mousei]}
        proj_ = project(wfs_, pca_model, embedding_model)
        ctrlpts_proj[mousei] = proj_

    # --- Build interpolated + extrapolated trajectory ---
    traj = interpolate_trajectory(ctrlpts_proj)
    traj = extrapolate_trajectory(traj)

    return traj, ctrlpts_proj, ctrlpts_labels 

