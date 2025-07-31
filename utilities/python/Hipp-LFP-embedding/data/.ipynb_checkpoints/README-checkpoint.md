# trajectory_points.npz

This file contains mean LFP waveforms used to compute the anatomical trajectory in the embedding space.

## Contents

`trajectory_points.npz` is a NumPy `.npz` archive containing the following keys:

- **'theta'** : ndarray of shape `(6, 13, 189)`  
  Mean theta waveforms from 6 mice  
  13 control points per mouse (including anatomical layers and interpolated positions)  
  189 time samples per waveform, aligned to the descending zero-crossing of the pyramidal-layer theta cycle

- **'sw'** : ndarray of shape `(6, 13, 629)`  
  Mean sharp wave waveforms from the same mice and control points  
  629 time samples per waveform, aligned to the peak of ripple-band amplitude

- **'layers'** : ndarray of shape `(13,)`  
  String labels for the anatomical layers and interpolated positions  
  Example: `'pyr'`, `'rad'`, `'rad-lm_interp1'`, etc.

## Notes

- LFPs were sampled at **1250 Hz**
