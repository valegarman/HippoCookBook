# Hipp-LFP-embedding

This code was developed by **Vítor Lopes dos Santos**.  
If you have questions, suggestions, or criticism, feel free to open an issue or reach out by email.

📧 [vtlsantos@gmail.com](mailto:vtlsantos@gmail.com)

This repository provides code, data, and pre-trained models for constructing a low-dimensional embedding of hippocampal LFP activity and computing an anatomical trajectory across hippocampal layers. 

It accompanies the paper:

**[Full citation placeholder — to be added]**

---

## Folder Structure

```text
Hipp-LFP-embedding/
├── data/                                      # Input data used to generate the embedding
│   └── trajectory_points.npz                  # Contains theta/SW waveforms and layer labels
│   └── README.md                              # Describes the structure and content of trajectory_points.npz
├── hipp_embedding/                            # Core Python package with embedding and trajectory functions
│   ├── __init__.py                            # Minimal init file to make it a package
│   ├── embedding.py                           # Functions to fit and project onto the embedding
│   └── trajectory.py                          # Functions to define trajectory
├── models/                                    # Pre-trained PCA and Isomap models used in the example notebook
│   ├── pca                                    # Trained PCA models per LFP feature (from original publication)
│   └── iso                                    # Trained Isomap embedding model (from original publication)
├── notebooks/                                 # Example notebooks for using the package
│   └── visualising_embedding_trajectory.ipynb # Demonstrates projection and trajectory visualization
├── requirements.txt                           # Exact package versions used (for reproducibility)
└── README.md                                  # Project overview and usage instructions

```

## Environment

Tested with:

- Python 3.10.12
- numpy 1.23.5
- scipy 1.15.2
- scikit-learn 1.2.2
- joblib 1.4.2
- matplotlib 3.10.1

> ⚠️ **Note:** The saved models may not be compatible with newer versions of `scikit-learn`. Please install the versions specified above.

---

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/your-username/Hipp-LFP-embedding.git
cd Hipp-LFP-embedding

# Create and activate a virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate

# Install required packages
pip install -r requirements.txt
