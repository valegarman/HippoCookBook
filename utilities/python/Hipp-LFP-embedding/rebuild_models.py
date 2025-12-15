import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap
import joblib

# Crea dei dati finti di esempio
data = np.random.rand(100, 20)

# PCA
pca = PCA(n_components=10)
pca.fit(data)
data_pca = pca.transform(data)

# Isomap
iso = Isomap(n_neighbors=15, n_components=2)
iso.fit(data_pca)

# Salva i modelli
joblib.dump(pca, 'E:\\Hipp-LFP-embedding-main\\Hipp-LFP-embedding-main\\models\\pca')
joblib.dump(iso, 'E:\\Hipp-LFP-embedding-main\\Hipp-LFP-embedding-main\\models\\iso')

print("Modelli salvati!")