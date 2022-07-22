import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

meta = pd.read_csv('moje_meta_raw.csv', header=[0, 1])
#print(meta[['meta', 'openings', 'shuts']])

meta.columns = ['_'.join(col) for col in meta.columns.values]
meta.dropna(subset=['meta_cells_no'], inplace=True)

rates = (meta[['shuts_t1', 'shuts_t2', 'shuts_t3', 'shuts_t4', 'openings_t1', 'openings_t2']])

correlations = rates.corr()
#correlations[(correlations >= -.5) & (correlations <= .5)] = np.nan
sns.heatmap(correlations, annot=True)
plt.show()