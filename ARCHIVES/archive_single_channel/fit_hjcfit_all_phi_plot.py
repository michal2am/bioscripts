import numpy as np
import pandas as pd
import plotly.express as px
import statsmodels
import argparse
import re
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import math
import matplotlib.ticker as ticker

phis = pd.read_csv('all_phis.csv')
print(phis)
sns.barplot(x='residue', y='phi', data=phis, capsize=0.2, errwidth=1, yerr=phis['std_err'])
plt.show()