import pandas as pd
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt

'''
distances = pd.read_csv('metdist_all.csv')
print(distances)
fig = px.histogram(distances, marginal='rug', title='MET268 to etomidate distance', histnorm='percent')
fig.write_html('metdistances.html')
fig.show()
'''

all_rmsd = pd.DataFrame()

for datfile in ['eto1st_rmsd_ll232a.dat', 'eto2nd_rmsd_ll232a.dat', 'eto1st_rmsd_m286al232a.dat', 'eto2nd_rmsd_m286al232a.dat']:

    rmsd = pd.read_csv(datfile, sep=' ')
    all_rmsd[datfile] = rmsd.iloc[:, 1]

print(all_rmsd)
rmsd_plot = sns.histplot(all_rmsd, kde=True)
plt.show()
rmsd_plot.figure.savefig('eto_rmsd.png')
