import pandas as pd
import plotly.express as px

psis = pd.read_csv('phe_psis.csv', index_col=0)
#psis = pd.read_csv('dists.csv', index_col=0)
print(psis)
psis = psis.abs()
# print(psis)

fig1 = px.line(psis, template='presentation', width=1200, height=600,)
#fig1.show()
fig1.write_html('phe_dihe_plot.html')
fig2 = px.violin(psis, template='presentation', width=600, height=600,)
#fig2.show()