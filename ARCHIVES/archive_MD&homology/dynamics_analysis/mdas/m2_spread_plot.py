import pandas as pd
import plotly.express as px

# gmx trjconv -f step5_input.pdb -s step7_1.tpr -o step5_input_chains.pdb

dists = pd.read_csv('M2_spread.dat', delimiter='\t', header=None)
print(dists)

fig = px.scatter(dists, trendline="ols")
#fig.add_trace(px.scatter(dists, trendline='ols').data[1])
#fig = px.violin(dists)
fig.show()

