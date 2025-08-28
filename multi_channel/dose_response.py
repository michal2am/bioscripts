import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns


def ll4(x, b, c, d, e):
    """
    This function is basically a copy of the LL.4 function from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50
     """
    return c + (d - c)/(1 + np.exp(b * (np.log(x) - np.log(e))))


dr_data = pd.read_csv('e153_doseresponse_noE153C.csv')
print(dr_data)

compoundData = dr_data.groupby(['type'])
fitData = []
for name,group in compoundData:
    fitCoefs, covMatrix = opt.curve_fit(ll4, group.concentration, group.response)
    resids = group.response-group.concentration.apply(lambda x: ll4(x,*fitCoefs))
    curFit = dict(zip(['b','c','d','e'],fitCoefs))
    curFit['type']=name
    curFit['residuals']=sum(resids**2)
    fitData.append(curFit)
fitCompound = [ item['type'] for item in fitData]

fitTable = pd.DataFrame(fitData).set_index('type')
print(fitTable)

sns.set_style()
sns.set_context("paper")

refDose = np.logspace(-5, 5, base=10, num=2000)
g = sns.relplot(x='concentration', y='response', data=dr_data, hue='type', palette = sns.xkcd_palette(["pale red", 'windows blue', 'green']),
                height=2, aspect=1.5, hue_order=['WT', 'E153K', 'E153A'])
#g.despine(top=False, right=False)

for fit in fitData:
    # g.map(sns.lineplot, x=refDose, y=[ll4(j,*[fit[i] for i in ['b','c','d','e']]) for j in refDose])
    # this will be colored as previous be default:
    plt.plot(refDose, [ll4(j, *[fit[i] for i in ['b', 'c', 'd', 'e']]) for j in refDose])

# using plt interface:
# plt.yticks([0, 0.5, 1])
# plt.xscale('log')

for ax in g.axes.flat:
    ax.set_yticks([0, 0.5, 1])

g.set(xlabel='[mM GABA]', ylabel='relative amplitude', xscale='log')

plt.savefig('dose_response.png', dpi=300)
plt.show()

