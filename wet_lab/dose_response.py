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


compoundData = dr_data.groupby(['type'])
fitData = []
for name,group in compoundData:
    fitCoefs, covMatrix = opt.curve_fit(ll4, group.concentration, group.response)
    resids = group.response-group.concentration.apply(lambda x: ll4(x,*fitCoefs))
    curFit = dict(zip(['b','c','d','e'],fitCoefs))
    curFit['type']=name[0]
    curFit['residuals']=sum(resids**2)
    curFit['x_concs'] = np.logspace(-5, 5, base=10, num=2000)
    curFit['y_fited'] = [ll4(j,*[curFit[i] for i in ['b','c','d','e']]) for j in curFit['x_concs']]
    fitData.append(curFit)
fitCompound = [ item['type'] for item in fitData]

fit_data = pd.DataFrame(fitData)
all_dr = pd.merge(dr_data, fit_data, on="type")
print(all_dr)
all_dr.to_csv('dose_reponse_fromplot.csv')

sns.set_style()
sns.set_context("paper")

curves = all_dr.groupby("type").first().reset_index()
curves = curves[["type", "x_concs", "y_fited"]].explode(["x_concs", "y_fited"])
curves["x_concs"] = curves["x_concs"].astype(float)
curves["y_fited"] = curves["y_fited"].astype(float)

# --- Step 2: make scatter base ---
palette = sns.xkcd_palette(["pale red", "windows blue", "green"])
order = ["WT", "E153K", "E153A"]

fig, ax = plt.subplots(figsize=(3, 2))

sns.scatterplot(
    data=all_dr,
    x="concentration", y="response",
    hue="type", palette=palette,
    hue_order=order,
    ax=ax, s=40, zorder=2
)

# --- Step 3: add SEM error bars ---
for t, sub in all_dr.groupby("type"):
    ax.errorbar(
        sub["concentration"], sub["response"],
        yerr=sub["response_sem"],
        fmt="none", ecolor=palette[order.index(t)],
        elinewidth=1, capsize=2, zorder=1
    )

# --- Step 4: overlay fitted curves ---
for t, sub in curves.groupby("type"):
    ax.plot(
        sub["x_concs"], sub["y_fited"],
        color=palette[order.index(t)],
        label=f"{t} fit", linewidth=1.5
    )

ax.set_xscale("log")  # often for concentration-response curves
ax.set_yticks([0, 0.5, 1])
ax.set_xlabel('[mM GABA]')
ax.set_ylabel('relative amplitude')
ax.legend().remove()


sns.despine()
plt.tight_layout()
plt.savefig('dose_response.png', dpi=300)
plt.show()

