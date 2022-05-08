import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
#from kneed import KneeLocator
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import plotly.express as px

#rates = rates.melt(id_vars=['project', 'type'])


# correlations
'''
rates = pd.read_csv('published_rates_noPK_WTmean.csv')
correlations = rates.corr()
correlations[(correlations >= -.5) & (correlations <= .5)] = np.nan
sns.heatmap(correlations, annot=True)
plt.show()
'''

# wykres rate vs rate
'''
rates = pd.read_csv('published_rates_noPK_WTmean.csv')
all_pairs = sns.PairGrid(rates, hue='type')
all_pairs.map_offdiag(sns.scatterplot)
all_pairs.map_diag(sns.histplot)
plt.show()
'''

# hists
''''
rates = pd.read_csv('published_rates_noPK_WTmean.csv')
rates = rates.melt(id_vars=['project', 'type'])
hists = sns.FacetGrid(rates, col='variable')
hists.map(sns.kdeplot, 'value')
plt.show()
'''

# k-means

'''
rates = pd.read_csv('published_rates_noPK_WTmean.csv')
rates = rates.set_index('type')
print(rates)

selected_rates = ['d2', 'r2']
rates_names = ['gamma', 'delta', 'alpha2', 'beta2', 'alpha2p', 'beta2p', 'd2', 'r2', 'd2p', 'r2p']
to_drop = [rate for rate in rates_names if rate not in selected_rates]

rates.drop('project', axis=1, inplace=True)
rates.drop(to_drop, axis=1, inplace=True)

wcss = []
for i in range(1, 11):
    kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0)
    kmeans.fit(rates)
    wcss.append(kmeans.inertia_)
plt.plot(range(1, 11), wcss)
plt.title('Elbow Method')
plt.xlabel('Number of clusters')
plt.ylabel('WCSS')
plt.show()

kmeans = KMeans(n_clusters=3, init='k-means++', max_iter=300, n_init=10, random_state=0)
kmeans.fit(rates)
labels = kmeans.predict(rates)
rates['cluster'] = labels

print(rates)

ax = sns.scatterplot(data=rates, x=selected_rates[0], y=selected_rates[1], hue="cluster")
for line in range(0,rates.shape[0]):
     ax.text(rates[selected_rates[0]][line], rates[selected_rates[1]][line], rates.index[line], horizontalalignment='center', size='medium', color='black', weight='semibold')
plt.show()
'''

def kmeans_missing(X, n_clusters, max_iter=10):
    """Perform K-Means clustering on data with missing values.

    Args:
      X: An [n_samples, n_features] array of data to cluster.
      n_clusters: Number of clusters to form.
      max_iter: Maximum number of EM iterations to perform.

    Returns:
      labels: An [n_samples] vector of integer labels.
      centroids: An [n_clusters, n_features] array of cluster centroids.
      X_hat: Copy of X with the missing values filled in.
    """

    # Initialize missing values to their column means
    missing = ~np.isfinite(X)
    mu = np.nanmean(X, 0, keepdims=1)
    X_hat = np.where(missing, mu, X)

    for i in range(max_iter):
        if i > 0:
            # initialize KMeans with the previous set of centroids. this is much
            # faster and makes it easier to check convergence (since labels
            # won't be permuted on every iteration), but might be more prone to
            # getting stuck in local minima.
            cls = KMeans(n_clusters, init=prev_centroids)
        else:
            # do multiple random initializations in parallel
            cls = KMeans(n_clusters, n_jobs=-1)

        # perform clustering on the filled-in data
        labels = cls.fit_predict(X_hat)
        centroids = cls.cluster_centers_

        # fill in the missing values based on their cluster centroids
        X_hat[missing] = centroids[labels][missing]

        # when the labels have stopped changing then we have converged
        if i > 0 and np.all(labels == prev_labels):
            break

        prev_labels = labels
        prev_centroids = cls.cluster_centers_

    return labels, centroids, X_hat

# missing k-means
'''

rates = pd.read_csv('published_rates_noPK_WTmean.csv')
rates = rates.set_index('type')
rates.drop('project', axis=1, inplace=True)


print(rates)

labels, centroids, X_hat = kmeans_missing(rates, 6)
rates['cluster'] = labels
print(rates)
'''
selected_rates = ['delta', 'gamma']
selected_project = 'P273'

rates = pd.read_csv('published_rates_noPK_WTmean.csv')
selected = rates[rates['project'].isin(['WT', selected_project])]
selected['forward'] = np.log(selected[selected_rates[0]])
selected['equilibrium'] = np.log(selected[selected_rates[0]]/selected[selected_rates[1]])

print(selected)

project_cumuCells = px.scatter(selected,
                               x='equilibrium', y='forward',
                               title=selected_project,
                               color='type', template='presentation', width=800, height=800,
                               marginal_x='rug', marginal_y='rug',
                               color_discrete_sequence=px.colors.qualitative.Dark24,
                               )
project_cumuCells.add_trace(
    px.scatter(selected, x='equilibrium', y='forward',
               trendline='ols',
               color_discrete_sequence=px.colors.qualitative.Dark24,
               ).data[1])

project_cumuCells.write_html(selected_project + '_plot.html')

