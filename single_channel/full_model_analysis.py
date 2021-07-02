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

# clusters
rates = pd.read_csv('published_rates_noPK_WTmean.csv')
rates = rates.set_index('type')
print(rates)

rates.drop('project', axis=1, inplace=True)
rates.drop(['gamma', 'delta', 'alpha2p', 'beta2p', 'd2', 'r2', 'd2p', 'r2p'], axis=1, inplace=True)

print(rates)

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
labels, centroids, X_hat = kmeans_missing(rates, 4)
rates['cluster'] = labels
print(rates)
'''
# regular k-means

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

ax = sns.scatterplot(data=rates, x="beta2", y="alpha2", hue="cluster")
for line in range(0,rates.shape[0]):
     ax.text(rates.beta2[line], rates.alpha2[line], rates.index[line], horizontalalignment='center', size='medium', color='black', weight='semibold')
plt.show()