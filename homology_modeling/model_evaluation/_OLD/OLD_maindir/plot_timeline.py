import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sns.set_style()
sns.set_context("paper")

test = pd.read_csv('6i53.csv')


test.drop(axis=1, columns=['Unnamed: 0'], inplace=True)
test['position'] = test.position.astype(float)
test['model'] = test.model.astype(float)
test['rmsf'] = test.rmsf.astype(float)


aggregated = test.pivot_table(index=['template', 'subunit', 'position'], columns='type', values='rmsf', aggfunc=np.mean)
aggregated['relative gly'] = aggregated['gly']/aggregated['wt']
aggregated['relative lys'] = aggregated['lys']/aggregated['wt']
#aggregated['relative gly'] = (aggregated['relative gly'] - aggregated['relative gly'].min())/(aggregated['relative gly'].max() - aggregated['relative gly'].min())
#aggregated['relative lys'] = (aggregated['relative lys'] - aggregated['relative lys'].min())/(aggregated['relative lys'].max() - aggregated['relative lys'].min())

aggregated['relative gly'] = aggregated['relative gly'] - aggregated['relative gly'].mean()
aggregated['relative lys'] = aggregated['relative lys'] - aggregated['relative lys'].mean()

aggregated.to_csv('test.csv')
#print(aggregated.sort_values('relative lys', ascending=False))
print(aggregated)
print(aggregated.sort_values('relative lys', ascending=False))


'''
g = sns.FacetGrid(new, col='type',
                  aspect=1, height=1.7125, sharex=False)
g.map(sns.lineplot, "position", "rmsf", "subunit")
#g.map(sns.lineplot, "position", "relative gly", "subunit", legend=False,)

# g.map(sns.lineplot, "position", "template_ene", "subunit", palette=sns.color_palette("pastel", 2))
#g.set_titles("{col_name}")
plt.tight_layout()

g.savefig("test.png", dpi=600)

plt.show()



aggregated = test.pivot_table(index=['template', 'subunit', 'position'], columns='type', values='rmsf', aggfunc=np.mean)

aggregated['relative lys'] = aggregated['lys']/aggregated['wt']
aggregated['relative gly'] = aggregated['gly']/aggregated['wt']

print(aggregated[aggregated['position'] == 46])


#aggregated['relative'] = (aggregated['relative'] - aggregated['relative'].min())/(aggregated['relative'].max() - aggregated['relative'].min())
#print(aggregated.sort_values('relative', ascending=False))

#print(aggregated.loc['6huo'].sort_values('relative', ascending=False))
aggregated.reset_index(inplace=True)


newtest = test.pivot_table(index=['template', 'subunit', 'position', 'type'], columns='model', values='rmsf')

#print(newtest)




#nlys = lys.pivot(index='position', columns='model', values='rmsf')

#print(test)

#new_test = test.pivot(index='position', columns='subunit', values='rmsf')


#print(new_test)



g = sns.FacetGrid(aggregated, col='template', col_order=['6i53', '6huk', '6hup', '6huo', '6huj', '6hug'],
                  aspect=1, height=1.7125,  col_wrap=4, sharex=False)
g.map(sns.lineplot, "position", "relative lys", "subunit", legend=False,)
g.map(sns.lineplot, "position", "relative gly", "subunit", legend=False,)

# g.map(sns.lineplot, "position", "template_ene", "subunit", palette=sns.color_palette("pastel", 2))
g.set_titles("{col_name}")
plt.tight_layout()

g.savefig("test.png", dpi=600)

plt.show()
'''