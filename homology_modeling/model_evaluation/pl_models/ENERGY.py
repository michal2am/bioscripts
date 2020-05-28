import vmd
import pandas as pd
import plotly.express as px

# TODO: those are some code chunks for incorporation energy data into new model table data format

d_fs = ['6i53_doubleCys_all_fit.pdb',  '6i53_doubleCysPatch_all_fit.pdb',
        '6hug_doubleCys_all_fit.pdb', '6hug_doubleCysPatch_all_fit.pdb']

e_fs = ['6i53_doubleCys_localeval.csv', '6i53_doubleCys_localeval.csv',
        '6hug_doubleCys_localeval.csv', '6hug_doubleCys_localeval.csv']

distances = []

for d_f, e_f in zip(d_fs, e_fs):

    molid = vmd.molecule.load('pdb', d_f)

    dict = {'1st': vmd.measure.bond(atom1=19, atom2=2898, molid=molid, first=0, last=99),
            '2nd': vmd.measure.bond(atom1=8174, atom2=11053, molid=molid, first=0, last=99),}
            #'1stPair': vmd.measure.bond(atom1=3956, atom2=600, molid=molid, first=0, last=99)}
    distance = pd.DataFrame(dict)

    meta = d_f.split('_')
    distance['template'] = meta[0]
    distance['patch'] = meta[1]

    distance.reset_index(inplace=True)
    distance.rename(columns={'index': 'model'}, inplace=True)

    long_distance = pd.melt(distance, id_vars=['template', 'patch', 'model'], value_vars=['1st', '2nd'],
                            var_name='interface', value_name='distance')



    energies =pd.read_csv(e_f)
    energy1 = energies[((energies['subunit'] == 'a1_a') & (energies['position'].isin(['12', '13', '14', '15', '16', '17', '18', '19', '81', '84', '85']))) | ((energies['subunit'] == 'b2_b')) & (energies['position'].isin(['28', '29', '30', '31', '32', '33', '34', '26', '163']))]
    energy2 = energies[((energies['subunit'] == 'a1_d') & (energies['position'].isin(['12', '13', '14', '15', '16', '17', '18', '19', '81', '84', '85']))) | ((energies['subunit'] == 'b2_e')) & (energies['position'].isin(['28', '29', '30', '31', '32', '33', '34', '26', '163']))]

    energy1.drop(columns=['residue number', 'residue', 'subunit', 'position'], inplace=True)
    energy1 = energy1.transpose()
    energy1 = energy1.astype(float)
    energy1['1st_ene'] = energy1.sum(axis=1)
    energy1.reset_index(inplace=True, drop=True)

    sums1 = energy1['1st_ene'].copy()

    energy2.drop(columns=['residue number', 'residue', 'subunit', 'position'], inplace=True)
    energy2 = energy2.transpose()
    energy2 = energy2.astype(float)
    energy2['2nd_ene'] = energy2.sum(axis=1)
    energy2.reset_index(inplace=True, drop=True)

    sums2 = energy2['2nd_ene'].copy()

    sums1 = sums1.to_frame()
    sums2 = sums2.to_frame()
    sums1.reset_index(inplace=True)
    sums2.reset_index(inplace=True)

    sums1.rename(columns={'index': 'model', '1st_ene': '1st'}, inplace=True)
    long_sums1 = pd.melt(sums1, id_vars=['model'], value_vars=['1st'],
                            var_name='interface', value_name='energy')
    sums2.rename(columns={'index': 'model', '2nd_ene': '2nd'}, inplace=True)
    long_sums2 = pd.melt(sums2, id_vars=['model'], value_vars=['2nd'],
                            var_name='interface', value_name='energy')


    all_sums = pd.concat([long_sums1, long_sums2], axis=0)
    all_sums.reset_index(inplace=True, drop=True)
    all_sums.drop(columns=['model', 'interface'], inplace=True)    # merge elegantly instead


    complete = pd.concat([long_distance, all_sums], axis=1)




    distances.append(complete)


distances = pd.concat(distances)
print(distances)
'''
fig = px.scatter(distances[distances['interface'] == '2nd'], x='2nd_ene', y='distance', facet_col='template', color='2nd_ene',
               template="presentation", width=700, height=500
               )

'''
fig1 = px.strip(distances[distances['patch'] == 'doubleCys'], x='template', y='distance', color='interface',
               template="presentation", width=700, height=500
               )

fig1 = px.strip(distances[distances['patch'] == 'doubleCys'], x='energy', y='distance', facet_col='interface', facet_row='template',
               template="presentation", width=700, height=500
               )

fig2 = px.strip(distances[distances['patch'] == 'doubleCysPatch'], x='interface', y='distance', facet_col="template",
               template="presentation", width=700, height=500
               )

'''
fig.add_shape(
    type="line", line_color="salmon", line_width=3, opacity=1, line_dash="dot",
    x0=0, x1=1, xref="paper", y0=5, y1=5, yref="y"
    )

#fig.layout.yaxis1.update(matches=None)
'''

fig1.write_html('test1.html')
fig2.write_html('test2.html')