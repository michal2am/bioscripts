import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--lipid')
parser.add_argument('-s', '--system')
parser.add_argument('-c', '--chains', nargs='+')
parser.add_argument('-cp', '--chainspos', nargs='+')
args = parser.parse_args()

dataset = pd.read_csv('Interaction_{}/Dataset_{}/Dataset.csv'.format(args.lipid, args.lipid))      # pylipid generated dataset file
dataset['Residue Name'] = dataset['Residue'].apply(lambda name: name[-3:])
new_dataset_file = '{}_{}_chains_dataset.csv'.format(args.lipid, args.system)           # name for new dataset file
pdb_cg = 'step5_charmm2gmx.pdb'                                         # raw structure file generated with charmm-gui
prot_chains = ['PROA', 'PROB', 'PROC', 'PROD', 'PROE']                  # segnames for protein chains in structure file
chains_type = dict(zip(prot_chains, args.chains))
chains_typepos = dict(zip(prot_chains, args.chainspos))

#chains_type = {'PROA': 'b2', 'PROB': 'a1', 'PROC': 'b2', 'PROD': 'a1', 'PROE': 'g2'}
#chains_typepos = {'PROA': 'b2_1', 'PROB': 'a1_1', 'PROC': 'b2_2', 'PROD': 'a1_2', 'PROE': 'g2_1'}

residues = []
residueIds = []
residueChains = []


# this parser is hardcoded, should work with charmm-gui generated martini files
with open(pdb_cg) as f:
    prev_resnum = 0
    residueId = 0
    for line in f:
        if any(chain in line for chain in prot_chains):
            atom = line.strip().split()
            if atom[4] != prev_resnum:
                residues.append(atom[4]+atom[3])
                residueChains.append(atom[10])
                residueIds.append(residueId)
                prev_resnum = atom[4]
                residueId += 1

resid_data = pd.DataFrame(list(zip(residues, residueIds, residueChains)), columns=['Residue', 'Residue ID', 'Residue Chain'])
resid_data['Lipid'] = args.lipid
resid_data['Residue Num'] = resid_data['Residue'].apply(lambda num: int(num[0:-3]))
resid_data['Residue Name'] = resid_data['Residue'].apply(lambda name: name[-3:])
resid_data['System'] = args.system
resid_data['Residue Chain Type'] = resid_data['Residue Chain'].apply(lambda chain: chains_type[chain])
resid_data['Residue Chain TypePos'] = resid_data['Residue Chain'].apply(lambda chain: chains_typepos[chain])

dataset.drop(['Residue'], axis=1, inplace=True)             # in case of pdb input for pylipid messed up numbering
new_dataset = resid_data.merge(dataset, on=['Residue ID', 'Residue Name'])

new_dataset.to_csv(new_dataset_file, index=False)
print(new_dataset)
if len(dataset) != len(new_dataset):
    print('Problem with residue data merging, some residues omitted!')


# selected = new_dataset.iloc[:, 0:15]
# selected.drop(axis=1, inplace=True, labels=['Occupancy std', 'Lipid Count', 'Lipid Count std'])
# print(selected.sort_values(by=['Occupancy'], ascending=False).iloc[0:50, :])

# identify at which chain(s) binding site 0 is located
# print(new_dataset[new_dataset['Binding Site ID'] == 0])
