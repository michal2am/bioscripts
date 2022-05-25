import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--lipid")
args = parser.parse_args()

dataset = pd.read_csv('Interaction_{}/Dataset_{}/Dataset.csv'.format(args.lipid, args.lipid))      # pylipid generated dataset file
new_dataset_file = '{}_chains_dataset.csv'.format(args.lipid)           # name for new dataset file
pdb_cg = 'step5_charmm2gmx.pdb'                                      # raw structure file generated with charmm-gui
prot_chains = ['PROA', 'PROB', 'PROC', 'PROD', 'PROE']                  # segnames for protein chains in structure file

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

new_dataset = resid_data.merge(dataset)
new_dataset.to_csv(new_dataset_file)

# selected = new_dataset.iloc[:, 0:15]
# selected.drop(axis=1, inplace=True, labels=['Occupancy std', 'Lipid Count', 'Lipid Count std'])
# print(selected.sort_values(by=['Occupancy'], ascending=False).iloc[0:50, :])

# identify at which chain(s) binding site 0 is located
# print(new_dataset[new_dataset['Binding Site ID'] == 0])
