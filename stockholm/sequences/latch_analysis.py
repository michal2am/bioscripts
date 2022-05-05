import pandas as pd
from Bio import SeqIO

latch01_list = []

with open('latch01.fasta') as fasta_file:  # Will close handle cleanly
    identifiers = []
    lengths = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        seq_name = seq_record.id.split('|')[2].split('_')[0]
        resid = seq_record.seq
        latch01_list.append([seq_name[0:3], seq_name, resid])

latch01_data = pd.DataFrame(latch01_list, columns=['receptor', 'seq_name', 'latch01'])

latch02_list = []

with open('latch02.fasta') as fasta_file:  # Will close handle cleanly
    identifiers = []
    lengths = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        seq_name = seq_record.id.split('|')[2].split('_')[0]
        resid = seq_record.seq
        latch02_list.append([seq_name[0:3], seq_name, resid])

latch02_data = pd.DataFrame(latch02_list, columns=['receptor', 'seq_name', 'latch02'])

# print(latch01_data)
# print(latch02_data)

latch_data = pd.merge(left=latch01_data, right=latch02_data)
# print(latch_data)
latch01_resnames = [str(seq)[0] for seq in latch_data.latch01.unique()]
latch02_resnames = [str(seq)[0] for seq in latch_data.latch02.unique()]
print('Latch01 (M3) residue types in pLGICs family: {}'.format(latch01_resnames))
print('Latch02 (M1) residue types in pLGICS family: {}'.format(latch02_resnames))

interface_list = []

for receptor in latch_data['receptor'].unique():
    selected_receptor = latch_data[latch_data['receptor'] == receptor]
    print(selected_receptor.sort_values(by='seq_name'))
    for seq_name01 in selected_receptor['seq_name']:
        for seq_name02 in selected_receptor['seq_name']:
            interface = ('{} : {}'.format(seq_name01, seq_name02))
            latch_01 = selected_receptor[selected_receptor['seq_name'] == seq_name01]['latch01'].values[0]
            latch_02 = selected_receptor[selected_receptor['seq_name'] == seq_name02]['latch02'].values[0]
            interface_list.append([receptor, interface, latch_01, latch_02])

interface_data = pd.DataFrame(interface_list, columns=['receptor', 'interface', 'latch01', 'latch02'])
print(interface_data)
