

from biopandas.pdb import PandasPdb
import os
import pandas as pd

dirs = os.listdir('../..')

read_pdbs = []

for file_name in dirs:
    if file_name.endswith('.pdb'):

        raw_read = PandasPdb().read_pdb(file_name).df['ATOM']
        raw_read.loc[:, 'model'] = int(file_name[-11:-8])
        raw_read.loc[:, 'template'] = '6i53'
        raw_read.loc[:, 'type'] = 'LYS'
        raw_read.drop(columns=['record_name', 'alt_loc', 'insertion', 'occupancy', 'b_factor', 'segment_id',
                                  'element_symbol', 'charge', 'line_idx', 'blank_1', 'blank_2', 'blank_3', 'blank_4'],
                      inplace=True)

        print(raw_read)
        print(raw_read.loc[(raw_read['residue_number'] == 46) & (raw_read['chain_id'].isin(['A', 'D']))])























        #parsed_read = raw_read.df['ATOM'][(raw_read.df['ATOM']['residue_number'] == 46) & (raw_read.df['ATOM']['chain_id'].isin(['A', 'D']))].copy()

        #read_pdbs.append(parsed_read)

#all_parsed = pd.concat(read_pdbs, axis=0)
#print(all_parsed)

        #print(parsed_read)
#print(new[['atom_number', 'atom_name', 'residue_name', 'chain_id', 'residue_number', 'x_coord', 'y_coord', 'z_coord']])