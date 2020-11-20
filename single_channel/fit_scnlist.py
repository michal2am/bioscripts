import pandas as pd

# find . -name \*.SCN -printf '%p\n' > scn_list.csv
# find . -name \*.SCN -exec cp {} workspace_hjcfit \;
scn_list = pd.read_csv('scn_list.csv', sep='/', header=None)
print(scn_list)
com_scn = scn_list.iloc[:, -2:]
com_scn.to_csv('scn_list_cells.csv')