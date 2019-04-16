import pandas as pd
import numpy as np
from functools import reduce


class ReadPropkaLog():

    def __init__(self, file_names, pre_chains, post_chains):
        self.files = file_names
        self.chains_map = dict(zip(pre_chains, post_chains))
        self.read = self.read_files()

    def read_files(self):

        dfs = []

        for file in self.files:

            with open(file) as rf:
                parsed = []
                for row in rf:
                    if any(substring in row[0:3] for substring in ["ARG", "LYS", "ASP", "GLU", "HIS", "TYR"]):
                        line = row.split()
                        if (len(line)) == 22:

                            if '*' in line[3]:
                                newline = [line[0], line[1], self.chains_map[line[2]], np.NaN]
                            else:
                                newline = [line[0], line[1], self.chains_map[line[2]], float(line[3])]

                            fullname = "{}{}{}".format(newline[0], newline[1], newline[2])
                            newline.insert(0, fullname)
                            parsed.append(newline)

                            labels = ["fullname", "resname", "resnum", "chain", file[0:4]]
                            df = pd.DataFrame.from_records(parsed, columns=labels)

            dfs.append(df)

        df_merged = reduce(lambda left, right: pd.merge(left, right, on=["fullname", "resname", "resnum", "chain"], how='outer'), dfs)

        return df_merged


test = ReadPropkaLog(['6hug_alig_prot.pka', '6huk_alig_prot.pka', '6hup_alig_prot.pka', '6huj_alig_prot.pka',
                      '6huo_alig_prot.pka', '6i53_alig_prot.pka'], ["B", "A", "E", "D", "C"], ["B3_1", "A1_1", "B3_2", "A1_2", "Y2"])

result = test.read
result['mean'] = result.iloc[:, 4:].mean(axis=1)
result['center'] = result['mean'] - 7.2
result = result.sort_values('center', axis=0)

print(result)
result.to_csv('test.csv')


