import numpy as np
import pandas as pd
import plotly.express as px

meta = pd.read_csv('moje_meta_raw.csv', header=[0, 1])
models = pd.read_csv('moje_meta_models_raw.csv', header=[0, 1])
models.drop_duplicates(inplace=True, ignore_index=True)

merged = meta.merge(models, left_on=[('meta', 'file')], right_on=[('meta', 'file')], how='left')

merged.to_csv('moje_meta_merged_raw.csv')
print(merged)

def rates_vs_res():

    model_vs_res = merged[['rates', 'meta']].droplevel(0, axis=1)
    model_vs_res = model_vs_res[model_vs_res['min_res'].notna()]                                                        # only first row with full cell data
    model_vs_res['min_res'] = pd.to_numeric(model_vs_res['min_res'])                                                    # res as numeric

    print(model_vs_res)

    for rate in ['b2', 'a2', 'gamma', 'delta']:
        fig = px.scatter(model_vs_res, x='min_res', y=rate, color='residue_mut', symbol='type')
        fig.show()
        # fig = px.scatter(model_vs_res, x='fc', y=rate, color='residue_mut', symbol='type')
        # fig.show()


event_times = merged[['shuts', 'openings']]

def mean_shut(cell):
    return cell[('shuts', 't1')] / 2

merged[('shuts', 't_test')] = merged.apply(lambda cell: mean_shut(cell), axis=1)
print(merged['shuts'])


#event_times.columns = [' '.join(col).strip() for col in event_times.columns.values]
#print(event_times)