import numpy as np
import pandas as pd
import plotly.express as px
from scikit_posthocs import outliers_iqr


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


# Bambi's meta file with all the data
meta = pd.read_csv('moje_meta_raw.csv', header=[0, 1])
# my file with Bambi's modeling results, just a copy with rows corresponding to the cell file name and rates
models = pd.read_csv('moje_meta_models_raw.csv', header=[0, 1])
models.drop_duplicates(inplace=True, ignore_index=True)                                                                 # necessary for duplicated WT cells

merged = meta.merge(models, left_on=[('meta', 'file')], right_on=[('meta', 'file')], how='left')                        # basically all Bambi's spreadsheets in one dataframe
merged.to_csv('moje_meta_merged_raw.csv')


# PLAYGROUND BELOW
# TODO: will not work, make great again by merging dataframe from ploter_all?

rates_vs_res()