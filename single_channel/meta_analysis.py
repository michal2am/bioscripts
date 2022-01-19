import pandas as pd
import plotly.express as px


meta = pd.read_csv('moje_meta_raw.csv', header=[0, 1])
models = pd.read_csv('moje_meta_models_raw.csv', header=[0, 1])
merged = meta.merge(models)

merged.to_csv('moje_meta_merged_raw.csv')

def rates_vs_res():

    model_vs_res = merged[['rates', 'meta']].droplevel(0, axis=1)
    model_vs_res = model_vs_res[model_vs_res['min_res'].notna()]
    model_vs_res['min_res'] = pd.to_numeric(model_vs_res['min_res'])

    print(model_vs_res)

    for rate in ['b2', 'a2', 'gamma', 'delta']:
        fig = px.scatter(model_vs_res, x='min_res', y=rate, color='residue_mut', symbol='type')
        fig.show()
        # fig = px.scatter(model_vs_res, x='fc', y=rate, color='residue_mut', symbol='type')
        # fig.show()


print(merged['t_crit'])