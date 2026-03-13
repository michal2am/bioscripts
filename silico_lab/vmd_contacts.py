import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px


# is StartPres really a start?

contactsAB = pd.read_csv('contacts_A_sys2.csv', dtype={'PresenceVector': str})
contactsCD = pd.read_csv('contacts_C_sys2.csv', dtype={'PresenceVector': str})

contacts = pd.concat([contactsAB, contactsCD], ignore_index=True)

contacts["PresenceList"] = contacts["PresenceVector"].apply(lambda s: [int(c) for c in s])
#print(contacts.columns)

contacts['StartPres'] = contacts["PresenceVector"].str[0].astype(int)
contacts["EquiStab200"] = contacts["PresenceList"].apply(lambda lst: sum(lst[:200]))/200
contacts["EquiStab"] = contacts["PresenceList"].apply(lambda lst: sum(lst[:2000]))/2000
contacts["RunStab"] = contacts["PresenceList"].apply(lambda lst: sum(lst[2001:])/(len(lst)-2000))

#contacts.drop(["PresenceVector", "PresenceList"], axis=1, inplace=True)


pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

contacts_without_vector = contacts[["ChainA", "ResnameA", "ResidA", "ChainB", "ResnameB", "ResidB", "Type", "StartPres", "EquiStab200", "EquiStab", "RunStab"]]

print(contacts_without_vector[(contacts['StartPres'] == 1)
                              & (contacts['EquiStab200'] < 0.75)
                              & (contacts['EquiStab'] < 1.01)
                              & (contacts['RunStab'] < 1.01)
                              & (contacts['ResidA'] < 220)])

# here EquiStab should exclude EquiStab200 part ...
print(contacts_without_vector[(contacts["StartPres"] == 1)
                              & (contacts['ResidA'] < 220)
                              & (contacts['EquiStab'] < 0.75*contacts['EquiStab200'])])

print(contacts_without_vector[(contacts["StartPres"] == 1)
                              & (contacts['ResidA'] < 220)
                              & (contacts['RunStab'] < 0.75*contacts['EquiStab'])])


# plot

selected_contacts = contacts_without_vector[(contacts_without_vector['StartPres'] == 1)
                                            & (contacts_without_vector['EquiStab200'] < 1.00)
                                            & (contacts_without_vector['EquiStab'] < 1.01)
                                            & (contacts_without_vector['RunStab'] < 1.01)
                                            & ((contacts_without_vector['ResidA'].between(190,220) |
                                                contacts_without_vector['ResidB'].between(190,220)))
                                            & (contacts_without_vector['Type'] == 'hbond')].copy()

# Build x-axis labels: ResnameA_ResidA – ResnameB_ResidB
selected_contacts["label"] = ( selected_contacts["ResnameA"] + " " + selected_contacts["ResidA"].astype(str) + " – "
                               + selected_contacts["ResnameB"] + " " + selected_contacts["ResidB"].astype(str))

selected_contacts['BS'] = np.where(selected_contacts['ChainA'].isin(['A', 'B']), "1st", "2nd")

selected_contacts.sort_values(['ResidA', 'ResidB'], inplace=True)
print(selected_contacts)

# Melt the three stability columns into long format
selected_contacts_melted = selected_contacts.melt(
    id_vars=["label", "Type", "ResidA", "ResidB", "BS"],
    value_vars=["EquiStab200", "EquiStab", "RunStab"],
    var_name="Metric",
    value_name="Stability",
)

selected_contacts_melted.sort_values(['ResidA', 'ResidB'], inplace=True)
ordered_labels = selected_contacts_melted['label'].tolist()


fig = px.bar(
    selected_contacts_melted[selected_contacts_melted['Type'] == 'hbond'],
    x="label",
    y="Stability",
    color="Metric",
    facet_row="BS",
    facet_col="Type",
    barmode="group",
    title="Interaction Stability per Residue Pair",
    labels={"label": "Residue A – Residue B"},
    height=700,
    template="plotly_white",
    category_orders={"label": ordered_labels},
)

fig.update_xaxes(tickangle=-45)
fig.show()