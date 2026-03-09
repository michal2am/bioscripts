import pandas as pd

# is StartPres really a start?

contacts = pd.read_csv('contacts.csv')
contacts["PresenceList"] = contacts["PresenceVector"].apply(lambda s: [int(c) for c in s])
#print(contacts.columns)

contacts['StartPres'] = contacts["PresenceVector"].str[0].astype(int)
contacts["EquiStab200"] = contacts["PresenceList"].apply(lambda lst: sum(lst[:200]))/200
contacts["EquiStab"] = contacts["PresenceList"].apply(lambda lst: sum(lst[:2000]))/2000
contacts["RunStab"] = contacts["PresenceList"].apply(lambda lst: sum(lst[2001:])/(len(lst)-2000))

#contacts.drop(["PresenceVector", "PresenceList"], axis=1, inplace=True)


pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

to_print = contacts[["ChainA", "ResnameA","ResidA", "ChainB", "ResnameB", "ResidB", "Type", "StartPres", "EquiStab200", "EquiStab", "RunStab"]]

print(to_print[(contacts['StartPres'] == 1)
               & (contacts['EquiStab200'] < 0.9)
               & (contacts['EquiStab'] < 1)
               & (contacts['RunStab'] < 0.5)
               & (contacts['ResidA'] < 220)])

# here EquiStab should exclude EquiStab200 part ...
print(to_print[(contacts["StartPres"] == 1)
               & (contacts['ResidA'] < 200)
               & (contacts['EquiStab'] < 0.75*contacts['EquiStab200'])])

print(to_print[(contacts["StartPres"] == 1)
               & (contacts['ResidA'] < 200)
               & (contacts['RunStab'] < 0.75*contacts['EquiStab'])])