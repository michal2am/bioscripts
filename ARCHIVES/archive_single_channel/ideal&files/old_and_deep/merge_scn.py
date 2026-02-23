from dcpyps import dataset
from dcpyps import dcio
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--files', nargs='+')
parser.add_argument('--tres', type=float)
parser.add_argument('--tcrit', type=float)
args = parser.parse_args()

rec1 = dataset.SCRecord(args.files, 100e-9, args.tres, args.tcrit)
#print(rec1.itint, [int(ampli) for ampli in rec1.iampl/0.0018939346773549914], rec1.iprops)

print(rec1)
# dcio.scn_write(rec1.itint, [int(ampli) for ampli in rec1.iampl/0.0018939346773549914], rec1.iprops,
#               calfac=0.0018939346773549914, ffilt=1.428999900817871, Emem=100.0)


# rec2 = dataset.SCRecord(['new_saved.SCN'], 100e-9, args.tres, args.tcrit)
