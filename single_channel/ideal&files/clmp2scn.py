import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from dcpyps import dcio


#clampfit_dwell_11 = pd.read_csv('WT_25_C1_Clampfit_11.csv')
#g = sns.relplot(kind='line', x='Start' , y='Level', data=clampfit_dwell_11, drawstyle='steps-pre', color='grey')
#plt.show()


idealized_deep = pd.read_csv('result_89145k1_2kFilter_Red.csv')
print(idealized_deep)

# reverse cum sum below is useless, as intervals are of constant values

'''
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def C(cumulative):
    yield cumulative[0]
    for a, b in pairwise(cumulative):
        yield b-a


intervals = []

for interval in C(idealized_deep['Time']):
    intervals.append(interval*10e5)
'''

predictions = list(idealized_deep['prediction'] - 2)
#flags = len(predictions)*[0]
#intervals = len(predictions)*[0.06]

resolution = idealized_deep['Time'][1]
print(resolution*1000)

previous = predictions[0]
current_length = 1
current_prediction = predictions[0]

c_intervals = []
c_predictions = []

for i, pred in enumerate(predictions[1:]):
    if previous == pred:
        #print('same event')
        current_length += 1
    else:
        #print('new event')
        c_intervals.append(current_length * 60)
        c_predictions.append(current_prediction)
        current_length = 1
        current_prediction = pred
    previous = pred

    #print(i, params[0], params[1])

print(c_intervals)
print(c_predictions)
c_flags = len(c_intervals)*[0]
#print(c_flags)
print(len(c_intervals))

dcio.scn_write(c_intervals, c_predictions, c_flags, filename='f45.SCN')

dcio.scn_write([0.1, 2, 4, 0.5, 4.0], [0, 1, 0, 1, 0], 5*[0], filename='unique.SCN')
dcio.scn_write([0.1, 2, 4, 0.5, 4.0], [0, 1, 1, 1, 0], 5*[0], filename='notunique.SCN')

#ioffset, nint, calfac2, header = dcio.scn_read_header('f45.SCN')
#print(type(header['ioffset']))
#dcio.scn_read_data('f45.SCN', dcio.scn_read_header('f45.SCN'))

print('xdd')