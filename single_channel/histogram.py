import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.scale as mscale
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.ticker as ticker
import argparse as ap


class SCh_histogram:

    def __init__(self, fname):

        self.fname = fname
        self.data, self.calc, self.mntime, self.mxtime, self.meven = self.read()

    def read(self):
        with open(self.fname) as f:

            raw = f.readlines()
            dataid, calcid, seriid = [], [], []

            for (id, line) in enumerate(raw):
                if 'raw' in line: dataid.append(id)
                if 'calc' in line: calcid.append(id)
                if len(line.split()) == 1: seriid.append(id)

            data = [map(float, line.split()) for line in raw[seriid[0]+1:seriid[0]+int(raw[seriid[0]])+1]]

            calc = []
            for i in range(len(seriid) - 1):
                calc.append([map(float, line.split()) for line in raw[seriid[i+1]+1:seriid[i+1]+int(raw[seriid[i+1]])+1]])

            data = pd.DataFrame(data).set_index(keys=0)
            calc = [pd.DataFrame(series).set_index(keys=0) for series in calc]

            meven = max([series.max().iloc[0] for series in calc])
            mxtime = max([series.index[-1] for series in calc])
            mntime = sorted([series.index[0] for series in calc])[1]                  # magic sorted, because missed event may be longer

            return data, calc, mntime, mxtime, meven

    def last_bin(self):
        bins = sch_hist.data.index.values
        diffs = np.diff(bins)
        for dif in diffs:
            print(dif)
        mpl.pyplot.plot(diffs)
        plt.show()

    def plot(self):

        # general plot settings
        mscale.register_scale(SquareRootScale)
        sns.set_style("white")
        sns.set_style("ticks")

        # figure generation
        fig, ax = plt.subplots(figsize=(6.5, 4))
        ax.set_xscale('log', basex=10)
        ax.set_yscale('squareroot')

        ### figure variables ###

        # data
        self.data.plot(ax=ax, drawstyle='steps-post', color='black', lw=3)
        self.calc[0].plot(ax=ax, color='black', lw=3)
        self.calc[-1].plot(ax=ax, color='grey', lw=3, style='--')
        for series in sch_hist.calc[1:-1]:
            series.plot(ax=ax, color='black', lw=1)

        # tick ranges
        ymax = int(self.meven ** 0.5)                               # magic need 'last before' power function argument
        if ymax % 2 == 1: ymax += 1

        ys = [(x * 2) ** 2 for x in range(0, int(ymax / 2 + 2))]    # magick scale if +2 so w need /2 and +2 for additional entry
        yms = [x ** 2 for x in np.arange(0, ymax + 2, 0.5)]         # magick scale is not divided so only +2 for additional entry

        last = ys[-1] if len(ys) % 2 == 1 else 0                    # too many ys? reduce, but preserve last
        lastm = yms[-1] if len(yms) % 2 ==1 else 0

        if len(ys) > 10:
            ys = ys[0:-1:4]
            yms = yms[0:-1:4]

        if len(ys) > 5:
            ys = ys[0:-1:2]
            yms = yms[0:-1:2]

        if last !=0: ys.append(last)
        if lastm !=0: yms.append(lastm)

        ax.set_yticks(ys)
        ax.set_yticks(yms, minor=True)
        ax.tick_params(labelsize=14)

        xmax = int(np.log10(self.mxtime))    # magic again, 'last before' base 10 log
        xmin = int(np.log10(self.mntime))    # magic again, first base 10 log

        # axes ranges
        ax.set_ylim([0, (ymax + 2) ** 2])
        ax.set_xlim([10 ** xmin, 10 ** (xmax + 1)])

        # overwrite labels and legend
        ax.set_xlabel('apparent open time [ms] (log scale)', fontsize=14)
        ax.set_ylabel('frequency [ ] (square root scale)', fontsize=14)
        ax.legend(['observed', 'fit', 'missed event correction', 'components fits'], loc='best', fontsize=14)


        sns.despine()
        plt.tight_layout()

        # plt.show()
        fig.savefig(self.fname.split('.')[0] + ".png")


class SquareRootScale(mscale.ScaleBase):
    """
    ScaleBase class for generating square root scale.
    """

    name = 'squareroot'

    def __init__(self, axis, **kwargs):
        mscale.ScaleBase.__init__(self)

    def set_default_locators_and_formatters(self, axis):
        axis.set_major_locator(ticker.AutoLocator())
        axis.set_major_formatter(ticker.ScalarFormatter())
        axis.set_minor_locator(ticker.NullLocator())
        axis.set_minor_formatter(ticker.NullFormatter())

    def limit_range_for_scale(self, vmin, vmax, minpos):
        return  max(0., vmin), vmax

    class SquareRootTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def transform_non_affine(self, a):
            return np.array(a)**0.5

        def inverted(self):
            return SquareRootScale.InvertedSquareRootTransform()

    class InvertedSquareRootTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def transform(self, a):
            return np.array(a)**2

        def inverted(self):
            return SquareRootScale.SquareRootTransform()

    def get_transform(self):
        return self.SquareRootTransform()


parser = ap.ArgumentParser()
parser.add_argument("-f", "--files", nargs='+', help="hcjfit ASCII files")
args = parser.parse_args()

for file in args.files:
    sch_hist = SCh_histogram(file)
    sch_hist.last_bin()
    sch_hist.plot()


'''
bin sizes:

>>> 24.76216-20.43879
4.323370000000001
>>> 20.43879-16.87026
3.5685300000000026
>>> 16.87026-13.92478
2.945479999999998
>>> 13.92478-11.49358
2.4312000000000005
>>> 11.49358-9.48684
2.006739999999999
>>>
'''