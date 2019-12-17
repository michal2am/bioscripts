import matplotlib.pyplot as plt
import math
import numpy as np
import yaml
from scalcs import popen
from scalcs import scplotlib as scpl
from scalcs import cjumps


def ratio_maxPopen(mec1, mec2, tres):
    emaxPopen1 = popen.maxPopen(mec1, tres)
    emaxPopen2 = popen.maxPopen(mec2, tres)
    return emaxPopen1 / emaxPopen2


tres = 30e-6
filename = "Flip_a1GlyRWT_Burz2004.yaml"
stream = open(filename, 'r')
mec0 = yaml.load(stream)
stream.close()
mec0.printout()
stream = open(filename, 'r')
mec = yaml.load(stream)