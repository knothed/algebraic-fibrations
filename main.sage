if False: # optimization
    from fibering import *
else:
    load('fibering.pyx')
load('more.sage')
load('results/results.sage')