if False: # optimization
    from fibering import *
else:
    load('fibering.pyx')
load('thesis.sage')
load('results.sage')