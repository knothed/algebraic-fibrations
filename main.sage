if False: # optimization
    from fibering import *
else:
    load('fibering.pyx')
load('martelli.sage')
load('virtual_fibering.sage')