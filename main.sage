if False: # optimization
    from bridge import *
else:
    load('bridge.pyx')
load('martelli.sage')
load('virtual_fibering.sage')