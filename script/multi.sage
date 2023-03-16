import sys
import os
os.chdir("..")
load("main.sage")

# arg 1: result file

n = 10
args = "14:0 -c"
total = 11515151
geng = "~/nauty/geng"
#geng = "/var/tmp/sage-9.7-current/local/bin/geng"

geng_fibering(n, geng, args, num_queues=40, queue_capacity=5, threads_per_queue=2,
              hyp_check=True, min_isos=1, write_to_file=sys.argv[1], total=total)