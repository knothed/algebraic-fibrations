import sys
import os
os.chdir("..")
load("main.sage")

# arg 1: result file

n = 11
args = "18:0 -c"
total = 1000000000
geng = "~/nauty/geng"
#geng = "/var/tmp/sage-9.7-current/local/bin/geng"

analyze_geng_stream(n, geng, args, num_queues=30, queue_capacity=3, threads_per_queue=4,
                    hyp_check=True, min_isos=1, write_to_file=sys.argv[1], total=total)
