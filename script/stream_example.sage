import sys
import os
os.chdir("..")
load("main.sage")

# arg 1: result output file

n = 15
input = "~/15.txt"
total = 100

analyze_stream(n, f"cat {input}", num_queues=15, queue_capacity=1, threads_per_queue=10,
               hyp_check=True, min_isos=1, write_to_file=sys.argv[1], total=total)