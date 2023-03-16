import os
os.chdir("..")
load("main.sage")

g = löbell_graph(7)
res = all_legal_orbits(g, num_threads = 80)

print()
print("Legal orbits of Löb(7):")
print(res)
