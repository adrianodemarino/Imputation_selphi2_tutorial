"""How the test lists were made: people with no parent or child in 1000 Genomes, per super-population.
test_samples.txt   100 people (20 per super-population): array imputation, steps 2-4
lowcov_samples.txt  20 of them (4 per super-population): low-coverage reads, step 5
usage: pick_test_samples.py data/20130606_g1k_3202_samples_ped_population.txt"""
import random
import sys

rows = [l.split() for l in open(sys.argv[1])][1:]
parents = {r[2] for r in rows} | {r[3] for r in rows}
free = [r for r in rows if r[2] == "0" and r[3] == "0" and r[1] not in parents]
rng = random.Random(22)
pops = ["AFR", "AMR", "EAS", "EUR", "SAS"]
pool = {sp: sorted(r[1] for r in free if r[6] == sp) for sp in pops}
low = {sp: sorted(rng.sample(pool[sp], 4)) for sp in pops}
more = {sp: rng.sample([s for s in pool[sp] if s not in low[sp]], 16) for sp in pops}
with open("lowcov_samples.txt", "w") as f:
    f.writelines(f"{s}\n" for sp in pops for s in low[sp])
with open("test_samples.txt", "w") as f:
    f.writelines(f"{s}\n" for sp in pops for s in sorted(low[sp] + more[sp]))
