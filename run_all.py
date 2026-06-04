import subprocess
import os

BASE = os.path.dirname(os.path.abspath(__file__))
POTENTIAL_DIRS = ["TaC_SNAP", "TaC_ACE", "TaC_MEAM"]
NP = 64


def lmp(script, cwd):
    out = script.replace(".in", "") + ".out"
    with open(os.path.join(cwd, out), "w") as f:
        subprocess.run(["mpirun", "-np", str(NP), "lmp", "-in", script],
                       cwd=cwd, stdout=f, check=True)


def run(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


for pot in POTENTIAL_DIRS:
    def d(sub): return os.path.join(BASE, pot, sub)

    lmp("lat_min.in",                  d("LatticeConstant"))
    lmp("in.elastic",                  d("ElasticConstants"))
    run(["bash", "gsf.sh"],            d("GsfCurves"))
    run(["python", "plot_gsf.py"],    d("GsfCurves"))
    lmp("in.dislocation",              d("EdgeRelaxation"))
    run(["python", "disregistry.py"], d("EdgeRelaxation"))
    lmp("applystress.in",              d("PeierlsStress"))
    run(["python", "plot.py"],        d("PeierlsStress"))
