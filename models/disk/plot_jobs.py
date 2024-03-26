from subprocess import Popen
import time
import sys
from glob import glob

# check arguments
if len(sys.argv) < 2:
  print("ERROR: number of cores missing usage is")
  print("python %s NCORES" % sys.argv[0])
  sys.exit()

ncores = int(sys.argv[1])  # number of cores from argument
nruns = len(glob("output_ngas*.dat")) - 1  # number of total models to run
cmd = "python plot.py"  # executable to run
tdump = 10  # print info to screen every tdump seconds

is_running = []
ps = []

# store initial time for screen dumping
tref = time.time()

# loop until there are process running or to run
while len(ps) < nruns or any(is_running):
  is_running = [p.poll() is None for p in ps]
  if sum(is_running) < ncores and len(is_running) < nruns:
    iproc = len(ps)  # id of the last process
    print("starting run number", iproc)

    cmdi = cmd + " " + str(iproc+1)

    # run program
    ps.append(Popen(cmdi.split(), universal_newlines=True, stdout=open("stdout_%03d.out" % iproc, "w"),
            stderr=open("stderr_%03d.out" % iproc, "w")))

  # sleep before the next polling to avoid too many calls
  time.sleep(0.2)

  # print to screen every tdump
  if time.time() - tref > tdump:
    tref = time.time()
    print("running: %d, processes: %d, models: %d, cores: %d" % (sum(is_running), len(ps), nruns, ncores))

print("running: %d, processes: %d, models: %d, cores: %d" % (sum(is_running), len(ps), nruns, ncores))
print("done!")
