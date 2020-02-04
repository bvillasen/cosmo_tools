import subprocess

for nSnap in range(170):
  command = 'mpirun -n 8 --map-by ppr:1:node python get_phase_diagram_2048.py {0}'.format(nSnap)
  process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  for line in process.stdout:
    print line
  process.wait()
  print process.returncode