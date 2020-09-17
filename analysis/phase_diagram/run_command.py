import subprocess


for nSnap in range(170):
  command = 'mpirun -n 4 --map-by ppr:1:node --oversubscribe python get_phase_diagram_512.py {0}'.format(nSnap)
  print('Command: {0}'.format( command ))
  process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  for line in process.stdout:
    print(line)
  process.wait()
  print(process.returncode)
# for n_cosmo in range(1,4):
#   for nSnap in range(1,16):
#     if n_cosmo == 1 and nSnap < 13: continue
#     command = 'mpirun -n 16 --map-by ppr:4:node --oversubscribe python get_phase_diagram_2048_cosmo.py {0} {1}'.format(nSnap, n_cosmo)
#     print('Command: {0}'.format( command ))
#     process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
#     for line in process.stdout:
#       print(line)
#     process.wait()
#     print(process.returncode)