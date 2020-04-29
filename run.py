import os
import json

path = os.path.dirname(os.path.realpath(__file__))
cwd = os.getcwd()

parameters_str = ""

with open(cwd+"/input.json", 'r') as f:
  data = json.load(f)
  for key,val in data.items():
    if "np" not in key:
      parameters_str += " -" + key + " " + str(val)  


os.system("matlab -nodisplay -nosplash -nodesktop -r \"cd " + cwd + ";run('jsonFileScriptForCircle.m');exit;\"")

for np in range(1, data["np"]+1):
    os.system("mpirun -np " + str(np) + \
              " "+path+"/main" + parameters_str)
