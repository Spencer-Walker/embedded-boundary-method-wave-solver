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

print("matlab -nodisplay -nosplash -nodesktop -r \"cd " + cwd + ";run('jsonFileScriptForInvEllipse.m');exit;\"")

os.system("matlab -nodisplay -nosplash -nodesktop -r \"cd " + cwd + ";run('jsonFileScriptForInvEllipse.m');exit;\"")

print("mpirun -np " + str(data["np"]) + \
  " " + path+"/main" + parameters_str)

os.system("mpirun -np " + str(data["np"]) + \
  " "+path+"/main" + parameters_str)
