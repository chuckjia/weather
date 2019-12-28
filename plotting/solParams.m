
param = readParams(projFolderPath)
Nt = param.Nt;
tf = param.tf;
Dt = tf / Nt;
numcsv = param.num_csv;
default_timesteps = 0:(tf / numcsv):tf;