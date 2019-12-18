
%% Plot solution

closeall

projFolderPath = "/Users/chuckjia/Documents/Workspace/Git/weather/";

solnameArray = ["theta", "qv", "qc", "qp"];
<<<<<<< Updated upstream
zlimitsArray = {[200, 380], [0, 12], [0, 0.035], [0, 1e-3]};
=======
zlimitsArray = {[-50, 100], [0, 12], [0, 0.035], [0, 0.05]};
>>>>>>> Stashed changes
zlimitsMap = containers.Map(solnameArray, zlimitsArray);

solname = "qv";
zlimits = zlimitsMap(solname);
Dt = 0.5;
<<<<<<< Updated upstream
Nt = 100;
timesteps = 0:Dt:Nt * Dt;
=======
Nt = 200;
>>>>>>> Stashed changes


timesteps = 0:Dt:Nt*Dt;


plotSolutions( ...
    projFolderPath, ... 
    'solutions', solname, ... 
    'timesteps', timesteps, ...
    'showGhost', false, ...
<<<<<<< Updated upstream
    'zlimits', zlimits, ...
    'view2d', true, ...
=======
    'zlimits', false, ...
    'view2d', true, ...
    'movie', false, ...
    'KtoC', solname == "theta", ...
    'contour', false, ...
    'snow', true);

%% Plot snow

closeall

projFolderPath = "/Users/chuckjia/Documents/Workspace/Git/weather/";

Dt = 0.5;
Nt = 200;

timesteps = 0:Dt:Nt*Dt;

plotSnow( ...
    projFolderPath, ... 
    'timesteps', timesteps, ...
    'zlimits', false, ...
    'contour', true, ...
>>>>>>> Stashed changes
    'movie', false);

