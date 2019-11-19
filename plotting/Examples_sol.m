
%% Plot solution

closeall

projFolderPath = "/Users/chuckjia/Documents/Workspace/Git/weather/";

solnameArray = ["theta", "qv", "qc", "qp"];
zlimitsArray = {[200, 380], [0, 12], [0, 0.035], [0, 1e-3]};
zlimitsMap = containers.Map(solnameArray, zlimitsArray);

solname = "qv";
zlimits = zlimitsMap(solname);
Dt = 0.5;
Nt = 100;
timesteps = 0:Dt:Nt * Dt;

plotSolutions(projFolderPath, 'solutions', solname, 'timesteps', timesteps, ...
    'showGhost', false, ...
    'zlimits', zlimits, ...
    'view2d', true, ...
    'movie', false);
