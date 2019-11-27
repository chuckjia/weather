
%% Plot solution

closeall

projFolderPath = "/Users/chuckjia/Documents/Workspace/Git/weather/";

solnameArray = ["theta", "qv", "qc", "qp"];
zlimitsArray = {[200, 380], [0, 12], [0, 0.035], [0, 0.1]};
zlimitsMap = containers.Map(solnameArray, zlimitsArray);

solname = "theta";
zlimits = zlimitsMap(solname);
Dt = 1;
Nt = 100;
timesteps = Nt * Dt;

plotSolutions(projFolderPath, 'solutions', solname, 'timesteps', timesteps, ...
    'showGhost', false, ...
    'zlimits', 0, ...
    'view2d', false, ...
    'movie', false);
