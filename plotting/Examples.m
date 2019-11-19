%% Plot mesh

closeall

projFolderPath = "/Users/chuckjia/Documents/Workspace/Git/weather/";
resultFolderPath = projFolderPath + "results/";

plotMesh(resultFolderPath, 'plotCellCenter', true, 'markGhostCell', true)


%% Plot velocities

closeall

projFolderPath = "/Users/chuckjia/Documents/Workspace/Git/weather/";
resultFolderPath = projFolderPath + "results/";

plotVelocities(resultFolderPath)


%% Plot exact solution

closeall

projFolderPath = "/Users/chuckjia/Documents/Workspace/Git/weather/";
plotSolutions(projFolderPath, 'solutions', "theta", 'timesteps', 0.25, 'prefix', "exact_");

