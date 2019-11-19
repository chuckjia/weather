function plotSolutions(projFolderPath, varargin)

p = inputParser;
addRequired(p, 'projFolderPath', @(x) ischar(x) || isstring(x));
addOptional(p, 'solutions', ["theta", "qv", "qc", "qp"]);
addOptional(p, 'timesteps', 0);
addOptional(p, 'filenamePrefix', "")
addOptional(p, 'showGhost', true)
addOptional(p, 'zlimits', 0)
addOptional(p, 'view2d', false)
addOptional(p, 'movie', false)
addOptional(p, 'frameRate', 5)
parse(p, projFolderPath, varargin{:});

removeGhostCells = ~p.Results.showGhost;
zlimits = p.Results.zlimits;
set_zlimits = length(zlimits) == 2;
view2d = p.Results.view2d;
filenamePrefix = p.Results.filenamePrefix;
movieFilename = p.Results.movie;
saveMovie = ischar(movieFilename) | isstring(movieFilename);

solnFolderPath = fullfile(projFolderPath, "solutions/");
[centersX, centersZ] = getCenterCoordMatrices(projFolderPath);

if removeGhostCells
    centersX = centersX(2:end-1, 2:end-1);
    centersZ = centersZ(2:end-1, 2:end-1);
end

frameNo = 0;

for solutionName = p.Results.solutions
    for t = p.Results.timesteps
        filename = sprintf("%s%s_%1.4fs.csv", filenamePrefix, solutionName, t);
        % fprintf("Plotting %s\n", filename)
        filepath = fullfile(solnFolderPath, filename);
        frameNo = frameNo + 1;
        
        sol = readmatrix(filepath);
        
        if removeGhostCells
            sol = sol(2:end-1, 2:end-1);
        end
        
        fig = surf(centersX, centersZ, sol);
        
        title({solutionName, sprintf("t = %1.4fs", t)})
        xlabel("x-axis"); ylabel("z-axis"); zlabel("Solution Value");
        colorbar
        
        if set_zlimits
            zlim(zlimits)
            caxis(zlimits)
        end
        
        if view2d
            view(2)
        end
        
        if saveMovie
            MovieMatrix(frameNo) = getframe(gcf);
            
            if frameNo == 1
                MovieMatrix = repmat(MovieMatrix(1), [1, 1, length(p.Results.solutions)]);
            end
        end
        
        drawnow
    end
end

if saveMovie
    fprintf("Saving movie to file.\n")
    if ~exist("results", "dir")
        mkdir("results");
    end
    v = VideoWriter(fullfile("results", movieFilename));
    v.FrameRate = p.Results.frameRate;
    open(v)
    writeVideo(v, MovieMatrix)
    close(v)
    fprintf("Movie successfully saved as: %s\n", movieFilename)
end