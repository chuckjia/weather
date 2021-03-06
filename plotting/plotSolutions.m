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
addOptional(p, 'KtoC', false)
addOptional(p, 'contour', false)
addOptional(p, 'snow', false)
addOptional(p, 'zoom', 1)
parse(p, projFolderPath, varargin{:});

removeGhostCells = ~p.Results.showGhost;
zlimits = p.Results.zlimits;
set_zlimits = length(zlimits) == 2;
view2d = p.Results.view2d;
filenamePrefix = p.Results.filenamePrefix;
movieFilename = p.Results.movie;
saveMovie = ischar(movieFilename) | isstring(movieFilename);
KtoC = p.Results.KtoC;
contour = p.Results.contour;
snow = p.Results.snow;
zoomRate = p.Results.zoom;

solnFolderPath = fullfile(projFolderPath, "solutions/");
[centersX, centersZ] = getCenterCoordMatrices(projFolderPath);

if removeGhostCells
    centersX = centersX(2:end-1, 2:end-1);
    centersZ = centersZ(2:end-1, 2:end-1);
end

ncol = size(centersX, 2);
bott_j = 2;
top_j = round(ncol * zoomRate);
centersX = centersX(:, bott_j:top_j);
centersZ = centersZ(:, bott_j:top_j);

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
        
        if KtoC
            sol = sol - 273.15;
        end
        
        sol = sol(:, bott_j:top_j);
        if contour
            fig = contourf(centersX, centersZ, sol);
        elseif snow
            fig = surf(centersX, centersZ, 1.*(sol < 0));
            view(2)
        else
            fig = surf(centersX, centersZ, sol);
        end
        
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