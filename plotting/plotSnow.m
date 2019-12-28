function plotSnow(projFolderPath, varargin)

p = inputParser;
addRequired(p, 'projFolderPath', @(x) ischar(x) || isstring(x))
addOptional(p, 'timesteps', 0)
addOptional(p, 'zlimits', 0)
addOptional(p, 'movie', false)
addOptional(p, 'contour', false)
addOptional(p, 'frameRate', 5)
addOptional(p, 'zoom', 1)
parse(p, projFolderPath, varargin{:});

removeGhostCells = true;
zlimits = p.Results.zlimits;
set_zlimits = length(zlimits) == 2;
movieFilename = p.Results.movie;
saveMovie = ischar(movieFilename) | isstring(movieFilename);
contour = p.Results.contour;
zoomRate = p.Results.zoom;

solnFolderPath = fullfile(projFolderPath, "solutions/");
[centersX, centersZ] = getCenterCoordMatrices(projFolderPath);

if removeGhostCells
    centersX = centersX(2:end-1, 2:end-1);
    centersZ = centersZ(2:end-1, 2:end-1);
end

ncol = size(centersX, 2);
top_j = round(ncol * zoomRate);
centersX = centersX(:, 1:top_j);
centersZ = centersZ(:, 1:top_j);

frameNo = 0;

for t = p.Results.timesteps
    T = readSolution(solnFolderPath, "theta", t);
    qp = readSolution(solnFolderPath, "qp", t);
    
    qs = (1 - alpha_fcn(T)) .* qp;
    
    frameNo = frameNo + 1;
    
    if removeGhostCells
        qs = qs(2:end-1, 2:end-1);
    end
    
    qs = qs(:, 1:top_j);
    qs = qs .* (qs > 0);
    
    if contour
        contourf(centersX, centersZ, qs);
    else
        surf(centersX, centersZ, qs);
    end
    %     surf(centersX, centersZ, T);
    view(2)
    
    title({"Snow", sprintf("t = %1.4fs", t)})
    xlabel("x-axis"); ylabel("z-axis"); zlabel("Solution Value");
    colorbar
    
    if set_zlimits
        zlim(zlimits)
        caxis(zlimits)
    end
    
    if saveMovie
        if frameNo >= 2
            MovieMatrix(frameNo - 1) = getframe(gcf);
        end
        if frameNo == 2
            MovieMatrix = repmat(MovieMatrix(1), [1, 1, length(p.Results.timesteps) - 1]);
        end
    end
    
    drawnow
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

end

%%

function alpha = alpha_fcn(T)

T_w = 273.15;
T_i = 263.15;

alpha = (T > T_i) .* (T < T_w) .* (T - T_i) ./ (T_w - T_i) + (T >= T_w);

end

function sol = readSolution(solnFolderPath, solutionName, t)

filename = sprintf("%s_%1.4fs.csv", solutionName, t);
filepath = fullfile(solnFolderPath, filename);
sol = readmatrix(filepath);

end
