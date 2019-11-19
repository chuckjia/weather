function plotMesh(folderPath, varargin)

p = inputParser;
addRequired(p, 'folderPath', @(x) ischar(x) || isstring(x));
addOptional(p, 'markGhostCell', true, @islogical);
addOptional(p, 'plotCellCenter', true, @islogical);
parse(p, folderPath, varargin{:});

% Read matrices from saved CSV files
x_bl_m = readmatrix(folderPath + "x_bl_m.csv");
z_bl_m = readmatrix(folderPath + "z_bl_m.csv");

x_br_m = readmatrix(folderPath + "x_br_m.csv");
z_br_m = readmatrix(folderPath + "z_br_m.csv");

x_tl_m = readmatrix(folderPath + "x_tl_m.csv");
z_tl_m = readmatrix(folderPath + "z_tl_m.csv");

x_tr_m = readmatrix(folderPath + "x_tr_m.csv");
z_tr_m = readmatrix(folderPath + "z_tr_m.csv");

% Generate grid matrices
gridX = [x_bl_m x_tl_m(:, end)];
gridX = [gridX; [x_br_m(end, :) x_tr_m(end, end)]];

gridZ = [z_bl_m z_tl_m(:, end)];
gridZ = [gridZ; [z_br_m(end, :) z_tr_m(end, end)]];

mesh(gridX, gridZ, zeros(size(gridX)));
view([0, 0, 1])
hold on

if p.Results.markGhostCell
    % Highlight ghost cell
    ghostCellColor = [218, 227, 243] ./ 255;  % Color used for the ghost cells. Light blue = [218, 227, 243]
    [nrow, ncol] = size(gridX);
    % Shade all ghost cells along the four sides
    mesh(gridX(1:2, :), gridZ(1:2, :), zeros(2, ncol), 'FaceColor', ghostCellColor);  % Left ghost cells
    mesh(gridX((end - 1):end, :), gridZ((end - 1):end, :), zeros(2, ncol), 'FaceColor', ghostCellColor);  % Right ghost cells
    mesh(gridX(:, 1:2), gridZ(:, 1:2), zeros(nrow, 2), 'FaceColor', ghostCellColor);  % Bottom ghost cells
    mesh(gridX(:, (ncol - 1):ncol), gridZ(:, (ncol - 1):ncol), zeros(nrow, 2), 'FaceColor', ghostCellColor);  % Top ghost cells
end

if p.Results.plotCellCenter
    % Graph cell centers. Need to plot *after* highlighting ghost cells
    centersX = readmatrix(folderPath + "x_ct_m.csv");
    centersZ = readmatrix(folderPath + "z_ct_m.csv");
    plot3(centersX, centersZ, zeros(size(centersX)), '.', 'MarkerSize', 10);
end

end

