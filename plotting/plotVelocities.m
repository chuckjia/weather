function plotVelocities(folderPath, varargin)
%PLOTVELOCITIES Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addRequired(p, 'folderPath', @(x) ischar(x) || isstring(x));
addOptional(p, 'xrange', [0, 1]);
addOptional(p, 'zrange', [0, 1]);
addOptional(p, 'normalize', 1);
addOptional(p, 'showFullRange', true);
parse(p, folderPath, varargin{:});

xrange = p.Results.xrange;
zrange = p.Results.zrange;
normalize = p.Results.normalize;
showFullRange = p.Results.showFullRange;

x_bl_m = readmatrix(folderPath + "x_bl_m.csv");
z_bl_m = readmatrix(folderPath + "z_bl_m.csv");

x_br_m = readmatrix(folderPath + "x_br_m.csv");
z_br_m = readmatrix(folderPath + "z_br_m.csv");

% x_tl_m = readmatrix(folderPath + "x_tl.csv");
% z_tl_m = readmatrix(folderPath + "z_tl.csv");
% 
% x_tr_m = readmatrix(folderPath + "x_tr.csv");
% z_tr_m = readmatrix(folderPath + "z_tr.csv");

u_b_m = readmatrix(folderPath + "u_b_m.csv");
w_b_m = readmatrix(folderPath + "w_b_m.csv");
x_b_m = 0.5 * (x_bl_m + x_br_m);
z_b_m = 0.5 * (z_bl_m + z_br_m);

if showFullRange
    quiver(x_b_m, z_b_m, u_b_m, w_b_m)
    figure
end

[nrow, ncol] = size(x_b_m);
r0 = round(1 + (nrow - 1) * xrange(1));
r1 = round(1 + (nrow - 1) * xrange(2));
c0 = round(1 + (ncol - 1) * zrange(1));
c1 = round(1 + (ncol - 1) * zrange(2));

x_b_m = x_b_m(r0:r1, c0:c1);
z_b_m = z_b_m(r0:r1, c0:c1);
u_b_m = u_b_m(r0:r1, c0:c1);
w_b_m = w_b_m(r0:r1, c0:c1);

if normalize ~= 0 && normalize ~= false
    vel_size_m = (u_b_m .^ 2 + w_b_m .^ 2) .^ 0.5 ./ normalize;
    u_b_m = u_b_m ./ vel_size_m;
    w_b_m = w_b_m ./ vel_size_m;
end

quiver(x_b_m, z_b_m, u_b_m, w_b_m)



end

