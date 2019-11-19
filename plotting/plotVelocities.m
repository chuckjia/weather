function plotVelocities(folderPath)
%PLOTVELOCITIES Summary of this function goes here
%   Detailed explanation goes here

x_bl_m = readmatrix(folderPath + "x_bl.csv");
z_bl_m = readmatrix(folderPath + "z_bl.csv");

x_br_m = readmatrix(folderPath + "x_br.csv");
z_br_m = readmatrix(folderPath + "z_br.csv");

% x_tl_m = readmatrix(folderPath + "x_tl.csv");
% z_tl_m = readmatrix(folderPath + "z_tl.csv");
% 
% x_tr_m = readmatrix(folderPath + "x_tr.csv");
% z_tr_m = readmatrix(folderPath + "z_tr.csv");

u_b_m = readmatrix(folderPath + "u_b_m.csv");
w_b_m = readmatrix(folderPath + "w_b_m.csv");
x_b_m = 0.5 * (x_bl_m + x_br_m);
z_b_m = 0.5 * (z_bl_m + z_br_m);

quiver(x_b_m, z_b_m, u_b_m, w_b_m)


end

