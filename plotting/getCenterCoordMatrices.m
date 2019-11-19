function [centersX, centersZ] = getCenterCoordMatrices(proj_folder_path)
%getCenterCoordMatrices Summary of this function goes here
%   Detailed explanation goes here

resultsFolderPath = fullfile(proj_folder_path, "results");
centersX = readmatrix(fullfile(resultsFolderPath, "x_ct_m.csv"));
centersZ = readmatrix(fullfile(resultsFolderPath, "z_ct_m.csv"));

end

