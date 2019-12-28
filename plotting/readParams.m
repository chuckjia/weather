function param = readParams(projFolderPath)
%READPARAMS Summary of this function goes here
%   Detailed explanation goes here

filepath = fullfile(projFolderPath, "results", "param.csv");
tab = readtable(filepath);

params = table2cell(tab(:, 2));
[x0, xf, zf, Nx, Nz, tf, Nt, num_csv] = params{:};
if num_csv > Nt
    num_csv = Nt;
end
param = ParamClass(x0, xf, zf, Nx, Nz, tf, Nt, num_csv);

end
