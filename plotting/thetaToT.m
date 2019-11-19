function T = thetaToT(theta, z)
%THETATOT Summary of this function goes here
%   Detailed explanation goes here

p_oo = 1e5;
R_d = 287.058;
c_p = 1004.68506;
c = R_d / c_p;

T = theta .* (p_e_fcn(z) ./ p_oo) .^ c;

end


function p = p_e_fcn( z)

T_0 = 288.15;
g = 9.8;
c_p = 1004.68506;
R_0 = 8.31582991;
c_1 = g / (c_p * T_0);
M = 0.02896968;
c_2 = c_p * M / R_0;
p_0 = 101325.0;

p = p_0 * (1 - c_1 .* z) .^ c_2;  % Eq.201910152104

end
