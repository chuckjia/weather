closeall

T0 = 288.15;
p0 = 101325.0;
g = 9.8;
c_p = 1004.68506;
R0 = 8.31582991;

c1_p_e = g / (c_p * T0);
M_p_e = 0.02896968;
c2_p_e = c_p * M_p_e / R0;


z = 0:100:16e3;
p_e = p0 * (1 - c1_p_e * z) .^ c2_p_e;

min(p_e)

plot(z, p_e)