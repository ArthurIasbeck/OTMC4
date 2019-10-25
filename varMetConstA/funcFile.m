function f = funcFile(x)

global rp

fObj = 1/3*(x(1) + 1)^3 + x(2);

g1 = -x(1) + 1;
g2 = -x(2);

f = fObj + rp*(max(0,g1)^2 + max(0,g2))^2;
