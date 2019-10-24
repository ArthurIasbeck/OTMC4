function [g, h] = constFile(x)

g = [-x(1) + 1;
    -x(2)];

h = [];