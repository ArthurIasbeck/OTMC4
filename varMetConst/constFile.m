function [g, h] = constFile(x)

% Cada linha da matriz 'g' representa uma restrição de desigualdade
g = [-x(1) + 1;
	 -x(2)];

% Cada linha da matriz 'h' representa uma restrição de igualdade
h = [];