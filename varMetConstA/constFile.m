function P = constFile(x)

g = [-x(1) + 1;
	 -x(2)];
 
h = [];
 
Pg = sum(max(0,g).^2);
Ph = sum(h.^2);

P = Pg + Ph;