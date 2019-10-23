function df = grad(f, x, h)

% Definindo valores padrão
if nargin < 3
    h = 1e-10;
end

n = length(x);

df = zeros(n,1);
for i = 1:n
    dx = x;
    dx(i) = dx(i) + h;
    df(i) = (f(dx) - f(x))/h;
end