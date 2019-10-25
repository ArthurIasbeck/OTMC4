init;

% Inicializar o problema ==================================================
% Ordem do problema -------------------------------------------------------
n = 2; 

% Palpite inicial para a solução ------------------------------------------
x0 = [0, 0]';

% Função objetivo ---------------------------------------------------------
% Definida a partir de um arquivo 
f = @funcFile;

% Definida diretamente 
% f = @(x) 1/3*(x(1) + 1).^3 + x(2);

% Gradiente analítico -----------------------------------------------------
% Definido a partir de um arquivo 
% df = @gradFile;

% Definido diretamente
% df = @(x) [8500*x(1) - 1000*x(2) - 2500*x(3) - 1000
%            3000*x(2) - 1000*x(1) - 500*x(3) - 2000
%            11500*x(3) - 500*x(2) - 2500*x(1) - 3000];

% Resolvendo o problema de otimização =====================================

% Tolerância para o critério de parada (Norma Euclidiana)
% tol = 1e-5;
% Valor de theta a ser utilizado no MVM
% theta = 1;
% Variação infinitesimal de x para a computação numérica do gradiente
% h = 1e-10;

% [xOpt, fOpt, nVal, k, alfaValues] = varMet(f, x0, df, tol, theta, h);
% [xOpt, fOpt, nVal, k, alfaValues] = varMet(f, x0, df, tol, theta);
% [xOpt, fOpt, nVal, k, alfaValues] = varMet(f, x0, df, tol);
% [xOpt, fOpt, nVal, k, alfaValues] = varMet(f, x0, df);
[xOpt, fOpt, nVal, k, alfaValues] = varMet(f, x0);

% Apresentando os resultados ==============================================
for i = 1:n
    fprintf(['x',num2str(i),'* = %.4f\n'], xOpt(i))
end

fprintf('f(x*) = %.4f\n', fOpt)
fprintf('Número de avaliações da função objetivo: %d\n', nVal)
fprintf('Número de iterações: %d\n', k)