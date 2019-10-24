init;

% Inicialização do problema ===============================================
% Palpite inicial para a solução ------------------------------------------
% x0 = [0, 0]';
x0 = rand(2,1);

% Ordem do problema -------------------------------------------------------
n = length(x0); 

% Função objetivo ---------------------------------------------------------
f = @(x) 1/3*(x(1) + 1).^3 + x(2);

% Resolvendo o problema de otimização =====================================
[xOpt, fOpt, nVal, k, alfaValues] = varMetConst(f, x0);

% Apresentando os resultados ==============================================
for i = 1:n
    fprintf(['x',num2str(i),'* = %.4f\n'], xOpt(i))
end

fprintf('f(x*) = %.4f\n', fOpt)
fprintf('Número de avaliações da função objetivo: %d\n', nVal)
fprintf('Número de iterações: %d\n', k)