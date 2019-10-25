% Configurações prévias ===================================================
% Antes de executar a função principal 'main.m', é preciso definir a função
% objetivo, seu gradiente (caso haja), e as restrições (caso haja). Para
% tanto basta editar os arquivos 'funcFile.m', 'gradFile.m', e
% 'constFile.m'

init; % Inicialização do ambiente de trabalho

% Resolvendo o problema de otimização =====================================
[xOpt, fOpt, nVal, k, alfaValues] = varMet();

% Apresentando os resultados ==============================================
for i = 1:length(xOpt)
    fprintf(['x',num2str(i),'* = %.4f\n'], xOpt(i))
end

fprintf('f(x*) = %.4f\n', fOpt)
fprintf('Número de avaliações da função objetivo: %d\n', nVal)
fprintf('Número de iterações: %d\n', k)

% Notas ===================================================================
% A seguir é introduzida a forma completa da função 'varMet'
% [xOpt, fOpt, nVal, k, alfaValues] = varMet(x0, theta, rp0, rpInc, tol, h)

% Para os parâmetros não informados pelo usuário serão adotados os valores
% padrão
% x0 = [0 0 0 ... 0]
% theta = 1
% rp0 = 1
% rpInc = 10
% tol = 1e-5
% h = 1e-10