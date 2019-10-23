function [xOpt, fOpt, nVal, k] = varMet(x0, f, df, tol, theta)

% Verificação de erros na entrada
if nargin < 2
    disp('Insira ao menos x0 e f');
    return;
end

% Verificação da entrada do gradiente
if nargin < 3 || isempty(df)
    df = @(x) grad(f,x);
end