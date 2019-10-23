init;

addpath('..');

% Parâmetros pra execução do algorítmo
n = 3;
numGrad = 0; 
x0 = [0, 0, 0]';
tol = 1e-4;
theta = 0;

% Definição da função objetivo
f = @(x) 4250*x(1).^2 - 1000*x(1)*x(2) - 2500*x(1)*x(3) - 1000*x(1) + ...
    1500*x(2).^2 - 500*x(2)*x(3) - 2000*x(2) + 5750*x(3).^2 - 3000*x(3);

if ~numGrad
    % Definição analítica do gradiente
    df = @(x) [8500*x(1) - 1000*x(2) - 2500*x(3) - 1000
               3000*x(2) - 1000*x(1) - 500*x(3) - 2000
               11500*x(3) - 500*x(2) - 2500*x(1) - 3000];
else 
    % Definição numérica do gradiente
    h = 1e-10;
    df = @(x) [(f([x(1) + h, x(2), x(3)]) - f([x(1), x(2), x(3)]))/h
               (f([x(1), x(2) + h, x(3)]) - f([x(1), x(2), x(3)]))/h
               (f([x(1), x(2), x(3) + h]) - f([x(1), x(2), x(3)]))/h];
end

% Variáveis para controle de execução
alfaOptValues = zeros(1,1);
k = 1;
nVal = 0;
H = eye(n);

while 1 
    % Reduzir a dimensão do problema de otimização
    g = @(alfa) f(x0 - alfa*H*df(x0));
    
    % Resolver o problema de otimização uni-dimensional
    [alfaOpt,~,nVal1] = aureaSec(g,-1,1,1e-4);
    
    % Atualizar a solução ótima
    x = x0 - alfaOpt*H*df(x0);
    
    % Armazenar dados de execução 
    alfaOptValues(k) = alfaOpt;
    if ~numGrad 
        nVal = nVal + nVal1 + 1; 
    else
        nVal = nVal + nVal1 + 6;
    end
    
    % OBS : Lembre-se que é necessária a computação do gradiente para
    % atualização de x. No entanto, se estivermos utilizando a aproximação
    % numérica para o gradiente, a computação do mesmo levará, neste caso a
    % 6 avaliações da função objetivo.
   
    % Verificar a condição de parada
    cp = norm(x - x0);
    if cp < tol
        break;
    end
    
    % Atualização de H (aproximação para a inversa da Matriz Hessiana)
    p = x - x0;
    y = df(x) - df(x0);
    sigma = p'*y;
    tal = y'*H*y;
    D = ((sigma + theta*tal)/sigma^2)*(p*p') ...
        + ((theta - 1)/tal)*(H*y)*(H*y)' ...
        - (theta/sigma)*(H*y*p' + p*(H*y)');
    
    H = H + D;

    % Atualizar variáveis para a próxima iteração
    x0 = x;
    k = k + 1;
end

xOpt = x;
fOpt = f(xOpt);

for i = 1:length(x)
    fprintf(['x',num2str(i),'* = %.4f\n'], xOpt(i));
end

fprintf('f(x*) = %.4f\n', fOpt);
fprintf('Número de avaliações da função objetivo: %d\n', nVal);
fprintf('Número de iterações: %d\n', k);