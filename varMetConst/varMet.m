% =========================================================================
% =========================================================================
% Implementação do Método da Variável Métrica =============================
% =========================================================================
% =========================================================================
function [xOpt, fOpt, nVal, k, alfaValues] = varMet(x0, theta, rp0, ...
    rpInc, tol, h)

% Verificação do palpite inicial fornecido pelo usuário
n = getOrder;
if nargin < 1 || isempty(x0)
    x0 = zeros(n,1);
else
    % Verificando a consistência entre x0 e a função objetivo
    if n ~= length(x0)
        error('O tamanho de x0 deve ser igual à ordem do problema.');
    end
end

% Verificando entradas da função e definindo valores padrão
if nargin < 2 || isempty(theta)
    theta = 1;
end
if nargin < 3 || isempty(rp0)
    rp0 = 1;
end
if nargin < 4 || isempty(rpInc)
    rpInc = 10;
end
if nargin < 5 || isempty(tol)
    tol = 1e-5;
end
if nargin < 6 || isempty(h)
    h = 1e-10;
end

% Verificando se o gradiente analítico foi definido
f = @(x) fObjConst(x);
dfTest = gradFile(x0);
if isempty(dfTest)
    % Gradiente numérico
    if nargin < 4 || isempty(h)
        df = @(x) grad(f,x);
    else
        df = @(x) grad(f,x,h);
    end
    numGrad = 1;
else
    % Gradiente analítico
    df = @(x) gradFile(x);
    numGrad = 0;
end

% Definição das entradas padrão
if nargin < 2 || isempty(tol)
    tol = 1e-5;
end

if nargin < 3 || isempty(theta)
    theta = 1;
end

% Variáveis para controle de execução
n = length(x0);
alfaValues = zeros(1,1);
k = 1;
nVal = 0;
H = eye(n);
global rp
rp = rp0;

% Implementação do processo de otimização
while 1
    % Reduzir a dimensão do problema de otimização
    g = @(alfa) f(x0 - alfa*H*df(x0));
    
    % Resolver o problema de otimização uni-dimensional
    [alfaOpt,~,nVal1] = aureaSec(g,-1,1,1e-4);
    
    % Atualizar a solução ótima
    x = x0 - alfaOpt*H*df(x0);
    
    % Armazenar dados de execução
    alfaValues(k) = alfaOpt;
    if numGrad
        nVal = nVal + nVal1 + 1;
    else
        nVal = nVal + nVal1 + 2*n;
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
    rp = rp*rpInc;
end

xOpt = x;
fOpt = f(xOpt);

% =========================================================================
% =========================================================================
% Obtenção numérica do gradiente ==========================================
% =========================================================================
% =========================================================================
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

% =========================================================================
% =========================================================================
% Implementação do Método da Seção Áurea ==================================
% =========================================================================
% =========================================================================
function [xOpt, fOpt, k] = aureaSec(f,a,b,tol)
tal = 0.618;

if nargin < 4
    tol = 1e-8;
end

alfa = a + (1 - tal)*(b - a);
beta = a + tal*(b - a);
fAlfa = f(alfa);
fBeta = f(beta);

k = 1;

while abs(a-b) > tol
    if fBeta < fAlfa
        a = alfa;
        alfa = beta;
        fAlfa = fBeta;
        beta = a + tal*(b - a);
        fBeta = f(beta);
    elseif fAlfa <= fBeta
        b = beta;
        beta = alfa;
        fBeta = fAlfa;
        alfa = a + (1 - tal)*(b - a);
        fAlfa = f(alfa);
    end
    k = k + 1;
end

xOpt = (alfa+beta)/2;
fOpt = f(xOpt);

% =========================================================================
% =========================================================================
% Função para cálculo das penalidades =====================================
% =========================================================================
% =========================================================================
function P = penalty(x)

% Determinação dos valores atribuídos a cada uma das restrições
[g, h] = constFile(x);

% Determinação das penalidades
Pg = sum(max(0,g).^2);
Ph = sum(h.^2);

% Cálculo das penalidades associadas tanto
P = Pg + Ph;

% =========================================================================
% =========================================================================
% Função objetivo penalizada ==============================================
% =========================================================================
% =========================================================================
function F = fObjConst(x)

global rp

% Computação da função objetivo original
fObj = funcFile(x);

% Computação das penalidades
P = penalty(x);

% Computação da função objetivo penalizada
F = fObj + rp*P;

% =========================================================================
% =========================================================================
% Função para verificação da ordem do problema com base na função objetivo 
% =========================================================================
% =========================================================================
function n = getOrder

x = zeros(1,100);
f0 = funcFile(x);
for i = 1:100
    x(i) = 1;
    f = funcFile(x);
    if f0 == f
        break
    end
    x(i) = 0;
end
n = i - 1;