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
    
end