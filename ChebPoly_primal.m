%%
% Given a multivariate monomial m_k and a domain Omega, 
% solves the optimization program
%
% minimize_{p,c} c 
% s.to: p is a polynomial of degree strictly less than deg(m_k)
%  and: (c >= m_k - p and c >= p - m_k) on Omega
%
% This is done by solving a semidefinite program (SDP) truncated at level t
% Note: to be executed, ChebPoly_primal requires GloptiPoly 3 

% Inputs:
% k = (k_1,...,k_n): the multidegree of the monomial
% Omega: the domain, which can be:
%  -> H: the hypercube
%  -> B: the euclidean ball
%  -> C: the cross-polytope
%  -> S: the simplex 
% t: the truncation level

% Outputs:
% err:      the optimal value of the SDP 
%           (this is an upper bound for the true error of best approximation)
% chebpoly: the polynomial m_k - p, p being a miiminizer of the SDP 
%           (this is a candidate Chebyshev polynomial)

% Written in May 2024
% Send comments to simon.foucart@centraliens.net

%%
function [err,chebpoly] = ChebPoly_primal(k,Omega,t)

%% the parameters of the problem
d = length(k);  % number of variables
n = sum(k);     % (total) degree of the monomial

%% the optimization variables
c = sdpvar(1);
x = sdpvar(d,1);
% the approximant
[p, coef_p] = polynomial([x],n-1);
% the monomial
mon = x(1)^k(1);
for i = 2:d
    mon = mon * x(i)^k(i);
end

%% the objective function
obj = c;

%% the polynomial descriptors of the domain Omega
if Omega == 'H'
    for h = 1:d
        g{h} = 1 - x(h)^2;
    end
end
if Omega == 'B'
    g{1} = 1;
    for h = 1:d
        g{1} = g{1} - x(h)^2;
    end
end
if Omega == 'C'
    signs = 2*(dec2bin(2^d-1:-1:0)-'0')-1; % produce all sign vectors of lenghth d
    for h = 1:size(signs,1)
        g{h} = 1;
        for j = 1:d
            g{h} = g{h} - signs(h,j)*x(j);
        end
    end
end
if Omega == 'S'
    for h = 1:d
        g{h} = x(h);
    end
    g{d+1} = 1;
    for i=1:d
        g{d+1} = g{d+1}-x(i);
    end
end

%% the constraints (involving q_plus and q_minus)
% the SOS multipliers
const_p = [];
for h = 1:length(g)
    [aux1, aux2] = polynomial([x],t-2);
    q_p{h} = aux1;
    coef_q_p{h} = aux2;
    const_p = [const_p; sos(q_p{h})];
end
const_p_aux = c - mon + p;
for h = 1:length(g)
    const_p_aux = const_p_aux - q_p{h}*g{h};
end
const_p = [const_p; sos(const_p_aux)];
const_m = [];
for h = 1:length(g)
    [aux1, aux2] = polynomial([x],t-2);
    q_m{h} = aux1;
    coef_q_m{h} = aux2;
    const_m = [const_m; sos(q_m{h})];
end
const_m_aux = c + mon - p;
for h = 1:length(g)
    const_m_aux = const_m_aux - q_m{h}*g{h};
end
const_m = [const_m; sos(const_m_aux)];
const = [const_p ; const_m];

%%
coeff_array = [c; coef_p]; 
for h = 1:length(coef_q_p)
    coeff_array = [coeff_array; coef_q_p{h}];
end
for h = 1:length(coef_q_m)
    coeff_array = [coeff_array; coef_q_m{h}];
end

%% solver options
SDPsolver = 'sedumi';
options = getSolverParams(SDPsolver);
if (~isempty(options))
    mset(options);
end

%% solving the optimization program
[~,~,~,res] = solvesos(const,obj,options,coeff_array);
fprintf('Residuals %f \n',res);

%% return the outputs
err = double(obj);
aux = ( fliplr(flipud(genpow(d,n))) == repmat(k,nchoosek(n+d-1,d-1),1) );
coef_mon = ( sum(aux,2) == d );
coef_chebpoly = num2cell([double(coef_p); coef_mon]);
chebpoly = flipud([sdisplay(monolist(x,n)), coef_chebpoly]);

end