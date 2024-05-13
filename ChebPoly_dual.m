%%
% Given the monomial m_k, we aim at solving the optimization program
%
% Given a multivariate monomial m_k and a domain Omega, 
% solves the optimization program
%
% maximize_{y+,y-} y+_k - y-_k 
% s.to: ( y+_l - y-_l = 0 for |l|<n, y+_0 + y-_0 = 1)
%  and: Hank(y+), Hank(y-), Hank(G_h y+), Hank(G_h y-)  are all PSD
%
% This is done by solving a semidefinite program (SDP) truncated at level t
% Note: to be executed, ChebPoly_dual requires GloptiPoly 3 


% Inputs:
% k = (k_1,...,k_n): the multidegree of the monomial
% Omega: the domain, which can be:
%  -> H: the hypercube
%  -> B: the euclidean ball
%  -> C: the cross-polytope
%  -> S: the simplex 
% t: the truncation level

% Outputs:
% err:         the optimal value of the SDP 
%              (this is an lower bound for the true error of best approximation)
% signature_p: the signature points where (m_k-p)(x) = + ||m_k-p||
% signature_m: the signature points where (m_k-p)(x) = - ||m_k-p||

% Written in May 2024
% Send comments to simon.foucart@centraliens.net

%%
function [err,signature_p,signature_m] = ChebPoly_dual(k,Omega,t)

%% the parameters of the problem
d = length(k);  % number of variables
n = sum(k);     % (total) degree of the monomial

%% the optimization variables
mpol('x_p',d);
mpol('x_m',d);
% Assign variables to measures
mu_p = meas(x_p);  % mu 'plus'
mu_m = meas(x_m);  % mu 'minus'
% the monomial
mon_p = x_p(1)^k(1);
mon_m = x_m(1)^k(1);
for i = 2:d
    mon_p = mon_p * x_p(i)^k(i);
    mon_m = mon_m * x_m(i)^k(i);
end

%% the objective function
obj = mom(mon_p) - mom(mon_m);

%% the polynomial descriptors of the domain Omega
if Omega == 'H'
    for h = 1:d
        g_p{h} = 1 - x_p(h)^2;
        g_m{h} = 1 - x_m(h)^2;
    end
end
if Omega == 'B'
    g_p{1} = 1;
    g_m{1} = 1;
    for h = 1:d
        g_p{1} = g_p{1} - x_p(h)^2;
        g_m{1} = g_m{1} - x_m(h)^2;
    end
end
if Omega == 'C'
    signs = 2*(dec2bin(2^d-1:-1:0)-'0')-1; % produce all sign vectors of lenghth d
    for h = 1:size(signs,1)
        g_p{h} = 1;
        g_m{h} = 1;
        for j = 1:d
            g_p{h} = g_p{h} - signs(h,j)*x_p(j);
            g_m{h} = g_m{h} - signs(h,j)*x_m(j);
        end
    end
end
if Omega == 'S'
    for h = 1:d
        g_p{h} = x_p(h);
        g_m{h} = x_m(h);
    end
    g_p{d+1} = 1;
    g_m{d+1} = 1;
    for i=1:d
        g_p{d+1} = g_p{d+1} - x_p(i);
        g_m{d+1} = g_m{d+1} - x_m(i);
    end
end

%% Constraints
% test functions
v_p = mmon(x_p,n-1);
v_m = mmon(x_m,n-1);
% moment constraints 
mom_const = [mom(v_p) - mom(v_m) == 0; mass(mu_p)+mass(mu_m) - 1 == 0];

%% Localization constraints 
loc_const = [];
for h = 1:length(g_p)
    loc_const = [loc_const; 0 <= g_p{h}];
end
for h = 1:length(g_m)
    loc_const = [loc_const; 0 <= g_m{h}];
end

%% build the problem
P = msdp(max(obj),mom_const,loc_const,t);

%% solver options
SDPsolver = 'sedumi';
options = getSolverParams(SDPsolver);
if (~isempty(options))
    mset(options);
end

%% solving the optimization program
[status,obj,M,~] = msol(P);
if(status < 0)
    disp('SDP could not be solved')
    return
end

%% return the outputs
err = double(obj);
try
     points_m = double(M(1));
     nb_points_m = size(points_m,3);
     signature_m = reshape(points_m,[d,nb_points_m]);

catch
     disp('Warning: support cannot be extracted')
     signature_m = [];
end
try
     points_p = double(M(2));
     nb_points_p = size(points_p,3);
     signature_p = reshape(points_p,[d,nb_points_p]);
catch
     disp('Warning: support cannot be extracted')
     signature_p = [];
end

end