%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproducible MATLAB file accompanying the paper
%        OPTIMIZATION-AIDED CONSTRUCTION OF
%        MULTIVARIATE CHEBYSHEV POLYNOMIALS          
% by M. Dressler, S. Foucart, E, de Klerk, M. Joldes, 
%    J. B. Lasserre, and Y. Xu
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here, one verifies ChebPoly_dual and ChebPoly_dual
% by comparing with known Chebyshev polynomials 
% Written in May 2024
% Send comments to simon.foucart@centraliens.net



%% %%%%%%%%%%%%%%%
%  THE HYPERCUBE %
% %%%%%%%%%%%%%%%%

%% An instance for the 2d-hypercube
Omega = 'H';
k = [3,2];
d = length(k);  % number of variables
n = sum(k);     % degree of the monomial
err_theory = 2^(-n+d);
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 5;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree

%% We can see that tensor-product structure of the signature

scatter(signature_m(1,:),signature_m(2,:),'bo')
hold on
scatter(signature_p(1,:),signature_p(2,:),'r+')
xlim([-1 1])
ylim([-1 1])
xticks([-1 1])
yticks([-1 1])


%% An instance for the 3d-hypercube
Omega = 'H';
k = [2,1,1];
d = length(k);  % number of variables
n = sum(k);     % degree of the monomial
err_theory = 2^(-n+d);
s = 5;
err_primal = ChebPoly_primal(k,Omega,s);
t = 6;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree



%% %%%%%%%%%%%%%%%%%%%%
%  THE EUCLIDEAN BALL %
% %%%%%%%%%%%%%%%%%%%%$

%% An instance for the 2d-euclidean ball
Omega = 'B';
k = [3,2];
d = length(k);  % number of variables
n = sum(k);     % degree of the monomial
err_theory = 2^(-n+1);
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 5;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree

%% The 3d-euclidean ball: particular case 1
Omega = 'B';
k = [1,1,1];
err_theory = 3^(-3/2);
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 5;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree

%% The 3d-euclidean ball: particular case 2
Omega = 'B';
k = [2,1,1];
err_theory = (3-sqrt(8))/2;
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 5;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree

%% The 3d-euclidean ball: particular case 3
Omega = 'B';
k = [3,1,1];
r = roots([9 -29 24 -29 9]);
a = min( r(abs(imag(r))<1e-4) );
err_theory = (1-a)*(a^3/5)^(1/4)/5;
s = 7;
err_primal = ChebPoly_primal(k,Omega,s);
t = 6;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree

%% The 3d-euclidean ball: particular case 4 
Omega = 'B';
k = [2,2,2];
err_theory = 1/72;
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 6;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree

%% The 3d-euclidean ball: particular case 5
% checking only the primal, as the dual takes too long...
Omega = 'B';
k = [4,4,4];
err_theory = 1/21.8935834/27^2;
s = 15;
err_primal = ChebPoly_primal(k,Omega,s);
[err_theory err_primal]
% as expected, the two values are the same


%% %%%%%%%%%%%%%
%  THE SIMPLEX %
% %%%%%%%%%%%%%%

%% An instance for the 2d-simplex
Omega = 'S';
k = [2,2];
d = length(k);  % number of variables
n = sum(k);     % degree of the monomial
err_theory = 2^(-2*n+1);
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 6;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree

%% The 3d-simplex: particular case 1
Omega = 'S';
k = [1,1,1];
err_theory = 1/72;
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 6;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree

%% The 3d-simplex: particular case 2
Omega = 'S';
k = [2,2,2];
err_theory = 1/21.8935834/27^2;
s = 10;
err_primal = ChebPoly_primal(k,Omega,s);
t = 8;
err_dual = ChebPoly_dual(k,Omega,t);
[err_dual err_theory err_primal]
% as expected, all the values agree
