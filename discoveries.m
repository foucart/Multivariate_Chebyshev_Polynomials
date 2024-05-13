%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproducible MATLAB file accompanying the paper
%        OPTIMIZATION-AIDED CONSTRUCTION OF
%        MULTIVARIATE CHEBYSHEV POLYNOMIALS          
% by M. Dressler, S. Foucart, E, de Klerk, M. Joldes, 
%    J. B. Lasserre, and Y. Xu
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here, one exploits ChebPoly_dual and ChebPoly_dual
% to discover Chebyshev polynomials for d=3 and n<=6
% Written in May 2024
% Send comments to simon.foucart@centraliens.net

format long


%% %%%%%%%%%%%%%%%%%%%%%%%
%  THE 3D-EUCLIDEAN BALL %
% %%%%%%%%%%%%%%%%%%%%%%%%

%% Case 1: k=(2,2,1)
% note this case is solved explicitly in the paper
Omega = 'B';
k = [2,2,1];
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 5;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are equal, so we can "confidentaly" say that
% the error of best approximation to x1^2*x2^2*x3 
% by polynomials of degree <5 on the ball is 0.0363000

%% notice that the signature lies on the sphere x^2+y^2+z^2=1
sum(signature_m.^2)
sum(signature_p.^2)

%% Case 2: k=(4,1,1)
Omega = 'B';
k = [4,1,1];
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 5;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are equal, so we can "confidentaly" say that
% the error of best approximation to x1^4*x2*x3 
% by polynomials of degree <6 on the ball is 0.0192378

%% notice that the signature lies on the sphere x^2+y^2+z^2=1
sum(signature_m.^2)
sum(signature_p.^2)

%% Case 3: k=(3,2,1)
Omega = 'B';
k = [3,2,1];
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 5;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are equal, so we can "confidentaly" say that
% the error of best approximation to x1^3*x2^2*x3 
% by polynomials of degree <5 on the ball is 0.0165297

%% in this case, we notice that the signature lies on the sphere x^2+y^2+z^2=1
sum(signature_m.^2)
sum(signature_p.^2)



%% %%%%%%%%%%%%%%%%%%%%%%
% THE 3D-CROSS-POLYTOPE %
% %%%%%%%%%%%%%%%%%%%%%%%

%% Case 1: k=(1,1,1)
Omega = 'C';
k = [1,1,1];
s = 7;
err_primal = ChebPoly_primal(k,Omega,s);
t = 8;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are equal, so we can "confidentaly" say that
% the error of best approximation to x1*x2*x3 
% by polynomials of degree <3 on the cross-polytope is 0.037037

%% notice that the signature lies on the L1-sphere |x|+|y|+|z|=1
sum(abs(signature_m))
sum(abs(signature_p))

%% Case 2: k=(2,1,1)
Omega = 'C';
k = [2,1,1];
s = 8;
err_primal = ChebPoly_primal(k,Omega,s);
t = 9;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are equal, so we can "confidentaly" say that
% the error of best approximation to x1^2*x2*x3 
% by polynomials of degree <4 on the cross-polytope is 0.012736

%% notice that the signature lies on the L1-sphere |x|+|y|+|z|=1
sum(abs(signature_m))
sum(abs(signature_p))

%% Case 3: k=(3,1,1)
Omega = 'C';
k = [3,1,1];
s = 10;
err_primal = ChebPoly_primal(k,Omega,s);
t = 11;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are "about" equal --- we can say that
% the error of best approximation to x1^3*x2*x3 
% by polynomials of degree <5 on the cross-polytope 
% is between 4.76444e-3 and 4.76445e-3

%% notice that the signature lies on the L1-sphere |x|+|y|+|z|=1
sum(abs(signature_m))
sum(abs(signature_p))

%% Case 4: k=(2,2,1)
Omega = 'C';
k = [2,2,1];
s = 10;
err_primal = ChebPoly_primal(k,Omega,s);
t = 11;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are equal, so we can "confidentaly" say that
% the error of best approximatio to x1^2*x2^2*x3 
% by polynomials of degree <5 on the cross-polytope is 0.00339887

%% notice that the signature lies on the L1-sphere |x|+|y|+|z|=1
sum(abs(signature_m))
sum(abs(signature_p))

%% Case 5: k=(4,1,1)
Omega = 'C';
k = [4,1,1];
s = 10;
err_primal = ChebPoly_primal(k,Omega,s);
t = 11;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are "about" equal --- we can say that
% the error of best approximatio to x1^4*x2*x3 
% by polynomials of degree <6 on the cross-polytope 
% is between 1.85306e-3 and 1.85307e-3

% note: the signature could not be extracted...

%% Case 6: k=(3,2,1)
Omega = 'C';
k = [3,2,1];
s = 10;
err_primal = ChebPoly_primal(k,Omega,s);
t = 11;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are "about" equal --- we can say that
% the error of best approximation to x1^3*x2^2*x3 
% by polynomials of degree <6 on the cross-polytope 
% is between 1.087950e-3 and 1.087955e-3

%% the signature may not lie on the L1-sphere |x|+|y|+|z|=1 ...
sum(abs(signature_m))
sum(abs(signature_p))

%% Case 7: k=(2,2,2)
Omega = 'C';
k = [2,2,2];
s = 14;
err_primal = ChebPoly_primal(k,Omega,s);
t = 13;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are "about" equal --- we can say that
% the error of best approximation to x1^2*x2^2*x3^2 
% by polynomials of degree <6 on the cross-polytope 
% is between 6.613e-4 and 6.614e-4

% note: the signature could not be extracted...



%% %%%%%%%%%%%%%%%
% THE 3D-SIMPLEX %
% %%%%%%%%%%%%%%%%

%% Case 1: k=(2,1,1)
% note this case is solved explicitly in the paper
Omega = 'S';
k = [2,1,1];
s = 6;
err_primal = ChebPoly_primal(k,Omega,s);
t = 7;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are equal, so we can "confidentaly" say that
% the error of best approximation to x1^2*x2*x3 
% by polynomials of degree <4 on the simplex is 0.0026885

%% notice that the signature lies on the face x+y+z=1
sum(signature_m)
sum(signature_p)

%% Case 2: k=(3,1,1)
Omega = 'S';
k = [3,1,1];
s = 10;
err_primal = ChebPoly_primal(k,Omega,s);
t = 8;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are "about" equal --- we can say that
% the error of best approximation to x1^3*x2*x3 
% by polynomials of degree <5 on the simplex 
% is between 5.9841e-4 and 5.9842e-4

% note: the signature could not be extracted...

%% Case 3: k=(2,2,1)
Omega = 'S';
k = [2,2,1];
s = 10;
err_primal = ChebPoly_primal(k,Omega,s);
t = 6;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are equal, so we can "confidentaly" say that
% the error of best approximation to x1^2*x2^2*x3 
% by polynomials of degree <5 on the simplex is 4.6956e-4

%% notice that the signature lies on the face x+y+z=1
sum(signature_m)
sum(signature_p)

%% Case 4: k=(4,1,1)
Omega = 'S';
k = [4,1,1];
s = 12;
err_primal = ChebPoly_primal(k,Omega,s);
t = 10;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are "about" equal --- we can say that
% the error of best approximation to x1^4*x2*x3 
% by polynomials of degree <6 on the simplex 
% is between 1.4051e-4 and 1.4052e-4 

% note: the signature could not be extracted...

%% Case 5: k=(3,2,1)
Omega = 'S';
k = [3,2,1];
s = 10;
err_primal = ChebPoly_primal(k,Omega,s);
t = 8;
[err_dual,signature_p,signature_m] = ChebPoly_dual(k,Omega,t);
[err_dual err_primal]
% these two values are "about" equal --- we can say that
% the error of best approximation to x1^3*x2^2*x3 
% by polynomials of degree <6 on the simplex 
% is between 1.00083e-4  and 1.00089e-4

% note: the signature could not be extracted...
