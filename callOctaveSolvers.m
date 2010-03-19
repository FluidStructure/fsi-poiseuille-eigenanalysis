function [] = solveEigenvalues()

more off

% Set some parameters
nummods = 8;
opts.maxit = 20000;
opts.disp = 1;
%opts.tol = eps*1e-3;
% Set point sigma: NOTE If: Period = lamda/Uinf => Freq = Uinf/lamda => Omega = 2*pi*Freq;
sigst = 'LR';
%sigst = 0 - 1*i;

disp('Loading matrix A... into octave')
load A.mat
disp('Done...')

disp('Starting calculation of eigenvalues...')

% Load an existing solution for the initial guess (needs some work)
%if exist('v0.mat') == 2
%    load v0.mat v0
%    opts.v0 = v0
%end
[Veigs,Deigs,flagc] = eigs(A,nummods,sigst,opts);
%[Veigs,Deigs,flagc] = eigs(RHS,LHS,nummods,sigst,opts);
if flagc == 0
    disp('ALL EIGENVALUES CONVERGED :-)')
else
    error('NOT ALL EIGENVALUES CONVERGED :-(')
end
evals = diag(Deigs)

load VARS.mat R
fname = ["results/evals_R" num2str(R) ".mat"]
save('-v7',fname,'Veigs','evals')

% Save a vector for an initial guess for the next step
%v0 = Veigs(:,1);
%save('-v7','v0.mat','v0')
