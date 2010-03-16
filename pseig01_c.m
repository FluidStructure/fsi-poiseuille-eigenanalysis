function [] = pseig01_c()

more off

load LHSRHSf.mat

disp('Performing LHS/RHS...')
tic;
RHS02 = LHS\RHS;
t = toc;
disp(['Matrix division completed in t = ' num2str(t) ' seconds'])
clear LHS RHS
save RHS02.mat

disp('Starting calculation of eigenvalues...')
nummods = 8;
sigst = 'LR';
% If: Period = lamda/Uinf => Freq = Uinf/lamda => Omega = 2*pi*Freq;
%sigst = 0 - 1*i;
%sigst = 0 + 0*i;
%sigst = 0 + 0*i;
opts.maxit = 20000;
opts.disp = 1;
%opts.tol = eps*1e-3;
%if exist('v0.mat') == 2
%    load v0.mat v0
%    opts.v0 = v0
%end
[Veigs,Deigs,flagc] = eigs(RHS02,nummods,sigst,opts);
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
v0 = real(Veigs(:,1));
save('-v7','v0.mat','v0')

