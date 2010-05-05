%  The script file orrsom.m computes the eigenvalues of the Orr-Sommerfeld
%  equation using NxN Chebyshev differentiation matrices.

% S.C. Reddy, J.A.C. Weideman 1998.  Code modified to display output
% by JACW, May 2003.

N = 40;            % Number of collocation points in the wall-normal direction

delta = 1;          % Convert Reynolds number to a boundary layer thickness for numerical calculation
ymax = 2;           % For truncation of the computational domain in the wall-normal y-domain.  The maximum y-value

%logm = [3:0.02:5];
%Rm = 10.^logm;
Rm = [1000:100:10000];                                              % Reynolds number
%km = [0.01*(2*pi/delta):0.01*(2*pi/delta):0.1*(2*pi/delta)];	% Range of 'k' values
km = [0.2:0.02:2.0];
%Rm = 10^(3.5);
%km = 0.1*(2*pi/delta);
%km = 1.2;

%------------------
%[x,D2,D1,phip,phim]=cheb2bc(N,[1 0 0;0 1 0]);
%nn = size(x,1);

[x,DM] = chebdif(N+2,2);                           % Compute second derivative
D2 = DM(2:N+1,2:N+1,2);                            % Enforce Dirichlet BCs

[x,D4] = cheb4c(N+2);                              % Compute fourth derivative
I = eye(size(D4));                                 % Identity matrix
%------------------

emat = zeros(size(km,2),size(Rm,2));
ind1 = 1;ind2 = 1;
cntr = 1;
for R = Rm;
    for k = km;

        % Get the other parameters from the Reynolds number and "delta"
        % (Reynolds number is based on delta)
        
        % Varying mu (or nu) (Uinf = 1 AND density = 1)
        Uinf = 1;rho = 1;
        mu = Uinf*rho*delta/R;
        
        % Varying Uinf (mu = 1 AND density = 1)
        %rho = 1;mu = 1;
        %Uinf = R*mu/(delta*rho);
    

        %---------------------------
        % Get the mean flow velocity profile
        %[udivU,vort,ddy2udivU] = pohlddy2((x+1)./delta);               % Polhausen approximation to the Blasius profile
        udivU = (1 - x.^2);ddy2udivU = -2*(delta.^2).*ones(size(x));    % Plane Pousille flow.
        
        U = udivU.*Uinf;
        d2Udy2 = ddy2udivU.*(Uinf./(delta.^2));
        %---------------------------

        %A = (D4-2*D2+I)/R-2*i*I-i*diag(U)*(D2-I);              % Set up A and B matrices
        A = (D4 - 2*(k^2)*D2 + I*(k^4)).*(mu/(rho)) + (i*k)*diag(d2Udy2)*I - (i*k)*diag(U)*(D2-I*(k^2));
        B = D2 - (k^2)*I;


        e = eig(A,B);                                           % Compute eigenvalues
        %[V,D] = eig(A,B);
        %e = diag(D);

        [m,l] = max(real(e));                                   % Find eigenvalue of largest

        emat(ind1,ind2) = e(l);
        ind1 = ind1 + 1;
        
        cntr = cntr + 1;
    end
    [m,l] = max(real(emat(:,ind2)));
        
    disp(['Re=' num2str(R) ',k=' num2str(k) ',Eigenvalue with largest real part = ' num2str(emat(l,ind2)) ' :: Step ' num2str(cntr)])            % real part

    
    ind1 = 1;
    ind2 = ind2 + 1;
end

save emat.mat
cs=contour(Rm,km,real(emat),[-1:0.01:0]);
clabel(cs);grid
xlabel('Re (based on channel half-width)')
ylabel('Wavenumber (k)')
title('Marginal Stability Curve for Plane Poiseuille flow')


% %---------------------------
% % Get the real-world coordinates
%
% b = 1 + (2*sL)/ymax;
% yn = sL.*(1+zeta)./(b-zeta);
% ynbc = sL.*(1+x)./(b-x);
%
% %---------------------------
% % Do a mapping of the derivation matrices
%
% % In this case, zeta = (y*b - L)/(y + L)
% % so to transform the differentiation matrixes, we need to know the value
% % of d(zeta)/dy and d^2(eta)/dy^2.
% % Therefore:
% zetad = (b./(ynbc+sL))-(((ynbc.*b)-sL)./((ynbc+sL).^2));
% zetadd = (-2.*b./((ynbc+sL).^2)) + (2.*((ynbc.*b)-sL)./((ynbc+sL).^3));
%
% % Get the first and second order differentiation matricies in the y-direction
% % In reality only the second order diff matrix is required.
% D1 = diag(zetad,0)*D1t;
% D2 = diag(zetadd,0)*D1t + diag((zetad.^2),0)*D2t;
