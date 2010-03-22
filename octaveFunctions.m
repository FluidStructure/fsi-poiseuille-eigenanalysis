function [] = setupMatrices()

%
% This matlab code tries to extract the eigenvalues and eigenmodes from a
% full pouseulle flow with one-channel wall being a flexible wall system.  
%
% The equations for the perturbations within the boundary layer are the
% linearised Navier-Stokes equations for pertuburbations to a mean flow, in
% streamfunction form.  This is the equation from which the
% Orr-Sommerfeld equation is derived by assuming perturbation shape for the
% streamfunction.
%
%
% NOTE: The solution eigenvector will be:
%       X = [T1 T2 T3 ... Tn]
% where:
%       Tn = Nth node of free-stream perturbation streamfunction.  The 1st
%       node will be the lower-left node at (x=0,y=0) and the numbers will
%       increase from left to right, and bottom to top respectively.
%
% 

more off

%-----------------------------------
% Set some parameters for the system:
R = 5000;		    % Reynolds number (based on displacement thickness of boundary layer) 
                    % instability from Orr-Sommerfeld is at 5800 (corresponding wavenumber = 1.025)

chebN = 16;         % Number of Chebyshev collocation points in the wall-normal direction
Nx = 200;           % Number of collocation points (panels) in the x-direction (for the flow elements)
NWxb1 = 90;         % Number of collocation points (panels) in the x-direction (for the solid wall elements) on the first (upstream) part of the lower wall
NCx = 20;           % Number of collocation points (panels) in the x-direction (for the compliant wall elements)

% REMEMBER THAT THE CHANNEL HAS A WIDTH OF 2!
%LT = 20;            % Length of the channel section
LT = 60;

% WALL PARAMETERS
%rhow = 2600;       % Density of the wall material
rhow = 1000;
h = 0.005;          % Thickness of the wall
d = 0;              % Damping coefficient of the backing-damping (N/m)
%E = 58.85e9;       % Youngs modulus of the wall material
E = 2.592e7;
K = 0;              % Stiffness of the backing springs

B = E*(h^3)/12;     % Flexural Rigidity of the wall

%%%%%%%%%%%%%%%%%%%%%%%%%%
% START THE SOLUTION
NWxb2 = Nx - NWxb1 - NCx;   % Number of collocation points (panels) in the x-direction (for the solid wall elements) on the
                            % second (downstream) part of the lower wall
NWxu = Nx;          	    % Number of collocation points (panels) in the x-direction (for the solid wall elements) on the upper surface

FS = LT*(NWxb1/Nx); 	    % Starting point of the flexible wall section
L = LT*(NCx/Nx);    	    % Length of the wall compliant (m)

delta = 1;                  % Convert Reynolds number to a boundary layer thickness for numerical calculation (delta = half channel width)
ymax = 2;                   % For truncation of the computational domain in the wall-normal y-domain.  The maximum y-value

Uinf = 1;                   % Mean flow velocity (m/s)
rho = 1;                    % Flow density (m/s)
mu = Uinf*rho*delta/R;      % Flow viscosity

nu = mu./rho;               % Flow kinematic viscosity

%-----------------------------------
% Establish the grid and differentiation matrix in the wall-normal
% direction

% Get the chebyshev grid points and 1st and 2nd order differentiation matrices over [-1,1]
% (but we really only need the 2nd order differentiation matrix in the
% wall-normal direction)
[ynbc,d2dy,d1dy,phip,phim]=cheb2bc(chebN,[0 1 0;0 1 0]);

% Plot the grid points
#fig1 = figure;
#plot(ynbc,'rs-');axis tight;grid

%---------------------------
% Get the height of each cell

% Get the node points of each cell and cell widths
yn = zeros(size(ynbc,1)+1,1);
yn(1,1) = ynbc(1);
for s = 1:(size(ynbc,1)-1)
    yn(s+1,1) = (ynbc(s+1) + ynbc(s))./2;
end
yn(s+2,1) = ynbc(size(ynbc,1),1);

% Get width and centres of each cell
dy = zeros(size(yn,1)-1,1);
for s = 1:(size(yn,1)-1)
    dy(s,1) = (yn(s+1) - yn(s));
end
dy = abs(dy);

% Get the vector of cell widths that correspond to each element
dhv = sparse(kron(eye(Nx),diag(dy)));

% Get the position of the flow elements (collocation points)
pcy = ynbc;
pcy(1,1) = (yn(1) + yn(2))./2;
pcy(size(pcy,1),1) = (yn(size(yn,1),1) + yn(size(yn,1)-1,1))./2;

% Plot the position of cell node points and centres
%figure(fig1)
%hold on;
%for s = 1:size(yn,1);plot([1 size(yn,1)],[yn(s,1) yn(s,1)],'r-');end
%for s = 1:size(pcy,1);plot([1 size(pcy,1)],[pcy(s,1) pcy(s,1)],'g-');end
%hold off

%---------------------------
% Get the x position of each element

% x coordinates of the flow elements
dxf = LT/Nx;
pcx = linspace(0.5*dxf,LT-(0.5*dxf),Nx);

% x-coordinates of the upper wall elements
dxu = (LT/NWxu);
cwu = linspace(0.5*dxu,LT-(0.5*dxu),NWxu);

% x-coordinates of the first part of the lower rigid wall (before the compliant section)
dxb1 = (FS/NWxb1);
cwb1 = linspace(0.5*dxb1,FS-(0.5*dxb1),NWxb1);

% x-coordinates of teh second part of the lower rigid wall (after the compliant section)
dxb2 = (LT-(FS+L))/NWxb2;
cwb2 = linspace(FS+L+(0.5*dxb2),LT-(0.5*dxb2),NWxb2);

% x-coordinates of the compliant wall section
dxc = (L/NCx);
xn = linspace(FS,FS+L,(NCx+1));
xc = (xn(1:(size(xn,2)-1)) + xn(2:size(xn,2)))./2;

%----------------------------
% Do a meshgrid for all the collocation points in the fluid domain
[pcxx,pcyy] = meshgrid(pcx,pcy);

% Reshape the square matrices of coordinates into vectors for calculating influence coefficients
n1 = size(pcxx,1);n2 = size(pcxx,2);
pcxxv = reshape(pcxx,(n1*n2),1);
pcyyv = reshape(pcyy,(n1*n2),1);

if exist('VARS.mat') == 2
    eVARS = load('VARS.mat');
end

save -v7 VARS.mat

%-------------------------------------
% Get the finite difference matrix for gradients in the x and y directions

% First derivative matrix in the x-direction (might need to use something other than a centred differencing here (Upwinding - to not get zeros on the diagonal)
%[dndx] = dndx1(pcx,dx,1,1);
%dndx(1,2) = 1*((1/(2*dx))); 	% Set the inlet boundary condition so that eta(-1)==0. Not Deta(1)/Dx==0 as is at the trailing edge. 
%OR:
[dndx] = dndx1uw(pcx',dxf,1,1);	% Use an upwinding method (shown to be more stable)
dndx = sparse(dndx);

% Second derivative matrix in the x-direction
[d2ndx2] = dndx2(pcx',dxf,1,1);
d2ndx2(1,2) = 1*((1/(dxf^2))); % Set the inlet boundary condition so that eta(-1)==0. Not Deta(1)/Dx==0 as is at the trailing edge.
d2ndx2 = sparse(d2ndx2);

% Apply the finite difference grids over the entire fluid domain (pcxxv,pcyyv) (kron does a tensor-product)
%dndxK = zeros(n1*n2);d2ndx2K = dndxK;
dndxK = sparse(kron(dndx,eye(length(pcy))));
d2ndx2K = sparse(kron(d2ndx2,eye(length(pcy))));

% y-differential (chebyshev) over the entire domain
d2ndy2K = sparse(kron(eye(length(pcx)),d2dy));

%---------------------------------
% A matrix that averages the nodes either side of a panel (for getting
% panel velocity from node velocity)
avmat = zeros(length(xc),length(xn));
for s = 1:length(xc)
    avmat(s,s) = 0.5;
    avmat(s,s+1) = 0.5;
end

%---------------------------------
% Save the finite-difference matrices
save -v7 FDMATS.mat


%--------------------------------
% Calculate the influence coefficients (if necessary)
remakeICs = 1;
if (exist('eVARS') == 1) & (exist('ICs.mat') == 2)
    if (eVARS.Nx==Nx) & (eVARS.chebN==chebN) & (eVARS.LT==LT)
        disp('Geometry seems to be the same. Using existing Influence Coefficient matrices...')
        remakeICs = 0;
    end
end

if remakeICs == 1
    F = system('python writeICs.py');

    load ICs.mat INww
    disp('Inverting wall-wall influence matrix...')
    invINww = inv(INww);
    disp('done')
    save invINww.mat invINww
end

endfunction


function [] = evaluateMatrices()

more off

load ICs.mat
load invINww.mat
load VARS.mat
load FDMATS.mat

%---------------------------
% Get the mean flow velocity profile
udivU = (1 - pcy.^2);ddy2udivU = -2*(delta.^2).*ones(size(pcy));    % Plane Pousille flow.
U = udivU.*Uinf;
d2Udy2 = ddy2udivU.*(Uinf./(delta.^2));

%---------------------------

% Matrix "A" - the coefficient that gives wall element strengths when
% pre-multiplied with the vector of flow element strengths and wall
% velocities
MA = [((-1*invINww*INfw)*dhv) zeros(length(INww),length(xn)) invINww(:,[1:length(xc)])*avmat];



% Matrix "B" - the coefficient that gives vertical velocity at the flow
% elements "vp" when pre-multiplied with the vector of flow element
% strengths and wall velocities
MB = [(INff*dhv) zeros(length(INff),length(xn)) zeros(length(INff),length(xn))];
MB = MB + [INwf*MA];



% Matrix "C" - Right hand side of fluid transport equation
MC = [diag(-1*Uinf*(1 - pcyyv.^2))*dndxK zeros(length(dndxK),length(xn)) zeros(length(dndxK),length(xn))];
MC = MC - MB*(2*Uinf);             % This is the second term in the "for" loop above
% nu = 0.0;
MC = MC + nu*[(d2ndx2K+d2ndy2K) zeros(length(dndxK),length(xn)) zeros(length(dndxK),length(xn))];



% Matrix "P" - the coefficient that gives tangential velocities at the wall
% element surfaces.
MP = [(ITfw*dhv) zeros(size(ITfw,1),length(xn)) zeros(size(ITfw,1),length(xn))];
MP = MP + (ITww*MA);

% Matrix "D" - Left hand side of fluid transport equation
MD = sparse([eye(size(dndxK)) zeros(length(dndxK),length(xn)) zeros(length(dndxK),length(xn))]);
% Add in the terms for the injection of vorticity
if 0    % To add in these terms or not (for debugging)
    cc = 0;
    s = 0;
    while s < size(ITfw,1)
        disp(['Injecting vorticity at ' num2str(s+1) ' and ' num2str(s+2) ' of ' num2str(size(ITfw,1))])
        cc = cc + 1;s = s + 1;
        MD(cc,:) = MD(cc,:) - MP(s,:)./dhv(cc,cc);      % This is probably not right because of the order of strengths in ITfw
        cc = cc + (length(pcy)-1);s = s + 1;
        MD(cc,:) = MD(cc,:) + MP(s,:)./dhv(cc,cc);
    end
end
%MD = MD + kron((-1*MP),diag([1;zeros(length(pcy)-2,1);1]));       % Check this

% Matrix "E" - Left hand side of equation for "v1_dot = v2"
ME = sparse([zeros(length(xn),length(dndxK)) eye(length(xn)) zeros(length(xn))]);

% Matrix "F" - Right hand side of equation for "v1_dot = v2"
MF = sparse([zeros(length(xn),length(dndxK)) zeros(length(xn)) eye(length(xn))]);

% Matrix "G" - Left hand side of wall equation
MG = [zeros(length(xn),length(dndxK)) zeros(length(xn)) (rhow*h)*eye(length(xn))];
% Add the pressure forcing term
%PF = (MP([1:length(xc)-1],:)+MP([2:length(xc)],:));  % This is not right
PF = dxc.*cumsum(MP([1:length(xc)-1],:),1);       % Integral of the tangential velocity at the wall (excluding the last node)
MG = MG + [zeros(1,(length(dndxK)+2*length(xn)));PF;zeros(1,(length(dndxK)+2*length(xn)))];
% If you turn off the second line above then this turns off presure forcing
% on the wall!!  (Might be useful for validation or very stiff walls).

% Matrix "H" - Right hand side of wall equation
dn4dx = dndx4(xn',dxc,2,2);     % (the 2's set the hinged-hinged BC's)
MH = [zeros(length(xn),length(dndxK)) -1*(B*dn4dx + K*eye(length(xn))) -1*d*eye(length(xn))];

% Construct the left and right hand matrices
LHS = sparse([MD;ME;MG]);
RHS = full([MC;MF;MH]);

% Remove the end-nodes (which don't move) from the set of equations
nn = [(length(pcyyv)+1) (length(pcyyv)+length(xn)) (length(pcyyv)+length(xn)+1) (length(pcyyv)+(2*length(xn)))]';
LHS(nn,:) = [];
LHS(:,nn) = [];
RHS(nn,:) = [];
RHS(:,nn) = [];

save LHSRHS.mat LHS RHS
%close all;clear all;

% Get the LHS and RHS for the fluid only
f = [1:length(dndxK)];
LHS = LHS(f,f);
RHS = RHS(f,f);
save LHSRHSf.mat LHS RHS

endfunction


function [] = calculateEigenmodes()

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

endfunction

function [] = plotResults()

more off

load VARS.mat
load evals.mat

% 20 gives nice looking results
evals
s = ss      % this is the eigenvalue to plot (from evals)

cmax = 0.008;
cmin = -0.008;

%T = 2*pi/s;
%t = linspace(0,3*T,60);

sn = length(pcyyv);
nn = length(xn)-2;  % -2 for the end-nodes that are removed

if length(Veigs) > sn
    % Get the wall modes
    wm = Veigs([sn+1:sn+nn],s);

    figure
    plot(real(wm));
    hold on;plot(imag(wm),'r');hold off
    title('Wall Eigenvalues (real=blue, red=imag)')
end

% Get the flow modes
fm = Veigs([1:sn],s);
fm = reshape(fm,length(pcy),length(pcx));

% Crop results in the x-direction (if you want)
%minx = 24;ind=find(pcx<minx);pcx(ind)=[];fm(:,ind)=[];
%maxx = 36;ind=find(pcx>maxx);pcx(ind)=[];fm(:,ind)=[];

figure
%contour(pcx,pcy,fm);
%contour(pcx,pcy,fm,[cmin:(cmax-cmin)/200:cmax])
%grid;colorbar;
pcolor(pcx,pcy,real(fm));caxis([cmin cmax]);shading interp;colorbar
print -dpng vortField.png

endfunction

function [] = plotFluidMotion()

load VARS.mat
load evals.mat

s = 7;      % this is the eigenvalue to plot (from evals)

cmax = 0.002;
cmin = -0.002;

ss = evals(s);
T = 2*pi/imag(ss);
t = linspace(0,2*T,80);

sn = length(pcyyv);
nn = length(xn)-2;  % -2 for the end-nodes that are removed

if length(Veigs) > sn
    % Get the flow modes
    fm = Veigs([1:sn],s);
    fm = reshape(fm,length(pcy),length(pcx));
    
    % Get the wall modes
    wm = Veigs([sn+1:sn+nn],s);
else
    % Get the flow modes
    fm = Veigs(:,s);
    fm = reshape(fm,length(pcy),length(pcx));
end

% Crop results in the x-direction (if you want)
minx = 24;ind=find(pcx<minx);pcx(ind)=[];fm(:,ind)=[];
maxx = 36;ind=find(pcx>maxx);pcx(ind)=[];fm(:,ind)=[];

fig1 = figure;
set(fig1,'Position',[100 100 800 200],'PaperPosition',[0.634517 6.34517 20 5])
if length(Veigs) > sn
    fig2 = figure;
    set(fig2,'Position',[100 100 800 200],'PaperPosition',[0.634517 6.34517 20 5])
end
c = 1;
for tt = t
    % DUMP THE FLUID 
    fmm = real(fm.*exp(-1*imag(ss)*i*tt));    % This is only the oscillatory part
    fmm = fmm*diag(exp(real(ss)*pcx));
    fmm = fmm./diag(exp(real(ss)*pcx(1)));    % Take into account convective instability
    figure(fig1);
    %contour(pcx,pcy,fm);
    contour(pcx,pcy,fmm,[cmin:(cmax-cmin)/40:cmax]);caxis([cmin cmax])
    %pcolor(pcx,pcy,fmm);caxis([cmin cmax]);shading interp
    grid;colorbar;
    title(['Frame ' num2str(c) ' of ' num2str(length(t))])
    
    fname = mvarname('jpgs/fluid','frame',c);
    fname = strcat([fname '.jpg']);
    print(fig1,'-djpeg100','-r150',fname);
    pause(0.000001)
    
    if length(Veigs) > sn
        % DUMP THE WALL
        wmm = real(wm.*exp(-1*imag(ss)*i*tt));
        wmm = wmm./max(abs(wm));
        wmm = [0;wmm;0];
        figure(fig2);
        plot(xn,wmm,'b');grid;
        axis([min(xn) max(xn) -1 1])

        fname = mvarname('jpgs/wall','frame',c);
        fname = strcat([fname '.jpg']);
        print(fig2,'-djpeg100','-r150',fname);
        pause(0.000001)
    end
    
    % STEP THE COUNTER
    c = c+1;
end

endfunction


function [xt,D2t,D1t,phip,phim]=cheb2bc(N,g)

% Program for computing first and second derivative matrices and
% and boundary condition functions for 2 point boundary conditions
%
%  a_1 u(1)  + b_1 u'(1)  = c_1
%  a_N u(-1) + b_N u'(-1) = c_N
%
%
% INPUT 
% N        =  number of Chebyshev points in [-1,1]
% g        =  boundary condition matrix = [a_1 b_1 c_1; a_N b_N c_N]
% 
% OUTPUT  
% xt       =  Chebyshev points corresponding to rows and columns
%             of D1t and D2t
% D1t      =  1st derivative matrix incorporating bc
% D2t      =  2nd derivative matrix incorporating bc
% phip     =  1st and 2nd derivative of bc function at x=1
%             (array with 2 columns)
% phim     =  1st and 2nd derivative of bc function at x=-1 
%             (array with 2 columns)

% S.C. Reddy, J.A.C. Weideman  1998


% Get differentiation matrices

   [x,DM]=chebdif(N,2);
   D0=eye(N,N);
   D1=DM(:,:,1);
   D2=DM(:,:,2);

% extract boundary condition coefficients

   a1=g(1,1); b1=g(1,2); c1=g(1,3);
   aN=g(2,1); bN=g(2,2); cN=g(2,3);

% Case 0: Invalid boundary condition information

if ((a1==0 & b1==0) | (aN==0 & bN==0)),

   fprintf('Invalid boundary condition information (no output) \n');


elseif (b1==0 & bN==0)          % Dirichlet/Dirichlet 

   J=2:N-1;
   K=(2:N-1)';
   D1t=D1(J,K);
   D2t=D2(J,K);
   phip=c1*[D1(K,1) D2(K,1)]/a1;          % phi_+
   phim=cN*[D1(K,N) D2(K,N)]/aN;          % phi_- 
   xt=x(K);                               % node vector 

elseif (b1~=0 & bN==0),         % Dirichlet x=-1, Robin x=1

   J=2:N-1; 
   K=(1:N-1)';
   xjrow=2*sin((J-1)*pi/2/(N-1)).^2;      % 1-x_j, using trig identity
   xkcol=2*sin((K-1)*pi/2/(N-1)).^2;      % 1-x_k, using trig identity
   oner=ones(size(xkcol));                % column of ones

   fac0 = oner*(1./xjrow);                % matrix -1/(1-x_j)
   fac1 = xkcol*(1./xjrow);               % matrix (1-x_k)/(1-x_j)
   D1t = fac1.*D1(K,J)-fac0.*D0(K,J);
   D2t = fac1.*D2(K,J)-2*fac0.*D1(K,J); 

   cfac = D1(1,1)+a1/b1;                  % compute phi'_1, phi''_1
   fcol1 = -cfac*D0(K,1)+(1+cfac*xkcol).*D1(K,1);
   fcol2 = -2*cfac*D1(K,1)+(1+cfac*xkcol).*D2(K,1);
   D1t  = [fcol1 D1t];                    
   D2t  = [fcol2 D2t];                    

   phim = xkcol.*D1(K,N)/2-D0(K,N)/2;     % phi'_-, phi''_- 
   phim = cN*[phim xkcol.*D2(K,N)/2-D1(K,N)]/aN;

   
   phip= -xkcol.*D1(K,1)+D0(K,1);         % phi'_+, phi''_+ 
   phip= c1*[phip -xkcol.*D2(K,1)+2*D1(K,1)]/b1;

   xt = x(K);                             % node vector

elseif (b1==0 & bN~=0),

% Case 3: Dirichlet at x=1 and Neumann or Robin boundary x=-1.

   J=2:N-1; 
   K=(2:N)';
   xjrow=2*cos((J-1)*pi/2/(N-1)).^2;      % 1+x_j, using trig identity
   xkcol=2*cos((K-1)*pi/2/(N-1)).^2;      % 1+x_k, using trig identity
   oner=ones(size(xkcol));                % column of ones

   fac0 = oner*(1./xjrow);                % matrix 1/(1+x_j)
   fac1 = xkcol*(1./xjrow);               % matrix (1+x_k)/(1+x_j)
   D1t = fac1.*D1(K,J)+fac0.*D0(K,J);
   D2t = fac1.*D2(K,J)+2*fac0.*D1(K,J); 

   cfac = D1(N,N)+aN/bN;                  % compute phi'_N, phi''_N
   lcol1 = -cfac*D0(K,N)+(1-cfac*xkcol).*D1(K,N);
   lcol2 = -2*cfac*D1(K,N)+(1-cfac*xkcol).*D2(K,N);
   D1t  = [D1t lcol1];                    
   D2t  = [D2t lcol2];                

   phip= xkcol.*D1(K,1)/2+D0(K,1);        % compute phi'_+,phi''_+
   phip= c1*[phip xkcol.*D2(K,1)/2+D1(K,1)]/a1;

   phim= xkcol.*D1(K,N)+D0(K,N);          % compute phi'_-,phi''_-
   phim= cN*[phim xkcol.*D2(K,N)+2*D1(K,N)]/bN;

   xt = x(K);                             % node vector

elseif (b1~=0 & bN~=0),

% Case 4: Neumann or Robin boundary conditions at both endpoints. 

   J=2:N-1;
   K=(1:N)';
   xkcol0=sin((K-1)*pi/(N-1)).^2;             % 1-x_k^2 using trig identity
   xkcol1=-2*x(K);                            % -2*x_k 
   xkcol2=-2*ones(size(xkcol0));              % -2
   xjrow=1./(sin((J-1)*pi/(N-1)).^2);         % 1-x_j^2 using trig identity

   fac0=xkcol0*xjrow;
   fac1=xkcol1*xjrow;
   fac2=xkcol2*xjrow;

   D1t=fac0.*D1(K,J)+fac1.*D0(K,J);
   D2t=fac0.*D2(K,J)+2*fac1.*D1(K,J)+fac2.*D0(K,J);

   omx=sin((K-1)*pi/2/(N-1)).^2;              % (1-x_k)/2 
   opx=cos((K-1)*pi/2/(N-1)).^2;              % (1+x_k)/2

   r0=opx+(0.5+D1(1,1)+a1/b1)*xkcol0/2;       % compute phi'_1, phi''_1
   r1=0.5-(0.5+D1(1,1)+a1/b1)*x;
   r2=-0.5-D1(1,1)-a1/b1;
   rcol1=r0.*D1(K,1)+r1.*D0(K,1);
   rcol2=r0.*D2(K,1)+2*r1.*D1(K,1)+r2.*D0(K,1);

   l0=omx+(0.5-D1(N,N)-aN/bN)*xkcol0/2;       % compute phi'_N, phi''_N
   l1=-0.5+(D1(N,N)+aN/bN-0.5)*x;
   l2=D1(N,N)+aN/bN-0.5;
   lcol1=l0.*D1(K,N)+l1.*D0(K,N);
   lcol2=l0.*D2(K,N)+2*l1.*D1(K,N)+l2.*D0(K,N);

   D1t=[rcol1 D1t lcol1];
   D2t=[rcol2 D2t lcol2];

   phim1=(xkcol0.*D1(K,N)+xkcol1.*D0(K,N))/2; 
   phim2=(xkcol0.*D2(K,N)+2*xkcol1.*D1(K,N)+xkcol2.*D0(K,N))/2;
   phim=cN*[phim1 phim2]/bN;                 % compute phi'_-, phi''_-

   phip1=(-xkcol0.*D1(K,1)-xkcol1.*D0(K,1))/2;
   phip2=(-xkcol0.*D2(K,1)-2*xkcol1.*D1(K,1)-xkcol2.*D0(K,1))/2;
   phip=c1*[phip1 phip2]/b1;                 % compute phi'_+, phi''_+

   xt=x(K);                                  % node vector

end

endfunction



function [uic,vic] = pvtfinicxy(elemat,xi,yi)

%
% This function works out the velocity at some
% evaluation points by the POINT vortex sheet 
% elements with properties given in elemat.
%
%
% FORMAT:
%  [uic,vic] = pvtfinic(elemat,xi,yi)
%
% WHERE:
%  elemat = an P*5 matrix, where "P" is the number
%           of vortex-sheet panels, WHERE:
%    elemat(:,1) = x-coordinate of panel
%    elemat(:,2) = y-coordinate of panel
%    elemat(:,3) = element angle
%    elemat(:,4) = element length
%    elemat(:,5) = vortex sheet strength
%
%  xi = A COLUMN VECTOR giving the x-coordinate of evaluation points
%  yi = A COLUMN VECTOR giving the y-coordinate of evaluation points
%
%  [uic,vic] = a MATRIX giving the influence coefficients of the velcity
%              at the evaluation points.  Note, this must be multiplied
%              by the vortex strength column vector to give velocities at
%              the evaluation points
%
% NOTE: AS WITH THE POINT VORTICES, THE VELOCITY FIELD IS 
%       ANTICLOCKWISE-POSITIVE
%

% Get information from "elemat"
xj = elemat(:,1);
yj = elemat(:,2);
eleangj = elemat(:,3);
elewidj = elemat(:,4);

% Put the vectors into a form where bulk addition and subtraction can occur
[xj,xi] = meshgrid(xj,xi);
[yj,yi] = meshgrid(yj,yi);
eleangj = meshgrid(eleangj,ones(size(xj,1),1));
elewidj = meshgrid(elewidj,ones(size(xj,1),1));

% Get the absolute distance between point "i" and "j"
%r = sqrt((yi - yj).^2 + (xi - xj).^2);

% Calculate the relative x and y distances from element "j"
rir = (((yi - yj)).*(cos(eleangj))-((xi - xj)).*(sin(eleangj)));
rur = (((xi - xj)).*(cos(eleangj))+((yi - yj)).*(sin(eleangj)));

% Calculate the relative x-component of velocity at "i" due to source sink
% MUST ENSURE THAT WE DO NOT DIVIDE BY ZERO OR LOG OF ZERO!
num = (rur + (elewidj/2)).^2 + rir.^2;
tmp1 = (num==0);
num = num + tmp1;
tmp1 = 1-tmp1;
den = (rur - (elewidj/2)).^2 + rir.^2;
tmp2 = (den==0);
den = den + tmp2;
tmp2 = 1-tmp2;
uic = (1/(4*pi))*(log(num./den)).*(tmp1).*(tmp2);

% Calculate the relative y-component of velocity at "i" due to source sink
% MUST ENSURE THAT WE DO NOT DIVIDE BY ZERO!
tmp1 = (rir==0);
rir = rir + tmp1;
tmp1 = 1 - tmp1;
vic = (1/(2*pi))*(-1*atan((rur - (elewidj/2))./rir) + atan((rur + (elewidj/2))./rir)).*(tmp1);

% Calculate the relative x-component of velocity at "i" due to vortex sheet
uicr = (-1)*vic;
% Calculate the relative y-component of velocity at "i" due to vortex sheet
vicr = uic;

% Calculate the actual global x-component of velocity at "i" due to vortex sheet
uic = uicr.*cos(eleangj) - vicr.*sin(eleangj);
% Calculate the actual global y-component of velocity at "i" due to vortex sheet
vic = vicr.*cos(eleangj) + uicr.*sin(eleangj);
endfunction


function [uic,vic] = pssfinicxy(elemat,xi,yi)

%
% This function works out the velocity at some
% evaluation points by the POINT Source/Sink sheet 
% elements with properties given in elemat.
% 
% FORMAT:
%  [uic,vic] = pssfinsh(elemat,xi,yi)
%
% WHERE:
%  elemat = an P*5 matrix, where "P" is the number
%           of vortex-sheet panels, WHERE:
%    elemat(:,1) = x-coordinate of panel
%    elemat(:,2) = y-coordinate of panel
%    elemat(:,3) = element angle
%    elemat(:,4) = element length
%    elemat(:,5) = Source/Sink sheet strength
%
%  xi = A COLUMN VECTOR giving the x-coordinate of evaluation points
%  yi = A COLUMN VECTOR giving the y-coordinate of evaluation points
%
%  [uic,vic] = COLUMN vectors giving the x and y velocities at the 
%             evaluation points.
%
% NOTE: AS WITH THE POINT VORTICES, THE VELOCITY FIELD IS 
%       ANTICLOCKWISE-POSITIVE
%

% Get information from "elemat"
xj = elemat(:,1);
yj = elemat(:,2);
eleangj = elemat(:,3);
elewidj = elemat(:,4);

% Put the vectors into a form where bulk addition and subtraction can occur
[xj,xi] = meshgrid(xj,xi);
[yj,yi] = meshgrid(yj,yi);
eleangj = meshgrid(eleangj,ones(size(xj,1),1));
elewidj = meshgrid(elewidj,ones(size(xj,1),1));

% Get the absolute distance between point "i" and "j"
%r = sqrt((yi - yj).^2 + (xi - xj).^2);

% Calculate the relative x and y distances from element "j"
rir = (((yi - yj)).*(cos(eleangj))-((xi - xj)).*(sin(eleangj)));
rur = (((xi - xj)).*(cos(eleangj))+((yi - yj)).*(sin(eleangj)));

% Calculate the relative x-component of velocity at "i" due to source sink
% MUST ENSURE THAT WE DO NOT DIVIDE BY ZERO OR LOG OF ZERO!
num = (rur + (elewidj/2)).^2 + rir.^2;
tmp1 = (num==0);
num = num + tmp1;
tmp1 = 1-tmp1;
den = (rur - (elewidj/2)).^2 + rir.^2;
tmp2 = (den==0);
den = den + tmp2;
tmp2 = 1-tmp2;
uicr = (1/(4*pi))*(log(num./den)).*(tmp1).*(tmp2);

% Calculate the relative y-component of velocity at "i" due to source sink
% MUST ENSURE THAT WE DO NOT DIVIDE BY ZERO!
tmp1 = (rir==0);
rir = rir + tmp1;
tmp1 = 1 - tmp1;
vicr = (1/(2*pi))*(-1*atan((rur - (elewidj/2))./rir) + atan((rur + (elewidj/2))./rir)).*(tmp1);

% Calculate the actual global x-component of velocity at "i"
uic = uicr.*cos(eleangj) - vicr.*sin(eleangj);
% Calculate the actual global y-component of velocity at "i"
vic = vicr.*cos(eleangj) + uicr.*sin(eleangj);
endfunction


function varname = mvarname(path,prefix,x)

%
% This is a function that outputs a variable name given the
% input of a "path" to the output directory, a "prefix" to all
% of the file names and a stepping variable "x"
%
%

y = num2str(x);

if x < 10
    varname = strcat(prefix,'000',y);
elseif x < 100 & x >= 10
    varname = strcat(prefix,'00',y);
elseif x < 1000 & x >= 100
    varname = strcat(prefix,'0',y);
elseif x < 10000 & x >= 1000
    varname = strcat(prefix,y);
else
    error('Stepping varible is greater than limit of 10000')
end

varname = strcat([path '/' varname]);
endfunction



function [dndx] = dndx4(nin,dx,econ1,econ2)

%
% Calculates the value of dndx^2 at each node on a plate
% defined by vertical positions "nin"
% 
% The output is a square matrix "dndx" which must
% be post-multiplied by the vector "n" in order to
% get the vector of derivatives.
%
%

x = size(nin,1); % Get teh size of vector "nin"

% Calculate the derivative across the bulk of the
% points
dndx = diag((ones(x,1)*(6)),0) + diag(ones(x-1,1)*(-4),1) + diag(ones(x-1,1)*(-4),-1) + diag(ones(x-2,1)*(1),2) + diag(ones(x-2,1)*(1),-2);

% Calculate the end conditions

% LHS END
if econ1 == 1
   % Built-in end conditions (dndx1 = 0)
   dndx(1,:) = zeros(1,size(dndx,2));
   dndx(2,2) = 7;
elseif econ1 == 2
   % Hinged end conditions (dndx2 = 0)
   dndx(1,:) = zeros(1,size(dndx,2));
   dndx(2,2) = 5;
end

% RHS End
k = size(dndx,1);
if econ2 == 1
   % Built-in end conditions (dndx1 = 0)
   dndx(k,:) = zeros(1,size(dndx,2));
   dndx(k-1,k-1) = 7;
elseif econ2 == 2
   % Hinged end conditions (dndx2 = 0)
   dndx(k,:) = zeros(1,size(dndx,2));;
   dndx(k-1,k-1) = 5;
end

dndx = dndx*((1/(dx^4)));

endfunction



function [dndx] = dndx2(nin,dx,econ1,econ2)

%
% Calculates the value of dndx^2 at each node on a plate
% defined by vertical positions "nin"
% 
% The output is a square matrix "dndx" which must
% be post-multiplied by the vector "n" in order to
% get the vector of derivatives.
%
%

x = size(nin,1); % Get teh size of vector "nin"

% Calculate the derivative across the bulk of the
% points

dndx = diag((ones(x,1)*(-2)),0) + diag(ones(x-1,1)*1,1) + diag(ones(x-1,1)*1,-1);


% Calculate the end conditions

% LHS END
if econ1 == 1
   % Built-in end conditions (dndx1 = 0)
   dndx(1,1) = -2;
   dndx(1,2) = 2;
elseif econ1 == 2
   % Hinged end conditions (dndx2 = 0)
   dndx(1,:) = zeros(1,size(dndx,2));  
end

% RHS End
k = size(dndx,1);
if econ2 == 1
   % Built-in end conditions (dndx1 = 0)
   dndx(k,k) = 2;
   dndx(k,k-1) = -2;
elseif econ2 == 2
   % Hinged end conditions (dndx2 = 0)
   dndx(k,:) = zeros(1,size(dndx,2));
end

dndx = dndx*((1/(dx^2)));

endfunction


function [dndx] = dndx1uw(nin,dx,econ1,econ2)

%
% Calculates the value of dndx at each node on a plate
% defined by vertical positions "nin"
% 
% The output is a square matrix "dndx" which must
% be post-multiplied by the vector "n" in order to
% get the vector of derivatives.
%
%

% Get teh size of vector "nin"
x = size(nin,1);

% Calculate the derivative across the bulk of the
% points

dndx = diag(ones(x,1)*1,0) + diag(ones(x-1,1)*-1,-1);


% Calculate the end conditions

% END CONDITIONS NOT FULLY IMPLEMENTED YET!

dndx = dndx*((1/dx));

endfunction


function [dndx] = dndx1(nin,dx,econ1,econ2)

%
% Calculates the value of dndx at each node on a plate
% defined by vertical positions "nin"
% 
% The output is a square matrix "dndx" which must
% be post-multiplied by the vector "n" in order to
% get the vector of derivatives.
%
%

% Get teh size of vector "nin"
x = size(nin,1);

% Calculate the derivative across the bulk of the
% points

dndx = diag(ones(x-1,1)*1,1) + diag(ones(x-1,1)*-1,-1);


% Calculate the end conditions

% LHS END
if econ1 == 1
   % Built-in end conditions (dndx1 = 0)
   dndx(1,:) = zeros(1,size(dndx,2));
elseif econ1 == 2
   % Hinged end conditions (dndx2 = 0)
   dndx(1,2) = 2;  
end

% RHS End
k = size(dndx,1);
if econ2 == 1
   % Built-in end conditions (dndx1 = 0)
   dndx(k,:) = zeros(1,size(dndx,2));
elseif econ2 == 2
   % Hinged end conditions (dndx2 = 0)
   dndx(k,k-1) = -2;
end

dndx = dndx*((1/(2*dx)));

endfunction



function [x, DM] = chebdif(N, M)

%  The function [x, DM] =  chebdif(N,M) computes the differentiation 
%  matrices D1, D2, ..., DM on Chebyshev nodes. 
% 
%  Input:
%  N:        Size of differentiation matrix.        
%  M:        Number of derivatives required (integer).
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced 
%  accuracy suggested by W. Don and S. Solomonoff in 
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric 
%  identities to avoid the computation of differences 
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%  Note added May 2003:  It may, in fact, be slightly better not to
%  implement the strategies (a) and (b).   Please consult the following
%  paper for details:   "Spectral Differencing with a Twist", by
%  R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp. 

%  J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by 
%  JACW, May 2003.

     I = eye(N);                          % Identity matrix.     
     L = logical(I);                      % Logical identity matrix.

    n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.

     k = [0:N-1]';                        % Compute theta vector.
    th = k*pi/(N-1);

     x = sin(pi*[N-1:-2:1-N]'/(2*(N-1))); % Compute Chebyshev points.

     T = repmat(th/2,1,N);                
    DX = 2*sin(T'+T).*sin(T'-T);          % Trigonometric identity. 
    DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];   % Flipping trick. 
 DX(L) = ones(N,1);                       % Put 1's on the main diagonal of DX.

     C = toeplitz((-1).^k);               % C is the matrix with 
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;     % entries c(k)/c(j)
C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;

     Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))  
  Z(L) = zeros(N,1);                      % with zeros on the diagonal.

     D = eye(N);                          % D contains diff. matrices.
                                          
for ell = 1:M
          D = ell*Z.*(C.*repmat(diag(D),1,N) - D); % Off-diagonals
       D(L) = -sum(D');                            % Correct main diagonal of D
DM(:,:,ell) = D;                                   % Store current D in DM
end

endfunction


function [x, D4] = cheb4c(N)

%  The function [x, D4] =  cheb4c(N) computes the fourth 
%  derivative matrix on Chebyshev interior points, incorporating 
%  the clamped boundary conditions u(1)=u'(1)=u(-1)=u'(-1)=0.
%
%  Input:
%  N:     N-2 = Order of differentiation matrix.  
%               (The interpolant has degree N+1.)
%
%  Output:
%  x:      Interior Chebyshev points (vector of length N-2)
%  D4:     Fourth derivative matrix  (size (N-2)x(N-2))
%
%  The code implements two strategies for enhanced 
%  accuracy suggested by W. Don and S. Solomonoff in 
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric 
%  identities to avoid the computation of differences 
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
   
%  J.A.C. Weideman, S.C. Reddy 1998.

    I = eye(N-2);                   % Identity matrix.
    L = logical(I);                 % Logical identity matrix.

   n1 = floor(N/2-1);               % n1, n2 are indices used 
   n2 = ceil(N/2-1);                % for the flipping trick.

    k = [1:N-2]';                   % Compute theta vector.
   th = k*pi/(N-1);                 

    x = sin(pi*[N-3:-2:3-N]'/(2*(N-1))); % Compute interior Chebyshev points.

    s = [sin(th(1:n1)); flipud(sin(th(1:n2)))];   % s = sin(theta)
                               
alpha = s.^4;                       % Compute weight function
beta1 = -4*s.^2.*x./alpha;          % and its derivatives.
beta2 =  4*(3*x.^2-1)./alpha;   
beta3 = 24*x./alpha;
beta4 = 24./alpha;
    B = [beta1'; beta2'; beta3'; beta4'];

    T = repmat(th/2,1,N-2);                
   DX = 2*sin(T'+T).*sin(T'-T);     % Trigonometric identity 
   DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];   % Flipping trick. 
DX(L) = ones(N-2,1);                % Put 1's on the main diagonal of DX.

   ss = s.^2.*(-1).^k;              % Compute the matrix with entries
    S = ss(:,ones(1,N-2));          % c(k)/c(j)
    C = S./S';                      

    Z = 1./DX;                      % Z contains entries 1/(x(k)-x(j)).
 Z(L) = zeros(size(x));             % with zeros on the diagonal.

    X = Z';                         % X is same as Z', but with 
 X(L) = [];                         % diagonal entries removed.
    X = reshape(X,N-3,N-2);

    Y = ones(N-3,N-2);              % Initialize Y and D vectors.
    D = eye(N-2);                   % Y contains matrix of cumulative sums,
                                    % D scaled differentiation matrices.
for ell = 1:4
          Y = cumsum([B(ell,:); ell*Y(1:N-3,:).*X]); % Recursion for diagonals
          D = ell*Z.*(C.*repmat(diag(D),1,N-2)-D);   % Off-diagonal
       D(L) = Y(N-2,:);                              % Correct the diagonal
DM(:,:,ell) = D;                                     % Store D in DM
end

   D4 = DM(:,:,4);                  % Extract fourth derivative matrix

endfunction



function [] = calcIww()

%------------------------------------
% Get the influence of every wall (source/sink) element on each other
% Vector of elements = {compliant;top wall;bottom1;bottom2}
INww = zeros(length(xc)+length(cwu)+length(cwb1)+length(cwb2));
ITww = INww;
xi = [cwu,cwb1,xc,cwb2]';
yi = [ones(size(cwu)),-1*ones(size(cwb1)),-1*ones(size(xc)),-1*ones(size(cwb2))]';
for c = 1:size(INww,1)
    disp(['Calculating influence (wall->wall) for element ' num2str(c) ' of ' num2str(size(INww,1))]);
    cc = 1;
    for s = 1:length(xc)
        if cc == c
            uic = 0;vic = 0.5;
        else
            [uic,vic] = pssfinicxy([xc(s) -1 0 dxc],xi(c),yi(c));
        end
        INww(c,cc) = vic;ITww(c,cc) = uic;
        cc = cc + 1;
    end
    for s = 1:length(cwu)
        if cc == c
            uic = 0;vic = -0.5;
        else
            [uic,vic] = pssfinicxy([cwu(s) 1 0 dxu],xi(c),yi(c));
        end
        INww(c,cc) = vic;ITww(c,cc) = uic;
        cc = cc + 1;
    end
    for s = 1:length(cwb1)
        if cc == c
            uic = 0;vic = 0.5;
        else
            [uic,vic] = pssfinicxy([cwb1(s) -1 0 dxb1],xi(c),yi(c));
        end
        INww(c,cc) = vic;ITww(c,cc) = uic;
        cc = cc + 1;
    end
    for s = 1:length(cwb2)
        if cc == c
            uic = 0;vic = 0.5;
        else
            [uic,vic] = pssfinicxy([cwb2(s) -1 0 dxb2],xi(c),yi(c));
        end
        INww(c,cc) = vic;ITww(c,cc) = uic;
        cc = cc + 1;
    end
end

endfunction


function [] = calcIwf()

%------------------------------------
% Get the influence of every wall (source/sink) element on the fluid
% Vector of elements = {compliant;top wall;bottom1;bottom2}
INwf = zeros(length(pcxxv),length(xc)+length(cwu)+length(cwb1)+length(cwb2));
ITwf = INwf;
xi = pcxxv;
yi = pcyyv;
for c = 1:size(INwf,1)
    disp(['Calculating influence (wall->fluid) for element ' num2str(c) ' of ' num2str(size(INwf,1))]);
    cc = 1;
    for s = 1:length(xc)
        [uic,vic] = pssfinicxy([xc(s) -1 0 dxc],xi(c),yi(c));
        INwf(c,cc) = vic;ITwf(c,cc) = uic;
        cc = cc + 1;
    end
    for s = 1:length(cwu)
        [uic,vic] = pssfinicxy([cwu(s) 1 0 dxu],xi(c),yi(c));
        INwf(c,cc) = vic;ITwf(c,cc) = uic;
        cc = cc + 1;
    end
    for s = 1:length(cwb1)
        [uic,vic] = pssfinicxy([cwb1(s) -1 0 dxb1],xi(c),yi(c));
        INwf(c,cc) = vic;ITwf(c,cc) = uic;
        cc = cc + 1;
    end
    for s = 1:length(cwb2)
        [uic,vic] = pssfinicxy([cwb2(s) -1 0 dxb2],xi(c),yi(c));
        INwf(c,cc) = vic;ITwf(c,cc) = uic;
        cc = cc + 1;
    end
end

endfunction


function [] = calcIfw()

%------------------------------------
% Get the influence of every fluid element on the wall elements
% Vector of elements = {from top-left -> top-to-bottom -> to bottom right}
INfw = zeros(length(xc)+length(cwu)+length(cwb1)+length(cwb2),length(pcxxv));
ITfw = INfw;
xi = [xc,cwu,cwb1,cwb2]';
yi = [ones(size(cwu)),-1*ones(size(cwb1)),-1*ones(size(xc)),-1*ones(size(cwb2))]';
for c = 1:size(INfw,1)
    disp(['Calculating influence (fluid->wall) for element ' num2str(c) ' of ' num2str(size(INfw,1))]);
    cc = 1;
    for s = 1:length(pcxxv)
        [uic,vic] = pssfinicxy([pcxxv(s) pcyyv(s) 0 dxf],xi(c),yi(c));
        INfw(c,cc) = vic;ITfw(c,cc) = uic;
        cc = cc + 1;
    end
end

endfunction


function [] = calcIff()

%------------------------------------
% Get the influence of every fluid element on themselves
% Vector of elements = {from top-left -> top-to-bottom -> to bottom right}
INff = zeros(length(pcxxv));
ITff = INff;
xi = pcxxv;
yi = pcyyv;
for c = 1:size(INff,1)
    disp(['Calculating influence (fluid->fluid) for element ' num2str(c) ' of ' num2str(size(INff,1))]);
    cc = 1;
    for s = 1:length(pcxxv)
        if cc == c
            uic = 0;vic = 0;
        else
            [uic,vic] = pvtfinicxy([pcxxv(s) pcyyv(s) 0 dxf],xi(c),yi(c));
        end
        INff(c,cc) = vic;ITff(c,cc) = uic;
        cc = cc + 1;
    end
end

endfunction




