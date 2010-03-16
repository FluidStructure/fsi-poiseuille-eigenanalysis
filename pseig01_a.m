function [] = pseig01_a()

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

