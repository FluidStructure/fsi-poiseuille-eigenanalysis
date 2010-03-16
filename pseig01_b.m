function [] = pseig01_b()

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


