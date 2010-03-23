function [] = callOctaveSolvers()

more off

% Load the parameters
load VARS.mat

if strcmp(method,'Naive')
    load A.mat
else
    A = []
end

if strcmp(solver,'eigs')

    % Set some parameters
    nummods = 8;
    opts.maxit = 20000;
    opts.disp = 1;
    %opts.tol = eps*1e-3;
    % Set point sigma: NOTE If: Period = lamda/Uinf => Freq = Uinf/lamda => Omega = 2*pi*Freq;
    sigst = 'LR';
    %sigst = 0 - 1*i;

    disp('Starting calculation of eigenvalues...')

    % Load an existing solution for the initial guess (needs some work)
    %if exist('v0.mat') == 2
    %    load v0.mat v0
    %    opts.v0 = v0
    %end

    if isempty(A) == 0
        disp('Loading matrix A... into octave')
        load A.mat;
        disp('Done...')
        [Veigs,Deigs,flagc] = eigs(A,nummods,sigst,opts);
    else
        disp('Calling eigenvalue solver with FMM function...')
        [Veigs,Deigs,flagc] = eigs(@callPythonEIG(x,sigma),nummods,sigst,opts);
    end
    
    if flagc == 0
        disp('ALL EIGENVALUES CONVERGED :-)')
    else
        error('NOT ALL EIGENVALUES CONVERGED :-(')
    end
    evals = diag(Deigs)

    fname = ["results/evals_R" num2str(R) ".mat"]
    save('-v7',fname,'Veigs','evals')

    % Save a vector for an initial guess for the next step
    %v0 = Veigs(:,1);
    %save('-v7','v0.mat','v0')

elseif strcmp(solver,'ode45')

    % Time steps to dump out
    TSPAN = linspace(0,20,200);
    if deterministicBCs == 'True'
        Ny = chebN-2
    else
        Ny = chebN;
    end
    N = Ny*Nx;
    Y0 = zeros(N,1);

    nn = round(Nx/4)*Ny + round(Ny/2);
    Y0(nn,1) = 1e-2;	% Initialise with a spot disturbance (generally not good)

    %nny = round(Ny/2)
    %dx = LT/Nx;
    %pcx = linspace(0.5*dx,LT-0.5*dx,Nx);
    %for nn = 1:Nx
    %	Y0(((nn-1)*Ny)+nny) = 1 - cos(pcx(nn));		% Initialise with a sin-wave disturbance along centre-line
    %end

    vopt = odeset("RelTol", 1e-3, "AbsTol", 1e-3, "NormControl", "on", "MaxStep", 1e-3, "OutputFcn", @odeSaveVars);
    ode45 (@callPythonODE, TSPAN, Y0, vopt, A);
    %[TOUT,YOUT] = ode45(@callPythonODE,TSPAN,Y0);
    %save RESODE45.mat

end

endfunction



function [x] = callPythonODE(t,v,A)

    if isempty(A)
        disp(['CALLING THE FUNCTION AGAIN!!, T = ' num2str(t)])
        save -v7 v.mat v
        F = system('python psEigs.py',0);
        load -v7 x.mat;
%        x = xRe + i*xIm;
        delete('v.mat')
        delete('x.mat')
    else
        disp(['Solving IVP using saved matrix A (Naive method), T = ' num2str(t)])
        x = A*v;
    end

endfunction




function [x] = callPythonEIG(v,sigma)

    global c

    % Evaluate python function
    save -v7 v.mat v sigma c
    system('python psEigs.py');
    load -v7 x.mat;
    x = xRe + i*xIm;

    if exist("pcy.mat") == 0
        save -v7 pcy.mat pcy
    end

    %keyboard

    delete('v.mat')
    delete('x.mat')

    c = c + 1;

endfunction



function [varargout] = odeSaveVars (vt, vy, vflag, varargin)

    %# No input argument check is done for a higher processing speed
    %persistent vfigure;
    %persistent vtold; 
    %persistent vyold;
    persistent vcounter;

    if (strcmp (vflag, 'init')) 
        %# Nothing to return, vt is either the time slot [tstart tstop]
        %# or [t0, t1, ..., tn], vy is the inital value vector 'vinit'
        %vfigure = figure;
        %vtold = vt(1,1); vyold = vy(:,1); 
        vcounter = 0;

        t = vt(1,1);y = vy(:,1);
        fwnam = mvarname('results','TStep',vcounter);
        fwnam = strcat([fwnam '.mat']);
        disp(['Saving result... to: '  fwnam])
        save('-v7',fwnam,'t','y')

    elseif (isempty (vflag))
        %# Return something in varargout{1}, either false for 'not stopping
        %# the integration' or true for 'stopping the integration'
        vcounter = vcounter + 1; 
        %figure (vfigure);
        %vtold(vcounter,1) = vt(1,1);vyold(:,vcounter) = vy(:,1);

        t = vt(1,1);y = vy(:,1);
        fwnam = mvarname('results','TStep',vcounter);
        fwnam = strcat([fwnam '.mat']);
        disp(['Saving result... to: '  fwnam])  
        save('-v7',fwnam,'t','y')

        %plot (vtold, vyold, '-o', 'markersize', 1); drawnow;
        varargout{1} = false;

    elseif (strcmp (vflag, 'done')) 
        %# Cleanup has to be done, clear the persistent variables because
        %# we don't need them anymore
        %clear ('vfigure', 'vtold', 'vyold', 'vcounter');
        clear ('vcounter');
    end

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

