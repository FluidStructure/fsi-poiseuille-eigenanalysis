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
