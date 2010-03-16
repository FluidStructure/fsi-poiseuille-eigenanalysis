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