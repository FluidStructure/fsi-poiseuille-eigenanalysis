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
