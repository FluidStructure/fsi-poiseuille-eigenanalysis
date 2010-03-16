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