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