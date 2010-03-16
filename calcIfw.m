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