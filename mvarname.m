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
