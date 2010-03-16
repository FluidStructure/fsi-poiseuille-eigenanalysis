function [dndx] = dndx1uw(nin,dx,econ1,econ2)

%
% Calculates the value of dndx at each node on a plate
% defined by vertical positions "nin"
% 
% The output is a square matrix "dndx" which must
% be post-multiplied by the vector "n" in order to
% get the vector of derivatives.
%
%

% Get teh size of vector "nin"
x = size(nin,1);

% Calculate the derivative across the bulk of the
% points

dndx = diag(ones(x,1)*1,0) + diag(ones(x-1,1)*-1,-1);


% Calculate the end conditions

% END CONDITIONS NOT FULLY IMPLEMENTED YET!

dndx = dndx*((1/dx));
