function [uic,vic] = pssfinicxy(elemat,xi,yi)

%
% This function works out the velocity at some
% evaluation points by the POINT Source/Sink sheet 
% elements with properties given in elemat.
% 
% FORMAT:
%  [uic,vic] = pssfinsh(elemat,xi,yi)
%
% WHERE:
%  elemat = an P*5 matrix, where "P" is the number
%           of vortex-sheet panels, WHERE:
%    elemat(:,1) = x-coordinate of panel
%    elemat(:,2) = y-coordinate of panel
%    elemat(:,3) = element angle
%    elemat(:,4) = element length
%    elemat(:,5) = Source/Sink sheet strength
%
%  xi = A COLUMN VECTOR giving the x-coordinate of evaluation points
%  yi = A COLUMN VECTOR giving the y-coordinate of evaluation points
%
%  [uic,vic] = COLUMN vectors giving the x and y velocities at the 
%             evaluation points.
%
% NOTE: AS WITH THE POINT VORTICES, THE VELOCITY FIELD IS 
%       ANTICLOCKWISE-POSITIVE
%

% Get information from "elemat"
xj = elemat(:,1);
yj = elemat(:,2);
eleangj = elemat(:,3);
elewidj = elemat(:,4);

% Put the vectors into a form where bulk addition and subtraction can occur
[xj,xi] = meshgrid(xj,xi);
[yj,yi] = meshgrid(yj,yi);
eleangj = meshgrid(eleangj,ones(size(xj,1),1));
elewidj = meshgrid(elewidj,ones(size(xj,1),1));

% Get the absolute distance between point "i" and "j"
%r = sqrt((yi - yj).^2 + (xi - xj).^2);

% Calculate the relative x and y distances from element "j"
rir = (((yi - yj)).*(cos(eleangj))-((xi - xj)).*(sin(eleangj)));
rur = (((xi - xj)).*(cos(eleangj))+((yi - yj)).*(sin(eleangj)));

% Calculate the relative x-component of velocity at "i" due to source sink
% MUST ENSURE THAT WE DO NOT DIVIDE BY ZERO OR LOG OF ZERO!
num = (rur + (elewidj/2)).^2 + rir.^2;
tmp1 = (num==0);
num = num + tmp1;
tmp1 = 1-tmp1;
den = (rur - (elewidj/2)).^2 + rir.^2;
tmp2 = (den==0);
den = den + tmp2;
tmp2 = 1-tmp2;
uicr = (1/(4*pi))*(log(num./den)).*(tmp1).*(tmp2);

% Calculate the relative y-component of velocity at "i" due to source sink
% MUST ENSURE THAT WE DO NOT DIVIDE BY ZERO!
tmp1 = (rir==0);
rir = rir + tmp1;
tmp1 = 1 - tmp1;
vicr = (1/(2*pi))*(-1*atan((rur - (elewidj/2))./rir) + atan((rur + (elewidj/2))./rir)).*(tmp1);

% Calculate the actual global x-component of velocity at "i"
uic = uicr.*cos(eleangj) - vicr.*sin(eleangj);
% Calculate the actual global y-component of velocity at "i"
vic = vicr.*cos(eleangj) + uicr.*sin(eleangj);
