function [dndx] = dndx2(nin,dx,econ1,econ2)

%
% Calculates the value of dndx^2 at each node on a plate
% defined by vertical positions "nin"
% 
% The output is a square matrix "dndx" which must
% be post-multiplied by the vector "n" in order to
% get the vector of derivatives.
%
%

x = size(nin,1); % Get teh size of vector "nin"

% Calculate the derivative across the bulk of the
% points

dndx = diag((ones(x,1)*(-2)),0) + diag(ones(x-1,1)*1,1) + diag(ones(x-1,1)*1,-1);


% Calculate the end conditions

% LHS END
if econ1 == 1
   % Built-in end conditions (dndx1 = 0)
   dndx(1,1) = -2;
   dndx(1,2) = 2;
elseif econ1 == 2
   % Hinged end conditions (dndx2 = 0)
   dndx(1,:) = zeros(1,size(dndx,2));  
end

% RHS End
k = size(dndx,1);
if econ2 == 1
   % Built-in end conditions (dndx1 = 0)
   dndx(k,k) = 2;
   dndx(k,k-1) = -2;
elseif econ2 == 2
   % Hinged end conditions (dndx2 = 0)
   dndx(k,:) = zeros(1,size(dndx,2));
end

dndx = dndx*((1/(dx^2)));
