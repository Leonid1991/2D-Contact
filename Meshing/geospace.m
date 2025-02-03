function [Y] = geospace(X1,X2,N,G)
% GEOSPACE Geometrically spaced vector.
%
%    Create a geometrically spaced vector. Similar to linspace 
%    but each step is G times bigger than the previous one.
%
%    GEOSPACE(X1, X2) generates a row vector of 100 linearly
%    equally spaced points between X1 and X2.
% 
%    GEOSPACE(X1, X2, N) generates N linearly spaced points 
%    between X1 and X2. If N<2 then X2 is returned.
%
%    GEOSPACE(X1, X2, N, G) generates N points between X1 and X2
%    with a spacing defined by G. Each step is G times bigger 
%    than the previous step.
% 
%   Example:
%       % Create a 2D mesh with geometric cell spacing.
%       xx = geospace(20,2000,10,1.2); 
%       [xx,yy] = meshgrid(xx,xx'); 
%       zz = diff(xx(2:end,:),[],2).*diff(yy(:,2:end),[],1); 
%       surf(xx,yy,ones(size(xx)),zz); 
%
%    See also linspace, logspace, :.
%

% check number of arguments
if nargin<3,   %% nargin-Number of function input arguments
  N = 100;
end;
if nargin<4,
  G = 1;
end;

if N < 2,  % too few points specified
  Y = X2;
else
  if (G==1),
    Y = X1:(X2-X1)/(N-1):X2;   % Prevent divide by zero
  else
    % Calculate geometrically spaced points
    Y = X1 + (cumprod(G*ones([1 N]))-G)/(G^N-G)*(X2-X1);  %% cumprod-Cumulative product
  end;
end;