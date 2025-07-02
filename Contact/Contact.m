function [Fc,Kc]= Contact(Body1,Body2,penalty,h,approach)

% Initialize the global contact forces
Fc = zeros(Body1.nx + Body2.nx,1);
Gap = NaN;
% Initialize the stiffnesses
Kc = zeros(Body1.nx + Body2.nx, Body1.nx + Body2.nx);

if approach ~=0  
    
    [Fc,Kc] = ContactForceVariation(Body1,Body2,penalty,h,approach);    

end


