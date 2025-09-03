function [Fc,Kc,GapNab]= Contact(Body1,Body2,penalty,approach)

% Initialize the global contact forces
Fc = zeros(Body1.nx + Body2.nx,1);

% Initialize the stiffnesses
Kc = zeros(Body1.nx + Body2.nx, Body1.nx + Body2.nx);

% Initialize the gap
Gap = 0;

% Initialize the gap variation
GapNab = zeros(Body1.nx + Body2.nx,1);

if approach ~=0  
    
    [Fc,Kc,GapNab] = ContactVariation(Body1,Body2,penalty,approach);    

end


