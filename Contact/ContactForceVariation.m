function [Fc,Kc] = ContactForceVariation(Body1,Body2,penalty,approach)
    

    % Backup original coordinates
    u1_backup = Body1.u;
    u2_backup = Body2.u;

    sqrtEps = sqrt(eps);
    nx = Body1.nx + Body2.nx;
    Kc = zeros(nx,nx);

    Fc = ContactForce(Body1,Body2,penalty,approach); % Body1 forces from the projection of Body2

    I_vec=zeros(nx,1);
    for ii = 1:nx
        I_vec(ii)=1;
        h = 2*sqrtEps;

        % this split is to distribute coord. between bodies 
        Body1.u = u1_backup - h*I_vec(1:Body1.nx); 
        Body2.u = u2_backup - h*I_vec(1+Body1.nx:end);
        Fch = ContactForce(Body1,Body2,penalty,approach); % force due to variation

        Kc(:,ii) = (Fc - Fch) / h; 
       
        % returning all back to normal
        I_vec(ii)=0;
        
    end    
    
    Body1.u = u1_backup;
    Body2.u = u2_backup;