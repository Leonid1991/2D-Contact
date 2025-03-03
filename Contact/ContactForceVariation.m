function [Fc,Kc] = ContactForceVariation(Body1,Body2,penalty,h,approach)

    nx = Body1.nx + Body2.nx;

    Kc = zeros(nx,nx);
    Fc = ContactForce(Body1,Body2,penalty,approach); % Body1 forces from the projection of Body2

    % variation of the variables
    I_vec=zeros(nx,1);
    for ii = 1:nx
        I_vec(ii)=1;
        
        % this split is to distribute coord. between bodies 
        Body1.u = Body1.u - h*I_vec(1:Body1.nx); 
        Body2.u = Body2.u - h*I_vec(1+Body1.nx:Body1.nx + Body2.nx);

        Fch = ContactForce(Body1,Body2,penalty,approach); % force due to variation
        
        Kc(:,ii) = (Fc - Fch) / h; 
    
        % returning all back to normal
        Body1.u = Body1.u + h*I_vec(1:Body1.nx);
        Body2.u = Body2.u + h*I_vec(1+Body1.nx:Body1.nx + Body2.nx);
        I_vec(ii)=0;
    end    