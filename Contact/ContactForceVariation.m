function [Fc,Kc] = ContactForceVariation(Body1,Body2,penalty,approach)
     
    sqrtEps = sqrt(eps);
    % Backup original coordinates
    u1_backup = Body1.u;
    u2_backup = Body2.u;
    u_backup = [u1_backup; u2_backup];


    nx = Body1.nx + Body2.nx;

    Kc = zeros(nx,nx);

    % disp('calcuation');
    Fc = ContactForce2(Body1,Body2,penalty,approach); % Body1 forces from the projection of Body2

    % variation of the variables
    I_vec=zeros(nx,1);

    % disp('variation');
    for ii = 1:nx
        I_vec(ii)=1;
        h = max(sqrtEps * abs(u_backup(ii)) , 2*sqrtEps);

        % this split is to distribute coord. between bodies 
        Body1.u = u1_backup - h*I_vec(1:Body1.nx); 
        Body2.u = u2_backup - h*I_vec(1+Body1.nx:Body1.nx + Body2.nx);
        Fch = ContactForce2(Body1,Body2,penalty,approach); % force due to variation

        Kc(:,ii) = (Fc - Fch) / h; 

        I_vec(ii)=0;
    end    
    
       
    % for ii = 1:nx
    %     I_vec(ii)=1;
    %     h = max(sqrtEps * abs(u_backup(ii)) , 2*sqrtEps); 
    % 
    %     % this split is to distribute coord. between bodies 
    %     Body1.u = u1_backup - h*I_vec(1:Body1.nx); 
    %     Body2.u = u2_backup - h*I_vec(1+Body1.nx:Body1.nx + Body2.nx);
    %     Fch = ContactForce2(Body1,Body2,penalty,approach); % force due to variation
    % 
    %     Body1.u = u1_backup + h*I_vec(1:Body1.nx); 
    %     Body2.u = u2_backup + h*I_vec(1+Body1.nx:Body1.nx + Body2.nx);
    %     Fch2 = ContactForce2(Body1,Body2,penalty,approach); % force due to variation
    % 
    %     Kc(:,ii) = (Fch2 - Fch) / (2*h); 
    % 
    %     I_vec(ii)=0;
    % end    
    
    % returning all back to normal
    Body1.u = u1_backup;
    Body2.u = u2_backup;