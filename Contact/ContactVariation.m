function [Fc,Kc,GapNab,GapDOFs,Gap] = ContactVariation(Body1,Body2,penalty,approach,ContactPointfunc,Gapfunc)
    
    sqrtEps = sqrt(eps);
    h = 2*sqrtEps;
     
    % Current (!!total!!) gap calculation 
    Gap = Gapfunc(Body1,Body2);
    
    % Current gap calculation with distribution over all DOFs
    GapDOFs = GapPointDistribution(Body1,Body2,ContactPointfunc);

    % Backup original coordinates
    u1_backup = Body1.u;
    u2_backup = Body2.u;

    nx = Body1.nx + Body2.nx;
    
    % Initiation 
    Kc = zeros(nx,nx);    
    GapNab = zeros(nx,1);
    I_vec=zeros(nx,1);    
    
    if (approach ~= 5) && (approach ~= 8) % all, but Lagrange multiplier
        Fc= ContactForce(Body1,Body2,penalty,approach,ContactPointfunc); % Body1 forces from the projection of Body2        
    else % Lagrange multiplier   
        Fc = zeros(nx,1); % we aren't interested 
        GapNab = GapDOFs;
    end
        
    for ii = 1:nx
        % start variationing 
        I_vec(ii)=1;
        
        % this split is to distribute coord. between bodies 
        Body1.u = u1_backup - h*I_vec(1:Body1.nx); 
        Body2.u = u2_backup - h*I_vec(1+Body1.nx:end);

        if (approach == 7) || (approach == 8) % Augumented Lagrange & Detailed Lagrange multiplier
            
            GapDOFsh = GapPointDistribution(Body1,Body2,ContactPointfunc);
            Kc(:,ii) = (GapDOFs - GapDOFsh) / h; 
            
        elseif (approach ~= 5) &&  (approach ~= 6) % Penalty- & Nitsche-based approaches
            
            Fch = ContactForce(Body1,Body2,penalty,approach,ContactPointfunc); % force due to variation
            Kc(:,ii) = (Fc - Fch) / h; 
        
        else % Lagrange multiplier & very simplified Penalty

            Gaph = Gapfunc(Body1,Body2);
            GapNab(ii) = (Gap - Gaph)/h;

        end

        % end variationing 
        I_vec(ii)=0;      
    end    

    % returning all back to normal
    Body1.u = u1_backup;
    Body2.u = u2_backup;