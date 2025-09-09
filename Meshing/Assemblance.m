function [ff_bc, K_bc, deltaf] = Assemblance(Body1, Body2, Fc,Kc,GapNab,approach,penalty)

    
    bc = [Body1.bc Body2.bc]; % total logical vector of constrains
    Fext = [Body1.Fext.vec; Body2.Fext.vec]; % Assemblance of external forces

    % Assemblance of elastic stiffness matrices
    Ke = [            Body1.Fint.K zeros(Body1.nx,Body2.nx);
          zeros(Body2.nx,Body1.nx)            Body2.Fint.K];

    Fe = [Body1.Fint.vec; Body2.Fint.vec]; % Assemblance of elastic forces
    
    ff =  Fe - Fext + Fc; % Residuals

    % Remuval of the fixed dofs
    Ke_bc = Ke(bc,bc);  
    Kc_bc = Kc(bc,bc);
    K_bc = Ke_bc + Kc_bc; % stiffness matrices assemblance
    ff_bc = ff(bc);
    deltaf = ff_bc/norm(Fext(bc)); 
    GapNab_bc = GapNab(bc);
    
    if approach == 5
        K_bc = [    Ke_bc GapNab_bc;
                GapNab_bc'       0];   
        ff_bc = [ff_bc; 0];
    elseif approach == 6

        Kc_bc = penalty * (GapNab_bc * GapNab_bc');

    end

    
    

    
