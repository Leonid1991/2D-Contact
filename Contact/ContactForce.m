function Fc = ContactForce(ContactBody,TargetBody,penalty,approach,ContactPointfunc)
    
    % Target Body is a body points are projected
    % Contact Body is a body points are taken for projection
    Fcont = zeros(ContactBody.nx,1);
    Ftarg = zeros(TargetBody.nx,1);

    [ContactPoints, ContactPointsElements] = ContactPointfunc(ContactBody);
    [count, ContactGeometry, TargetGeometry, Gaps, Normals] = Projection(ContactPoints,ContactPointsElements,ContactBody,TargetBody); 

    if count~=0 % we have contact
        
        for i = 1:count % loop all over contact points

            xi_cont = ContactGeometry.CoordsXi(:,i) ;    
            xi_targ = TargetGeometry.CoordsXi(:,i);
            
            Gap = Gaps(i);
            Normal = Normals(:,i);
            
            Normal_cont = -Normal;
            Normal_targ =  Normal; 
            
            % penalty approach
            if (approach == 1) || (approach == 6) || (approach == 7) || (approach == 8)     
            
               % calculation of the forces applied to the nodes of contact elemnet 
               Fcont_loc = penalty * Gap * Normal_cont;                                                                              
               Ftarg_loc = penalty * Gap * Normal_targ; 

            % Nitsche approaches   
            elseif (approach > 1) && (approach < 5)
                
               % contact 
               X_cont = ContactGeometry.Coords(:,i);
               U_cont = ContactGeometry.Disp(:,i);                
               Sigma_cont = Sigma_2412(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont(1),xi_cont(2));
               
               % Traget
               X_targ = ContactGeometry.Coords(:,i);
               U_targ = ContactGeometry.Disp(:,i);
               Sigma_targ = Sigma_2412(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ(1),xi_targ(2));
                          
                         
               % Normal force difference 
               Sigma_n = Normal_cont' * Sigma_cont * Normal_cont - Normal_targ' * Sigma_targ * Normal_targ;                          
               lambda = Gap * norm(Sigma_n);
               
               d_lambda_targ = norm(Sigma_n)*Normal_targ;
               d_lambda_cont = norm(Sigma_n)*Normal_cont; 
                          
               if approach >2 % adding nonlinear of gap    
                    
                   nabla_sigma_targ = nabla_sigma_2412(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ(1),xi_targ(2));                
                   Nabla_Sigma_n_targ = NablaMultiplication(nabla_sigma_targ,Normal_targ);
    
                   nabla_sigma_cont = nabla_sigma_2412(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont(1),xi_cont(2));
                   Nabla_Sigma_n_cont = NablaMultiplication(nabla_sigma_cont,Normal_cont); 
    
                   Nabla_Sigma_n = Normal_cont' * Nabla_Sigma_n_cont * Normal_cont + Normal_targ' * Nabla_Sigma_n_targ * Normal_targ; 
                   d_lambda_targ = d_lambda_targ + Gap * Nabla_Sigma_n * Normal_targ;
                   d_lambda_cont = d_lambda_cont + Gap * Nabla_Sigma_n * Normal_cont; 
               end 
                      
               if approach >3 % adiing all terms
                  Nabla_n = (eye(2) - Normal_cont * Normal_cont')./Gap;  
                  Sigma = Sigma_targ + Sigma_cont;
                  d_lambda_targ = d_lambda_targ + Gap * ( Normal_targ' * Nabla_n'* Sigma * Normal_targ + Normal_targ' * Sigma * Nabla_n * Normal_targ + Normal_targ' * Sigma * Normal_targ * Nabla_n) * Normal_targ;
                  d_lambda_cont = d_lambda_cont + Gap * ( Normal_cont' * Nabla_n'* Sigma * Normal_cont + Normal_cont' * Sigma * Nabla_n * Normal_cont + Normal_cont' * Sigma * Normal_cont * Nabla_n) * Normal_cont;
               end
    
               % calculation  & redistribution over the nodes of target element
               Ftarg_loc = penalty * Gap * Normal_targ + lambda * Normal_targ + Gap * d_lambda_targ;
               Fcont_loc = penalty * Gap * Normal_cont + lambda * Normal_cont + Gap * d_lambda_cont;  
            end   

                        
            % redistribution over the nodes
            DOFpositions_cont = ContactGeometry.Dofs(:,i);
            DOFpositions_targ = TargetGeometry.Dofs(:,i);

            Fcont(DOFpositions_cont) = Fcont(DOFpositions_cont) + Nm_2412(xi_cont(1),xi_cont(2))' * Fcont_loc;
            Ftarg(DOFpositions_targ) = Ftarg(DOFpositions_targ) + Nm_2412(xi_targ(1),xi_targ(2))' * Ftarg_loc;   

        end 
    end  
    
    % Assemblace
    Fc = [Fcont;Ftarg]; 

    
