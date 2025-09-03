function Fc =  ContactForce(ContactBody,TargetBody,penalty,approach)

% Target Body is a body points are projected
% Contact Body is a body points are taken for projection

Fcont = zeros(ContactBody.nx,1);
Ftarg = zeros(TargetBody.nx,1);

DofsAtNode_cont = ContactBody.DofsAtNode;
ContactNode_cont = ContactBody.contact.nodalid;

% current position of the "possible contact" nodes of the Contact Body
ContactPoints_X = ContactBody.q(xlocChosen(DofsAtNode_cont, ContactNode_cont,1)) + ...
                  ContactBody.u(xlocChosen(DofsAtNode_cont, ContactNode_cont,1));   ... % coords on X axis;

ContactPoints_Y =  ContactBody.q(xlocChosen(DofsAtNode_cont,ContactNode_cont,2)) + ... % coords on Y axis
                   ContactBody.u(xlocChosen(DofsAtNode_cont,ContactNode_cont,2));

ContactPoints = [ContactPoints_X ContactPoints_Y]; % nodes of the contact surfaces of the contact body 

for ii = 1:size(ContactPoints,1) % loop over all contact points
  
    ContactPoint = ContactPoints(ii,:);  
    Outcome = FindTargetPoint_fast(TargetBody,ContactPoint);
    
    % Checking the condition of the penalty approach
    if Outcome.Gap < 0 % we have meaningful outcome from the search

        % negative because after penetration direction from contact point
        % to its projection is opposite to outwards normal

        % DOFs of the contact element under consideration
        DOFpositions = xlocChosen(DofsAtNode_cont, ContactNode_cont(ii),1:2);

        % Now we need to find forces applied to the target body due to the interation
        % index - element under consideration
        [X,U] = GetCoorDisp(Outcome.Index,TargetBody.nloc,TargetBody.P0,TargetBody.u);

        [xi,eta] = FindIsoCoord(X,U,Outcome.Position); % finding isoparametric coodinates of the point
                
        Normal = Outcome.Normal; % outwards normal (from the targeted body to contact one)
        Gap = abs(Outcome.Gap);
        Index = Outcome.Index;

        Normal_cont = -Normal;
        Normal_targ =  Normal;
        
        % penalty approach
        if approach == 1
          
           % calculation of the forces applied to the nodes of contact elemnet 
           Fcont_loc = penalty * Gap * Normal_cont;                                                                              
           Ftarg_loc = penalty * Gap * Normal_targ;
            
        % Nitsche approach   
        elseif approach > 1
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % finding the isocoordinates and element numbers of the corrseponded
           % contact element
           Node = ContactBody.contact.nodalid(ii);
           ElemenNumber = find(any(ContactBody.nloc == Node, 2), 1, 'first');
           [X2,U2] = GetCoorDisp(ElemenNumber,ContactBody.nloc,ContactBody.P0,ContactBody.u);
           [xi2,eta2] = FindIsoCoord(X2,U2,ContactPoint'); 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           % Traget
           Sigma_targ = Sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);
                      
           % contact
           Sigma_cont = Sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
                      
           % Normal force difference 
           Sigma_n = Normal_cont' * Sigma_cont * Normal_cont - Normal_targ' * Sigma_targ * Normal_targ;           
           lambda = Gap * norm(Sigma_n);
           
           d_lambda_targ = norm(Sigma_n)*Normal_targ;
           d_lambda_cont = norm(Sigma_n)*Normal_cont; 
                      
           if approach >2 % adding nonlinear of gap    
                
               nabla_sigma_targ = nabla_sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);                
               Nabla_Sigma_n_targ = NablaMultiplication(nabla_sigma_targ,Normal_targ);

               nabla_sigma_cont = nabla_sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
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


        Ftarg(TargetBody.xloc(Index,:)) = Ftarg(TargetBody.xloc(Index,:)) + Nm_2412(xi,eta)' * Ftarg_loc; % redistribution over the nodes of target element 
        Fcont(DOFpositions) = Fcont(DOFpositions) + Fcont_loc;
    end       
end 

% Assemblace
Fc = [Fcont;Ftarg]; 
