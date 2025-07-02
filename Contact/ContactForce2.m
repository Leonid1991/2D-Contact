function Fc =  ContactForce2(ContactBody,TargetBody,penalty,approach)
% Target Body is a body points are projected
% Contact Body is a body points are taken for projection

Fcont = zeros(ContactBody.nx,1);
Ftarg = zeros(TargetBody.nx,1);


% current position of the "possible contact" nodes of the Contact Body
ContactPoints_X = ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1)) + ...
                  ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1));   ... % coords on X axis;

ContactPoints_Y =  ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2)) + ... % coords on Y axis
                   ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2));

ContactPoints = [ContactPoints_X ContactPoints_Y]; % nodes of the contact body of the contact surfaces


for ii = 1:size(ContactPoints,1) % loop over all contact points
  
    ContactPoint = ContactPoints(ii,:);

    % Searchin the attributes of the corresponding point on the target surface 
    Outcome = FindTargetPoint(TargetBody,ContactPoint);

    % Checking the condition of the penalty approach
    if Outcome.Gap < 0 % we have meaningful outcome from the search
        % negative because after penetration diraction from contact point
        % to its projection is opposite to outwards normal
        % because of it, all signs are opposite to the paper

        % DOFs of the contact element under consideration
        DOFpositions = xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid(ii),1:2);

        % Now we need to find forces applied to the target body due to the interation
        % index - element under consideration
        [X,U] = GetCoorDisp(Outcome.Index,TargetBody.nloc,TargetBody.P0,TargetBody.u);
        [xi,eta] = FindIsoCoord2(X,U,Outcome.Position); % finding isoparametric coodinates of the point
        
        Normal = Outcome.Normal; % outwards normal (from the targeted body)
        Gap = Outcome.Gap;
        Index = Outcome.Index;

        % penalty approach
        if approach == 1
          
           % calculation of the forces applied to the nodes of contact elemnet 
           Fcont(DOFpositions) = Fcont(DOFpositions) - penalty * Gap * Normal;% we put minus here due to gap negativity and normal is from the target body                                                                         
           % calculation  & redistribution over the nodes of target element 
           Ftarg_loc = penalty * Gap * Nm_2412(xi,eta)'*Normal;
           Ftarg(TargetBody.xloc(Index,:)) = Ftarg(TargetBody.xloc(Index,:)) + Ftarg_loc; % storing contact forces to target 

        % Nitsche approach   
        elseif approach == 2
                      
           % Traget
           Sigma = Sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);
           nabla_sigma = nabla_sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);       
           nabla_sigma_n = nabla_sigma * Normal; % the first n 
           nabla_sigma_n_tensor = [nabla_sigma_n(1) nabla_sigma_n(3);
                                   nabla_sigma_n(3) nabla_sigma_n(2)];
           d_lambda = nabla_sigma_n_tensor * Normal; % the second n 
           

           % contact

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % just to find element DOFs
           eta2 = -1;
           if ii == 1
               xi2 = -1;
               a = ContactBody.contact.nodalid(ii);
               b = ContactBody.contact.nodalid(ii+1);
           else
               xi2 = 1;
               a = ContactBody.contact.nodalid(ii-1);
               b = ContactBody.contact.nodalid(ii);
           end 
           ElemenNumber = find(any(ContactBody.nloc == a, 2) & any(ContactBody.nloc == b, 2)); 
           [X2,U2] = GetCoorDisp(ElemenNumber,ContactBody.nloc,ContactBody.P0,ContactBody.u); % position of nodes of the element 
           
            
           % just in case, becaus e we have a point at the edge  
           if norm(ContactPoint' - Nm_2412(xi2,eta2)*(X2 + U2))> 1e-6
               error('Point')
           end    
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

           Sigma2 = Sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
           nabla_sigma2 = nabla_sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
           nabla_sigma_n2 = nabla_sigma2 * Normal; % the first n 
           nabla_sigma_n_tensor2 = [nabla_sigma_n2(1) nabla_sigma_n2(3);
                                    nabla_sigma_n2(3) nabla_sigma_n2(2)];
           d_lambda2 = nabla_sigma_n_tensor2 * Normal; % the second n 


           lambda_1 =  ((Sigma2 + Sigma) *Normal)' * Normal;
           lambda = abs(lambda_1); 


           d_lambda2 = sign(lambda_1)* d_lambda2;
           d_lambda  = sign(lambda_1)* d_lambda; 

           % calculation  & redistribution over the nodes of target elemen
           Ftarg_loc = penalty * Gap * Nm_2412(xi,eta)'*Normal  + lambda*Nm_2412(xi,eta)'*Normal + Gap * Nm_2412(xi,eta)'*d_lambda;
           Ftarg(TargetBody.xloc(Index,:)) = Ftarg(TargetBody.xloc(Index,:)) + Ftarg_loc; % storing contact forces to target 
 
           
           Fcont(DOFpositions) = Fcont(DOFpositions) - penalty * Gap * Normal - lambda * Normal - Gap * d_lambda2;   % storing contact forces, we put minus here due to gap negativity                                                                                    

        end 

        % Gap= Gap + abs(Outcome.Gap);
    end   
end 

% Assemblace
Fc = [Fcont;Ftarg]; 

