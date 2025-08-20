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
    Outcome = FindTargetPoint(TargetBody,ContactPoint);

    % Searchin the attributes of the corresponding point on the target surface 
    % Outcome = FindTargetPoint(TargetBody,ContactPoint);

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
           Fcont_loc = penalty * Gap * Normal;                                                                              
           Ftarg_loc = penalty * Gap * Nm_2412(xi,eta)'*Normal; % calculation  & redistribution over the nodes of target element 
            
        % Nitsche approach   
        elseif approach == 2
                
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Node = ContactBody.contact.nodalid(ii);
           ElemenNumber = find(any(ContactBody.nloc == Node, 2), 1, 'first');
           [X2,U2] = GetCoorDisp(ElemenNumber,ContactBody.nloc,ContactBody.P0,ContactBody.u);
           [xi2,eta2] = FindIsoCoord2(X2,U2,ContactPoint'); 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % 
           % % Traget
           % Sigma = Sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);
           % nabla_sigma = nabla_sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);                
           % d_lambda = NablaMultiplication(nabla_sigma,Normal,Normal);
           % 
           % % contact
           % Sigma2 = Sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
           % nabla_sigma2 = nabla_sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
           % d_lambda2 = NablaMultiplication(nabla_sigma2,Normal,Normal); 
           % 
           % % assemblance
           % lambda_1 =  Normal' * (Sigma2 - Sigma) * Normal;
           % lambda = norm(lambda_1); 
           % 
           % d_lambda  = sign(lambda_1)*d_lambda; 
           % d_lambda2 = sign(lambda_1)*d_lambda2;
           % 
           % % calculation  & redistribution over the nodes of target elemen
           % Ftarg_loc = penalty * Gap * Nm_2412(xi,eta)'*Normal  + lambda*Nm_2412(xi,eta)'*Normal + Gap * Nm_2412(xi,eta)'*d_lambda;
           % Fcont_loc = penalty * Gap * Normal + lambda * Normal + Gap * d_lambda2;   % storing contact forces, we put minus here due to gap negativity 

           % vector Gap
           vec = Gap^2 * Normal; % Gap^2 = abs(rb-ra), 

           % Traget
           Sigma = Sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);
           nabla_sigma = nabla_sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);                
           d_lambda = NablaMultiplication(nabla_sigma,Normal) * vec;

           % contact
           Sigma2 = Sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
           nabla_sigma2 = nabla_sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
           d_lambda2 = NablaMultiplication(nabla_sigma2,Normal) * vec; 

           % assemblance
           lambda = vec' * (Sigma + Sigma2) * Normal; 

           % calculation  & redistribution over the nodes of target elemen
           Ftarg_loc = penalty * Gap * Nm_2412(xi,eta)'*Normal  + lambda*Nm_2412(xi,eta)'*Normal + Gap * Nm_2412(xi,eta)'*d_lambda;
           Fcont_loc = penalty * Gap * Normal + lambda * Normal + Gap * d_lambda2;   % storing contact forces, we put minus here due to gap negativity 
        end 

        Ftarg(TargetBody.xloc(Index,:)) = Ftarg(TargetBody.xloc(Index,:)) + Ftarg_loc; % storing contact forces to target
        Fcont(DOFpositions) = Fcont(DOFpositions) - Fcont_loc; % we put minus here due to gap negativity and normal is from the target body
    end   
end 

% Assemblace
Fc = [Fcont;Ftarg]; 

