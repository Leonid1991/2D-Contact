function Fc =  ContactForce(ContactBody,TargetBody,penalty,approach)

Fcont = zeros(ContactBody.nx,1);
Ftarg = zeros(TargetBody.nx,1);

% current position of the "possible contact" nodes of the Contact Body
ContactPoints = [ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1)) + ...
                 ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1))   ... % coords on X axis
                 ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2)) + ... % coords on Y axis
                 ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2))];

for ii = 1:length(ContactBody.contact.nodalid) % loop over all contact nodes
  
    ContactPoint = ContactPoints(ii,:);

    % Searchin the attributes of the corresponding point on the target surface 
    Outcome = FindTargetPoint(TargetBody,ContactPoint);

    % Checking the condition of the penalty approach
    if Outcome.Gap < 0 % we have meaningful outcome from the search

       % DOFs of the contact elemnet under consideration
       DOFpositions = xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid(ii),1:2);

       % Now we need to find forces applied to the target body due to the interation
       % index - element under consideration
       [X,U] = GetCoorDisp(Outcome.Index,TargetBody.nloc,TargetBody.P0,TargetBody.u);
       [xi,eta] = FindIsoCoord2(X,U,Outcome.Position); % finding isoparametric coodinates of the point
           
       % penalty approach
       if approach == 1
          Contact_c = penalty*Outcome.Gap;

       % Nitsche approach   
       % elseif approach == 2 
       %    penalty_coef = penalty*Outcome.Gap;
       %    Nitshe_coef = (nabla_u_2412(U,X,xi,eta)*normal)'*vec;
       %    Contact_c = penalty_coef + Nitshe_coef;

       end    
            
       % calculation of the forces applied to the nodes of contact elemnet 
       Fcont(DOFpositions) = Fcont(DOFpositions) - Contact_c*Outcome.Normal;   % storing contact forces
                                                                       % we put minus here due to gap negativity  
                                                                       
       % calculation  & redistribution over the nodes of target element
       Ftarg_loc = Contact_c*Nm_2412(xi,eta)'*Outcome.Normal;
       Ftarg(TargetBody.xloc(Outcome.Index,:)) = Ftarg(TargetBody.xloc(Outcome.Index,:)) + Ftarg_loc; % storing contact forces to target        
       
    end   
end 

% Assemblace
Fc = [Fcont;Ftarg]; 

