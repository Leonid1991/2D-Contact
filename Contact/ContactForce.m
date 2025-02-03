function Fc =  ContactForce(ContactBody,TargetBody,penalty,approach)

Fcont = zeros(ContactBody.nx,1);
Ftarg = zeros(TargetBody.nx,1);

% current position of the "possible contact" nodes of the Contact Body
ContactPoints = [ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1)) + ...
                 ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1))   ...
                 ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2)) + ...
                 ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2))];


for ii = 1:length(ContactBody.contact.nodalid) % loop over all contact nodes
  
    % Searchin the attributes of the corresponding point on the target surface 
    [GlobPosition,index,normal,boundaryCheck] = FindTargetPoint(TargetBody,ContactPoints(ii,:));

    % The vector from projection to the point
    vec = ContactPoints(ii,:)' - GlobPosition';
    gap = normal' * vec;  % gap calculation
    
    % Checking the condition of the penalty approach
    if (gap < 0) && boundaryCheck>0 && boundaryCheck<1   % additional (boundary) conditions are the safe mechanism
                                                         % because projection outside element gives also values 0 and 1         
       % DOFs of the contact elemnet under consideration
       DOFpositions = xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid(ii),1:2);
       % Now we need to find forces applied to the target body due to interation
       % index - element under consideration
       [X,U] = GetCoorDisp(index,TargetBody.nloc,TargetBody.P0,TargetBody.u);
       [xi,eta] = FindIsoCoord(X,U,GlobPosition); % finding isoparametric coodinates of the point
           
       % penalty approahc
       if approach == 1
          penalty_coef = penalty*gap;
          Contact_c = penalty_coef;
       % Nitsche approach   
       elseif approach == 2 
          penalty_coef = penalty*gap;
          Nitshe_coef = (nabla_u_2412(U,X,xi,eta)*normal)'*vec;
          Contact_c = penalty_coef + Nitshe_coef;
       end    
            
       % calculation of the forces applied to the nodes of contact elemnet 
       Fcont(DOFpositions) = Fcont(DOFpositions) - Contact_c*normal;   % storing contact forces
       % calculation  & redistribution over the nodes of target element
       Ftarg_loc = Contact_c*Nm_2412(xi,eta)'*normal;
       Ftarg(TargetBody.xloc(index,:)) = Ftarg(TargetBody.xloc(index,:)) + Ftarg_loc; % storing contact forces to target           
    end   
end 

% Assemblace
Fc = [Fcont;Ftarg]; 

