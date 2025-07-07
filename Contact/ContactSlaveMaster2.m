function [Fcont,Ftarg] =  ContactSlaveMaster2(ContactBody,TargetBody,penalty,approach)
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
    Outcome = FindTargetPoint2(TargetBody,ContactPoint);
    
    
    % Checking the condition of the penalty approach
    if Outcome.Gap < 0 % we have meaningful outcome from the search
        % negative because after penetration direction from contact point
        % to its projection is opposite to outwards normal

        % DOFs of the contact element under consideration
        DOFpositions = xlocChosen(DofsAtNode_cont, ContactNode_cont(ii),1:2);

        % Now we need to find forces applied to the target body due to the interation
        % index - element under consideration
        [X,U] = GetCoorDisp(Outcome.Index,TargetBody.nloc,TargetBody.P0,TargetBody.u);

        [xi,eta] = FindIsoCoord2(X,U,Outcome.Position); % finding isoparametric coodinates of the point
                
        Normal = Outcome.Normal; % outwards normal (from the targeted body)
        Gap = abs(Outcome.Gap);
        Index = Outcome.Index;

        Normal_cont = -Normal;
        Normal_targ =  Normal;
        
        % penalty approach
        % calculation of the forces applied to the nodes of contact elemnet 
        Fcont_loc =  penalty * Gap * Normal_cont;                                                                              
        Ftarg_loc =  penalty * Gap * Normal_targ;
                   
        Ftarg(TargetBody.xloc(Index,:)) = Ftarg(TargetBody.xloc(Index,:)) + Nm_2412(xi,eta)' * Ftarg_loc; % redistribution over the nodes of target element 
        Fcont(DOFpositions) = Fcont(DOFpositions) + Fcont_loc;
    end    
end 
