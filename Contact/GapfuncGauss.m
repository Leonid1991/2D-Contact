function Gap = GapfuncGauss(ContactBody,TargetBody,n)

Gap = 0;

% current position of the "possible contact" nodes of the Contact Body
ContactPoints_X = ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1)) + ...
                  ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1));   ... % coords on X axis;

ContactPoints_Y =  ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2)) + ... % coords on Y axis
                   ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2));

ContactPoints = [ContactPoints_X ContactPoints_Y]; % nodes of the contact body of the contact surfaces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition to Gauss points
ContactPoints2 = []; 
for ii = 2:size(ContactPoints,1)
    % taking two consecutive nodes
    a = ContactPoints(ii-1,:);
    b = ContactPoints(ii,:);    
    % idea that on the edge, two nodes are uniquely belong to one element only 
    [xx,~] = gauleg2(a(1),b(1),n); % split in x- axis
    [yy,~] = gauleg2(a(2),b(2),n); % split in y- axis

    % Following the Intercept theorem (Thales's theorem according Russian naming):
    % the function devide x- and y- axes in the same propotions 
    ContactPoints2 = [ContactPoints2; xx yy];
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:size(ContactPoints2,1) % loop over all contact points
  
    ContactPoint = ContactPoints2(ii,:);  
    Outcome = FindTargetPoint_fast(TargetBody,ContactPoint);
    % Checking the condition of the penalty approach
    if Outcome.Gap < 0 % we have meaningful outcome from the search

       Gap= Gap + abs(Outcome.Gap);
        
    end 
end 


