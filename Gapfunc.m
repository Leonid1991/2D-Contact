function Gap = Gapfunc(ContactBody,TargetBody)


Gap = 0;

% current position of the "possible contact" nodes of the Contact Body
ContactPoints_X = ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1)) + ...
                  ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,1));   ... % coords on X axis;

ContactPoints_Y =  ContactBody.q(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2)) + ... % coords on Y axis
                   ContactBody.u(xlocChosen(ContactBody.DofsAtNode,ContactBody.contact.nodalid,2));

ContactPoints = [ContactPoints_X ContactPoints_Y]; % nodes of the contact body of the contact surfaces
TagetPoints = [];
for ii = 1:size(ContactPoints,1) % loop over all contact points
  
    ContactPoint = ContactPoints(ii,:);

    % Searchin the attributes of the corresponding point on the target surface 
    Outcome = FindTargetPoint(TargetBody,ContactPoint);
    % TagetPoints = [TagetPoints; Outcome.Position'];

    % Outcome.Normal

    % Checking the condition of the penalty approach
    if Outcome.Gap < 0 % we have meaningful outcome from the search
        Gap= Gap + abs(Outcome.Gap);
    end   
end
% disp('here eeeeee');
% figure();
% plot(ContactPoints(:,1), ContactPoints(:,2),'.-k');
% hold on 
% plot(TagetPoints(:,1), TagetPoints(:,2),'.-b');

