function Result = FindTargetPoint2(TargetBody,ContactPoint)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % we always have something to work with (just for saving)
    Result.Gap = 0; 
    Result.Index = 1;
    Result.Normal = [0;1]; % normal on the surf.
    Result.Position = ContactPoint'; % point projection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nloc = TargetBody.nloc; 
    DofsAtNode = TargetBody.DofsAtNode;
    % current position of the contact line
    ContactPoints_X = TargetBody.q(xlocChosen(DofsAtNode,TargetBody.contact.nodalid,1)) + ...
                      TargetBody.u(xlocChosen(DofsAtNode,TargetBody.contact.nodalid,1));   ... % coords on X axis;
    
    ContactPoints_Y =  TargetBody.q(xlocChosen(DofsAtNode,TargetBody.contact.nodalid,2)) + ... % coords on Y axis
                       TargetBody.u(xlocChosen(DofsAtNode,TargetBody.contact.nodalid,2));
    
    ContactLine= [ContactPoints_X ContactPoints_Y]; % nodes of the contact body of the contact surfaces
        
    [xy, distance, t_a]   = distance2curve(ContactLine,ContactPoint,'linear');
                      
    % take two consecutive (the closest) nodes of the target body
    dists = sqrt(sum((ContactLine - xy).^2, 2)); % differences
    [~, sortedIdx] = sort(dists);
    twoClosestIdx = sortedIdx(1:2);                     % indices of closest points  
    a = TargetBody.contact.nodalid(twoClosestIdx(1));   % indice a
    b = TargetBody.contact.nodalid(twoClosestIdx(2));   % indices b
    position_a =  ContactLine(twoClosestIdx(1), :)';    % coordinates of a
    position_b =  ContactLine(twoClosestIdx(2), :)';    % coordinates of b
    
    info = [];
    tol = 1e-6; % Tolerance for error margin
    if ~( (t_a < tol || abs(t_a - 1) < tol) && distance > tol )  % sanity check that the point isn't outside.  
                                                                 % the point doesn't have t_a = 0, 1 and having distance > 0 at the same time
         
        ElemenNumber = find(any(nloc == a, 2) & any(nloc == b, 2));        % idea that on the edge, two nodes are uniquely belong to one element only
        [X,U] = GetCoorDisp(ElemenNumber,nloc,TargetBody.P0,TargetBody.u); % position of the element nodes
               
        central = Nm_2412(0,0)*(X+U); % Position of the central point of the chosen element
        
        % Finding external normal to the element (a central point helps identify the outward direction)
        normal = Normal3points(central,position_a,position_b); % algorithm doesn't depend on the order a and b  
        info = [xy distance normal ElemenNumber];
    end

    if isempty(info) == false % we actually have connection

        minDistance = min(info(:,3)); % minimal distance finding (the nearest point to the line)        
        index = min(find( minDistance == info(:,3))); % "min" is to remove possible ambiguity, such as the node can be dropped between elements   
    
        % info extraction
        TargPosition = info(index,1:2);
        Normal = info(index,4:5)';

        % Saving
        Result.Index = info(index,6); % element on which point projected       
        Result.Gap = (ContactPoint - TargPosition) * Normal;  % gap calculation
        Result.Normal = Normal; % normal on the surf.
        Result.Position = TargPosition'; % point projection
   end

    