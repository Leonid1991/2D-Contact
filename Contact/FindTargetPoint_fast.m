function Result = FindTargetPoint_fast(TargetBody,ContactPoint)


    DofsAtNode_targ = TargetBody.DofsAtNode;
    ContactNode_targ = TargetBody.contact.nodalid;
    nloc = TargetBody.nloc;
     
    % current position of the "possible contact" nodes of the Contact Body
    TargetPoints_X = TargetBody.q(xlocChosen(DofsAtNode_targ, ContactNode_targ,1)) + ...
                     TargetBody.u(xlocChosen(DofsAtNode_targ, ContactNode_targ,1));   ... % coords on X axis;
    
    TargetPoints_Y = TargetBody.q(xlocChosen(DofsAtNode_targ,ContactNode_targ,2)) + ... % coords on Y axis
                     TargetBody.u(xlocChosen(DofsAtNode_targ,ContactNode_targ,2));
    
    TargetPoints = [TargetPoints_X TargetPoints_Y]; % nodes of the contact surfaces of the contact body 
    
    distances = vecnorm(TargetPoints - ContactPoint, 2, 2); % distances between all target contact node and the point
    
    % sorting and choosing two closest
    [~, sortedIndices] = sort(distances); 
    a = ContactNode_targ(sortedIndices(1));
    b = ContactNode_targ(sortedIndices(2));

    position_a = TargetPoints(sortedIndices(1),:);
    position_b = TargetPoints(sortedIndices(2),:);

    [xy, distance, t_a]  = distance2curve([position_a; position_b],ContactPoint,'linear');
    tol = 1e-6; % Tolerance for error margin    
     
    info = []; % array, where we will store info of the contact point projection 
    if ~( (t_a < tol || abs(t_a - 1) < tol) && distance > tol )  % sanity check that the point isn't outside  (t_a ~= 0, 1)
                                                                 % and having distance > 0 at the same time                
        % To what element these nodes belong
        % idea that on the edge, two nodes are uniquely belong to one element only 
        ElemenNumber = find(any(nloc == a, 2) & any(nloc == b, 2)); 
        if ~isempty(ElemenNumber)
            [X,U] = GetCoorDisp(ElemenNumber,nloc,TargetBody.P0,TargetBody.u); % position of nodes of the element     
            % Position of the central point of the chosen element
            central = Nm_2412(0,0)*(X+U); 
        
            % Finding external normal to the element (a central point helps identify the outward direction)
            Normal = Normal3points(central,position_a',position_b'); % outwards normal of element (algorithm doesn't depend on the order a and b)            
            info = [info; xy distance Normal ElemenNumber];

        else
            Result.Index = 0;
        end             
    end 

    Result.Gap = 0; % we always have something to work with;
    

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

    