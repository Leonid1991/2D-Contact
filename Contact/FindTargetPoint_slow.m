function Result = FindTargetPoint(TargetBody,ContactPoint)


DofsAtNode_targ = TargetBody.DofsAtNode;
ContactNode_targ = TargetBody.contact.nodalid;
nloc = TargetBody.nloc;

Result.Gap = 0; % we always have something to work with;
ArrayLength = length(ContactNode_targ)-1; % element number 
info=[]; % array, where we will store info of the contact point projection 

%% TODO: instead of all points, takes three closest 

    for jj = 1:ArrayLength % dealing with only element's line, because the whole contact line is a piece-wise
                           % its simultanious approximation migh result in error
        
        % take two consecutive nodes of the target body  
        a = ContactNode_targ(jj); 
        b = ContactNode_targ(jj+1);

        % Nodes position in global coordinate system
        position_a = TargetBody.q(xlocChosen(DofsAtNode_targ,a,1:2)) + TargetBody.u(xlocChosen(DofsAtNode_targ,a,1:2));
        position_b = TargetBody.q(xlocChosen(DofsAtNode_targ,b,1:2)) + TargetBody.u(xlocChosen(DofsAtNode_targ,b,1:2));  
        

        % Projection of the current contact node on the line: target nodes (ii) -- (ii+1)
        % The line is border of the element (two points is enough, because we use elements with a linear interpolation)

        [xy, distance, t_a]  = distance2curve([position_a'; position_b'],ContactPoint,'linear');
        tol = 1e-6; % Tolerance for error margin    

        if ~( (t_a < tol || abs(t_a - 1) < tol) && distance > tol )  % sanity check that the point isn't outside  (t_a ~= 0, 1)
                                                                     % and having distance > 0 at the same time                
            % To what element these nodes belong
            % idea that on the edge, two nodes are uniquely belong to one element only 
            ElemenNumber = find(any(nloc == a, 2) & any(nloc == b, 2)); 
            [X,U] = GetCoorDisp(ElemenNumber,nloc,TargetBody.P0,TargetBody.u); % position of nodes of the element 
    
            % Position of the central point of the chosen element
            central = Nm_2412(0,0)*(X+U); 
    
            % Finding external normal to the element (a central point helps identify the outward direction)
            Normal = Normal3points(central,position_a,position_b); % outwards normal of element (algorithm doesn't depend on the order a and b)            
            info = [info; xy distance Normal ElemenNumber];       
        end           

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

    