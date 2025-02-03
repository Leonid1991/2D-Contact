function [GlobPosition,index,normal,boundaryCheck] = FindTargetPoint(TargetBody,ContactPoint)
nloc = TargetBody.nloc;
ArrayLength = length(TargetBody.contact.nodalid)-1; % number of all target elements
distanceinfo=zeros(ArrayLength,4);
    for jj = 1:ArrayLength
        
        % take two consecutive nodes of the target body  
        a = TargetBody.contact.nodalid(jj); 
        b = TargetBody.contact.nodalid(jj+1);

        % what element they (these nodes) belong to
        % idea that on the edge, two nodes are uniquely belong to one element only 
        ElemenNumber = find(any(nloc == a, 2) & any(nloc == b, 2)); 
        [X,U] = GetCoorDisp(ElemenNumber,nloc,TargetBody.P0,TargetBody.u); % position of nodes of the element 
        
        % Nodes position in global coordinate system
        position_a = TargetBody.q(xlocChosen(TargetBody.DofsAtNode,a,1:2)) + TargetBody.u(xlocChosen(TargetBody.DofsAtNode,a,1:2));
        position_b = TargetBody.q(xlocChosen(TargetBody.DofsAtNode,b,1:2)) + TargetBody.u(xlocChosen(TargetBody.DofsAtNode,b,1:2));  
        % Position of the central point of the chosen element
        central = Nm_2412(0,0)*(X+U); 

        % Finding external normal to the element (central element helps identify the outward direction)
        normal = Normal3points(central,position_a,position_b); % algorithm doesn't depend on the order a and b  

        % Projection of the current contact node on the line: target nodes (ii) -- (ii+1)
        % The line is border of the element (two points is enough, because we use elements with a linear interpolation)
        [distanceinfo(jj,1:2), distanceinfo(jj,3), distanceinfo(jj,4)]   = distance2curve([position_a'; position_b'],ContactPoint,'linear');
    end

    minDistance = min(distanceinfo(:,3)); % minimL distance finding 
    index = min(find( minDistance == distanceinfo(:,3))); % "min" is to remove possible ambiguity because it can be dropped on the node between elements
    
    GlobPosition = distanceinfo(index,1:2);

    boundaryCheck=distanceinfo(index,4);