
function Body = CreateFext(ii,steps,Body)
    
    DofsAtNode = Body.DofsAtNode;
    Fext = zeros(Body.nx,1); % Initialize vector of ext forces
    
    NodalIDAt = FindGlobNodalID(Body.P0,Body.Fext.loc,Body.shift);
    NumberOfIDs = length(NodalIDAt);


    Load_x = (Body.Fext.x*ii/steps)/(NumberOfIDs-1);
    Load_y = (Body.Fext.y*ii/steps)/(NumberOfIDs-1);
    
    NodalIDAtEdge1 = FindGlobNodalID(Body.P0,Body.edge1.loc,Body.shift);
    NodalIDAtEdge2 = FindGlobNodalID(Body.P0,Body.edge2.loc,Body.shift);

    nodeEdge1 = intersect(NodalIDAt, NodalIDAtEdge1);
    nodeEdge2 = intersect(NodalIDAt, NodalIDAtEdge2);

    Fext(xlocChosen(DofsAtNode,NodalIDAt,1))= Load_x;       
    Fext(xlocChosen(DofsAtNode,NodalIDAt,2))= Load_y;
        
    Fext(xlocChosen(DofsAtNode,nodeEdge1,1))= Load_x/2;
    Fext(xlocChosen(DofsAtNode,nodeEdge1,2))= Load_y/2;

    Fext(xlocChosen(DofsAtNode,nodeEdge2,1))= Load_x/2;
    Fext(xlocChosen(DofsAtNode,nodeEdge2,2))= Load_y/2;

    
    Body.Fext.vec = Fext;
    Body.Fext.nodalid = NodalIDAt;