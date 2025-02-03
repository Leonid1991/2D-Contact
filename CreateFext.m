
function Body = CreateFext(ii,steps,Body)
    
    DofsAtNode = Body.DofsAtNode;
    Fext = zeros(Body.nx,1); % Initialize vector of ext forces
    
    NodalIDAt = FindGlobNodalID(Body.P0,Body.Fext.loc,Body.shift);
    NumberOfIDs = length(NodalIDAt);

    Fext(xlocChosen(DofsAtNode,NodalIDAt,1))= Body.Fext.x*ii/(steps*NumberOfIDs);       
    Fext(xlocChosen(DofsAtNode,NodalIDAt,2))= Body.Fext.y*ii/(steps*NumberOfIDs);

    Body.Fext.vec = Fext;
    Body.Fext.nodalid = NodalIDAt;