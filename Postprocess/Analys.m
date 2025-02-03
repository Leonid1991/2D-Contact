function Body = Analys(Body)

NumberOfIDs=length(Body.Fext.nodalid);
DofsAtNode = Body.DofsAtNode;
Results = [];




ux = Body.u(xlocChosen(DofsAtNode,Body.Fext.nodalid,1));
uy = Body.u(xlocChosen(DofsAtNode,Body.Fext.nodalid,2));
% averaged displacements.
uxavg=sum(ux)/(NumberOfIDs);
uyavg=sum(uy)/(NumberOfIDs);
Results = [Results; Body.nElems.x Body.nElems.y Body.ndof uxavg uyavg]; 

Body.results = Results;