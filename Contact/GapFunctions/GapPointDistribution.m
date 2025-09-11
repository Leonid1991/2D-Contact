function GapDOFs = GapPointDistribution(ContactBody,TargetBody,ContactPointfunc)
    
    % Target Body is a body points are projected
    % Contact Body is a body points are taken for projection
    GapDOFs_cont = zeros(ContactBody.nx,1);
    GapDOFs_targ = zeros(TargetBody.nx,1);

    [GapPoints, GapPointsElements] = ContactPointfunc(ContactBody);
    [count, ContactGeometry, TargetGeometry, Gaps, Normals] = Projection(GapPoints, GapPointsElements,ContactBody,TargetBody); 

    if count~=0 % we have contact
        
        for i = 1:count % loop all over contact points

            xi_cont = ContactGeometry.CoordsXi(:,i) ;    
            xi_targ = TargetGeometry.CoordsXi(:,i);
            
            Gap = Gaps(i);
            Normal = Normals(:,i);
            
            Normal_cont = -Normal;
            Normal_targ =  Normal; 

                        
            % redistribution over the nodes
            DOFpositions_cont = ContactGeometry.Dofs(:,i);
            DOFpositions_targ = TargetGeometry.Dofs(:,i);

            GapDOFs_cont(DOFpositions_cont) = GapDOFs_cont(DOFpositions_cont) + Gap * Nm_2412(xi_cont(1),xi_cont(2))' * Normal_cont;
            GapDOFs_targ(DOFpositions_targ) = GapDOFs_targ(DOFpositions_targ) + Gap * Nm_2412(xi_targ(1),xi_targ(2))' * Normal_targ;   
            
        end 
    end  
    
    % Assemblace
    GapDOFs = [GapDOFs_cont; GapDOFs_targ];
    
