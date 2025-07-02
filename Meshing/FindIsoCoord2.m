function [xi,eta] = FindIsoCoord2(X,U,point)

    r = Nm_2412(0,0)*(X + U);        
    Xi = [0;0];
    
    
    while (norm(r - point) > 1e-6) 

          Jac = [Nm_2412_xi(Xi(1),Xi(2))*(X+U) Nm_2412_eta(Xi(1),Xi(2))*(X+U)];   

          Xi = Xi - Jac^-1 * (r - point);  
        
          % removing numerical overflow near the edge
          Xi(1)  = max(min(Xi(1), 1), -1);
          Xi(2)  = max(min(Xi(2), 1), -1); 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          r = Nm_2412(Xi(1),Xi(2))*(X + U);

    end    
    
    xi = Xi(1);
    eta = Xi(2); 