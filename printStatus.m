function status = printStatus(deltaf, uu_bc, tol, ii, jj, imax, steps, titertot, Gap)
    
     persistent previousForce previousDisp previousDisp_2

     if jj == 1
        previousForce = 0; 
        previousDisp = 0; 
        previousDisp_2 = 0;
     end 
        
     if all(abs(deltaf)<tol) || ...
         all(abs(uu_bc)<tol^2) ||  abs(norm(uu_bc) - previousDisp)<tol^2 || ... % stop when the change is too small    
         (abs(norm(uu_bc) - previousDisp)<tol && abs(norm(uu_bc) - previousDisp_2)<tol) || ... % stop when small changes are repeating
         ( all(abs(deltaf) - previousForce)<tol && abs(norm(uu_bc) - previousDisp)<tol ) % additional condition for exit  
         if ~isnan(Gap)
             fprintf('Convergence: %10.4f, Displacements norm: %10.4f, Total gap: %10.7f\n', norm(abs(deltaf)), norm(uu_bc), Gap);            
         else
             fprintf('Convergence: %10.4f, Displacements norm: %10.4f\n', norm(abs(deltaf)), norm(u_bc));            
         end
         fprintf('Solution for %d / %d step  is found on %d iteration, Total CPU-time: %f\n', ii, steps, jj, titertot);
         status = true;
     elseif jj==imax 
         fprintf('The solution for %d step is not found. The maximum number of iterations is reached. Total CPU-time: %.3f\n', ii, titertot);
         status = false;   
     else     
         if ~isnan(Gap)
            fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f, Total gap: %10.7f\n', jj, norm(abs(deltaf)), norm(uu_bc), Gap);
         else
            fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f\n', jj, norm(abs(deltaf)), norm(uu_bc));
         end
         status = false;
     end 

  
     % Update previous steps' meanings
     if (jj > 1) && (jj < imax)
        previousForce = deltaf;
        previousDisp_2 = previousDisp;
        previousDisp = norm(uu_bc);
        
     end