function status = printStatus(deltaf, uu_bc, tol, ii, jj, imax, steps, titertot, Gap)
    
     persistent previousDisp previousDisp_2

     if jj == 1
        previousDisp = 0; 
        previousDisp_2 = 0;
     end 
        
     if all(abs(deltaf)<tol) || all(abs(uu_bc)<tol) || (abs(2*norm(uu_bc) - previousDisp - previousDisp_2)<tol) 
         if ~isnan(Gap)
             fprintf('Convergence: %10.4f, Displacements norm: %10.4f, Total gap: %10.7f\n', norm(abs(deltaf)), norm(uu_bc), Gap);            
         else
             fprintf('Convergence: %10.4f, Displacements norm: %10.4f\n', norm(abs(deltaf)), norm(u_bc));            
         end
         fprintf('Solution for %d / %d step  is found on %d iteration, Total CPU-time: %f\n', ii, steps, jj, titertot);
         status = true;
     elseif jj==imax 
         fprintf('The solution for %d step is not found. The maximum number of iterations is reached. Total CPU-time: %d\n', ii, jj);
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
        previousDisp_2 = previousDisp;
        previousDisp = norm(uu_bc);
        
     end