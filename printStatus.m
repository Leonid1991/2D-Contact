function status = printStatus(deltaf, uu_bc, tol, ii, jj, imax, steps, titertot, Gap)
    
     if  nargin < 9 % Gap is optional
         Gap= NaN;
     end   
        
     if all(abs(deltaf)<tol) || all(abs(uu_bc)<tol^4)
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