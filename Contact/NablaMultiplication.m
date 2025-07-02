function nabla_sigma_n_tensor = NablaMultiplication(NablaSigma,vec1)

        % NablaVector = zeros(2,1);
        % 
        % i = 1;
        % NablaSigma_xi = [NablaSigma(1,i) NablaSigma(3,i);
        %                  NablaSigma(3,i) NablaSigma(2,i)];
        % 
        % i = 2;
        % NablaSigma_eta = [NablaSigma(1,i) NablaSigma(3,i);
        %                   NablaSigma(3,i) NablaSigma(2,i)];
        % 
        % 
        % NablaVector(1) = vec1' * NablaSigma_xi * vec2;
        % NablaVector(2) = vec1' * NablaSigma_eta * vec2;
        
        nabla_sigma_n = NablaSigma * vec1; % the first n (over the components of given by nabla) 
        nabla_sigma_n_tensor = [nabla_sigma_n(1) nabla_sigma_n(3);
                                nabla_sigma_n(3) nabla_sigma_n(2)];
        % NablaVector = nabla_sigma_n_tensor * vec2; % the second n 

