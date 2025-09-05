function Fc =  ContactForceGauss(ContactBody,TargetBody,penalty,approach,n)
% It is the same as 'ContactForce.m', but instead of nodes we want to
% consider Gauss points, according to Thesis Gary R. (1993) - see Mendeley
% Section 3.3 (Mendeley page - 41 p.)

% Target Body is a body points are projected
% Contact Body is a body points are taken for projection

Fcont = zeros(ContactBody.nx,1);
Ftarg = zeros(TargetBody.nx,1);

DofsAtNode_cont = ContactBody.DofsAtNode;
ContactNode_cont = ContactBody.contact.nodalid;

% current position of the "possible contact" nodes of the Contact Body
ContactPoints_X = ContactBody.q(xlocChosen(DofsAtNode_cont, ContactNode_cont,1)) + ...
                  ContactBody.u(xlocChosen(DofsAtNode_cont, ContactNode_cont,1));   ... % coords on X axis;

ContactPoints_Y =  ContactBody.q(xlocChosen(DofsAtNode_cont,ContactNode_cont,2)) + ... % coords on Y axis
                   ContactBody.u(xlocChosen(DofsAtNode_cont,ContactNode_cont,2));

ContactPoints = [ContactPoints_X ContactPoints_Y]; % nodes of the contact surfaces of the contact body 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition to Gauss points
nloc_cont = ContactBody.nloc;
ContactPoints2 = []; 
ContactPointsElements = [];
for ii = 2:size(ContactPoints,1)
    % taking two consecutive nodes
    a = ContactPoints(ii-1,:);
    b = ContactPoints(ii,:);    
    % idea that on the edge, two nodes are uniquely belong to one element only 
    ElemenNumber = find(any(nloc_cont == ii-1, 2) & any(nloc_cont == ii, 2));     
    % [~,w] = gauleg2(1,-1,n); % 
    [xx,~] = gauleg2(a(1),b(1),n); % split in x- axis
    [yy,~] = gauleg2(a(2),b(2),n); % split in y- axis

    % just a test, that it correlates with "ContactForceNode.m" 
    % xx = ContactPoints(ii,1); % split in x- axis
    % yy = ContactPoints(ii,2); % split in x- axis
        

    % Following the Intercept theorem (Thales's theorem according Russian naming):
    % the function devide x- and y- axes in the same propotions 
    ContactPoints2 = [ContactPoints2; xx yy];
    ContactPointsElements = [ContactPointsElements; ElemenNumber*ones(n,1)]; % multiplication to correlate with points
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:size(ContactPoints2,1) % loop over all contact points
  
    ContactPoint = ContactPoints2(ii,:);  
    Outcome = FindTargetPoint_fast(TargetBody,ContactPoint);
    nloc_targ = TargetBody.nloc;
    % Checking the condition of the penalty approach
    if Outcome.Gap < 0 % we have meaningful outcome from the search

        % negative because after penetration direction from contact point
        % to its projection is opposite to outwards normal

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ElemenNumber = ContactPointsElements(ii);                
        DOFpositions_cont = ContactBody.xloc(ElemenNumber,:); % DOFs of the contact element under consideration
        [X2,U2] = GetCoorDisp(ElemenNumber,nloc_cont,ContactBody.P0,ContactBody.u);
        [xi2,eta2] = FindIsoCoord(X2,U2,ContactPoint');  % finding the isocoordinates of the contact point

        % index - element under consideration
        Index = Outcome.Index;
        DOFpositions_targ = TargetBody.xloc(Index,:);
        [X,U] = GetCoorDisp(Index,nloc_targ,TargetBody.P0,TargetBody.u);
        [xi,eta] = FindIsoCoord(X,U,Outcome.Position); % finding isoparametric coodinates of the target point    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Normal = Outcome.Normal; % outwards normal (from the targeted body to contact one)
        Gap = abs(Outcome.Gap);
        
        Normal_cont = -Normal;
        Normal_targ =  Normal;
        
        % penalty approach
        if approach == 1
          
           % calculation of the forces applied to the nodes of contact elemnet 
           Fcont_loc = penalty * Gap * Normal_cont;                                                                              
           Ftarg_loc = penalty * Gap * Normal_targ;
            
        % Nitsche approach   
        elseif (approach > 1) && (approach < 5)
                     
           % Traget
           Sigma_targ = Sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);
                      
           % contact
           Sigma_cont = Sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
                      
           % Normal force difference 
           Sigma_n = Normal_cont' * Sigma_cont * Normal_cont - Normal_targ' * Sigma_targ * Normal_targ;           
           lambda = Gap * norm(Sigma_n);
           
           d_lambda_targ = norm(Sigma_n)*Normal_targ;
           d_lambda_cont = norm(Sigma_n)*Normal_cont; 
                      
           if approach >2 % adding nonlinear of gap    
                
               nabla_sigma_targ = nabla_sigma_2412(TargetBody.E,TargetBody.nu,U,X,xi,eta);                
               Nabla_Sigma_n_targ = NablaMultiplication(nabla_sigma_targ,Normal_targ);

               nabla_sigma_cont = nabla_sigma_2412(ContactBody.E,ContactBody.nu,U2,X2,xi2,eta2);
               Nabla_Sigma_n_cont = NablaMultiplication(nabla_sigma_cont,Normal_cont); 

               Nabla_Sigma_n = Normal_cont' * Nabla_Sigma_n_cont * Normal_cont + Normal_targ' * Nabla_Sigma_n_targ * Normal_targ; 
               d_lambda_targ = d_lambda_targ + Gap * Nabla_Sigma_n * Normal_targ;
               d_lambda_cont = d_lambda_cont + Gap * Nabla_Sigma_n * Normal_cont; 
           end 
                  
           if approach >3 % adiing all terms
              Nabla_n = (eye(2) - Normal_cont * Normal_cont')./Gap;  
              Sigma = Sigma_targ + Sigma_cont;
              d_lambda_targ = d_lambda_targ + Gap * ( Normal_targ' * Nabla_n'* Sigma * Normal_targ + Normal_targ' * Sigma * Nabla_n * Normal_targ + Normal_targ' * Sigma * Normal_targ * Nabla_n) * Normal_targ;
              d_lambda_cont = d_lambda_cont + Gap * ( Normal_cont' * Nabla_n'* Sigma * Normal_cont + Normal_cont' * Sigma * Nabla_n * Normal_cont + Normal_cont' * Sigma * Normal_cont * Nabla_n) * Normal_cont;
           end

           % calculation  & redistribution over the nodes of target element
           Ftarg_loc = penalty * Gap * Normal_targ + lambda * Normal_targ + Gap * d_lambda_targ;
           Fcont_loc = penalty * Gap * Normal_cont + lambda * Normal_cont + Gap * d_lambda_cont;  
        end   

        % redistribution over the nodes 
        Ftarg(DOFpositions_targ) = Ftarg(DOFpositions_targ) + Nm_2412(xi,eta)'   * Ftarg_loc; 
        Fcont(DOFpositions_cont) = Fcont(DOFpositions_cont) + Nm_2412(xi2,eta2)' * Fcont_loc;
    end       
end 

% Assemblace
Fc = [Fcont;Ftarg]; 
