function [ContactPoints, ContactPointsElements] = ContactPointGauss(ContactBody,n)     
    
    % we want to consider Gauss points, according to Thesis Gary R. (1993) - see Mendeley
    % Section 3.3 (Mendeley page - 41 p.)

    nloc_cont = ContactBody.nloc;
    DofsAtNode_cont = ContactBody.DofsAtNode;
    ContactNode_cont = ContactBody.contact.nodalid;
    
    % current position of the "possible contact" nodes of the Contact Body
    ContactNodesPosition_X = ContactBody.q(xlocChosen(DofsAtNode_cont, ContactNode_cont,1)) + ...
                             ContactBody.u(xlocChosen(DofsAtNode_cont, ContactNode_cont,1));   ... % coords on X axis;
    
    ContactNodesPosition_Y =  ContactBody.q(xlocChosen(DofsAtNode_cont,ContactNode_cont,2)) + ... % coords on Y axis
                              ContactBody.u(xlocChosen(DofsAtNode_cont,ContactNode_cont,2));
    
    ContactNodesPosition = [ContactNodesPosition_X ContactNodesPosition_Y]; % nodes of the contact surfaces of the contact body     
    
    ContactPoints = []; % taking the first node position, otherwise later it will be omitted 
    ContactPointsElements = []; % taking the first node number
    Contactweight = [];
    % transition to Linspace points  
    for ii = 2:size(ContactNodesPosition,1)
        % taking two consecutive nodes
        a = ContactNodesPosition(ii-1,:);
        b = ContactNodesPosition(ii,:);    
        % idea that on the edge, two nodes are uniquely belong to one element only 
        ElemenNumber = find(any(nloc_cont == ii-1, 2) & any(nloc_cont == ii, 2));     
        [xx,ww] = gauleg2(a(1),b(1),n); % split in x- axis
        [yy,~] = gauleg2(a(2),b(2),n); % split in y- axis
       
        % just a test, that it correlates with Node positioning 
        % xx = ContactPoints(ii,1); % split in x- axis
        % yy = ContactPoints(ii,2); % split in x- axis
        
        % Following the Intercept theorem (Thales's theorem according Russian naming):
        % the function devide x- and y- axes in the same propotions 
        ContactPoints = [ContactPoints; xx yy];
        ContactPointsElements = [ContactPointsElements; ElemenNumber*ones(n,1)]; % multiplication to correlate with points
        Contactweight = [Contactweight; ww];
    end    