function ContactForceByPoints = ContacPointSetting(ContactPoints)

    n = ContactPoints.n;  
    name = ContactPoints.name;
    % Sanity check
    if name == "nodes"
        ContactForce = @ContactForceNode;
        n = 1; % we aren't interested 
    elseif (name == "Gauss") || (name == "LinSpace")
         
        if n < 2 
           warning("Number of points is not enough, it is set to 2 ")
           n = 2; 
        end  
         
        if  name == "Gauss" 
            ContactForce  = @ContactForceGauss;                
        elseif  name == "LinSpace"   
            ContactForce  = @ContactForceLinSpace;
        end

    else
        warning("ContactForceByPoints is set to node-based calculations")
        ContactForce  = @ContactForceNode;
        n = 1; % we aren't interested 
    end

    ContactForceByPoints.n = n;
    ContactForceByPoints.name = ContactForce;