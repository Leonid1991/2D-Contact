function ContactForce = ContactForceSetting(ContactPoints,n)
 
    % Sanity check
    if ContactPoints == "nodes"

        ContactForce = @ContactForceNode;

    elseif (ContactPoints == "Gauss") || (ContactPoints == "LinSpace")
         
        if n < 2 
           warning("Number of points is not enough, it is set to 2 ")
           n = 2; 
        end  
         
        if  ContactPoints == "Gauss" 

            ContactForce  = @(ContactBody,TargetBody,penalty,approach) ContactForceGauss(ContactBody,TargetBody,penalty,approach,n);                
        
        elseif  ContactPoints == "LinSpace"

            ContactForce  = @(ContactBody,TargetBody,penalty,approach) ContactForceLinSpace(ContactBody,TargetBody,penalty,approach,n);
        
        end

    else

        warning("ContactForceByPoints is set to node-based calculations")
        ContactForce  = @ContactForceNode;

    end
