function ContactPointfunc = ContactPointSetting(ContactPoints,n)
 
    addpath("Contact\ContactPointFunctions")
    
    % Sanity check
    if ContactPoints == "nodes"

        ContactPointfunc = @ContactPointNode;

    elseif (ContactPoints == "Gauss") || (ContactPoints == "LinSpace")
         
        
         
        if  ContactPoints == "Gauss" 

            ContactPointfunc  = @(ContactBody) ContactPointGauss(ContactBody,n);                
        
        elseif  ContactPoints == "LinSpace"

            if n < 2 
               warning("Number of points is not enough, it is set to 2 ")
               n = 2; 
            end  
            
            ContactPointfunc  = @(ContactBody) ContactPointLinSpace(ContactBody,n);
        
        end

    else

        warning("Contact Points are nodes")
        ContactPointfunc  = @ContactPointNode;

    end
