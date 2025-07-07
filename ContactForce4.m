function Fc =  ContactForce4(Body1,Body2,penalty,approach)  

         [Fc1,FcTarg1] = ContactSlaveMaster2(Body1, Body2, penalty,approach); % first is a slave, the second is a master 
         [Fc2,FcTarg2] = ContactSlaveMaster2(Body2, Body1, penalty,approach);

         % Fc = [Fc1; Fc2];
         % Fc = [FcTarg2; FcTarg1];
         Fc = [Fc1+FcTarg2; Fc2+FcTarg1];