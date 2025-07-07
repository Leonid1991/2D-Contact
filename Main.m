clc, clear, close all
format long;
addpath("InnerForce")
addpath("Meshing")
addpath("Postprocess")
addpath("Contact")
%########## Reads element's data ###############################
ElementData;   
%########## Reads problem's data ###############################
ProblemData;
%########## Element positioning (from (0.0) coord. ) ###########
Body1.shift.x = 0;
Body1.shift.y = 0;

Body2.shift.x = 0;
Body2.shift.y = -Body2.Ly;
%#################### Mesh #########################################
dx = 5;
dy = 1;

Body1.nElems.x = dx;
Body1.nElems.y = dy;

Body2.nElems.x = dx;
Body2.nElems.y = dy;

Body1 = CreateFEMesh(DofsAtNode,Body1);
Body2 = CreateFEMesh(DofsAtNode,Body2);
% % %#################### BC  ###########################################
Body1.loc.x = 0; 
Body1.loc.y = 'all';  % Number (Location of nodes along the axis) or 'all' can be an option

Body2.loc.x = 0; 
Body2.loc.y = 'all'; 

Body1 = CreateBC(Body1);
Body2 = CreateBC(Body2); 

bc = [Body1.bc Body2.bc];

%##################### Loadings ######################
% local positions (assuming all bodies in (0,0) )
Body1.Fext.x = 0; 

Body1.Fext.y =-62.5*10^(6);

Body1.Fext.loc.x = Body1.Lx;
Body1.Fext.loc.y = 'all';

Body2.Fext.y = 0*62.5*10^(6); 
Body2.Fext.x = 0;
Body2.Fext.loc.x = Body2.Lx;
Body2.Fext.loc.y = 'all';

%##################### Egde nodes #########################
Body1.edge1.loc.x = Body1.Lx;
Body1.edge1.loc.y = 0;

Body1.edge2.loc.x = Body1.Lx;
Body1.edge2.loc.y = Body1.Ly;

Body2.edge1.loc.x = Body2.Lx;
Body2.edge1.loc.y = 0;

Body2.edge2.loc.x = Body2.Lx;
Body2.edge2.loc.y = Body2.Ly;

%##################### Contact ############################
approach = 1; % 0 - none; 1- penalty, 2- Nitsche (linear of gap), 3- Nitsche (nonlinear of gap), 4 - all items    
pn = 1e12;
penalty = pn;
 
% Identification of possble contact surfaces
% local positions (assuming all bodies in (0,0) )
Body1.contact.loc.x = 'all';
Body1.contact.loc.y = 0;   
Body1.contact.nodalid = FindGlobNodalID(Body1.P0,Body1.contact.loc,Body1.shift);

Body2.contact.loc.x = 'all';
Body2.contact.loc.y = Body2.Ly;  
Body2.contact.nodalid = FindGlobNodalID(Body2.P0,Body2.contact.loc,Body2.shift);% 

%##################### Newton iter. parameters ######################
imax=15;
tol=1e-4;         
steps= 20;
total_steps = 0;
titertot=0;  
% %#################### Processing ######################
for ii = 1:steps

    % Update forces, supported loading types: linear, exponential, quadratic, cubic;
    type = "linear";
    Body1 = CreateFext(ii,steps,Body1,type);
    Body2 = CreateFext(ii,steps,Body2,type);
    
    % contact convergence
    for jj = 1:imax
        tic;

        total_steps = total_steps + 1;
        % interacation of two bodies
        [Fc,Kc] = Contact(Body1,Body2,penalty,approach);

        Body1 = Elastic(Body1);
        Body2 = Elastic(Body2);
         
        % Assemblance of stiffnesses
        Ke = [            Body1.Fint.K zeros(Body1.nx,Body2.nx);
              zeros(Body2.nx,Body1.nx)            Body2.Fint.K];
        K = Ke + Kc;
    
        % Assemblance of forces
        Fe = [Body1.Fint.vec; Body2.Fint.vec];
        Fext = [Body1.Fext.vec; Body2.Fext.vec];
        ff =  Fe - Fext + Fc;
        
        % Calculations
        K_bc = K(bc,bc); 
        ff_bc = ff(bc);
        deltaf = ff_bc/norm(Fext(bc)); 
        uu_bc = -K_bc\ff_bc;      
    
        % Separation
        Body1.u(Body1.bc) = Body1.u(Body1.bc) + uu_bc(1:Body1.ndof);
        Body2.u(Body2.bc) = Body2.u(Body2.bc) + uu_bc(Body1.ndof + 1:end);
         
        Body1 = Analys(Body1);
        Body2 = Analys(Body2);
        
        Gap = Gapfunc(Body1,Body2);
        titer=toc;
        titertot=titertot+titer;

        if printStatus(deltaf, uu_bc, tol*sum(bc), ii, jj, imax, steps, titertot,Gap)
            break;  
        end 

    end
    
end
% %##################### Post-Processing ######################
fig_number = 1; 
ShowNodeNumbers = false;
disp('Static test')
PostProcessing(Body1,fig_number,'b',ShowNodeNumbers)
PostProcessing(Body2,fig_number,'r',ShowNodeNumbers);


Gap = Gapfunc(Body1,Body2);
fprintf('total steps is %d  \n', total_steps );
fprintf('total gap is %10.6f  \n', Gap )


