 
function PostProcessing(Body,fig_number,color,ShowNodeNumbers)
    
    FontSize = 10;
    DofsAtNode = Body.DofsAtNode;
    Results = Body.results;
    disp(' Mesh & DOFs & ux & uy \\)')
    for k=1:size(Results,1) 
        disp(sprintf('%dx%d & %4d  & %10.6f & %10.6f & %10.6f  \\\\',Results(k,:)))
    end   
    figure(fig_number)
    DrawingLimits=[0 2.5 -1.25 0.5 -0.1 0.1]; 
    % DrawMesh2412(Body.P0,Body.nloc, zeros(Body.nx,1),DofsAtNode,DrawingLimits,'b')
    DrawMesh2412(Body.P0,Body.nloc,Body.u,DofsAtNode,DrawingLimits,color)  
    set(gca, 'FontSize', [FontSize], 'FontName','Times New Roman');
    set(text, 'FontSize', [FontSize], 'FontName','Times New Roman');
    xlabel('{\it{X}} [m]','FontName','Times New Roman','FontSize',[FontSize]),ylabel('{\it{Y}} [m]','FontName','Times New Roman','FontSize',[FontSize]),zlabel('Z [m]','FontName','Times New Roman','FontSize',[FontSize]);
    % title(['Displacements at free end: {\it{u_x}} = ' ,num2str(Results(end,4)) ,'  ,{\it{u_y}} = ', num2str(Results(end,5))],'FontName','Times New Roman','FontSize',[FontSize]);     
    %% Writes nodal indeces
    if ShowNodeNumbers 
        nn= Body.nx/DofsAtNode;
        [k1,k2] = size(Body.P0);
        P = zeros(k1,k2); 
        for ii=1:nn
            P(ii,:)= Body.q((ii-1)*DofsAtNode+1:(ii-1)*DofsAtNode+DofsAtNode) +...
                     Body.u((ii-1)*DofsAtNode+1:(ii-1)*DofsAtNode+DofsAtNode);    
            text(P(ii,1),P(ii,2),0, int2str(ii))
        end
    end
