 
function Visualization(Body,fig_number,vis, ShowNodeNumbers)
    
    FontSize = 10;
    DofsAtNode = Body.DofsAtNode;
    figure(fig_number)

    node = reshape(Body.q + Body.u, 2, []).';
    element = Body.nloc;
    Displacementing = reshape(Body.u, 2, []).';

    if vis == 'x'
        patching = Displacementing(:,1);
        % fprintf('presenting deformation along %s\n', vis);
    elseif vis == 'y'
        patching = Displacementing(:,2);
         % fprintf('presenting deformation along %s\n', vis);
        
    elseif vis == 'total' 
        patching = sqrt(Displacementing(:,1).^2 + Displacementing(:,2).^2);
       %  fprintf('presenting deformation along %s\n', vis);
    else
        warning('***Choose the correct axis for deformaion visualization\n')

        patching = Displacementing(:,2);
        fprintf('presenting deformation along %s\n', vis);
    end    
    h = patch('vertices', node, 'faces', element, 'facevertexcdata', patching, 'facecolor', 'interp');
    colormap(jet); 
    cb = colorbar;
    ylabel(cb, "Displacement: "+ vis, 'FontSize', [FontSize]);    % name it
    axis equal
    
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
