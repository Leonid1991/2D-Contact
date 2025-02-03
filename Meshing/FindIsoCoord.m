function [xi,eta] = FindIsoCoord(X,U,point)


    r_eta1=Nm_2412(0,-1)*(X+U);
    r_eta2=Nm_2412(0, 1)*(X+U);
    r_xi0_line=[r_eta1';
                r_eta2'];
    [~,~,eta_01] = distance2curve(r_xi0_line,point,'linear');
    eta = -1 + 2 * eta_01;

    r_xi1=Nm_2412(-1,0)*(X+U);
    r_xi2=Nm_2412( 1,0)*(X+U);
    r_eta0_line=[r_xi1';
                 r_xi2'];
    [~,~,xi_01] = distance2curve(r_eta0_line,point,'linear');
    xi = -1 + 2 * xi_01;
  