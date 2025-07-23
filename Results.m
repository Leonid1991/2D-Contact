clc, clear, close all;

Point = [3 5 9 17]; % 33

Nitsche6 = [0.007781 0.013337 0.014509 0.025585]./Point;
Nitsche8 = [0.007589 0.012836 0.013964 0.023988]./Point;
Nitsche10 = [0.002289 0.002728 0.002800 0.001727]./Point;
Nitsche12 = [0.000032 0.000034 0.000034 0.000041]./Point; % 0.000024 


Penalty6 = [0.768880 1.289705 2.248321 4.064437]./Point; % 7.361179 
Penalty8 = [0.264694 0.327033 0.373850 0.398130]./Point; % 0.400238 
Penalty10 = [0.003247 0.003487 0.003501  0.001923]./Point; % 0.002920
Penalty12 = [0.000033 0.000034 0.000035  0.0000738]./Point;

% Draw
h1=loglog(Point,Nitsche6,'-sk',...
          Point,Nitsche8,'--k',...
          Point,Nitsche10,'-vk',...
          Point,Nitsche12,'-ok',...
          Point,Penalty6,'-sb',...
          Point,Penalty8,'--b',...
          Point,Penalty10,'-vb',...
          Point,Penalty12,'-ob');
set(gca,'Fontsize',14);
set(gca,'FontName','Times New Roman');
grid on
grid minor
ylabel('$\frac{g}{n_{cont}}$', 'Interpreter', 'latex')
xlabel('$n_{cont}$', 'Interpreter', 'latex')
legend('Nitsche 1e6','Nitsche 1e8','Nitsche 1e10','Nitsche 1e12','Penalty 1e6','Penalty 1e8','Penalty 1e10','Penalty 1e12')
