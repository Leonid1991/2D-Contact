function [X,U] = GetCoorDisp(ElemID,nloc,P0,uvec)
% Function finds initial coordinates X and displacement coordinates
% U for element ElemID 
% Initiate to zero
X_temp = zeros(4,2);
U_temp = zeros(4,2);
for ii = 1:4
    nodeno =  nloc(ElemID,ii);
    X_temp(ii,1) = P0(nodeno,1);
    X_temp(ii,2) = P0(nodeno,2);
    U_temp(ii,1) = uvec(2*nodeno-1);
    U_temp(ii,2) = uvec(2*nodeno);
end

X(1:2:7)=X_temp(:,1);
X(2:2:8)=X_temp(:,2);
U(1:2:7)=U_temp(:,1);
U(2:2:8)=U_temp(:,2);
X=X(:);
U=U(:);

