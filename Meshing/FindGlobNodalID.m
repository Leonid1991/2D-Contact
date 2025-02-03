function NodalID = FindGlobNodalID(P,loc,shift)

x = loc.x;
y = loc.y;
shiftX = shift.x;
shiftY = shift.y;

% Returns global nodal ID 
NodalID=[];
tol=2*sqrt(eps);
jj=1;
for ii=1:size(P,1) % all over all points
  
    if isnumeric(x) && isnumeric(y)
       if (abs(P(ii,1)-x-shiftX) < tol) && (abs(P(ii,2)-y-shiftY) < tol)
          NodalID(jj)=ii;
          jj=jj+1;   
       end
    elseif isnumeric(x) && strcmp(y, 'all')
        if (abs(P(ii,1)-x-shiftX) < tol)
            NodalID(jj)=ii;
            jj=jj+1;   
        end
    elseif isnumeric(y) && strcmp(x, 'all')
        if (abs(P(ii,2)-y-shiftY) < tol)
            NodalID(jj)=ii; 
            jj=jj+1;   
        end     
    elseif strcmp(y, 'all') && strcmp(x, 'all')
        NodalID(jj)=ii; 
        jj=jj+1;   
    end          
end

