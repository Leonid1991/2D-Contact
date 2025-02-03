function xlocAll= xlocAll2412(nloc)

[nl,m] = size(nloc);
xlocAll = zeros(nl,2*4);
for k = 1:nl
  n1 = [nloc(k,1)*2-1:1:nloc(k,1)*2];
  n2 = [nloc(k,2)*2-1:1:nloc(k,2)*2];
  n3 = [nloc(k,3)*2-1:1:nloc(k,3)*2];
  n4 = [nloc(k,4)*2-1:1:nloc(k,4)*2];
  xlocAll(k,:) = [n1 n2 n3 n4];
end
