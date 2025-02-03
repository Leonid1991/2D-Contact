function loc = xlocChosen(DofsAtNode,nodes,comps)

loc = [];
for n=1:length(nodes)
  nn = nodes(n);
  loc = [loc DofsAtNode*(nn-1)+comps];
end
  