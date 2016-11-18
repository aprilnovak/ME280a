function [plot] = mesh_plots(coordinates, num_nodes, ro, LM)

% mesh of the semi-circle, with node numbering
x = coordinates(:,1);
y = coordinates(:,2);
scatter(x,y)
a = [1:num_nodes]';
b = num2str(a);
c = cellstr(b);
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
text(x+dx, y+dy, c);
xlim([-ro-2*dx, ro+5*dx])
ylim([-ro/2, ro + ro/2])
%saveas(gcf, 'Mesh', 'jpeg')

% output the LM for report
for e = 1:length(LM)
   fprintf('%i & %i & %i & %i \\\\\n', LM(e, 1), LM(e,2), LM(e,3), LM(e,4))
end

plot = 1;