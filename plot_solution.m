clear all;

load 'variables.m';
load 'cells.m';

numberReports = variables(6);

for report = 0:numberReports 
   
   filename = strcat('solution', sprintf('%03d',[report]), '.m');

   u = load(filename);   
   
   figure(report+1);
   trimesh(cells,u(:,1),u(:,2),u(:,3));
   axis square;
   axis tight;
   xlabel('x');
   ylabel('y');
   zlabel('u');

 end
