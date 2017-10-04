clear all;

load 'variables.m';
load 'cells.m';

t_init = variables(4);
output_t = variables(5);
numberReports = variables(6);
dt = (output_t-t_init)/numberReports;

for report = 0:numberReports 
   
   filename = strcat('solution', sprintf('%03d',[report]), '.m');

   u = load(filename);
   
   time = t_init + dt*report
   
   figure(report+1);
   trimesh(cells,u(:,1),u(:,2),u(:,3));
   axis square;
   axis tight;
   xlabel('x');
   ylabel('y');
   zlabel('u');
   title(strcat('Porous Medium Equation Solution at t=', num2str(time)))

 end
