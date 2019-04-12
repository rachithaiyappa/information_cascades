% fileID = fopen('phi_star_data.txt','r');
% A = fscanf(fileID,'%f\t%f')

infile = 'phi_star_data.txt';
data = load(infile,'-ascii');

x = data(:,1);
y = data(:,2);
z = data(:,3);

scatter(x,y,11,z,'filled')