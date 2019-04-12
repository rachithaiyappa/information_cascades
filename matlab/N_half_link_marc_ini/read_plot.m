% fileID = fopen('phi_star_data.txt','r');
% A = fscanf(fileID,'%f\t%f')

% infile = 'phi_star_data.txt';
% data = load(infile,'-ascii');
% 
% x = data(:,1);
% y = data(:,2);
% z = data(:,3);
% 
% scatter(x,y,200,z,'filled')

infile1 = 'phi_star_data_1.txt';
data = load(infile1,'-ascii');

x1 = data(:,1);
y1 = data(:,2);
z1 = data(:,3);

infile2 = 'phi_star_data_2.txt';
data = load(infile2,'-ascii');

x2 = data(:,1);
y2 = data(:,2);
z2 = data(:,3);

infile3 = 'phi_star_data_3.txt';
data = load(infile3,'-ascii');

x3 = data(:,1);
y3 = data(:,2);
z3 = data(:,3);

infile4 = 'phi_star_data_4.txt';
data = load(infile4,'-ascii');

x4 = data(:,1);
y4 = data(:,2);
z4 = data(:,3);n

infile5 = 'phi_star_data_5.txt';
data = load(infile5,'-ascii');

x5 = data(:,1);
y5 = data(:,2);
z5 = data(:,3);

x = (x1+x2+x3+x4+x5)/5;
y = (y1+y2+y3+y4+y5)/5;
z = (z1+z2+z3+z4+z5)/5;

scatter(x,y,200,z,'filled')