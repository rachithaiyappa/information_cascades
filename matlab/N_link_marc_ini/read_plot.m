% fileID = fopen('phi_star_data_corrections.txt','r');
% A = fscanf(fileID,'%f\t%f')

infile = 'phi_star_data_corrections.txt';
data = load(infile,'-ascii');

x = data(:,1);
y = data(:,2);
z = data(:,3);

scatter(x,y,200,z,'filled')

% infile1 = 'phi_star_data_N_1.txt';
% data = load(infile1,'-ascii');
% 
% x1 = data(:,1);
% y1 = data(:,2);
% z1 = data(:,3);
% 
% infile2 = 'phi_star_data_N_2.txt';
% data = load(infile2,'-ascii');
% 
% x2 = data(:,1);
% y2 = data(:,2);
% z2 = data(:,3);
% 
% infile3 = 'phi_star_data_N_3.txt';
% data = load(infile3,'-ascii');
% 
% x3 = data(:,1);
% y3 = data(:,2);
% z3 = data(:,3);
% 
% infile4 = 'phi_star_data_N_4.txt';
% data = load(infile4,'-ascii');
% 
% x4 = data(:,1);
% y4 = data(:,2);
% z4 = data(:,3);
% 
% infile5 = 'phi_star_data_N_5.txt';
% data = load(infile5,'-ascii');
% 
% x5 = data(:,1);
% y5 = data(:,2);
% z5 = data(:,3);
% 
% infile6 = 'phi_star_data_N_6.txt';
% data = load(infile6,'-ascii');
% 
% x6 = data(:,1);
% y6 = data(:,2);
% z6 = data(:,3);
% 
% infile7 = 'phi_star_data_N_7.txt';
% data = load(infile7,'-ascii');
% 
% x7 = data(:,1);
% y7 = data(:,2);
% z7 = data(:,3);
% 
% infile8 = 'phi_star_data_N_8.txt';
% data = load(infile8,'-ascii');
% 
% x8 = data(:,1);
% y8 = data(:,2);
% z8 = data(:,3);
% 
% infile9 = 'phi_star_data_N_9.txt';
% data = load(infile9,'-ascii');
% 
% x9 = data(:,1);
% y9 = data(:,2);
% z9 = data(:,3);
% 
% infile10 = 'phi_star_data_N_10.txt';
% data = load(infile10,'-ascii');
% 
% x10 = data(:,1);
% y10 = data(:,2);
% z10 = data(:,3);
% 
% x = (x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10;
% y = (y1+y2+y3+y4+y5+y6+y7+y8+y9+y10)/10;
% z = (z1+z2+z3+z4+z5+z6+z7+z8+z9+z10)/10;
% 
% scatter(x,y,200,z,'filled')
