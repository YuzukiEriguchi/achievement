%clear
SFS_start;
%%
% データの個数
endA = length(oxyzn(:,1))+1;
A = endA - 1;
%A = 1090 - 1;

% oとnをoxyznから表にそれぞれ格納
o = zeros(A, 1);
n = zeros(A, 1);

for i = 1 : A
    o(i) = oxyzn(i, 1);
    n(i) = oxyzn(i, 5);
end

% パワーの計算
pow = zeros(A, 1);
power1 = zeros(A, 1);
power = zeros(A, 1);

for i = 1:A
    for e = 1 : n(i)
        pow(i,e) = waveforms(i, e) .* waveforms(i, e); 
    end
end

pow(isnan(pow)) = 0;

for i = 1:A
    power1(i) = sum(pow(i, :)); %powerの正確なデータ

    power(i) = power1(i)^0.5;
end

powerd = 20 * log10(power);
waveforms(isnan(waveforms)) = 0;
%% 近接4点法の座標の計算
c = 343.9397;
fs = 48000;

dx = 0.05108;
dy = 0.05108;
dz = 0.05056;

x = zeros(A, 1);
y = zeros(A, 1);
z = zeros(A, 1);
d = zeros(A, 1);
e = zeros(A, 1);
dis = zeros(A, 1);

for i = 1:A
    to = o(i) ;
    tx = to + oxyzn(i, 2) / 16;
    ty = to + oxyzn(i, 3) / 16;
    tz = to + oxyzn(i, 4) / 16;
    
    
    ro = c * to / fs;
    rx = c * tx / fs;
    ry = c * ty / fs;
    rz = c * tz / fs;
    
    
    x(i) = (dx^2 + ro^2 - rx^2) / (2*dx);
    y(i) = (dy^2 + ro^2 - ry^2) / (2*dy);
    z(i) = (dz^2 + ro^2 - rz^2) / (2*dz);
end

%時間差o_sort(:)
for i = 1:length(x)
o_sort(i) = o(i)/48000; %oを格納
end

xyz = [x,y,z];
%%
%3次元
%clear
conf = SFS_config;
conf.resolution = 100;
fs = 44100;
c = 340;
f = 2000;
coord = [x, y, z]; % come from plane-wave indent direction
y_theta = zeros(length(x));
y_phi = zeros(length(x));
y_theta_from = zeros(length(x));
y_phi_from = zeros(length(x));
y_r = zeros(length(x));
y_come_from_sphy_unit = zeros(length(x),3);
for i = 1:length(x)
[y_theta(i), y_phi(i), ~] = cartesianToSpherical(coord(i,1),coord(i,2),coord(i,3));
[y_theta_from(i), y_phi_from(i), y_r(i)] = cartesianToSpherical(coord(i,1),coord(i,2),coord(i,3));
%y_sphy_unit = [sin(y_theta)*cos(y_phi) sin(y_theta)*sin(y_phi) cos(y_theta)];
y_come_from_sphy_unit(i,:) = [sin(y_phi_from(i))*cos(y_theta_from(i)) sin(y_theta_from(i))*sin(y_phi_from(i)) cos(y_phi_from(i))];
end
conf.plot.loudspeakers = false;

% set and grid creation
X = [-1 1];
Y = [-1 1];
Z = [-1 1];
[xx,yy,zz] = xyz_grid(X,Y,Z,conf);
[~,x1] = xyz_axes_selection(xx,yy,zz); % get first non-singleton axis
[theta, phi, r] = cartesianToSpherical(xx,yy,zz);

conf.secondary_sources.geometry = 'sphere'; % or 'circular'
conf.secondary_sources.number = 64;
x0 = secondary_source_positions(conf);
conf.plot.colormap = 'parula';
conf.plot.usenormalisation = 'max';

yl = x0(:,1:3);
[yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
[yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
%%
%2次元
%clear
conf = SFS_config;
conf.resolution = 160;
fs = 44100;
c = 340;
f = 2000;
coord = [x, y, z]; % come from plane-wave indent direction
y_theta = zeros(length(x));
y_phi = zeros(length(x));
y_theta_from = zeros(length(x));
y_phi_from = zeros(length(x));
y_r = zeros(length(x));
y_come_from_sphy_unit = zeros(length(x),3);
for i = 1:length(x)
[y_theta(i), y_phi(i), ~] = cartesianToSpherical(coord(i,1),coord(i,2),coord(i,3));
[y_theta_from(i), y_phi_from(i), y_r(i)] = cartesianToSpherical(coord(i,1),coord(i,2),coord(i,3));
%y_sphy_unit = [sin(y_theta)*cos(y_phi) sin(y_theta)*sin(y_phi) cos(y_theta)];
y_come_from_sphy_unit(i,:) = [sin(y_phi_from(i))*cos(y_theta_from(i)) sin(y_theta_from(i))*sin(y_phi_from(i)) cos(y_phi_from(i))];
end
conf.plot.loudspeakers = false;

% set and grid creation
X = [-1.6 1.6];
Y = [-1.6 1.6];
Z = 0;
[xx,yy,zz] = xyz_grid(X,Y,Z,conf);
[~,x1] = xyz_axes_selection(xx,yy,zz); % get first non-singleton axis
[theta, phi, r] = cartesianToSpherical(xx,yy,zz);

conf.secondary_sources.geometry = 'sphere'; % or 'circular'
conf.secondary_sources.number = 64;
x0 = secondary_source_positions(conf);
conf.plot.colormap = 'parula';
conf.plot.usenormalisation = 'max';

yl = x0(:,1:3);
[yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
[yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);


%%
%2次元(波面見る用)
%clear
conf = SFS_config;
conf.resolution = 160;
fs = 44100;
c = 340;
%f = 2000;
coord = [x, y, z]; % come from plane-wave indent direction
y_theta = zeros(length(x));
y_phi = zeros(length(x));
y_theta_from = zeros(length(x));
y_phi_from = zeros(length(x));
y_r = zeros(length(x));
y_come_from_sphy_unit = zeros(length(x),3);
for i = 1:length(x)
[y_theta(i), y_phi(i), ~] = cartesianToSpherical(coord(i,1),coord(i,2),coord(i,3));
[y_theta_from(i), y_phi_from(i), y_r(i)] = cartesianToSpherical(coord(i,1),coord(i,2),coord(i,3));
%y_sphy_unit = [sin(y_theta)*cos(y_phi) sin(y_theta)*sin(y_phi) cos(y_theta)];
y_come_from_sphy_unit(i,:) = [sin(y_phi_from(i))*cos(y_theta_from(i)) sin(y_theta_from(i))*sin(y_phi_from(i)) cos(y_phi_from(i))];
end
conf.plot.loudspeakers = false;

% set and grid creation
X = [-1.6 1.6];
Y = [-1.6 1.6];
Z = 0;
[xx,yy,zz] = xyz_grid(X,Y,Z,conf);
[~,x1] = xyz_axes_selection(xx,yy,zz); % get first non-singleton axis
[theta, phi, r] = cartesianToSpherical(xx,yy,zz);

conf.secondary_sources.geometry = 'sphere'; % or 'circular'
conf.secondary_sources.number = 64;
x0 = secondary_source_positions(conf);
conf.plot.colormap = 'parula';
conf.plot.usenormalisation = 'max';

yl = x0(:,1:3);
[yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
[yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
