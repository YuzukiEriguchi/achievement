% データの個数
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
%% 座標をグラフ化
figure;
% 色などの設定
colors = zeros(A, 3);
colorsd = zeros(A, 3); % 初始化为白色
delta_to = zeros(A, 1);
delta_to1 = zeros(A, 1);

colors(1, :) = [1, 0, 0];
colorsd(1, :) = [1, 0, 0];

color = {'r', 'g', 'm', 'c', 'w'};
description = {'0ms', '50ms', '100ms', '200ms', '400ms'};

%インパルス応答の初めの時間から
for i = 1:A
    %delta_to(i) = (o(i, 1) - o(i-1, 1)) * 250 / 48000;
    delta_to(i) = o(i) / 48000;
    if delta_to(i) < 0.05 
        colors(i, :) = [1, 0, 0]; % 红色
    elseif delta_to(i) < 0.1
        colors(i, :) = [0, 1, 0]; % 绿色
    elseif delta_to(i) < 0.2
        colors(i, :) = [1, 0, 1]; % 洋红色
    elseif delta_to(i) < 0.4
        colors(i, :) = [0, 1, 1]; % 青色
    else
        colors(i, :) = [0, 0, 0]; % 黒色
    end
end

%直接音が入ってきた時間から
for i = 1:A
    %delta_to(i) = (o(i, 1) - o(i-1, 1)) * 250 / 48000;
    delta_to1(i) = o(i) / 48000 - o(1)/ 48000;

    if delta_to1(i) < 0.05 
        colorsd(i, :) = [1, 0, 0]; % 红色
    elseif delta_to1(i) < 0.1
        colorsd(i, :) = [0, 1, 0]; % 绿色
    elseif delta_to1(i) < 0.2
        colorsd(i, :) = [1, 0, 1]; % 洋红色
    elseif delta_to1(i) < 0.4
        colorsd(i, :) = [0, 1, 1]; % 青色
    else
        colorsd(i, :) = [0, 0, 0]; % 黒色
    end
end

% 3次元のグラフとして表示
scatter3(x, y, z, 7, colors, 'filled'); 
xlabel('x-axis (m)','Color','k');
ylabel('y-axis (m)','Color','k');
zlabel('z-axis (m)','Color','k');
fontsize(40, "points");
%title('dvnn','Color','k');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k'); 
set(gcf, 'Color', 'w'); 

xlim([-1000 1000]);
ylim([-800 800]);
zlim([-600 600]);

grid on; % 开启网格
axis equal; % 各轴刻度间隔保持一致
view(3); 

% パワーを表示
hold on;
theta = linspace(0, 2*pi, 100); % 定义圆的参数方程的角度

for i = 1:length(x)
    r = power(i)/7000; % powerの円の半径
    % xy平面上での表示
    xc_xy = r * cos(theta) + x(i);
    yc_xy = r * sin(theta) + y(i);
    zc_xy = z(i) * ones(size(theta));
    plot3(xc_xy, yc_xy, zc_xy, 'Color', colors(i, :), 'LineWidth', 0.3);

    % yz平面上での表示
    xc_yz = x(i) * ones(size(theta));
    yc_yz = r * cos(theta) + y(i);
    zc_yz = r * sin(theta) + z(i);
    plot3(xc_yz, yc_yz, zc_yz, 'Color', colors(i, :), 'LineWidth', 0.3);

    % xz平面上での表示
    xc_xz = r * cos(theta) + x(i);
    yc_xz = y(i) * ones(size(theta));
    zc_xz = r * sin(theta) + z(i);
    plot3(xc_xz, yc_xz, zc_xz, 'Color', colors(i, :), 'LineWidth', 0.3);
end

hold off;

% 凡例のカラーバーの表示
hold on; 
h(1) = patch([0 1 1 0], [0 0 1 1], 'r', 'EdgeColor', 'none'); 
h(2) = patch([2 3 3 2], [0 0 1 1], 'g', 'EdgeColor', 'none'); 
h(3) = patch([4 5 5 4], [0 0 1 1], 'm', 'EdgeColor', 'none'); 
h(4) = patch([6 7 7 6], [0 0 1 1], 'c', 'EdgeColor', 'none'); 
h(5) = patch([8 9 9 8], [0 0 1 1], 'k', 'EdgeColor', 'none');

hold off; 

% 凡例の表示
lgd = legend(h, '0 ms', '50 ms', '100 ms', '200 ms', '400 ms', 'Location', 'best');

set(lgd, 'TextColor', 'k');

patches = findobj(lgd, 'type', 'patch');

for p = patches'
    pos = get(p, 'Vertices');
    center = mean(pos, 1);
    new_pos = pos * 0.05 + center * 0.05;
    set(p, 'Vertices', new_pos);
end

%% powerの方向と強さを表示
figure;
hold on;
for i = 1:A
    dir = [x(i), y(i), z(i)];
    len = norm(dir); 
    if len ~= 0 
        dir = dir / len; 
    end
    
    scaled_dir = dir * (powerd(i));
   
    plot3([0, scaled_dir(1)], [0, scaled_dir(2)], [0, scaled_dir(3)], 'Color', colorsd(i, :));
end
xlabel('x-axis [m]','Color','k');
ylabel('y-axis [m]','Color','k');
zlabel('z-axis [m]','Color','k');
%title('direct','Color','k');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k'); 
set(gcf, 'Color', 'w'); 


%
xlim([-100 100]);
ylim([-100 100]);
zlim([-100 100]);
%

%{
xticks(0:20:30:40:50);
yticks(0:20:30:40:50);
zticks(0:20:30:40:50);
%}

grid on;
axis equal;
view(3);
hold off;

% 凡例
hold on; 
h(1) = patch([0 1 1 0], [0 0 1 1], 'r', 'EdgeColor', 'none'); 
h(2) = patch([2 3 3 2], [0 0 1 1], 'g', 'EdgeColor', 'none'); 
h(3) = patch([4 5 5 4], [0 0 1 1], 'm', 'EdgeColor', 'none'); 
h(4) = patch([6 7 7 6], [0 0 1 1], 'c', 'EdgeColor', 'none'); 
h(5) = patch([8 9 9 8], [0 0 1 1], 'k', 'EdgeColor', 'none');
hold off; 

lgd = legend(h, '0 ms', '50 ms', '100 ms', '200 ms', '400 ms', 'Location', 'best');

set(lgd, 'TextColor', 'k');

patches = findobj(lgd, 'type', 'patch');

for p = patches'
    pos = get(p, 'Vertices');
    center = mean(pos, 1);
    new_pos = pos * 0.05 + center * 0.05;
    set(p, 'Vertices', new_pos);
end
