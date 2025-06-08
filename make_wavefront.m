polar_coords = zeros(3,3);
o_sort = zeros(20,1); %時間差
%A = zeros(20,1); %power
for i = 1:length(x)

%fprintf('Powerが最大の仮想音源の座標は: (x, y, z) = (%.2f, %.2f, %.2f)\n', x(i), y(i), z(i));
%fprintf('最大のPowerの値は: %.2f\n', powerd(i));
%fprintf('この時のidxの値は：%d\n', i);
%A(i, :) = powerd(i); %powerを格納

% r, theta, phi の計算
r = sqrt(x(i)^2 + y(i)^2 + z(i)^2);
theta = atan2(y(i), x(i));  % atan2 は y/x の角度
phi = acos(z(i) / r); 
polar_coords(i, :) = [r, theta, phi];
r1(i,:) = r;
thetax(i, :) = theta;
phix(i, :) = phi;
o_sort(i, :) = o(i); %oを格納
end
disp('条件を満たした座標の極座標 (r, theta, phi) 一覧:');
disp('      r          theta       phi');
disp(polar_coords);

%%
%波面を表示させる用(高速化)
% Define the spatial grid
x1 = -9:c/(fs*2):9;
y1 = -9:c/(fs*2):9;
z1 = -9:c/(fs*2):9;
[X, Y, Z] = meshgrid(x1, y1, z1);
t1 = 0:1/6000:0.1;
L = length(x1);
L2 = 4800;
c = 343.9397; %音速
fs = 6000;
%Fs = c/; %サンプリング波数?

x_range = find(x1 >= -0.3 & x1 <= 0.3);
y_range = find(y1 >= -0.3 & y1 <= 0.3);
z_range = find(z1 >= -0.3 & z1 <= 0.3);

P1 = zeros(length(x_range), length(y_range), length(z_range), length(t1));
P1_sum = zeros(length(x_range), length(y_range), length(z_range),1, length(t1));

%P1 = zeros(length(x1), length(y1), length(z1), length(t1));
P_sub = zeros(length(x_range), length(y_range), length(z_range));

total_iterations =  L/2;
count = 0;

j = 1;
%for j = 1:2
    
P_freq = zeros(length(x_range), length(y_range), length(z_range), L);
%P_freq = zeros(length(x1), length(y1), length(z1), (L-1)/2);

for freq_idx = 1:L/2
    %freq_idx = 6;
    %freq = (freq_idx-1)*(Fs/L)*(c/(2*pi));
    freq = (freq_idx-1) * fs / L2;

P = zeros(length(x1), length(y1), length(z1));
    
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c;

kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));

% X, Y, Z は x1, y1, z1 のメッシュグリッドから生成された行列
kx_x = kx * x1; % 1D ベクトル
ky_y = ky * y1; % 1D ベクトル
kz_z = kz * z1; % 1D ベクトル
    
% 行列演算でPを計算
P = spec1(j, freq_idx) * exp(-1i * (kx_x' + ky_y + reshape(kz_z, 1, 1, [])));
%P = power(j)* exp(-1i * (kx_x' + ky_y + reshape(kz_z, 1, 1, [])));

    P_sub = P(x_range(1):x_range(length(x_range)), y_range(1):y_range(length(y_range)), z_range(1):z_range(length(z_range)));

    P_freq(:,:,:,freq_idx)=P_sub;
    
    %for ii4 = 1:length(t1)
    %P_freq(:,:,:,freq_idx,ii4)=P_sub*exp(1i*w*(t1(ii4)-o_sort(j)/fs));
    %end
    
    % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end

    
    
end

%P_sum = sum(P_freq, 4);
P1 = sum(P_freq, 4);

%P1_sum = P1_sum + P1;

% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
%end


%%
%波面を表示させる用(完成版)
fs = 6000; %サンプリング周波数
% Define the spatial grid
x1 = -9:c/(fs*2):9;
y1 = -9:c/(fs*2):9;
z1 = -9:c/(fs*2):9;
[X, Y, Z] = meshgrid(x1, y1, z1);
t1 = 0:1/6000:0.1;
L = length(x1);
%L2 = 4800;
c = 343.9397; %音速
ks = 2*pi*fs/c; %サンプリング波数?

x_range = find(x1 >= -0.3 & x1 <= 0.3);
y_range = find(y1 >= -0.3 & y1 <= 0.3);
z_range = find(z1 >= -0.3 & z1 <= 0.3);

P1 = zeros(length(x_range), length(y_range), length(z_range), length(t1));
P1_sum = zeros(length(x_range), length(y_range), length(z_range),1, length(t1));

%P1 = zeros(length(x1), length(y1), length(z1), length(t1));
P_sub = zeros(length(x_range), length(y_range), length(z_range));

total_iterations = 46 * L;
count = 0;

%j = 1;
for j = 1:46
    
P_freq = zeros(length(x_range), length(y_range), length(z_range), L);
%P_freq = zeros(length(x1), length(y1), length(z1), (L-1)/2);

for freq_idx = 1:L/2
    %freq = (freq_idx-1)*(Fs/L)*(c/(2*pi));
    freq = (freq_idx-1) * fs / L;

P = zeros(length(x1), length(y1), length(z1));
    
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c;

kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));

% X, Y, Z は x1, y1, z1 のメッシュグリッドから生成された行列
kx_x = kx * x1; % 1D ベクトル
ky_y = ky * y1; % 1D ベクトル
kz_z = kz * z1; % 1D ベクトル
    
% 行列演算でPを計算
P = spec1(j, freq_idx) * exp(-1i * (kx_x' + ky_y + reshape(kz_z, 1, 1, [])));
%P = power(j)* exp(-1i * (kx_x' + ky_y + reshape(kz_z, 1, 1, [])));

    P_sub = P(x_range(1):x_range(length(x_range)), y_range(1):y_range(length(y_range)), z_range(1):z_range(length(z_range)));

    %P_freq(:,:,:,freq_idx)=P_sub;
    for ii4 = 1:length(t1)
    P_freq(:,:,:,freq_idx,ii4)=P_sub*exp(1i*w*(t1(ii4)-o_sort(j)/48000));
    end
    % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end

    
    
end

%P_sum = sum(P_freq, 4);
P1 = sum(P_freq, 4);

P1_sum = P1_sum + P1;

% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end
%%
% 動画ファイルの名前と形式の指定
video_name = 'P1_sum_2.mp4';  % 保存する動画のファイル名
video = VideoWriter(video_name, 'MPEG-4');  % 動画の形式（ここではMP4）
video.FrameRate = 20;  % フレームレートの指定 (1秒あたりのフレーム数)

% 動画の書き込みを開始
open(video);

figure;
for ii4 = 1:length(t1)
slice(x2, y2, z2, real(-2*P1_sum(:,:,:,:,ii4)), 0, 0, 0);
xlabel('x-axis [m]');
ylabel('y-axis [m]');
zlabel('z-axis [m]');
xlim([-0.3 0.3]);
ylim([-0.3 0.3]);
zlim([-0.3 0.3]);
title(sprintf('t=%.4f [s]', t1(ii4)));
shading interp; % スムーズなシェーディング
colorbar; % カラーバーを表示
fontsize(14, "points");
clim([-30000 30000]);
axis square;
% 現在の図をフレームとしてキャプチャ
frame = getframe(gcf);  % 現在のfigureウィンドウをキャプチャ
writeVideo(video, frame);  % 動画にフレームを書き込む

%pause(1/20);  % フレームごとに少し停止して次に進む
end

% 動画ファイルを閉じて保存
close(video);

disp('動画が保存されました！');

%%
fs=6000;
x2 = -0.3:c/(fs*2):0.3;
y2 = -0.3:c/(fs*2):0.3;
z2 = -0.3:c/(fs*2):0.3;

%[X2, Y2, Z2] = meshgrid(x2, y2, z2);
figure;
for ii4 = 1:length(t1)
slice(x2, y2, z2, real(-2*P1_sum(:,:,:,:,ii4)), 0, 0, 0);
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
xlim([-0.3 0.3]);
ylim([-0.3 0.3]);
zlim([-0.3 0.3]);
title(sprintf('3D Pressure Field\n t=%.4f', t1(ii4)));
shading interp; % スムーズなシェーディング
colorbar; % カラーバーを表示
clim([-3*10^4 3*10^4]);
axis square;
pause(1/50);  % フレームごとに少し停止して次に進む
fontsize(20, "points");
end
