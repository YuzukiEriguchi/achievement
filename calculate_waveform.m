%%

%波形求める用(原点) パルス信号
% Define the spatial grid
x1 = -8:0.1:8;
y1 = -8:0.1:8;
z1 = -8:0.1:8;
[X, Y, Z] = meshgrid(x1, y1, z1);
Fs = 10; %空間サンプリング波数(1mあたり何回サンプリングしているか)
L = length(x1); %データの数
fs=48000;
t2 = 0:1/48000:0.1;
L2 = length(t2)-1;
P5 = zeros(1, length(t2));
total_iterations = length(y1) * length(z1)  *46* L2;
count = 0;
j=1;
%for j = 1:46
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 0:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算

c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数


kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (power(j)/L)*exp(-1i*(kx*x1(81)))*exp(-1i*(ky*y1(81)))*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end
    
    
    
    for ii4 = 1:length(t2)
        P5(ii4) = P5(ii4) + P(81)*exp(1i*w*(t2(ii4)-o_sort(j)/fs));
    end

    
    
    
    
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
%end



%%
%0~0.1s
%波形求める用(原点)
% Define the spatial grid
x1 = -8:0.1:8;
y1 = -8:0.1:8;
z1 = -8:0.1:8;
[X, Y, Z] = meshgrid(x1, y1, z1);
Fs = 10; %空間サンプリング波数(1mあたり何回サンプリングしているか)
L = length(x1); %データの数
fs = 48000;

t2 = 0:1/48000:0.1;
L2 = length(t2);
P2 = zeros(1, length(t2));
indices2=1:46;
total_iterations = length(z1)*length(indices2)* L2;
count = 0;
j=1;
%for j = indices2
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P2(ii4) = P2(ii4) + P(81)*exp(1i*w*(t2(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
%end

%%
%0.1~0.2s
t3 = 0.1:1/48000:0.2;
P3 = zeros(1, length(t3));
indices3=[47:127,1158:1161];
total_iterations = length(z1)*length(indices3)* L2;
count = 0;

for j = indices3
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P3(ii4) = P3(ii4) + P(81)*exp(1i*w*(t3(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%0.2~0.3s
t4 = 0.2:1/48000:0.3;
P4 = zeros(1, length(t4));
indices4=[128:220,1162:1168];
total_iterations = length(z1)*length(indices4)* L2;
count = 0;

for j = indices4
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P4(ii4) = P4(ii4) + P(81)*exp(1i*w*(t4(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%0.3~0.4s
t5 = 0.3:1/48000:0.4;
P5 = zeros(1, length(t5));
indices5=[221:308,1169:1177];
total_iterations = length(z1)*length(indices5)* L2;
count = 0;

for j = indices5
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P5(ii4) = P5(ii4) + P(81)*exp(1i*w*(t5(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    %fprintf('%d回目の計算が終わりました\n', j);
end

%%
%0.4~0.5s
t6 = 0.4:1/48000:0.5;
P6 = zeros(1, length(t6));
indices6=[309:393,1178:1179];
total_iterations = length(z1)*length(indices6)* L2;
count = 0;

for j = indices6
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P6(ii4) = P6(ii4) + P(81)*exp(1i*w*(t6(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%0.5~0.6s
t7 = 0.5:1/48000:0.6;
P7 = zeros(1, length(t7));
indices7=[394:462,1180:1181];
total_iterations = length(z1)*length(indices7)* L2;
count = 0;

for j = indices7
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P7(ii4) = P7(ii4) + P(81)*exp(1i*w*(t7(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%0.6~0.7s
t8 = 0.6:1/48000:0.7;
P8 = zeros(1, length(t8));
indices8=[463:547,1182:1190];
total_iterations = length(z1)*length(indices8)* L2;
count = 0;

for j = indices8
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P8(ii4) = P8(ii4) + P(81)*exp(1i*w*(t8(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%0.7~0.8s
t9 = 0.7:1/48000:0.8;
P9 = zeros(1, length(t9));
indices9=[548:633,1191:1194];
total_iterations = length(z1)*length(indices9)* L2;
count = 0;

for j = indices9
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P9(ii4) = P9(ii4) + P(81)*exp(1i*w*(t9(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%0.8~0.9s
t10 = 0.8:1/48000:0.9;
P10 = zeros(1, length(t10));
indices10=[634:715,1195:1197];
total_iterations = length(z1)*length(indices10)* L2;
count = 0;

for j = indices10
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P10(ii4) = P10(ii4) + P(81)*exp(1i*w*(t10(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%0.9~1.0s
t11 = 0.9:1/48000:1.0;
P11 = zeros(1, length(t11));
indices11=[716:800,1198:1203];
total_iterations = length(z1)*length(indices11)* L2;
count = 0;

for j = indices11
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P11(ii4) = P11(ii4) + P(81)*exp(1i*w*(t11(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%1.0~1.1s
t12 = 1.0:1/48000:1.1;
P12 = zeros(1, length(t12));
indices12=[801:902,1204:1208];
total_iterations = length(z1)*length(indices12)* L2;
count = 0;

for j = indices12
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P12(ii4) = P12(ii4) + P(81)*exp(1i*w*(t12(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%1.1~1.2s
t13 = 1.1:1/48000:1.2;
P13 = zeros(1, length(t13));
indices13=[903:993,1209:1220];
total_iterations = length(z1)*length(indices13)* L2;
count = 0;

for j = indices13
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P13(ii4) = P13(ii4) + P(81)*exp(1i*w*(t13(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%1.2~1.3s
t14 = 1.2:1/48000:1.3;
P14 = zeros(1, length(t14));
indices14=[994:1092,1221:1227];
total_iterations = length(z1)*length(indices14)* L2;
count = 0;

for j = indices14
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P14(ii4) = P14(ii4) + P(81)*exp(1i*w*(t14(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end

%%
%1.3~1.4s
t15 = 1.3:1/48000:1.4;
P15 = zeros(1, length(t15));
indices15=[1093:1157,1228:1234];
total_iterations = length(z1)*length(indices15)* L2;
count = 0;

for j = indices15
    P = zeros(length(z1));
    spec_indices = 0;
for freq = 1:fs/L2:((L2-1)/L2)*fs
spec_indices = spec_indices+1; % 全てのインデックスを丸めて計算
c = 343.9397; %音速
w = 2*pi*freq; %角周波数
k = w/c; %波数
kx = k*cos(thetax(j))*sin(phix(j));
ky = k*sin(thetax(j))*sin(phix(j));
kz = k*cos(phix(j));
%for ii1 = 1:length(x1)
    %for ii2 = 1:length(y1)
        for ii3 = 1:length(z1)
            P(ii3) = (spec2(j, spec_indices))*exp(-1i*(kx*x1(81))).*exp(-1i*(ky*y1(81)))'.*exp(-1i*(kz*z1(ii3)));
            %P(ii1, ii2, ii3) = (spec(j, freq*L2/fs)/L)*exp(-1i*(kx*x1(ii1)))*exp(-1i*(ky*y1(ii2)))*exp(-1i*(kz*z1(ii3)));
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('進捗: %.2f%% 完了\n', (count / total_iterations)*100);
            end
            

        end         
    %end
%end    
    for ii4 = 1:length(t2)
        P15(ii4) = P15(ii4) + P(81)*exp(1i*w*(t15(ii4)-o_sort(j)/fs));
    end
end
% n回目の終了時点でメッセージを表示
    fprintf('%d回目の計算が終わりました\n', j);
end


%%
P3a=P3(2:end);
P4a=P4(2:end);
P5a=P5(2:end);
P6a=P6(2:end);
P7a=P7(2:end);
P8a=P8(2:end);
P9a=P9(2:end);
P10a=P10(2:end);
P11a=P11(2:end);
P12a=P12(2:end);
P13a=P13(2:end);
P14a=P14(2:end);
P15a=P15(2:end);

%%
P_all = [P2, P3a, P4a, P5a, P6a, P7a, P8a, P9a, P10a, P11a, P12a, P13a, P14a, P15a];
