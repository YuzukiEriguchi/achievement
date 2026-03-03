%%
%波面用
%スピーカーアレイ 1~5
f = 1;
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=1;
end_wavn=5;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
t1 = 0.03:1/5000:0.1;
for order=min_order:max_order
    min_alpha=10;
    max_alpha=10;
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        count = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
       
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_1_5 = sum(T_sum, 3);

    end  
end


% スピーカーアレイ 6~10
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=6;
end_wavn=10;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        count = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_6_10 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 11~15
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=11;
end_wavn=15;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_11_15 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 16~20
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=16;
end_wavn=20;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
       
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_16_20 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 21~25
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=21;
end_wavn=25;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_21_25 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 26~30
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=26;
end_wavn=30;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_26_30 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 31~35
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=31;
end_wavn=35;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_31_35 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 36~40
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=36;
end_wavn=40;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_36_40 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 41~46
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=41;
end_wavn=46;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_41_46 = sum(T_sum, 3);  

        
    end  
end

%%
%波形用
%スピーカーアレイ 1~5
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=1;
end_wavn=5;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
t1 = 0.03:1/48000:0.1;
for order=min_order:max_order
    min_alpha=10;
    max_alpha=10;
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        count = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
       
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_1_5 = sum(T_sum, 3);

    end  
end


% スピーカーアレイ 6~10
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=6;
end_wavn=10;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        count = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_6_10 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 11~15
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=11;
end_wavn=15;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_11_15 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 16~20
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=16;
end_wavn=20;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
       
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_16_20 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 21~25
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=21;
end_wavn=25;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_21_25 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 26~30
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=26;
end_wavn=30;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_26_30 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 31~35
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=31;
end_wavn=35;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_31_35 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 36~40
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=36;
end_wavn=40;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_36_40 = sum(T_sum, 3);     
    end  
end

%スピーカーアレイ 41~45
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=41;
end_wavn=45;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_41_45 = sum(T_sum, 3);  

        
    end  
end

%スピーカーアレイ 46
N = round(2*pi*f*yl_r(1)/c,TieBreaker="plusinf");
min_order=7;
max_order=7;
L=50;
count=0;
start_wavn=46;
end_wavn=46;
wavn = end_wavn - start_wavn + 1;
fs = 48000;
Fs = 4800;
for order=min_order:max_order
    
    error_Lame=zeros(1,max_alpha-min_alpha+1);
    total_iterations = wavn*L*(max_alpha-min_alpha+1);
    for alpha = min_alpha:max_alpha
        conf = SFS_config;
        conf.secondary_sources.geometry = 'sphere';
        conf.secondary_sources.number = 256;
        x0 = secondary_source_positions(conf);
        conf.secondary_sources.geometry = 'custom';
        conf.plot.colormap = 'parula';
        conf.plot.usenormalisation = 'max';
        yl = x0(:,1:3);
        [yl_x,yl_y,yl_z] = xyz_grid(yl(:,1),yl(:,2),yl(:,3),conf);
        [yl_theta, yl_phi, yl_r] = cartesianToSpherical(yl_x,yl_y,yl_z);
        s = yl_r;
        % yl_r = abs(round(s.*(cos(yl_theta).^alpha.*sin(yl_phi).^alpha + sin(yl_theta).^alpha.*sin(yl_phi).^alpha+cos(yl_phi).^alpha).^(-1/alpha),3));
        yl_r = round(s.*(abs(cos(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha + abs(sin(yl_theta)).^alpha.*abs(sin(yl_phi)).^alpha+abs(cos(yl_phi)).^alpha).^(-1/alpha),3);
        x0_3d = yl_r.*[cos(yl_theta).*sin(yl_phi) sin(yl_theta).*sin(yl_phi) cos(yl_phi)];
        x0(:,1:3)=x0_3d;
   
        wav_idx = 0;
        
        T_ft = zeros(length(xx(1,:)), length(yy(1,:)),wavn, L, length(t1));
        for i = start_wavn:end_wavn
        wav_idx = wav_idx + 1;
        % Culculate speaker coeff
        
        %clear a_filter;
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        P = zeros(order^2+2*order+1,conf.secondary_sources.number);
        b = zeros(order^2+2*order+1,1);
        warning('off','all')
        row = 0;
        for n = 0:order
            Plm = legendre(n,cos(yl_phi),'norm');
            Plm_b = legendre(n,cos(y_phi_from(i,1)),'norm');
            rn = Rn(n,k,yl_r');
            for m = -n:n
                row = row+1;
                Pnm = Plm(abs(m)+1,:);
                Pnm_b = Plm_b(abs(m)+1,:);
                em = exp(-1i*m.*yl_theta');
                em_b = exp(-1i*m*y_theta_from(i,1));
                P(row,:) = rn.*(Pnm.*em);
                b(row,:) = Pnm_b*em_b;
            end
        end
        a_filter(:,fidx) = lsqr(P,b,4e-2,70);
        end
        %Calculate sound pressure
        
        for fidx = 1:50
        freq = (fs/Fs)*fidx;
        w = 2*pi*freq;
        k = 2*pi*freq/c;
        kx = k*cos(y_theta(i,1));
        ky = sqrt(k^2-kx^2);
        T = zeros(order+1,length(xx(1,:)),length(xx(1,:)));
        for n = 0:order
            xn = Xn(n,k,r);
            rn = Rn(n,k,yl_r);
            p = zeros(2*n+1,length(xx(1,:)),length(xx(1,:)));
            index = 1;
            for m=-n:n
                ynm = Ynm(n,m,theta,phi);
                ynm_conj = Ynm_conj(n,m,yl_theta,yl_phi);
                pe = a_filter(:,fidx).*(ynm_conj.').*rn;
                p_sp= sum(pe);
                const = 4*pi*p_sp;
                loc = ynm;
                p(index,:,:) = const.*loc;
                index = index +1;
            end
            if n~=0
                p_sum = sum(p);
                p_m = squeeze(p_sum);
            else
                p_m = squeeze(p);
            end
            T(n+1,:,:)= p_m.*xn;
        end
        if max_order ~= 0
            T = squeeze(sum(T));
        else
            T = squeeze(T);
        end
        for t = 1:length(t1)
            T_ft(:,:,wav_idx,fidx,t) = spec2(i,fidx)*T*exp(1i*w*(t1(t)-o_sort(i)));
            %T_ft(:,:,i,fidx,t) = T*exp(1i*w*(t1(t)-o_sort(i)/48000));
        end
            % 進捗の計算と表示
            count = count + 1;
            if mod(count, floor(total_iterations/100)) == 0
                fprintf('[%s] %.2f%% 完了\n', string(datetime('now','Format','HH:mm:ss')), (count / total_iterations)*100);
            end
        end
        T_sum = sum(T_ft, 4);
        % n回目の終了時点でメッセージを表示
        fprintf('%d回目の計算が終わりました\n', i);
        end
        T_sum_all_46 = sum(T_sum, 3);  

        
    end  
end
%%
%波面用
T_sum_all_1_46_Lame10 = T_sum_all_1_5+T_sum_all_6_10+T_sum_all_11_15+T_sum_all_16_20+T_sum_all_21_25+T_sum_all_26_30+T_sum_all_31_35+T_sum_all_36_40+T_sum_all_41_46;
%%
%波形用
T_sum_all_1_46_Lame10 = T_sum_all_1_5+T_sum_all_6_10+T_sum_all_11_15+T_sum_all_16_20+T_sum_all_21_25+T_sum_all_26_30+T_sum_all_31_35+T_sum_all_36_40+T_sum_all_41_45+T_sum_all_46;
%%
T_lame10_center = T_sum_all_1_46_Lame10(81,81,1,1,:);
T_lame10_center = squeeze(T_sum_center);
t2 = 0.03:1/5000:0.1;

%%
T_sum_center_lame10 = T_sum_all_1_46_Lame10(32,32,1,1,:);
T_lame10_center = squeeze(T_sum_center_lame10);
%%
figure;
subplot(2,1,1);
plot(t2, P2);
xlim([0.03 0.1]);
ylim([-1500 1000]);
%title('Ideal waveform');
xlabel('時間[s]');
ylabel('振幅');
fontsize(30, "points");

subplot(2,1,2);
plot(t2, T_sum_center);
xlim([0.03 0.1]);
ylim([-1500 1000]);
%title('Simulated waveform');
xlabel('時間[s]');
ylabel('振幅');
fontsize(30, "points");

%%
figure;
plot(t2, P2);
xlim([0.03 0.1]);
ylim([-1100 1100]);
%title('Ideal waveform');
xlabel('時間[s]');
ylabel('振幅');
fontsize(30, "points");

hold on;
plot(t2, T_sum_center, 'r');
xlim([0.03 0.1]);
ylim([-1100 1100]);
%title('Simulated waveform');
xlabel('時間[s]');
ylabel('振幅');
fontsize(30, "points");
hold off;

legend('Ideal waveform', 'Simulated waveform')

%%
% 波面全体の確認

xx1 = xx(1, :);
yy1 = yy(:, 1);
figure;
for t = 1:length(t1)
surf(xx1, yy1, real(T_sum_all_1_5(:,:,:,:,t)));
xlim([-1.6 1.6]);
ylim([-1.6 1.6]);
%clim([-1 500]);
clim([min(real(T_sum_all_1_5(81,81,1,1,:))) max(real(T_sum_all_1_5(81,81,1,1,:)))]);
title(sprintf('Pressure Field\n t=%.4f s', t1(t)));
view(2);
shading interp;
axis square;
colorbar;
pause(1/10);
end

%%
figure;
scatter3(yl_x, yl_y, yl_z, 50, 'filled'); % 50はマーカーサイズ、'filled'で塗りつぶし
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('3D Speaker Array Positions');
view(2);
axis equal;       % 各軸のスケールを同じにする（球形の確認に便利）
grid on;          % グリッドを表示

%%
m=10;
theta = linspace(0, 2 * pi, 100);
r = 1.5*(cos(theta).^m + sin(theta).^m).^(-1/m);
x_lame = r .* cos(theta);
y_lame = r .* sin(theta);
xx1 = xx(1, :);
yy1 = yy(:, 1);
figure;
for t = 1:length(t1)
    surf(xx1, yy1, real(T_sum_all_1_46_Lame10(:,:,:,:,t)));
    xlim([-1.6 1.6]);
    ylim([-1.6 1.6]);
    clim([min(real(T_sum_all_1_46_Lame10(81,81,1,1,:))) max(real(T_sum_all_1_46_Lame10(81,81,1,1,:)))]);
    title(sprintf('Pressure Field\n t=%.4f s', t1(t)));
    view(2);
    shading interp;
    axis square;
    colorbar;
    
    hold on;

    %{
    % 半径1.5mの点線円を原点中心に表示
    r = 1.5;
    theta = linspace(0, 2*pi, 300);
    x_circle = r * cos(theta);
    y_circle = r * sin(theta);
    plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);  % 黒の点線
    %}



    plot(x_lame, y_lame, 'LineWidth', 2);
    %axis equal;
    grid on;

    hold off;

    pause(1/10);
end

%% %正規化後
m=10;
theta = linspace(0, 2 * pi, 100);
r = 1.5*(cos(theta).^m + sin(theta).^m).^(-1/m);
x_lame = r .* cos(theta);
y_lame = r .* sin(theta);
xx1 = xx(1, :);
yy1 = yy(:, 1);
figure;
for t = 1:length(t1)
    surf(xx1, yy1, real(T_lame10_norm(:,:,t)));
    xlim([-1.6 1.6]);
    ylim([-1.6 1.6]);
    %clim([min(real(T_sum_all_1_46_Lame10(81,81,1,1,:))) max(real(T_sum_all_1_46_Lame10(81,81,1,1,:)))]);
    %clim([min(abs(real(T_lame10_norm(:)))) max(abs(real(T_lame10_norm(:))))]);
    clim([-1 1])
    title(sprintf('Pressure Field\n t=%.4f s', t1(t)));
    view(2);
    shading interp;
    axis square;
    colorbar;
    
    hold on;

    %{
    % 半径1.5mの点線円を原点中心に表示
    r = 1.5;
    theta = linspace(0, 2*pi, 300);
    x_circle = r * cos(theta);
    y_circle = r * sin(theta);
    plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);  % 黒の点線
    %}



    plot(x_lame, y_lame, 'LineWidth', 2);
    %axis equal;
    grid on;

    hold off;

    pause(1/10);
end
%%
%録画
t1 = 0.03:1/5000:0.1;
videoname = 'T_Lame10_2D.mp4';
video = VideoWriter(videoname, 'MPEG-4');
video.FrameRate = 12;
m=10;
theta = linspace(0, 2 * pi, 100);
r = 1.5*(cos(theta).^m + sin(theta).^m).^(-1/m);
x_lame = r .* cos(theta);
y_lame = r .* sin(theta);
z = 1500*ones(size(theta));
open(video);
figure;
for t = 1:length(t1)
surf(xx1, yy1, real(T_sum_all_1_46_Lame10(:,:,:,:,t)));
xlim([-1.6 1.6]);
ylim([-1.6 1.6]);
%clim([-1 500]);
clim([min(real(T_sum_all_1_46_Lame10(81,81,1,1,:))) max(real(T_sum_all_1_46_Lame10(81,81,1,1,:)))]);
title(sprintf('Pressure Field\n t=%.4f s', t1(t)));
view(2);
shading interp;
axis square;
colorbar;

hold on;
plot3(x_lame, y_lame, z, 'k--', 'LineWidth', 2);
hold off;

pause(1/5);

frame = getframe(gcf);
writeVideo(video, frame);
end

close(video);

disp('動画が保存されました！');

%%
T_Lame10_squ = squeeze(T_sum_all_1_46_Lame10);