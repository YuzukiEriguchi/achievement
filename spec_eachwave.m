%% ゼロ詰め
waveforms1 = zeros(length(x), 629);
waveforms1(1:length(x), 1:300) = waveforms;

waveforms2 = zeros(length(x), 4800);
waveforms2(1:length(x), 1:300) = waveforms;


%% フーリエ変換(4800)
L3 = 4800;

fidx = (0:(L3-1))*fs/L3;
fft_waveforms2 = fft(waveforms2, [], 2);
spec2 = fft_waveforms2/L3;
phase = angle(fft_waveforms2);

%% フーリエ変換(300)
L4 = 629;
fidx1 = (0:(L4-1))*fs/L4;
fft_waveforms1 = fft(waveforms1, [], 2);
spec1 = fft_waveforms1/L4;
%%
fidx = (0:(L4-1))*(Fs/L4)*(c/(2*pi));
figure;
plot(fidx, 20*log10(spec(1, :)));
%xlim([0, fs/2]);
