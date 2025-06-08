%% verファイルの読み込み
filename = 'tu11t.ver';  
fileID = fopen(filename, 'r');
if fileID == -1
    error('Cannot open file: %s', filename);
end

data = fread(fileID, 'int16');  %short

fclose(fileID);  

% txtファイルに出力
txtFilename = 'tu11t.txt';
txtFileID = fopen(txtFilename, 'w');
if txtFileID == -1
    error('Cannot open file: %s', txtFilename);
end
fprintf(txtFileID, '%d\n', data);
fclose(txtFileID);

%% txtファイルの読み込み
fileID = fopen(txtFilename, 'r');
if fileID == -1
    error('Cannot open file: %s', txtFilename);
end
data = fscanf(fileID, '%f');
fclose(fileID);

%% txtファイルからoxyznとwaveformsの表を作成
a = 1;
i = 1;
oxyzn = [];  
waveforms = [];  

% oxyzn
while a < 2000
    %while i <= length(data)
    %txtデータの構成は, "oxyzn→波形データ"のループ
    o = data(i); i = i + 1;
    if o < 0
        o = o + 65536;
    end
    x = data(i); i = i + 1;
    y = data(i); i = i + 1;
    z = data(i); i = i + 1;
    N = data(i); i = i + 1;
   % waveform = data(i:i+N-1);
   if N == 0
       break
   end

    i = i + N;

    a = a + 1;
  
    oxyzn = [oxyzn; o, x, y, z, N];

   % waveforms = [waveforms; waveform', NaN(1, size(waveforms, 2) - N)];

   % end
end

maxN = max(oxyzn(:,5));
endA = a;

a = 1;
i = 6;

% waveforms
while a < endA
    waveform = data(i:i+oxyzn(a,5)-1);
    waveforms = [waveforms; waveform', NaN(1, maxN - (oxyzn(a,5)))];
    if oxyzn(a,5) == 0
        break
    end
    i = i + oxyzn(a,5) + 5;
    a = a + 1;
end