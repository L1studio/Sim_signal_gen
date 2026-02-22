function [sigJ] = Gen_jssig(N,dem_flag,Fs,jsId)
global Jsobj

x = find(jsId == Jsobj.jsId_V);
type = Jsobj.type_V(x);
transmitPower=Jsobj.transmitPower_V(x);

LFSI_band = Jsobj.LFSI_band_V(x);
LFSI_sweep_rate= Jsobj.LFSI_sweep_rate_V(x);
WBN_band = Jsobj.WBN_band_V(x);
if type == "STI"
    %% 生成 STI 单音干扰
    init_phase = rand()*2*pi;
    sigJ = ones(1,N)*exp(init_phase);%载波信号 sqrt(2*Pj)
     sigJ = sigJ/sqrt(mean(abs(sigJ).^2));
    sigJ=sigJ*sqrt(10^(transmitPower/10));
    if dem_flag == 1
        figure()
        plot(10*log10(abs(fftshift(fft(sigJ)/length(sigJ))).^2));
        title('STI')
    end

elseif type == "LFSI"
    %% LFSI线性扫频干扰
    init_phase = rand()*2*pi;
    sigJ = generateChirpSignal(LFSI_band, LFSI_sweep_rate, init_phase, 10^(transmitPower/10), Fs, N);
    if dem_flag == 1
        fft_mat = con_fft_new((sigJ),512,512,0);  %
        figure()
        imagesc(10*log10(fft_mat)) % 绘制瀑布图
    end

elseif type =="NBI"
    %% 窄带干扰
    [~,sigJ] = Gen_basesig(N,Fs,Jsobj.jsId_V(x),0);
    if dem_flag ==1
        figure();
        plot(10*log10(abs(fftshift(fft(sigJ)/length(sigJ))).^2));
        title('NBI');
    end

elseif type =="WBN"
    WBN = overlap_retention(randn(1,N)+1i*randn(1,N),Fs,WBN_band);%产生复高斯白噪
    WBN = WBN/sqrt(mean(abs(WBN).^2));
    sigJ=WBN*sqrt(10^(transmitPower/10));%将信号功率转为指定的功率
    if dem_flag ==1
        figure();
        plot(10*log10(abs(fftshift(fft(sigJ)/length(sigJ))).^2));
        title('NBI');
    end
end


end
function signal = generateChirpSignal(band, sweep_rate, init_phase, power, fs, datalen)
% 参数说明：
% f_start: 初始频率 (Hz)
% f_end: 结束频率 (Hz)
% sweep_rate: 扫频速度 (Hz/s)
% init_phase: 初始相位 (弧度)
% power: 信号功率
% fs: 采样率 (Hz)
% duration: 信号持续时间 (秒)

% 时间向量
t = 0:1/fs:(datalen-1)/fs;
f_start = -band/2;
f_end = band/2;
% 计算瞬时频率
% 计算瞬时频率
instant_freq = f_start + sweep_rate * t;

% 处理超过扫频范围的情况
sweep_range = f_end - f_start;
instant_freq = mod(instant_freq - f_start, sweep_range) + f_start;

% 计算相位
phase = zeros(size(t));
cumulative_phase = init_phase;

for i = 2:length(t)
    % 计算当前频率下的相位增量
    delta_phase = 2 * pi * instant_freq(i-1) / fs;
    cumulative_phase = cumulative_phase + delta_phase;

    % 检查是否超过扫频范围
    if instant_freq(i) < instant_freq(i-1)
        % 计算频率重置后的相位调整
        cumulative_phase = cumulative_phase + 2 * pi * sweep_range / fs;
    end

    phase(i) = cumulative_phase;
end

% 生成复数形式的信号
signal = sqrt(power) * exp(1i * phase);

end
