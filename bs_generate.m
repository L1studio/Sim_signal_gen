clc;clear;clear global;
close all;
tic
global Txobj Rxobj
for sh = 1:70
    rng("shuffle");
    waterfall_raws = 34;%瀑布图行数
    Get_Tx('Tx_para.json');Get_Rx('Rx_para.json');
    Get_Ch('Ch_para.json');Get_Js('Js_para.json');
    if(1)
        shuffled_signals = rand_json(Rxobj.freq_rf,Rxobj.sample_band);
        Get_Tx('Tx_para.json');
    end
    %%
    % 信号生成
    Rxobj.fs = Rxobj.sample_band;
    Rxobj.F_ = Rxobj.fs/(Rxobj.nfft-1);
    Rxobj.offset = round(Rxobj.nfft/4);
    Rxobj.T_ = Rxobj.offset/Rxobj.fs;
    Rxobj.ts = Rxobj.T_* waterfall_raws;
    dataLen = round(Rxobj.ts*Rxobj.fs);
    waterfall_raws = floor((dataLen-Rxobj.nfft)/round(Rxobj.nfft/4))+1;
    % --------- 预分配（用普通数组承接 parfor 输出） ---------

    numTx   = Txobj.Num;
    txIds   = Txobj.txId_V;           % 先取出需要的向量，避免在 parfor 里索引 Txobj
    band    = zeros(1, numTx);        % 原来是 Txobj.band，用数组替代
    BaseSig = zeros(numTx, dataLen);
    fs = Rxobj.fs;
    % 如果你想在 parfor 内写 shuffled_signals(x).number，
    % 必须保证该字段已存在（先统一创建）
    if ~isfield(shuffled_signals, 'number')
        [shuffled_signals(1:numTx).number] = deal(0);
    end

    %% 接收信号生成  噪声
    %% 产生信号源
    txobj = Txobj;
  
    parfor x = 1:numTx
        [bx, sigx] = Gen_basesig(dataLen, fs, txIds(x), 1,txobj);
        band(x)    = bx;              % 只写普通数组
        BaseSig(x,:) = sigx;
    end
 
    for x = 1:numTx
        shuffled_signals(x).number = x;
    end
    Txobj.band = band;  

    for ss = 1:10
        tic
        stringx = "第"+num2str(sh)+"轮"+"，"+"第"+num2str(ss)+"次";
        disp([stringx sh]);
        close all;
        rfSig = zeros(numTx,dataLen);
        shuffled_signals = rand_except_basejson(Rxobj.freq_rf,Rxobj.sample_band,shuffled_signals);
        Get_Tx('Tx_para.json');
        freqC_V = Txobj.freqC_V;
        freq_rf = Rxobj.freq_rf;
        idx = [shuffled_signals.number];   
        BaseSig_perm = BaseSig(idx, :);
        parfor x = 1:numTx
            rfSig(x,:)  = Fc_change(BaseSig_perm(x,:) ,(freqC_V(x) - freq_rf)/fs);
        end
        f = (0:Rxobj.F_:Rxobj.fs)-Rxobj.fs/2+Rxobj.freq_rf;% 频率向量
        t = 0:Rxobj.T_:Rxobj.ts-Rxobj.T_;
        [~,pink_noise] = Gen_pinknoise(f(1),f(end),length(rfSig));
        gaussian_noise = Gen_gaussiannoise(length(rfSig));
        noise = gaussian_noise*0.5+pink_noise;
        spec = fftshift(abs(fft(noise,length(rfSig))/length(rfSig)).^2);

        f_= linspace(-Rxobj.fs/2, Rxobj.fs/2, length(rfSig))+Rxobj.freq_rf;     % 基带频率 [-fs/2, +fs/2]
        
        for i = 1:length(shuffled_signals)
            idx = find(f_ >= shuffled_signals(i).start_freq & f_ <= shuffled_signals(i).end_freq);
            P_noise(i) = sum((spec(idx))) ;  % 带内噪声功率
            P_sig = sum(abs(rfSig(i,:)).^2) / length(rfSig(i,:));     % 带内信号功率
            target_snr = 10^(shuffled_signals(i).receive_snr / 10);
            target_Psig = target_snr * P_noise(i);
            scale = sqrt(target_Psig / P_sig);
            rfSig(i,:) = rfSig(i,:) * scale;
        end
        
        RxSig = sum(rfSig,1)+noise;
        fft_mat = con_fft_new((RxSig),Rxobj.nfft,Rxobj.nfft,round(Rxobj.nfft/4));  %
        fft_mat = single(fft_mat);
        % figure()
        % imagesc(f, t,10*log10(fft_mat)) % 绘制瀑布图

        % figure
        spec = sum(fft_mat(1:30,:));spec(end-150:end) =spec(end-300:end-150) ;
        spec = single(spec);
        shuffled_signals = snr_annotate(spec,f,shuffled_signals);

        % plot(f,10*log10(spec))
        % hold on

        % colors = lines(length(shuffled_signals));  % 使用 lines colormap 生成不同颜色
        % for i = 1:length(shuffled_signals)
        %     outline = round((shuffled_signals(i).start_freq-f(1))/Rxobj.F_):round((shuffled_signals(i).end_freq-f(1))/Rxobj.F_);
        %     plot(f(outline),10*log10(spec(outline)),'-r','LineWidth',1);
        %     hold on
        % end

        % cut_point = 4096;
        % for i = 1:floor(length(spec)/cut_point)
        %     line([f(cut_point*(i)), f(cut_point*(i))], [10*log10(min(spec)), 10*log10(max(spec))], 'Color', 'r', 'LineStyle', '--') ;
        %     hold on
        % end
        % for i = 1:floor(length(spec)/cut_point)
        %     figure()
        %     fft_son_mat = fft_mat(:,cut_point*(i-1)+1:cut_point*(i));
        %     imagesc(f(cut_point*(i-1)+1:cut_point*(i)), t,10*log10(fft_son_mat)) % 绘制瀑布图
        % end
        save_path = "E:\m_pro\博二\宋钰\第三章\data";  % 确保 save_path 是 string 类型
        timestamp = datestr(now, 'yyyymmdd');  % 获取当前日期（年月日格式）
        save_fold = save_path + "\" + "simu_" + timestamp;  % 使用 + 拼接路径，确保是 string 类型

        % 检查文件夹是否存在，不存在则创建
        if ~exist(save_fold, 'dir')  % 使用 'dir' 参数，明确检查文件夹
            mkdir(save_fold);  % 创建文件夹
        end
        count = length(dir(save_fold))-2+1;
        count = sprintf('%05d', count);  % 格式化为5位数字，不足前面补零

        filename = "ori_matrix_(dBW)_chapter4"+"_"+count;
        image = 10*log10(spec);
        % fft_mat = 10*log10((fft_mat));

        maxVal = max(image);
        minVal = min(image);

        F_start = -Rxobj.fs/2+Rxobj.freq_rf;
        F_resolution = Rxobj.F_;
        annotation =[];
        save(fullfile(save_fold, filename), "fft_mat","Rxobj","shuffled_signals","maxVal","minVal","F_start","F_resolution","image","annotation");
        toc
    end
end