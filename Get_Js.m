function Get_Js(filePath)
global Jsobj
fileID = fopen(filePath, 'r');
rawData = fread(fileID, inf, 'uint8');  %二进制数据
fclose(fileID);
jsonData = jsondecode(char(rawData'));  %json解码
nodeNum = numel(jsonData);              %干扰源个数
%% 分配空间
Jsobj.Num=nodeNum;
Jsobj.jsId_V = strings(nodeNum, 1);%干扰源ID
Jsobj.type_V = strings(nodeNum, 1);%干扰源类型
Jsobj.transmitPower_V = zeros(nodeNum, 1);
Jsobj.freq_c_V = zeros(nodeNum, 1);

Jsobj.LFSI_band_V=zeros(nodeNum,1);
Jsobj.LFSI_sweep_rate_V=zeros(nodeNum,1);

Jsobj.NBI_modType_V = strings(nodeNum, 1);
Jsobj.NBI_multiplexingType_V = strings(nodeNum, 1);
Jsobj.NBI_symbolRate_V = zeros(nodeNum, 1);
Jsobj.NBI_shapingType_V = strings(nodeNum, 1);
Jsobj.NBI_arfa_V = zeros(nodeNum, 1);
Jsobj.NBI_modDepth_V = zeros(nodeNum, 1);
Jsobj.NBI_contPhase_V = strings(nodeNum, 1);
Jsobj.NBI_systemType_V=strings(nodeNum,1);

Jsobj.WBN_band_V=zeros(nodeNum,1);

%% 赋值
for i = 1:nodeNum
    Jsobj.jsId_V(i) = jsonData(i).js_id;
    Jsobj.type_V(i) = jsonData(i).type;
    Jsobj.transmitPower_V(i) = jsonData(i).transmit_power;
    Jsobj.freq_c_V(i) = jsonData(i).freq_c;

    Jsobj.LFSI_band_V(i)=jsonData(i).LFSI_band;
    Jsobj.LFSI_sweep_rate_V(i)=jsonData(i).LFSI_sweep_rate;

    Jsobj.NBI_modType_V (i)= jsonData(i).NBI_mod_type;
    Jsobj.NBI_multiplexingType_V(i) = jsonData(i).NBI_multiplexing_type;
    Jsobj.NBI_symbolRate_V(i) = jsonData(i).NBI_symbol_rate;
    Jsobj.NBI_shapingType_V(i) = jsonData(i).NBI_shaping_type;
    Jsobj.NBI_arfa_V(i) = jsonData(i).NBI_arfa;
    Jsobj.NBI_modDepth_V(i) = jsonData(i).NBI_mod_depth;
    Jsobj.NBI_contPhase_V(i) = jsonData(i).NBI_cont_Phase;
    Jsobj.NBI_systemType_V(i)=jsonData(i).NBI_system_type;
    Jsobj.WBN_band_V(i) = jsonData(i).WBN_band;
end
end
