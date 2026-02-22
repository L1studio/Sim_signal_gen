function Get_Ch(filePath)
% function Get_Ch(~)
% 打开并读取文件内容
global Chobj
fileID = fopen(filePath, 'r');
rawData = fread(fileID, inf, 'uint8');
fclose(fileID);

% 解码 JSON 数据
jsonData = jsondecode(char(rawData'));
nodeNum = numel(jsonData);
% 初始化输出变量
Chobj.Num = nodeNum;

Chobj.noise_power_dbw = zeros(nodeNum, 1);
Chobj.chId_V = strings(nodeNum, 1);%干扰源ID

Chobj.noise_power_dbw = zeros(nodeNum, 1);
for i = 1:nodeNum
   Chobj.chId_V(i) = jsonData(i).ch_id;

    Chobj.noise_power_dbw(i) = jsonData(i).noise_power_dbw;
end

end
