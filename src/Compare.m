function [lowGap, highGap] = Compare(Input, z)
% return the gap between DVC and optimal value
% à¯êîzÇÃÉxÉNÉgÉãÇÃÇ›

% need search each structure's size

[~, strNum] = size(Input);
voxelNum = 0;
endNum = 0;
initNum = 1;

% initialize
lp = -1;
ld = -1;
up = 110;
ud = 1000;
lpIndex = -1;
upIndex = 110;

lowGap = [];
highGap = [];

for str = 1:strNum
    endNum = endNum + Input{str}.size;
    voxelNum = Input{str}.size;
    strZ = sort(z(initNum:endNum));
    
    for ind = 1:max(size(Input{str}.lp))
        lp = Input{str}.lp(ind);
        ld = Input{str}.ld(ind);
        lpIndex = ceil(voxelNum*lp); % index of alpha rate of the DVC
        lowGap = [lowGap; (ld - strZ(voxelNum - lpIndex + 1))]; % if satisfies lower DVC then lowgap is necative
        
        fprintf('===Structure{%d}===\n', str);
        fprintf('lp={%.2f}, ld={%.1f}, gap={%.2f}\n', lp, ld, ld-strZ(voxelNum - lpIndex + 1));
        
        
    end
    
    for ind = 1:max(size(Input{str}.up))
        up = Input{str}.up(ind);
        ud = Input{str}.ud(ind);
        upIndex = ceil(voxelNum*up); % index of alpha rate of the DVC
        highGap = [highGap; (strZ(voxelNum - upIndex + 1) - ud)]; % if satisfies upper DVC then lowgap is necative
        fprintf('up={%.2f}, ud={%.1f}, gap={%.2f}\n', up, ud, strZ(voxelNum - upIndex + 1)-ud);
    end
    initNum = initNum + Input{str}.size;
end

end