function [z, IDs] = setTG119CShape()
global planC
indexS = planC{end}

IMC = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,1);
IMT = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,2);
%WC = getVoxelWeight(planC{indexS.IM}(end).IMDosimetry,1,planC);
%WT = getVoxelWeight(planC{indexS.IM}(end).IMDosimetry,2,planC);
IDs = [1,2];


A=[]; %->Outer Targetのデータ
sz = size(IMT); %return [r, c] in matrix IMT
sz = sz(1); %return first element => num of rows of IMT
A.size = sz;
A.mat = IMT;
A.val = ones(sz,1)*52; 
A.lp = [0.95];
A.ld = [  50];
A.lparam = [1]; %(5d)のunderline{P}
A.up = [0.10];
A.ud = [  55];
A.uparam = [1];
%A.weight = WT;
A.or_mat = A.mat;
A.or_val = A.val;

B=[]; %->Coreのデータ
sz = size(IMC);
sz = sz(1);
B.size = sz;
B.mat = IMC;
B.val = ones(sz,1) * 0;
B.lp = [];
B.ld = [];
B.up = [0.10];
B.ud = [  25]; %harder->10  easier->25
B.uparam = [1];
% B.uparam = [10];
%B.weight = WC;
B.or_mat = B.mat;
B.or_val = B.val;

z = [{A} {B}];

%Aにouter target、Bにcoreの情報を入れてそれらをCshapeに入れて返り値としている。