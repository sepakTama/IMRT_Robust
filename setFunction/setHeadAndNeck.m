function [z, IDs] = setHeadAndNeck()
global planC
indexS = planC{end}

IMPTV = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,4);
IMCord = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,5);
IMLtParotid = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,6);
IMRtParotid = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,7);
IDs = [4,5,6,7];

A=[];
[sz tmp] = size(IMPTV);
A.size = sz;
A.mat = IMPTV;
A.val = ones(sz, 1) * 52;
A.lp = [ 0.90 0.99];
A.ld = [ 50 46.5];
A.lparam = [1 1];
A.up = [0.20];
A.ud = [  55];
A.uparam = [1];
A.or_mat = A.mat;
A.or_val = A.val;

B=[];
[sz tmp] = size(IMCord);
B.size = sz;
B.mat = IMCord;
B.val = ones(sz, 1) * 0;
B.lp = [];
B.ld = [];
B.up = [1e-15];
B.ud = [  40];
B.uparam = [1];
B.or_mat = B.mat;
B.or_val = B.val;

C=[];
[sz tmp] = size(IMLtParotid);
C.size = sz;
C.mat = IMLtParotid;
C.val = ones(sz, 1) * 0;
C.lp = [];
C.ld = [];
C.up = [0.50];
C.ud = [  20];
C.uparam = [1];
C.or_mat = C.mat;
C.or_val = C.val;

D=[];
[sz tmp] = size(IMRtParotid);
D.size = sz;
D.mat = IMRtParotid;
D.val = ones(sz, 1) * 0;
D.lp = [];
D.ld = [];
D.up = [0.50];
D.ud = [  20];
D.uparam = [1];
D.or_mat = D.mat;
D.or_val = D.val;

z = [{A} {B} {C} {D}];