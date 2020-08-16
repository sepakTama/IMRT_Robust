function [z,IDs] = setProstate()

global planC
indexS = planC{end}

IMPro = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,14);
IMRec = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,17);
IMBla = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,16);
IDs = [14,17,16];

A=[];
[sz tmp] = size(IMPro);
A.size = sz;
A.mat = IMPro;
A.val = ones(sz,1)*52;
A.lp = [  0.95];
A.ld = [  75.6];
A.lparam = [1];
A.up = [0.05];
A.ud = [  83];
A.uparam = [1];
A.or_mat = A.mat;
A.or_val = A.val;

B=[];
[sz tmp] = size(IMRec);
B.size = sz;
B.mat = IMRec;
B.val = ones(sz,1) * 0;
B.lp = [];
B.ld = [];
B.lparam = [];
B.up = [0.30 0.10];
B.ud = [  70   75];
B.uparam = [1 1];
B.or_mat = B.mat;
B.or_val = B.val;

C=[];
[sz tmp] = size(IMBla);
C.size = sz;
C.mat = IMBla;
C.val = ones(sz,1) * 0;
C.lp = [];
C.ld = [];
C.lparam = [];
C.up = [0.30 0.10];
C.ud = [  70 75];
C.uparam = [1 1];
C.or_mat = C.mat;
C.or_val = C.val;

z = [{A} {C} {B}] ;

