function [z, IDs] = setMultiTarget() 
% function z = setMultiTarget(IMCen, IMSup, IMInf) 

global planC
indexS = planC{end}

IMCen = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,10);
IMSup = getInfluenceM(planC{indexS.IM}(end).IMDosimetry,11);
IMInf= getInfluenceM(planC{indexS.IM}(end).IMDosimetry,12);
IDs = [10, 11, 12]

A=[];
[sz tmp] = size(IMCen);
A.size = sz;
A.mat = IMCen;
A.val = ones(sz,1)*52;
A.lp = [0.99];
A.ld = [  50];
A.lparam = [1];
A.up = [0.10];
A.ud = [  53];
A.uparam = [1];
A.or_mat = A.mat;
A.or_val = A.val;

B=[];
[sz tmp] = size(IMSup);
B.size = sz;
B.mat = IMSup;
B.val = ones(sz,1) * 0;
B.lp = [0.99];
B.ld = [25];
B.lparam = [1];
B.up = [0.10];
B.ud = [  35];
B.uparam = [1];
B.or_mat = B.mat;
B.or_val = B.val;

C=[];
[sz tmp] = size(IMInf);
C.size = sz;
C.mat = IMInf;
C.val = ones(sz,1) * 0;
C.lp = [0.99];
C.ld = [12.5];
C.lparam = [1];
C.up = [0.10];
C.ud = [  25];
C.uparam = [1];
C.or_mat = C.mat;
C.or_val = C.val;

z = [{A} {B} {C}] ;

