function [input, x, t_history, time_history, total_time] = KishimotoIterateA(input, alpha, iterationNum, option)
% variables = [beamlet intensity, voxel intensity, and more]
% x means current good data.
% alpha means percentage of using voxels.
% in this case, now, there are no beta.

% this is version 4. and improvement of version 2.

% WARNING!
% this version is pending now.
% there is no code at lower_bound.


start_total = tic; 

if nargin <= 3 
    options = 'zzz'; % This option will be set in a future for something
end
    

t_history = [];
time_history = [];

EPS = 1e-2;
INF = 120;
boundAdjust = 0;

[tmp beamletNum] = size(input{1}.mat);
[tmp strNum] = size(input);

if alpha<0.99 
    for i=1:strNum
     [input{i}.mat input{i}.val] = chooseMatrix(input{i}.mat, input{i}.val, alpha);
     [input{i}.size tmp] = size(input{i}.mat);
    end
end
if alpha > 1.01 
   for i=1:strNum
       cur = min(1, alpha / input{i}.size); 
       [input{i}.mat input{i}.val] = chooseMatrix(input{i}.mat, input{i}.val, cur);
       [input{i}.size tmp] = size(input{i}.mat);
   end
end

for str = 1:strNum 
   fprintf('===Structure{%d}===', str);
   input{str}
end

% set equal constraints. 
Aeq = [];
voxelNum = 0;
exvariable = 0;


for str = 1:strNum
    Aeq = [Aeq; input{str}.mat]; 
    voxelNum = voxelNum + input{str}.size;
    exvariable = exvariable + (input{str}.size+1) * (max(size(input{str}.lp))+ max(size(input{str}.up))); 
end
% this is minimize variable.
exvariable = exvariable + 1;

Aeq = [Aeq -speye(voxelNum) ]; 
bcon = zeros(3,1);

beq = zeros(voxelNum, 1);

%ctmp = ones(beamletNum,1);

% set objective function
%c = [zeros(beamletNum, 1); ones(voxelNum, 1); zeros(exvariable, 1)];
c = [zeros(beamletNum, 1); zeros(voxelNum, 1); 1; zeros(exvariable, 1)];

% cntConst means "count constraints", cntSetVariable means "count length of
% spi" i.e. "the number of non-zero elements."
cntConst = 0;
cntSetVariable = 0;

prevZeta = [];
for str = 1:strNum
    strVoxelNum = input{str}.size;
    DVCNum = max(size(input{str}.lp)) + max(size(input{str}.up));
    
    prevZeta = [prevZeta; zeros(DVCNum, 1)];
    
    %penalty function
    %cntSetVariable = cntSetVariable + 2 * strVoxelNum * DVCNum;
    %cntConst = cntConst + strVoxelNum * DVCNum;
    
    %delta_j+w_j-eta_s
    cntSetVariable = cntSetVariable + 3 * strVoxelNum * DVCNum;
    cntConst = cntConst + strVoxelNum * DVCNum;
    
    % C-VaR constaraint 
    cntSetVariable = cntSetVariable + (1 + 1 + strVoxelNum) * DVCNum ;
    cntConst = cntConst + DVCNum;
end

% set initial constraints.

for str = 1:strNum
%   voxelDose = input{str}.mat * x;
%   maxim = max(voxelDose);
%   minim = min(voxelDose);
   minim = -1;
   maxim = 100;

    lsize = max(size(input{str}.lp));
   usize = max(size(input{str}.up));
   input{str}.lt = ones(lsize, 1) * minim;
   input{str}.ut = ones(usize, 1) * maxim;   
end

% this is answer variable.
finalx = [zeros(beamletNum + voxelNum, 1) ; 1000];
for z=1:iterationNum
    
    spi = zeros(cntSetVariable, 1); spj = zeros(cntSetVariable, 1); sps = zeros(cntSetVariable, 1);
    CVaRIndex = [];
% original ub and lb value.
    uborigin = [ones(beamletNum, 1)*INF; ones(voxelNum, 1) * 100];  
    lborigin = zeros(beamletNum + voxelNum, 1);  
    
    fixIndex = beamletNum;
    curConst = 1;
    curPos = 1;
    fixIndex = beamletNum;
    
    wIndex = beamletNum + voxelNum + 2;
%    param : minimize this parameter.
    paramIndex = wIndex-1;
    
    etaIndexes = [];
    exitflag = 1;
    
    coefValues = [];
    
    x0 = [finalx(1:(beamletNum + voxelNum+1)) ;zeros(exvariable, 1)  ];
    prevZetaIndex = 1;
    
    for str=1:strNum
        input{str}.prevDose = input{str}.mat * finalx(1:beamletNum);
        for ind = 1:max(size(input{str}.lp))
            if input{str}.lp(ind)>1-(1e-3)
               ld = input{str}.ld(ind);
               for voxel=1:input{str}.size
                  lborigin(fixIndex+voxel) = max(lborigin(fixIndex+voxel), ld);
               end
               continue; 
            end
            
            useCnt = 0;
            unuseIndex = zeros(input{str}.size, 1);
            for voxel = 1:input{str}.size
                
                checker = input{str}.lt(ind) - EPS;
                if input{str}.prevDose(voxel) < checker
                    unuseIndex(voxel) = 1;
                else
                    useCnt = useCnt+1;
                end
            end
            
            etaIndex = wIndex + useCnt ;
            etaIndexes = [etaIndexes; etaIndex];
            initIndex = wIndex-1;
            for voxel = 1:input{str}.size
                if unuseIndex(voxel)==1 
                    continue;
                end
                spi(curPos:curPos+2) = [curConst curConst curConst];
                spj(curPos:curPos+2) = [fixIndex+voxel wIndex etaIndex];
                sps(curPos:curPos+2) = [-1 -1 1];
                bcon(curConst) = 0;
                
                x0(wIndex) = max(0, prevZeta(prevZetaIndex) - finalx(fixIndex+voxel));
                
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 3;
            end
            unuseCnt = input{str}.size - useCnt;
            remainCnt = (1-input{str}.lp(ind)) * input{str}.size - unuseCnt;
            if remainCnt <=0
               exitflag = 0; 
            end
            coef = 1/remainCnt;
            
            cntDiff = 0;
            for voxel=1:input{str}.size
                if unuseIndex(voxel)==1
                    cntDiff = cntDiff + 1;
                    continue;
                end
                spi(curPos) = curConst;
                spj(curPos) = initIndex + voxel - cntDiff;
                sps(curPos) = coef;
                curPos = curPos + 1;
            end
            spi(curPos) = curConst;
            spj(curPos) = etaIndex;
            sps(curPos) = -1;
            
            %x0(wIndex) = prevZeta(prevZetaIndex);
            x0(wIndex) = input{str}.lt(ind);
            
            curPos = curPos + 1;
            wIndex = wIndex + 1;
            
            spi(curPos) = curConst;
            spj(curPos) = paramIndex;
            sps(curPos) = -input{str}.lparam(ind);
            
            ld = input{str}.ld(ind) + boundAdjust;
            bcon(curConst) = -ld;
            CVaRIndex = [CVaRIndex; curConst];
            
            curPos = curPos + 1;
            curConst = curConst + 1;
            
            coefValues = [coefValues; coef];
            if 0
            info = [unuseCnt remainCnt coef bcon(curConst-1)]
            end
            
            prevZetaIndex = prevZetaIndex + 1;
        end
        for ind = 1:max(size(input{str}.up))
            if input{str}.up(ind) < 1e-3
                ud = input{str}.ud(ind);
                for voxel=1:input{str}.size
                   uborigin(fixIndex+voxel) = min(uborigin(fixIndex+voxel), ud);
                   
                   spi(curPos:curPos+1) = [curConst curConst ];
                   spj(curPos:curPos+1) = [paramIndex fixIndex + voxel];
                   sps(curPos:curPos+1) = [-1 1];
                   
                   bcon(curConst) = ud;
                   curConst = curConst + 1;
                   curPos = curPos + 2;
                   
                end
                
               continue; 
            end
            
            useCnt = 0;
            unuseIndex = zeros(input{str}.size, 1);
            for voxel = 1:input{str}.size
                
               checker = input{str}.ut(ind) + EPS;
               if input{str}.prevDose(voxel) > checker
                  unuseIndex(voxel) = 1;
               else
                   useCnt = useCnt + 1;
               end
            end
            
            etaIndex = wIndex + useCnt;
            etaIndexes = [etaIndexes; etaIndex];
            initIndex = wIndex-1;
            for voxel = 1:input{str}.size
                if unuseIndex(voxel) == 1
                    continue; 
                end
                spi(curPos:curPos+2) = [curConst curConst curConst];
                spj(curPos:curPos+2) = [fixIndex+voxel wIndex etaIndex];
                sps(curPos:curPos+2) = [1 -1 -1];
                bcon(curConst) = 0;
                
                x0(wIndex) = max(0, finalx(fixIndex+voxel) - prevZeta(prevZetaIndex));
                
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 3;
            end
            unuseCnt = input{str}.size - useCnt;
            % unuseCnt
            remainCnt = input{str}.up(ind) * input{str}.size - unuseCnt;
            if remainCnt <=0 
                exitflag = 0;
            end
            coef = 1/remainCnt;
            
            cntDiff = 0;
            for voxel = 1:input{str}.size
                if unuseIndex(voxel) == 1
                   cntDiff = cntDiff + 1;
                   continue; 
                end
                spi(curPos) = curConst;
                spj(curPos) = initIndex + voxel - cntDiff;
                sps(curPos) = coef;
                curPos = curPos + 1;
            end
            spi(curPos) = curConst;
            spj(curPos) = etaIndex;
            sps(curPos) = 1;
            
            % x0(wIndex) = prevZeta(prevZetaIndex);
            x0(wIndex) = input{str}.ut(ind);
            
            curPos = curPos + 1;
            wIndex = wIndex + 1;
            
            spi(curPos) = curConst;
            spj(curPos) = paramIndex;
            sps(curPos) = -input{str}.uparam(ind);
            
            ud = input{str}.ud(ind) - boundAdjust;
            bcon(curConst) = ud;
            
            CVaRIndex = [CVaRIndex; curConst];
            
            curPos = curPos + 1;
            curConst = curConst + 1;
            
            coefValues = [coefValues; coef];
            info = [unuseCnt remainCnt coef bcon(curConst-1)];
            
            prevZetaIndex = prevZetaIndex + 1;
        end
        fixIndex = fixIndex + input{str}.size;
    end
%     [beamletNum voxelNum exvariable]
%     [cntConst cntSetVariable]
%     [size(spi)]
%     [spi spj sps]
%     bcon
%     etaIndexes
    
    range = (1:(curPos-1));
    Acon = sparse(spi(range), spj(range), sps(range));
    bcon = bcon(1:(curConst-1));
  
    cCur = c(1:wIndex-1);
    AeqCur = [Aeq sparse(voxelNum, wIndex - 1 - beamletNum - voxelNum)];
    
    lb = [lborigin; sparse(wIndex - 1 - beamletNum - voxelNum, 1)];
    lb(paramIndex) = -INF;
    
    ub = [uborigin; ones( wIndex - 1 - beamletNum - voxelNum, 1) * INF];
    nnzCnt = [nnz(AeqCur) nnz(Acon)];
    if 0
    %size(cCur)
    size(Acon)
    size(AeqCur)
    %size(bcon)
    %size(beq)
    %size(lb)
    size(spi)
    end
    
    x0 = x0(1:wIndex-1);
    
    % thread number
    if isunix() == 1 % on Linux
        [tmp, NumThreadsStr] = system('nproc');
    else % on Windows
        NumThreadsStr = getenv('NUMBER_OF_PROCESSORS');
    end
    NumThreads = str2num(NumThreadsStr);
    
    options = cplexoptimset('cplex');
    options.threads = ceil(NumThreads/2); % half threads will be used
    options.lpmethod = 4;
    options.display = 'on';
    options.barrier.convergetol = 1e-4;
    options.output.clonelog = -1;
 %   options.barrier.algorithm=3;
    %options.barrier.ordering = 3;
%     options.parallel = 8;
%     options.threads = 8;
    
    %options = cplexoptimset(options, 'Display', 'iter');
    
    %options = optimset(options,'Algorithm','interior-point');
    tic
    fval = -10;
    if exitflag == 1
        cplex_start = tic;
        [newx, fval, exitflag, output, lambda] = cplexlp(cCur, Acon, bcon, AeqCur, beq, lb, ub, x0, options); 
        t_history = [t_history; fval];
        time_history = [time_history; toc(cplex_start)];
        total_time = toc(start_total);
        
        fprintf('**** %d iter ****\n', z);
        t_history
        time_history
        total_time
    end
    %[newx, fval, exitflag, output, lambda] = linprog(c, Acon, bcon, Aeq, beq,lb, ub, []);
    toc
    
    range = beamletNum:wIndex-1;
%    normVal = norm(x0(range) - newx(range))    
%    [x0(range) - newx(range) x0(range) newx(range)]
    
    prevZeta = newx(etaIndexes);
    if 0
        prevZeta
        etaIndexes
    end
        
    
    newParamValue = newx(paramIndex);
      
%     rangeVal1 = paramIndex+1:etaIndexes(1)-1;
%     rangeVal2 = etaIndexes(1)+1:etaIndexes(2)-1;
%     rangeVal3 = etaIndexes(2)+1:etaIndexes(3)-1;
%     
%     [51.5 - newParamValue newx(etaIndexes(1)) - sum(newx(rangeVal1))  * coefValues(1)]
%     [53.5 + newParamValue newx(etaIndexes(2)) + sum(newx(rangeVal2)) * coefValues(2)]
%     [23.5 + newParamValue newx(etaIndexes(3)) + sum(newx(rangeVal3)) * coefValues(3)]
%     
%     
%      max(newx(rangeVal1))
%      max(newx(rangeVal2))
%      max(newx(rangeVal3))
if 0
    CVaRIndex
    lambda.ineqlin(CVaRIndex)
    size(bcon)
    size(lambda.ineqlin)
    transpose(lambda.ineqlin) * bcon
end

%     for i=1:max(size(bcon))
%        if bcon(i) > 0
%           i0.
%        end
%     end
    
    DEBUG = [fval exitflag];
    for str=1:strNum
        for ind=1:max(size(input{str}.lp))
            input{str}.lt(ind) = input{str}.ld(ind) - input{str}.lparam(ind) * newParamValue;
        end
        for ind=1:max(size(input{str}.up))
            input{str}.ut(ind) = input{str}.ud(ind) + input{str}.uparam(ind) * newParamValue;
        end
    end
    finalx = newx;
end

%x = finalx;
x = finalx(1:beamletNum);

if strcmp(option, 'xzt')
    x = finalx(1:(beamletNum + voxelNum + 1));
end

% for str=1:strNum
%     input{str}
%     calcError(input{str}, x)
% end

appendix = [];
total_time = toc(start_total);

fprintf('=== Information for right hand side ===\n');
for str = 1:strNum
    input{str}.le = zeros(size(input{str}.lp));
    input{str}.ue = zeros(size(input{str}.up));
    for ind = 1 : max(size(input{str}.lp))
        fprintf('L_%d^%.2f = %.1f  => %.4f\n', str, input{str}.lp(ind), ...
            input{str}.ld(ind), input{str}.lt(ind));
    end
    for ind = 1 : max(size(input{str}.up))
        fprintf('U_%d^%.2f = %.1f  => %.4f\n', str, input{str}.up(ind), ...
            input{str}.ud(ind), input{str}.ut(ind));
    end
end

fprintf('== t_history == ');
fprintf('%.2f ', t_history);
fprintf('\n');
fprintf('== time_history ==');
fprintf('%.3f ', time_history);
fprintf('\n');
fprintf('Total time = %.3f\n', total_time);




