function [y] = dualIte(Input, alpha, iterationNum, option, lambdaMinus, lambdaPlus)

% Input : Information of structures
% alpha : If you lack of the memories, use this rate that we calcute how
% many voxels
% iterationNum : the number of repeat to calculate LP
% option : use to decide output opt, for example, 'onlyZ' is output dose
% lambda : each penalty, for example, [100;10]

start_total = tic; 

% time recorder
t_history = [];
time_history = [];

EPS = 1e-2;
INF = 120; % sufficient large integer
boundAdjust = 0; % 消してもいい？

[~, beamletNum] = size(Input{1}.mat);
[~, strNum] = size(Input);

lagrangian = [];

global init;

% chose randomly voxels to calculate
if alpha<1
    for i=1:strNum
     [Input{i}.mat, Input{i}.val] = chooseMatrix(Input{i}.mat, Input{i}.val, alpha);
     [Input{i}.size,  ~] = size(Input{i}.mat);
    end
elseif alpha > 1 
   for i=1:strNum
       cur = min(1, alpha / Input{i}.size); 
       [Input{i}.mat, Input{i}.val] = chooseMatrix(Input{i}.mat, Input{i}.val, cur);
       [Input{i}.size, ~] = size(Input{i}.mat);
   end
elseif alpha == 1 && isempty(init) == 0
    for i=1:strNum
         [Input{i}.mat, Input{i}.val] = chooseMatrix(Input{i}.mat, Input{i}.val, alpha);
         [Input{i}.size,  ~] = size(Input{i}.mat);
    end
end

for str = 1:strNum 
   fprintf('===Structure{%d}===', str);
   Input{str}
end

% set equal constraints. 
Aeq = [];
voxelNum = 0;
exvariable = 0;

for str = 1:strNum
    Aeq = [Aeq; Input{str}.mat];
    voxelNum = voxelNum + Input{str}.size;
    % exvariable : 
    %(\zeta - z)^+...voxelNum + (\theta - z)^+...voxelNum + \zeta...1
    % there are upper and lower constraints, so multiple 2
    exvariable = exvariable + 2 * (Input{str}.size+1) * (max(size(Input{str}.lp))+ max(size(Input{str}.up))); 
end

% this is minimize variable. -> おそらくいらない
% exvariable = exvariable + 1;

% Dx = z <=> [D -I][x , z]^T = 0 <=> Aeq * (variable) = beq
Aeq = [Aeq -speye(voxelNum) ]; 
beq = zeros(voxelNum, 1);

% to use constant in inequations
bcon = zeros(3,1);

% set objective function -> c = [x, z, t, exvariavles]
c = [zeros(beamletNum, 1); zeros(voxelNum, 1); 1; zeros(exvariable, 1)];

% cntSetVariable means "count length of spi" 
% i.e. "the number of non-zero elements."
% prevZeta holds previous Zeta

cntSetVariable = 0;
prevZeta = [];
for str = 1:strNum 
    strVoxelNum = Input{str}.size;
    DVCNum = max(size(Input{str}.lp)) + max(size(Input{str}.up));
    
    prevZeta = [prevZeta; zeros(DVCNum, 1)];
    
    % \zeta, z, \omega=()^+ -> use in DVC constraints with CVaR
    cntSetVariable = cntSetVariable + 3 * strVoxelNum * DVCNum;
    
    % z, \tau=(\theta-z)^+ -> use in objective function
    cntSetVariable = cntSetVariable + 2 * strVoxelNum * DVCNum;
    
    % \zeta, ()^+, t -> use in C-VaR constraint 
    cntSetVariable = cntSetVariable + (1 + 1 + strVoxelNum) * DVCNum ;
end

% set initial constraints.
for str = 1:strNum
   minim = -1;
   maxim = 100;

   lsize = max(size(Input{str}.lp));
   usize = max(size(Input{str}.up));
   Input{str}.lt = ones(lsize, 1) * minim;
   Input{str}.ut = ones(usize, 1) * maxim; 
   % -> use either be in the set R or not
end

% set the number of unusing voxel
CVaRIndex = [];

% this is answer variable.
finalx = [zeros(beamletNum + voxelNum, 1) ; 1000];

% main routine
for z=1:iterationNum
    % component of the matrix with nonzero -> (spi, spj) = sps
    spi = zeros(cntSetVariable, 1); 
    spj = zeros(cntSetVariable, 1); 
    sps = zeros(cntSetVariable, 1);
    
    % initialize the index of penalty \lambda 
    lambdaMinusIndex = 1;
    lambdaPlusIndex = 1;
    
    % original ub and lb value.
    uborigin = [ones(beamletNum, 1)*INF; ones(voxelNum, 1) * 100];  
    lborigin = zeros(beamletNum + voxelNum, 1);  
    
    % the index of component of the matrix with nonzero
    curConst = 1;
    curPos = 1;
    
    % variables = [beamletNum, voxelNum, 1(t), exvariable]
    fixIndex = beamletNum; 
    wIndex = beamletNum + voxelNum + 2;
    paramIndex = wIndex-1; 
    
    % initialize c -> initialize penalty \lambda
    c(wIndex:wIndex+exvariable-1) = 0;
    
    etaIndexes = []; % set the index of Zeta
    exitflag = 1; % use if occur a bug
    
    coefValues = []; % set the value of coefficient in DVC with CVaR
    
    % x0 is starting point (->初期点でそんなに時間が変わらない)
    x0 = [finalx(1:(beamletNum + voxelNum+1)) ;zeros(exvariable, 1)  ];
    prevZetaIndex = 1;
    
    % strLIndexはデータのバグで\alpha=0,1が二つ以上あるときに使用
    strLIndex = 0;
    strUIndex = 0;
    for str=1:strNum
        Input{str}.prevDose = Input{str}.mat * finalx(1:beamletNum);
        % l - Pt <= z
        for ind = 1:max(size(Input{str}.lp))
            if Input{str}.lp(ind) > 1-(1e-3) 
               ld = Input{str}.ld(ind);
               for voxel=1:Input{str}.size
                   lborigin(fixIndex+strLIndex+voxel) = max(lborigin(fixIndex+strLIndex+voxel), ld); 
                  
                   spi(curPos:curPos+1) = [curConst curConst]; 
                   spj(curPos:curPos+1) = [paramIndex fixIndex + voxel];
                   sps(curPos:curPos+1) = [-1 -1];
                   
                   bcon(curConst) = -ld;
                   curConst = curConst + 1;
                   curPos = curPos + 2;
               end
               strLIndex = strLIndex + Input{str}.size;
               % ここに新しくpenalty項を入れる -> wIndexの数を考慮する必要ある？
               strsize = Input{str}.size;
               for voxel = 1:Input{str}.size
                   spi(curPos:curPos+1) = [curConst, curConst];
                   spj(curPos:curPos+1) = [fixIndex+voxel, wIndex];
                   sps(curPos:curPos+1) = [-1, -1];
                   bcon(curConst) = -ld;
                   
                   x0(wIndex) = 0;
                   c(wIndex) = lambdaMinus(lambdaMinusIndex) / strsize;
                   
                   wIndex = wIndex + 1;
                   curConst = curConst + 1;
                   curPos = curPos + 2;
               end
               lambdaMinusIndex = lambdaMinusIndex + 1;
               continue; 
            end
            
            % count the number of the voxel in set R
            useCnt = 0;
            unuseIndex = zeros(Input{str}.size, 1);
            for voxel = 1:Input{str}.size
                checker = Input{str}.lt(ind) - EPS;
                % if z<L-Pt then this coresponding index turns from 0 to 1
                if Input{str}.prevDose(voxel) < checker
                    unuseIndex(voxel) = 1;
                else
                    useCnt = useCnt+1;
                end
            end
            
            % eta is used when zeta is needed in inequarity
            etaIndex = wIndex + useCnt ;
            etaIndexes = [etaIndexes; etaIndex];
            initIndex = wIndex-1;
            for voxel = 1:Input{str}.size
                if unuseIndex(voxel)==1 
                    continue;
                end
                % zeta(etaIndex) - z(fixIndex+voxel) - w(wIndex) <= 0(bcon)
                spi(curPos:curPos+2) = [curConst curConst curConst];
                spj(curPos:curPos+2) = [fixIndex+voxel wIndex etaIndex];
                sps(curPos:curPos+2) = [-1 -1 1];
                bcon(curConst) = 0;
                
                x0(wIndex) = max(0, prevZeta(prevZetaIndex) - finalx(fixIndex+voxel));
                
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 3;
            end
            unuseCnt = Input{str}.size - useCnt;
            % remainCnt = (1 - alpha)*v_s - |R|
            remainCnt = (1-Input{str}.lp(ind)) * Input{str}.size - unuseCnt;
            if remainCnt <=0
               exitflag = 0; 
            end
            coef = 1/remainCnt;
            
            cntDiff = 0;
            for voxel=1:Input{str}.size
                if unuseIndex(voxel)==1
                    cntDiff = cntDiff + 1;
                    continue;
                end
                % coef * w /// w = (zeta - z)^+
                % from this, spi does not change, to express "zeta - coef*w >= L - Pt"
                spi(curPos) = curConst;
                spj(curPos) = initIndex + voxel - cntDiff;
                sps(curPos) = coef;
                curPos = curPos + 1;
            end
            spi(curPos) = curConst;
            spj(curPos) = etaIndex;
            sps(curPos) = -1;
            
            x0(wIndex) = Input{str}.lt(ind);
            
            curPos = curPos + 1;
            wIndex = wIndex + 1;
            
            spi(curPos) = curConst;
            spj(curPos) = paramIndex;
            sps(curPos) = -Input{str}.lparam(ind);
            
            ld = Input{str}.ld(ind) + boundAdjust;
            bcon(curConst) = -ld;
            CVaRIndex = [CVaRIndex; unuseCnt];
            
            if z == iterationNum
            lagrangian = [lagrangian ; curConst]; %新規で追加
            end
            
            curPos = curPos + 1;
            curConst = curConst + 1;
            
            % coefの情報はあまり必要ない？|R|だけでいいのでは？
            coefValues = [coefValues; coef];
            
            % penalty term (theta - z)^+ /// - z - w <= - theta(ld)
            strsize = useCnt;         
            for voxel = 1:Input{str}.size
                if unuseIndex(voxel)==1
                    continue;
                end
                spi(curPos:curPos+1) = [curConst, curConst];
                spj(curPos:curPos+1) = [fixIndex+voxel, wIndex];
                sps(curPos:curPos+1) = [-1, -1];
                bcon(curConst) = -ld;
                
                x0(wIndex) = 0; % 初期点を入れるかは考え中
                % set penalty to objective function c(wIndex)
                c(wIndex) = lambdaMinus(lambdaMinusIndex) / strsize;
                
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 2;
            end
            
            if 0
            info = [unuseCnt remainCnt coef bcon(curConst-1)]
            end
            
            prevZetaIndex = prevZetaIndex + 1;
            lambdaMinusIndex = lambdaMinusIndex + 1;
        end
        
        % from this, Upper DVC constraints with CVaR
        for ind = 1:max(size(Input{str}.up))
            if Input{str}.up(ind) < 1e-3
                ud = Input{str}.ud(ind);
                for voxel=1:Input{str}.size
                   uborigin(fixIndex+strUIndex+voxel) = min(uborigin(fixIndex+strUIndex+voxel), ud);
                   
                   spi(curPos:curPos+1) = [curConst curConst ];
                   spj(curPos:curPos+1) = [paramIndex fixIndex + voxel];
                   sps(curPos:curPos+1) = [-1 1];
                   
                   bcon(curConst) = ud;
                   curConst = curConst + 1;
                   curPos = curPos + 2; 
                end
                strUIndex = strUIndex + Input{str}.size;
                % ここにpenalty項をいれる
                strsize = Input{str}.size;
                for voxel = 1:Input{str}.size
                    spi(curPos:curPos+1) = [curConst, curConst];
                    spj(curPos:curPos+1) = [fixIndex+voxel, wIndex];
                    sps(curPos:curPos+1) = [1, -1];
                    bcon(curConst) = ud;
                    
                    x0(wIndex) = 0; % 初期点
                    c(wIndex) = lambdaPlus(lambdaPlusIndex) / strsize;
                    
                    wIndex = wIndex + 1;
                    curConst = curConst + 1;
                    curPos = curPos + 2;
                end
                lambdaPlusIndex = lambdaPlusIndex + 1;
                continue; 
            end
            
            useCnt = 0;
            unuseIndex = zeros(Input{str}.size, 1);
            for voxel = 1:Input{str}.size
                
               checker = Input{str}.ut(ind) + EPS;
               if Input{str}.prevDose(voxel) > checker
                  unuseIndex(voxel) = 1;
               else
                   useCnt = useCnt + 1;
               end
            end
            
            etaIndex = wIndex + useCnt;
            etaIndexes = [etaIndexes; etaIndex];
            initIndex = wIndex-1;
            for voxel = 1:Input{str}.size
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
            unuseCnt = Input{str}.size - useCnt;
            remainCnt = Input{str}.up(ind) * Input{str}.size - unuseCnt;
            if remainCnt <=0 
                exitflag = 0;
            end
            coef = 1/remainCnt;
            
            cntDiff = 0;
            for voxel = 1:Input{str}.size
                if unuseIndex(voxel) == 1
                   cntDiff = cntDiff + 1;
                   continue; 
                end
                % zeta + coef*w <= U + Pt
                spi(curPos) = curConst;
                spj(curPos) = initIndex + voxel - cntDiff;
                sps(curPos) = coef;
                curPos = curPos + 1;
            end
            spi(curPos) = curConst;
            spj(curPos) = etaIndex;
            sps(curPos) = 1;

            x0(wIndex) = Input{str}.ut(ind);
            
            curPos = curPos + 1;
            wIndex = wIndex + 1;
            
            spi(curPos) = curConst;
            spj(curPos) = paramIndex;
            sps(curPos) = -Input{str}.uparam(ind);
            
            ud = Input{str}.ud(ind) - boundAdjust;
            bcon(curConst) = ud;
            
            CVaRIndex = [CVaRIndex; unuseCnt];
            
            if z == iterationNum
            lagrangian = [lagrangian ; curConst]; %新規で追加
            end
            
            curPos = curPos + 1;
            curConst = curConst + 1;
            
            % coefの情報はいらない？
            coefValues = [coefValues; coef];
            
            % from this, penalty term
            strsize = useCnt;
            for voxel = 1:Input{str}.size
                if unuseIndex(voxel) == 1
                    continue;
                end
                spi(curPos:curPos+1) = [curConst, curConst];
                spj(curPos:curPos+1) = [fixIndex+voxel, wIndex];
                sps(curPos:curPos+1) = [1, -1];
                bcon(curConst) = ud;
                
                x0(wIndex) = 0; % 変更すると速くなる？
                c(wIndex) = lambdaPlus(lambdaPlusIndex) / strsize;
                
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 2;
            end
            
            info = [unuseCnt remainCnt coef bcon(curConst-1)];
            
            prevZetaIndex = prevZetaIndex + 1;
            lambdaPlusIndex = lambdaPlusIndex + 1;
        end
        fixIndex = fixIndex + Input{str}.size;
    end
    % from this, construct the form to use CPLEX
    range = (1:(curPos-1));
    Acon = sparse(spi(range), spj(range), sps(range));
    bcon = bcon(1:(curConst-1));
  
    cCur = c(1:wIndex-1);
    AeqCur = [Aeq sparse(voxelNum, wIndex - 1 - beamletNum - voxelNum)];
    
    lb = [lborigin; sparse(wIndex - 1 - beamletNum - voxelNum, 1)];
    lb(paramIndex) = -INF;
    
    ub = [uborigin; ones( wIndex - 1 - beamletNum - voxelNum, 1) * INF];
    nnzCnt = [nnz(AeqCur) nnz(Acon)];

    x0 = x0(1:wIndex-1);
    
    % thread number
    if isunix() == 1 % on Linux
        [~, NumThreadsStr] = system('nproc');
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
%   options.barrier.ordering = 3;
%   options.parallel = 8;
%   options.threads = 8;
%   options = cplexoptimset(options, 'Display', 'iter');
%   options = optimset(options,'Algorithm','interior-point');

    tic
    fval = -10;
    if exitflag == 1
        cplex_start = tic;
        [newx, fval, exitflag, output, lambda] = cplexlp(cCur, Acon, bcon, AeqCur, beq, lb, ub, x0, options); 
        t_history = [t_history; newx(paramIndex)];
        time_history = [time_history; toc(cplex_start)];
        total_time = toc(start_total);
        fprintf('**** %d iter ****\n', z);
        t_history
        time_history
        total_time
    end
    toc
    
    range = beamletNum:wIndex-1;
    
    prevZeta = newx(etaIndexes);
    if 0
        prevZeta
        etaIndexes
    end
        
    
    newParamValue = newx(paramIndex);

    DEBUG = [fval exitflag];
    % calculate  L - Pt and U + Pt, and put this information to structure
    for str=1:strNum
        for ind=1:max(size(Input{str}.lp))
            Input{str}.lt(ind) = Input{str}.ld(ind) - Input{str}.lparam(ind) * newParamValue;
        end
        for ind=1:max(size(Input{str}.up))
            Input{str}.ut(ind) = Input{str}.ud(ind) + Input{str}.uparam(ind) * newParamValue;
        end
    end
    finalx = newx;
    
    % if the value of t is nonpositive then LP iteration is finished.
    % if newx(paramIndex) < 0
    %    break;
    % end
end

x = finalx(1:beamletNum);

% it is useful to use output of opt for changing option
% if you use the function "lambdaIteration.m", then use option "onlyZ"
if strcmp(option, 'all')
    opt = finalx;
elseif strcmp(option, 'onlyX')
    opt = finalx(1:beamletNum);
elseif strcmp(option, 'onlyZ')
    opt = finalx((beamletNum+1):(beamletNum+voxelNum)); 
elseif strcmp(option, 'xzt')
    opt = finalx(1:(beamletNum + voxelNum + 1));
end

appendix = [];
total_time = toc(start_total);

fprintf('=== Information for right hand side ===\n');
for str = 1:strNum
    Input{str}.le = zeros(size(Input{str}.lp));
    Input{str}.ue = zeros(size(Input{str}.up));
    for ind = 1 : max(size(Input{str}.lp))
        fprintf('L_%d^%.2f = %.1f  => %.4f\n', str, Input{str}.lp(ind), ...
            Input{str}.ld(ind), Input{str}.lt(ind));
    end
    for ind = 1 : max(size(Input{str}.up))
        fprintf('U_%d^%.2f = %.1f  => %.4f\n', str, Input{str}.up(ind), ...
            Input{str}.ud(ind), Input{str}.ut(ind));
    end
end

%CVaRの表示
fprintf('== CVaRIndex == ');
fprintf('%.2f ', CVaRIndex);
fprintf('\n');

% 変更点 lambdaの表示
fprintf('== lambdaPlus == ');
fprintf('%.3f ', lambdaPlus);
fprintf('\n');
fprintf('== lambdaMinus == ');
fprintf('%.3f ', lambdaMinus);
fprintf('\n');

fprintf('== t_history == ');
fprintf('%.2f ', t_history);
fprintf('\n');
fprintf('== time_history ==');
fprintf('%.3f ', time_history);
fprintf('\n');
fprintf('Total time = %.3f\n', total_time);

A_A = AeqCur;
assignin('base', 'A_A' , A_A); 
assignin('base', 'output', output);
assignin('base', 'lambda', lambda);

fprintf('== lagrangian == \n');
fprintf('%.3f \n', lambda.ineqlin(lagrangian));
y = 0;




