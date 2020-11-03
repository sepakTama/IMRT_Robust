function [Input, x, opt,  t_history, time_history, total_time] = UncertainIterate(Input, alpha, iterationNum, option, lambdaMinus, lambdaPlus, ratio, uncertainBound)
% Input : CShape. The infomation about DVC, InfluenceMatrix and so on.
% alpha : if alpha<1, consider the voxel at the ratio alpha.
%         if alpha>1, consider the alpha voxel.
%         if alpha=1, consider the all voxel.
% iterationNum : the number of iteration fot LP. In paper, k.
% option : get output opt, depend on option.
%          'all'   : return all output opt
%          'onlyX' : return x, that is the intence of beamlet
%          'onlyZ' : return z, that is the irradiated dose
%          'xzt'   : [x; z; t], whose t is the deviation out of DVC. in paper t
% lambdaMinus : the weight vector of lower DVC
% lambdaPlus  : the upper vector of lower DVC

% start timer
start_total = tic; 

% deviation in each iteration
t_history = [];

% sum time until each iteration
time_history = [];

EPS = 1e-2;
INF = 120;

% ?ｽp?ｽ?ｽ?ｽﾄゑｿｽ?ｽﾈゑｿｽ
boundAdjust = 0;

% beamletNum: number of beamlets, strNum: number of PTV and OAR
[~, beamletNum] = size(Input{1}.mat);
[~, strNum] = size(Input);

global init;

% decide how many we use depende on alpha
% set Influence Matrix with chooseMatrix
% Input{*}.val ?ｽﾍ用?ｽ?ｽ?ｽﾄゑｿｽ?ｽﾈゑｿｽ?ｽB
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

% add noise to Influence Matrix
[Input, uncertainBound] = mat2uncertain(Input, 'Normal', ratio, uncertainBound);

% show information of PTV/OAR in details
for str = 1:strNum 
   fprintf('===Structure{%d}===', str);
   Input{str}
end

% from here
% set coefficient matrix for substituting in cplex
Aeq = [];
voxelNum = 0;
exvariable = 0;


for str = 1:strNum
    Aeq = [Aeq; Input{str}.matMinus];
    Aeq = [Aeq; Input{str}.matPlus];
    voxelNum = voxelNum + 2 * Input{str}.size;
    % exvariable : 
    % (\zeta - z)^+...voxelNum
    % (\theta - z)^+...voxelNum
    % \zeta...1
    % there are upper and lower constraints
    % ?ｽC?ｽ?ｽ?ｽ?ｽ?ｽﾄみゑｿｽ
    exvariable = exvariable + (2*Input{str}.size + 1) * (max(size(Input{str}.lp))+max(size(Input{str}.up))); 
end

% this is minimize variable. -> ?ｽ?ｽ?ｽ?ｽ?ｽ轤ｭ?ｽ?ｽ?ｽ?ｽﾈゑｿｽ
% exvariable = exvariable + 1;

% ?ｽ?ｽ?ｽ?ｽ
% [D+;D-]x = z <=> [D+;D-, -I][x , z]^T = 0 <=> Aeq * (variable) = beq
Aeq = [Aeq -speye(voxelNum)]; 
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
    % ?ｽ?ｽ?ｽ?ｽ?ｽ轤ｭstrVoxelNum?ｽﾍゑｿｽ?ｽ?ｽﾈゑｿｽ
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
% ?ｽW?ｽ?ｽR?ｽﾉ含まゑｿｽ?ｽvoxel?ｽﾌ個撰ｿｽ
CVaRIndex = [];

% ?ｽo?ｽﾍゑｿｽopt?ｽﾌ鯉ｿｽ?ｽ^
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
    
    % z ?ｽﾌ包ｿｽ?ｽ?ｽ?ｽ?ｽ?ｽC?ｽ?ｽ?ｽK?ｽv? -> ?ｽC?ｽ?ｽ?ｽ?ｽ?ｽﾈゑｿｽ?ｽﾄゑｿｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ
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
    
    % x0 is starting point (->?ｽ?ｽ?ｽ?ｽ_?ｽﾅゑｿｽ?ｽ?ｽﾈに趣ｿｽ?ｽﾔゑｿｽ?ｽﾏゑｿｽ?ｽﾈゑｿｽ)
    x0 = [finalx(1:(beamletNum + voxelNum+1)) ;zeros(exvariable, 1)  ];
    prevZetaIndex = 1;
    
    % ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ
    for str=1:strNum
        Input{str}.prevDoseMinus = Input{str}.matMinus * finalx(1:beamletNum);
        Input{str}.prevDosePlus = Input{str}.matPlus * finalx(1:beamletNum);
        % DVC with no ratio alpha
        % l - Pt <= z
        for ind = 1:max(size(Input{str}.lp))
            if Input{str}.lp(ind) > 1-(1e-3) 
               ld = Input{str}.ld(ind);
               for voxel=1:Input{str}.size
                   % lower bound on variable
                   % each bound is for z- z+, respectively
                   lborigin(fixIndex+voxel) = max(lborigin(fixIndex+voxel), ld); % ?ｽ?ｽ?ｽ?ｽ?ｽC?ｽ?ｽ?ｽ?ｽ?ｽ?ｽB?ｽ?ｽﾂゑｿｽ?ｽ?ｽmax?ｽﾆゑｿｽ?ｽﾄゑｿｽ?ｽ?ｽB
                   lborigin(fixIndex+Input{str}.size+voxel) = max(lborigin(fixIndex+Input{str}.size+voxel), ld);
                   
                   % -(z-) - Pt <= -L
                   spi(curPos:curPos+1) = [curConst curConst]; 
                   spj(curPos:curPos+1) = [fixIndex+voxel paramIndex];
                   sps(curPos:curPos+1) = [-1 -1]; % P?ｽ?ｽ?ｽl?ｽ?ｽ?ｽ?ｽ?ｽ?ｽﾈゑｿｽﾎゑｿｽ?ｽ?ｽ?ｽ?ｽﾏ更
                   bcon(curConst) = -ld;
                   
                   % upload the index of sparse matrix
                   curConst = curConst + 1;
                   curPos = curPos + 2;
               end
               
               % ?ｽ?ｽ?ｽ?ｽ?ｽﾉ新?ｽ?ｽ?ｽ?ｽpenalty?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ -> wIndex?ｽﾌ撰ｿｽ?ｽ?ｽ?ｽl?ｽ?ｽ?ｽ?ｽ?ｽ?ｽK?ｽv?ｽ?ｽ?ｽ?ｽH
               % penalty terms ("nu" is the artificial variable)
               strsize = Input{str}.size;
               for voxel = 1:Input{str}.size
                   % -(z-) - nu <= - L
                   spi(curPos:curPos+1) = [curConst, curConst];
                   spj(curPos:curPos+1) = [fixIndex+voxel, wIndex];
                   sps(curPos:curPos+1) = [-1, -1];
                   bcon(curConst) = -ld;
                   % initial point
                   x0(wIndex) = 0;
                   % coefficient of objective function
                   c(wIndex) = lambdaMinus(lambdaMinusIndex) / strsize;
                   
                   % upload the index
                   wIndex = wIndex + 1;
                   curConst = curConst + 1;
                   curPos = curPos + 2;
               end
               lambdaMinusIndex = lambdaMinusIndex + 1;
               continue; 
            end
            
            % DVC with ratio alpha
            % count the number of the voxel in set R
            useCnt = 0;
            unuseIndex = zeros(Input{str}.size, 1);
            for voxel = 1:Input{str}.size
                checker = Input{str}.lt(ind) - EPS;
                % if z+ < L - Pt, then this coresponding index turns from 0 to 1
                % z- ?ｽﾅ撰ｿｽ?ｽ?ｽ?ｽ?ｽ?ｽﾌゑｿｽ?ｽH?ｽv?ｽ`?ｽF?ｽb?ｽN
                if Input{str}.prevDosePlus(voxel) < checker
                    unuseIndex(voxel) = 1;
                else
                    useCnt = useCnt+1;
                end
            end
            
            % eta is used when zeta is needed in inequarity
            % ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽeta?ｽﾍ論?ｽ?ｽ?ｽ?ｽzeta
            etaIndex = wIndex + useCnt ;
            etaIndexes = [etaIndexes; etaIndex];
            initIndex = wIndex-1;
            for voxel = 1:Input{str}.size
                % check voxel is included in the set R
                if unuseIndex(voxel)==1 
                    continue;
                end
                %  - z(fixIndex+voxel) - w(wIndex) + zeta(etaIndex) <= 0(bcon)
                spi(curPos:curPos+2) = [curConst curConst curConst];
                spj(curPos:curPos+2) = [fixIndex+voxel wIndex etaIndex];
                sps(curPos:curPos+2) = [-1 -1 1];
                bcon(curConst) = 0;
                
                % initial point
                x0(wIndex) = max(0, prevZeta(prevZetaIndex) - finalx(fixIndex+voxel));
                
                % upload index
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 3;
            end
            % from here
            % set the DVC with ratio alpha,\
            unuseCnt = Input{str}.size - useCnt;
            % remainCnt = (1 - alpha)*v_s - |R|
            remainCnt = (1-Input{str}.lp(ind)) * Input{str}.size - unuseCnt;
            if remainCnt <=0
               exitflag = 0; 
            end
            coef = 1/remainCnt;
            
            cntDiff = 0; 
            for voxel=1:Input{str}.size
                % if unuseIndex is true then skip this iteration
                % else, then set the sparse matrix about (zeta - z)^+
                % and then use cntDiff in order that exvariables are used successfully
                if unuseIndex(voxel)==1
                    cntDiff = cntDiff + 1;
                    continue;
                end
                % coef * w /// w = (zeta - z)^+
                % express "zeta - coef*w >= L - Pt"
                spi(curPos) = curConst;
                spj(curPos) = initIndex + voxel - cntDiff;
                sps(curPos) = coef;
                curPos = curPos + 1;
            end
            % zeta
            spi(curPos) = curConst;
            spj(curPos) = etaIndex;
            sps(curPos) = -1;
            % initial point
            x0(wIndex) = Input{str}.lt(ind);
            % upload index
            curPos = curPos + 1;
            wIndex = wIndex + 1;
            % t
            spi(curPos) = curConst;
            spj(curPos) = paramIndex;
            sps(curPos) = -Input{str}.lparam(ind);
            % constant
            ld = Input{str}.ld(ind) + boundAdjust;
            bcon(curConst) = -ld;
            CVaRIndex = [CVaRIndex; unuseCnt];
            % upload index
            curPos = curPos + 1;
            curConst = curConst + 1;
            % until here, about DVC with ratio alpha
            
            % coef?ｽﾌ擾ｿｽ?ｽﾍゑｿｽ?ｽﾜゑｿｽK?ｽv?ｽﾈゑｿｽ?ｽH|R|?ｽ?ｽ?ｽ?ｽ?ｽﾅゑｿｽ?ｽ?ｽ?ｽﾌでは？
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
                
                x0(wIndex) = 0; % ?ｽ?ｽ?ｽ?ｽ_?ｽ?ｽ?ｽ?ｽ驍ｩ?ｽﾍ考?ｽ?ｽ?ｽ?ｽ
                % set penalty to objective function c(wIndex)
                c(wIndex) = lambdaMinus(lambdaMinusIndex) / strsize;
                % upload index
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 2;
            end
            
            if 0
                info = [unuseCnt remainCnt coef bcon(curConst-1)]
            end
            % upload index of zeta and lambda
            prevZetaIndex = prevZetaIndex + 1;
            lambdaMinusIndex = lambdaMinusIndex + 1;
        end
        
        % from here, Upper DVC constraints
        % be careful that variable is [x, z-, z+, t, exvariable]
        for ind = 1:max(size(Input{str}.up))
            % DVC with no ratio alpha
            % z+ <= U + Pt
            strVoxelNum = Input{str}.size; % for using z+
            if Input{str}.up(ind) < 1e-3
                ud = Input{str}.ud(ind);
                for voxel=1:Input{str}.size
                   uborigin(fixIndex+strVoxelNum+voxel) = min(uborigin(fixIndex+strVoxelNum+voxel), ud);
                   
                   spi(curPos:curPos+1) = [curConst curConst];
                   spj(curPos:curPos+1) = [fixIndex+strVoxelNum+voxel paramIndex];
                   sps(curPos:curPos+1) = [1 -1];
                   bcon(curConst) = ud;
                   % upload index
                   curConst = curConst + 1;
                   curPos = curPos + 2; 
                end
                % ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽpenalty?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ?ｽ
                % penalty terms ("nu" is the artificial variable)
                % z+ - nu <= U
                strsize = Input{str}.size;
                for voxel = 1:Input{str}.size
                    spi(curPos:curPos+1) = [curConst, curConst];
                    spj(curPos:curPos+1) = [fixIndex+strVoxelNum+voxel, wIndex];
                    sps(curPos:curPos+1) = [1, -1];
                    bcon(curConst) = ud;
                    
                    % initial point
                    x0(wIndex) = 0;
                    
                    % coefficient of the objective function
                    c(wIndex) = lambdaPlus(lambdaPlusIndex) / strsize;
                    
                    % upload index
                    wIndex = wIndex + 1;
                    curConst = curConst + 1;
                    curPos = curPos + 2;
                end
                lambdaPlusIndex = lambdaPlusIndex + 1;
                continue; 
            end
            
            % DVC with ratio alpha
            useCnt = 0;
            unuseIndex = zeros(Input{str}.size, 1);
            % check the voxel is included in the set R
            for voxel = 1:Input{str}.size
               checker = Input{str}.ut(ind) + EPS;
               if Input{str}.prevDoseMinus(voxel) > checker
                  unuseIndex(voxel) = 1;
               else
                   useCnt = useCnt + 1;
               end
            end
            
            % this eta is zeta in paper
            etaIndex = wIndex + useCnt;
            etaIndexes = [etaIndexes; etaIndex];
            initIndex = wIndex-1;
            for voxel = 1:Input{str}.size
                if unuseIndex(voxel) == 1
                    continue; 
                end
                % (z+ - zeta)^+ -> z+ - nu - zeta <= 0
                spi(curPos:curPos+2) = [curConst curConst curConst];
                spj(curPos:curPos+2) = [fixIndex+strVoxelNum+voxel wIndex etaIndex];
                sps(curPos:curPos+2) = [1 -1 -1];
                bcon(curConst) = 0;
                
                % initial point
                x0(wIndex) = max(0, finalx(fixIndex+strVoxelNum+voxel) - prevZeta(prevZetaIndex));
                
                % upload index
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 3;
            end
            % remainCnt is (alpha * v - |R|)
            unuseCnt = Input{str}.size - useCnt;
            remainCnt = Input{str}.up(ind) * Input{str}.size - unuseCnt;
            if remainCnt <= 0 
                exitflag = 0;
            end
            coef = 1/remainCnt;
            % if unuseIndex is true then skip this iteration
            % else, then set the sparse matrix about (zeta - z)^+
            % and then use cntDiff in order that exvariables are used successfully
            cntDiff = 0;
            for voxel = 1:Input{str}.size
                if unuseIndex(voxel) == 1
                   cntDiff = cntDiff + 1;
                   continue; 
                end
                % zeta + coef * w <= U + Pt
                % w
                spi(curPos) = curConst;
                spj(curPos) = initIndex + voxel - cntDiff;
                sps(curPos) = coef;
                curPos = curPos + 1;
            end
            % zeta
            spi(curPos) = curConst;
            spj(curPos) = etaIndex;
            sps(curPos) = 1;
            
            % initial point
            x0(wIndex) = Input{str}.ut(ind);
            
            % upload index
            curPos = curPos + 1;
            wIndex = wIndex + 1;
            
            % t
            spi(curPos) = curConst;
            spj(curPos) = paramIndex;
            sps(curPos) = -Input{str}.uparam(ind);
            
            % constant
            ud = Input{str}.ud(ind) - boundAdjust;
            bcon(curConst) = ud;
            
            CVaRIndex = [CVaRIndex; unuseCnt];
            
            curPos = curPos + 1;
            curConst = curConst + 1;
            
            % coef?ｽﾌ擾ｿｽ?ｽﾍゑｿｽ?ｽ?ｽﾈゑｿｽ?ｽH
            coefValues = [coefValues; coef];
            
            % from here, penalty term
            strsize = useCnt;
            for voxel = 1:Input{str}.size
                if unuseIndex(voxel) == 1
                    continue;
                end
                % z+ - nu <= U
                spi(curPos:curPos+1) = [curConst, curConst];
                spj(curPos:curPos+1) = [fixIndex+strVoxelNum+voxel, wIndex];
                sps(curPos:curPos+1) = [1, -1];
                bcon(curConst) = ud;
                
                % initial point
                x0(wIndex) = 0; % ?ｽﾏ更?ｽ?ｽ?ｽ?ｽﾆ托ｿｽ?ｽ?ｽ?ｽﾈゑｿｽH
                
                % coefficient of the objective function
                c(wIndex) = lambdaPlus(lambdaPlusIndex) / strsize;
                
                % upload index
                wIndex = wIndex + 1;
                curConst = curConst + 1;
                curPos = curPos + 2;
            end
            
            info = [unuseCnt remainCnt coef bcon(curConst-1)];
            
            prevZetaIndex = prevZetaIndex + 1;
            lambdaPlusIndex = lambdaPlusIndex + 1;
        end
        fixIndex = fixIndex + 2*Input{str}.size; % z+ z-
    end
    % from here, construct the form to use CPLEX
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
    %opt = finalx((beamletNum+1):(beamletNum+voxelNum)); 
    for str = 1:strNum
        A0 = [];
        A0 = [A0; Input{str}.mat0];
    end
    opt = A0 * finalx(1:beamletNum);
elseif strcmp(option, 'xzt')
    % opt = finalx(1:(beamletNum + voxelNum + 1));
    A0 = [];
    for str = 1:strNum
        A0 = [A0; Input{str}.mat0];
    end
    opt_z = A0 * finalx(1:beamletNum);
    opt = [finalx(1:beamletNum); opt_z; finalx(beamlet+2*voxelNum+1)];
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

%CVaR?ｽﾌ表?ｽ?ｽ
fprintf('== CVaRIndex == ');
fprintf('%.2f ', CVaRIndex);
fprintf('\n');

% ?ｽﾏ更?ｽ_ lambda?ｽﾌ表?ｽ?ｽ
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



