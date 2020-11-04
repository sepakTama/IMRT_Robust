% register functions that has different paths
% we need to change paths depend on local machine
% addpath(genpath('/Applications/CERR-master'));
% addpath(genpath('/Applications/CPLEX_Studio128/cplex/matlab/x86-64_osx'));
addpath(genpath('../CERR'));
addpath(genpath('/opt/ibm/ILOG/CPLEX_Studio1262/cplex/matlab/x86-64_linux'));

% /setFunction has the informations of PTVs and OARs, and these functions are used after caluculating IMRTP
addpath('./setFunction')
addpath('./src')
addpath('./method')

% start CERR
CERR;
set(findobj('Tag', 'CERRStartupFig'), 'visible', 'off');
sliceCallBack('init');

fprintf('Cshape      =>  1 (Core),2 (Outer Target)\n');
fprintf('HeadAndNeck =>  4 (PTV), 5 (Cord), 6 (LtParotid), 7(RtParotid)\n');
fprintf('Prostate    => 14 (Prostate PTV), 17 (Rectum), 16 (Bladder)\n');
fprintf('MultiTarget => 10 (Center), 11 (Inferior), 12 (Superior)\n');


