function [matret, dret] = chooseMatrix(mat, d, alpha)
% It makes rmatrix that has only probability alpha.
% Rows is chosen randomly.

sz = size(mat);
width = sz(2);  
sz = sz(1); 


newsz = ceil(sz * alpha); 

vec = zeros(newsz,1); 
for i = 1:newsz
    vec(i) = i; %(1,2,3,...,newsz)
end

vec = [vec; zeros(sz-newsz, 1)];

for i = 1:sz
    j = randi(i);
    tmp = vec(i);
    vec(i) = vec(j);
    vec(j) = tmp;
end 

ref = zeros(newsz, 1); 
for i = 1:sz
   if vec(i)>=1
       ref(vec(i)) = i;
   end
end 
matret = []; dret = [];
tic 

%{
mat = transpose(mat);
d   = transpose(d); 
%for i = 1:newsz
    matret = mat(:, ref(1:newsz)); 
    dret   = d(:, ref(1:newsz));
%end

matret = transpose(matret);
dret   = transpose(dret);
%}
matret = mat(ref(1:newsz),:); 
dret   = d(ref(1:newsz),:);

toc 
ended = 'endChooseMatrixPhase.'



