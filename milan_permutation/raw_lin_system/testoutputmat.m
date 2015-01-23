clear all 
AA=[11 22  26   0   0;
    33 44  55  66  0;
    0  77  88  99  0;
    0  100 110 120 130;
    0  0   0   140 150];

% First get the column with the transpose matrix
% Note that this is actually the row of the transpose
% but MATLAB traverses matrices
%[rowAA,colAA,valAA]=find(AA');
[colAA,~,valAA]=find(AA');

colAA = colAA - 1;

nnzAA = nnz(AA);

% Get the row_pt with the normal matrix.
nnz_per_row = sum(AA~=0,2);

