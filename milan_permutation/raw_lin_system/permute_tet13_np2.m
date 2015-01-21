clear all;
%p=input('Input the prder of approximation (1/2): ');

%%% First we try this out on jac1_np2

% The files are jac1_np2r0 and jac1_np2r1

% First load r0
fname=('jac13_np2r0');
rawmat=load(fname);
rawmat(:,1)=rawmat(:,1)+1; rawmat(:,2)=rawmat(:,2)+1;
Amat0=spconvert(rawmat); % put it in Amatr0
clear rawmat;

% Then  load r1
fname=('jac13_np2r1');
rawmat=load(fname);
rawmat(:,1)=rawmat(:,1)+1; rawmat(:,2)=rawmat(:,2)+1;
Amat1=spconvert(rawmat); % put it in Amatr1
clear rawmat;

% Concatente these matrices into Amat, but remember the nrow_local to
% separate them.
nrow_local0=size(Amat0,1);
nrow_local1=size(Amat1,1);

Amat=[Amat0;
      Amat1];

clear Amat0 Amat1;


%%% Continue with the stuff as normal.
nrow=size(Amat,2); % 
fprintf('Problem size: %7i\n',nrow);
nnzA=nnz(Amat);
fprintf('Number of non-zero elements: %9i\n',nnzA);
sts=nnzA/nrow;
fprintf('Average stencil size: %6.2f\n',sts);

%tst=tic; % RAY To start time
%Pamd=amd(Amat); % A is the jacobian
%Aamd=Amat(Pamd,Pamd);
%texe=toc(tst); % RAY END TIME
%fprintf('Approximate minimum degree reordering time %6.2f seconds\n',texe);

%tst=tic; % RAY START NEW TIME
Prcm=symrcm(Amat);% RAY permutation vector for Approximate reverse Cuthill-McKee 
Arcm=Amat(Prcm,Prcm); % RAY perform the permutation
clear Amat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that we have the permutation vector, we: 
% 1) load all the residuals
% 2) Concatenate them
% 3) Permute them
% 4) Split them according to the nrow_local
% 5) Output them.

% 1) Only residual1_np2r0 and residual1_np2r1
fname=('residual13_np2r0'); % Change this!
rawres=load(fname);
res0=rawres(:,2); % Change this!
clear rawres;

fname=('residual13_np2r1'); % Change this!
rawres=load(fname);
res1=rawres(:,2); % Change this!
clear rawres;

% 2.1) concat res0 and res1 into res
% 2.2) clear res0 and res1
res = [res0;res1];
clear res0 res1;

% 3) Permute the residual
res=res(Prcm);


%%% Recall that we need to split this up.
cut_min = 1;
cut_max = nrow_local0;
res0=res(cut_min:cut_max); % change this

cut_min = cut_max+1;
cut_max = cut_max + nrow_local1;
res1=res(cut_min:cut_max); % change this

% clear res
clear res;

disp('Outputting residual');
fid=fopen('res13_np2r0','w+');
for i=1:nrow_local0
   fprintf(fid,'%12.15f\n',res0(i));
end
fclose(fid);
clear res0

fid=fopen('res13_np2r1','w+');
for i=1:nrow_local1
   fprintf(fid,'%12.15f\n',res1(i));
end
fclose(fid);
clear res1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Now we output the jacobian

% First we split the jacobian.
cut_min = 1;
cut_max = nrow_local0;
Arcm0 = Arcm(cut_min:cut_max,1:nrow);

cut_min = cut_max+1;
cut_max = cut_max + nrow_local1;
Arcm1 = Arcm(cut_min:cut_max,1:nrow);

clear Arcm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Now output the val, col and row_pt for r0
Arcm = Arcm0; % Set this
nrow_local = nrow_local0; % Set this
fname='jac13_np2r0';

% First get the column with the transpose matrix
% Note that this is actually the row of the transpose
% but MATLAB traverses matrices
%[rowAA,colAA,valAA]=find(AA');
[colArcm,~,valArcm]=find(Arcm');

colArcm = colArcm - 1;

nnzArcm = nnz(Arcm);

disp('Outputting val');
fnameval=strcat(fname,'_val');
fid=fopen(fnameval,'w+');
for i=1:nnzArcm
   fprintf(fid,'%12.15f\n',valArcm(i));
end
fclose(fid);

disp('Outputting col');
fnamecol=strcat(fname,'_col');
fid=fopen(fnamecol,'w+');
for i=1:nnzArcm
   fprintf(fid,'%7i\n',colArcm(i));
end
fclose(fid);

disp('Outputting rowpt');
fnamerow=strcat(fname,'_row');
fid=fopen(fnamerow,'w+');
fprintf(fid,'%7i\n',0);

% Get the row_pt with the normal matrix.
nnz_per_row = sum(Arcm~=0,2);
nnz_running = 0;
for i=1:nrow_local
   nnz_running = nnz_running + nnz_per_row(i);
   fprintf(fid,'%7i\n',nnz_running);
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Now output the val, col and row_pt for r1
%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Now output the val, col and row_pt for r0
Arcm = Arcm1; % Set this
nrow_local = nrow_local1; % Set this
fname='jac13_np2r1';

% First get the column with the transpose matrix
% Note that this is actually the row of the transpose
% but MATLAB traverses matrices
%[rowAA,colAA,valAA]=find(AA');
[colArcm,~,valArcm]=find(Arcm');

colArcm = colArcm - 1;

nnzArcm = nnz(Arcm);

disp('Outputting val');
fnameval=strcat(fname,'_val');
fid=fopen(fnameval,'w+');
for i=1:nnzArcm
   fprintf(fid,'%12.15f\n',valArcm(i));
end
fclose(fid);

disp('Outputting col');
fnamecol=strcat(fname,'_col');
fid=fopen(fnamecol,'w+');
for i=1:nnzArcm
   fprintf(fid,'%7i\n',colArcm(i));
end
fclose(fid);

disp('Outputting rowpt');
fnamerow=strcat(fname,'_row');
fid=fopen(fnamerow,'w+');
fprintf(fid,'%7i\n',0);

% Get the row_pt with the normal matrix.
nnz_per_row = sum(Arcm~=0,2);
nnz_running = 0;
for i=1:nrow_local
   nnz_running = nnz_running + nnz_per_row(i);
   fprintf(fid,'%7i\n',nnz_running);
end
fclose(fid);



















