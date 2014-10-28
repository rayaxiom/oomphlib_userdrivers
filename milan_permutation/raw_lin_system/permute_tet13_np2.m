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


[colArcm,~,valArcm]=find(Arcm);
nnzA = nnz(Arcm);

disp('Outputting val');
fid=fopen('jac13_np2r0_val','w+'); % CCC
for i=1:nnzA
   fprintf(fid,'%12.15f\n',valArcm(i));
end
fclose(fid);

disp('Outputting col');
fid=fopen('jac13_np2r0_col','w+'); % CCC
for i=1:nnzA
   fprintf(fid,'%7i\n',colArcm(i)-1);
end
fclose(fid);

disp('Outputting rowpt');
fid=fopen('jac13_np2r0_row','w+'); % CCC
fprintf(fid,'%7i\n',0);
nnz_per_row = sum(Arcm~=0,2);
nnz_running = 0;
for i=1:nrow_local
   nnz_running = nnz_running + nnz_per_row(i);
   fprintf(fid,'%7i\n',nnz_running);
end
fclose(fid);


%%% Now output the val, col and row_pt for r0
Arcm = Arcm1; % Set this
nrow_local = nrow_local1; % Set this

[colArcm,~,valArcm]=find(Arcm);
nnzA = nnz(Arcm);


disp('Outputting val');
fid=fopen('jac13_np2r1_val','w+'); % CCC
for i=1:nnzA
   fprintf(fid,'%12.15f\n',valArcm(i));
end
fclose(fid);

disp('Outputting col');
fid=fopen('jac13_np2r1_col','w+'); % CCC
for i=1:nnzA
   fprintf(fid,'%7i\n',colArcm(i)-1);
end
fclose(fid);

disp('Outputting rowpt');
fid=fopen('jac13_np2r1_row','w+'); % CCC
fprintf(fid,'%7i\n',0);
nnz_per_row = sum(Arcm~=0,2);
nnz_running = 0;
for i=1:nrow_local
   nnz_running = nnz_running + nnz_per_row(i);
   fprintf(fid,'%7i\n',nnz_running);
end
fclose(fid);



%texe=toc(tst); % RAY END TIME
%fprintf('Approximate reverse Cuthill-McKee reordering time %6.2f seconds\n',texe);


%ff=figure(1);
%set(ff,'Position',[0 420 1270 505]);
%subplot(1,3,1); spy(A); title('Default ordering');
%subplot(1,3,2); spy(Aamd); title ('Approximate minimum degree ordering');
%subplot(1,3,3); spy(Arcm); title ('Reverese Cuthill-McKee ordering');

%%
%[colAmat,rowAmat,valAmat]=find(Amat);
%[colAamd,rowAamd,valAamd]=find(Aamd); 
%%%%[colArcm,rowArcm,valArcm]=find(Arcm);

%
%fprintf('DEFAULT ORDERING:\n');
%bwA=0;
%for i=1:nnzA
%   bwA=max(bwA,abs(rowA(i)-colA(i))); 
%end
%fprintf('   max(i-j)=%7i\n',bwA);
%fid=fopen('J0.txt','w');
%fprintf(fid,'%7i\n',n);
%fprintf(fid,'%9i\n',nnzA);
%for i=1:nnzA
%fprintf(fid,'%7i %7i %12.6f\n',rowA(i)-1,colA(i)-1,valA(i));
%end
%fclose(fid);

%
%fprintf('APPROXIMATE MINIMUM DEGREE ORDERING:\n');
%bwAamd=0;
%for i=1:nnzA
%   bwAamd=max(bwAamd,abs(rowAamd(i)-colAamd(i))); 
%end
%fprintf('   max(i-j)=%7i\n',bwAamd);
%fid=fopen('J1.txt','w');
%fprintf(fid,'%7i\n',n);
%fprintf(fid,'%9i\n',nnzA);
%for i=1:nnzA
%   fprintf(fid,'%7i %7i %12.6f\n',rowAamd(i)-1,colAamd(i)-1,valAamd(i));
%end
%fclose(fid);

%
%fprintf('REVERSE CUTHILL-MCKEE ORDERING:\n');
%bwArcm=0;
%for i=1:nnzA
%   bwArcm=max(bwArcm,abs(rowArcm(i)-colArcm(i))); 
%end
%fprintf('   max(i-j)=%7i\n',bwArcm);

%%%%disp('Outputting val');
%%%%fid=fopen('jac12_np1r0_val','w+');
%%%%for i=1:nnzA
%%%%   fprintf(fid,'%12.15f\n',valArcm(i));
%%%%end
%%%%fclose(fid);

%%%%disp('Outputting col');
%%%%fid=fopen('jac12_np1r0_col','w+');
%%%%for i=1:nnzA
%%%%   fprintf(fid,'%7i\n',colArcm(i)-1);
%%%%end
%%%%fclose(fid);

%%%%disp('Outputting rowpt');
%%%%fid=fopen('jac12_np1r0_row','w+');
%%%%fprintf(fid,'%7i\n',0);
%%%%nnz_per_row = sum(Amat~=0,2);
%%%%nnz_running = 0;
%%%%for i=1:nrow
%%%%   nnz_running = nnz_running + nnz_per_row(i);
%%%%   fprintf(fid,'%7i\n',nnz_running);
%%%%end
%%%%fclose(fid);



%rowptrArcm(1:n+1)=0; rowptrArcm(1)=1; k=2; %% RAY rowptrArcm - row pointer array
%for i=1:nnzA-1
%   if(rowArcm(i+1)~=rowArcm(i))
%      rowptrArcm(k)=i+1;
%      k=k+1;
%   end
%end 
%rowptrArcm(k)=nnzA+1;
%if(k~=n+1)
%    fprintf('Error in converting Arcm to CRS format\n');
%end

%% n - the size of the matrix
%% nnzA - the number of non-zeros
%% rowptrArcm - row pointer array
%% colArcm - column array
%% valArcm - the elements
%}