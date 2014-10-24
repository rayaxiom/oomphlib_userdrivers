p=input('Input the prder of approximation (1/2): ');

% This switch is horrible - change this.
switch(p)
    case 1 
       fname=('jac_np1r0');
%       fname=('Jac_l');
    case 2 
       fname=('Jac_q'); 
    otherwise
       return;
end


X=load(fname);
X(:,1)=X(:,1)+1; X(:,2)=X(:,2)+1;
A=spconvert(X);
n=max(X(:,1));
fprintf('Problem size: %7i\n',n);
nnzA=nnz(A);
fprintf('Number of non-zero elements: %9i\n',nnzA);
sts=nnzA/n;
fprintf('Average stencil size: %6.2f\n',sts);
clear X;


tst=tic; % RAY To start time
Pamd=amd(A); % A is the jacobian
Aamd=A(Pamd,Pamd);
texe=toc(tst); % RAY END TIME
fprintf('Approximate minimum degree reordering time %6.2f seconds\n',texe);


tst=tic; % RAY START NEW TIME
Prcm=symrcm(A);% RAY permutation vector for Approximate reverse Cuthill-McKee 
Arcm=A(Prcm,Prcm); 
texe=toc(tst); % RAY END TIME
fprintf('Approximate reverse Cuthill-McKee reordering time %6.2f seconds\n',texe);


ff=figure(1);
set(ff,'Position',[0 420 1270 505]);
subplot(1,3,1); spy(A); title('Default ordering');
subplot(1,3,2); spy(Aamd); title ('Approximate minimum degree ordering');
subplot(1,3,3); spy(Arcm); title ('Reverese Cuthill-McKee ordering');

%%
[colA,rowA,valA]=find(A);
[colAamd,rowAamd,valAamd]=find(Aamd); 
[colArcm,rowArcm,valArcm]=find(Arcm); 

%
fprintf('DEFAULT ORDERING:\n');
bwA=0;
for i=1:nnzA
   bwA=max(bwA,abs(rowA(i)-colA(i))); 
end
fprintf('   max(i-j)=%7i\n',bwA);
fid=fopen('J0.txt','w');
fprintf(fid,'%7i\n',n);
fprintf(fid,'%9i\n',nnzA);
for i=1:nnzA
fprintf(fid,'%7i %7i %12.6f\n',rowA(i)-1,colA(i)-1,valA(i));
end
fclose(fid);

%
fprintf('APPROXIMATE MINIMUM DEGREE ORDERING:\n');
bwAamd=0;
for i=1:nnzA
   bwAamd=max(bwAamd,abs(rowAamd(i)-colAamd(i))); 
end
fprintf('   max(i-j)=%7i\n',bwAamd);
fid=fopen('J1.txt','w');
fprintf(fid,'%7i\n',n);
fprintf(fid,'%9i\n',nnzA);
for i=1:nnzA
   fprintf(fid,'%7i %7i %12.6f\n',rowAamd(i)-1,colAamd(i)-1,valAamd(i));
end
fclose(fid);

%
fprintf('REVERSE CUTHILL-MCKEE ORDERING:\n');
bwArcm=0;
for i=1:nnzA
   bwArcm=max(bwArcm,abs(rowArcm(i)-colArcm(i))); 
end
fprintf('   max(i-j)=%7i\n',bwArcm);
fid=fopen('J2.txt','w');
fprintf(fid,'%7i\n',n);
fprintf(fid,'%9i\n',nnzA);
for i=1:nnzA
   fprintf(fid,'%7i %7i %12.6f\n',rowArcm(i)-1,colArcm(i)-1,valArcm(i));
end
fclose(fid);


rowptrArcm(1:n+1)=0; rowptrArcm(1)=1; k=2; %% RAY rowptrArcm - row pointer array
for i=1:nnzA-1
   if(rowArcm(i+1)~=rowArcm(i))
      rowptrArcm(k)=i+1;
      k=k+1;
   end
end 
rowptrArcm(k)=nnzA+1;
if(k~=n+1)
    fprintf('Error in converting Arcm to CRS format\n');
end
%% n - the size of the matrix
%% nnzA - the number of non-zeros
%% rowptrArcm - row pointer array
%% colArcm - column array
%% valArcm - the elements
