

fname=('AAjac1_np2r1');
rawmat=load(fname);
rawmat(:,1)=rawmat(:,1)+1; rawmat(:,2)=rawmat(:,2)+1;
AAmat1=spconvert(rawmat); % put it in Amatr1

AAmat1full = full(AAmat1);
Arcm1full=full(Arcm1);

diffmat = abs(AAmat1full - Arcm1full);

Maxele=max(max(diffmat));

