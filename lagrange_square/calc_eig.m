clear all; close all;
%cd A30V0R0N4;
%cd A30V0R0N4
load j_0000;
load j_0001;
load j_0002;
load j_0003;
load j_0004;
load j_0005;

load j_0100;
load j_0101;
load j_0102;
load j_0103;
load j_0104;
load j_0105;

load j_0200;
load j_0201;
load j_0202;
load j_0203;
load j_0204;
load j_0205;

load j_0300;
load j_0301;
load j_0302;
load j_0303;
load j_0304;
load j_0305;

load j_0400;
load j_0401;
load j_0402;
load j_0403;
load j_0404;
load j_0405;


load j_0500;
load j_0501;
load j_0502;
load j_0503;
load j_0504;
load j_0505;

temp=j_0000;
j00 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0001;
j01 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0002;
j02 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0003;
j03 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0004;
j04 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0005;
j05 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);

temp=j_0100;
j10 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0101;
j11 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0102;
j12 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0103;
j13 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0104;
j14 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0105;
j15 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);

temp=j_0200;
j20 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0201;
j21 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0202;
j22 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0203;
j23 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0204;
j24 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0205;
j25 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);



temp=j_0300;
j30 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0301;
j31 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0302;
j32 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0303;
j33 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0304;
j34 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0305;
j35 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);


temp=j_0400;
j40 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0401;
j41 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0402;
j42 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0403;
j43 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0404;
j44 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0405;
j45 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);


temp=j_0500;
j50 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0501;
j51 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0502;
j52 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0503;
j53 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0504;
j54 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);
temp=j_0505;
j55 = spconvert([temp(1:end,1)+1,temp(1:end,2)+1,temp(1:end,3)]);

clear load j_0000 j_0001 j_0002 0003 j_0004 j_0005 j_0100 j_0101 j_0102...
j_0103 j_0104 j_0105 j_0200 j_0201 j_0202 j_0203 j_0204 j_0205 j_0300...
j_0301 j_0302 j_0303 j_0304 j_0305 j_0400 j_0401 j_0402 j_0403 j_0404...
j_0405 j_0500 j_0501 j_0502 j_0503 j_0504 j_0505;

%     b   b   c   c   p   L
J = [j00 j01 j02 j03 j04 j05
     j10 j11 j12 j13 j14 j15
     j20 j21 j22 j23 j24 j25
     j30 j31 j32 j33 j34 j35
     j40 j41 j42 j43 j44 j45
     j50 j51 j52 j53 j54 j55];
 
Mx = j52;
My = j53;
Mxt = j25;
Myt = j35;
z05 = zeros(size(j05,1),size(j05,2));
z15 = zeros(size(j15,1),size(j15,2));
z25 = zeros(size(j25,1),size(j25,2));
z35 = zeros(size(j35,1),size(j35,2));
z45 = zeros(size(j45,1),size(j45,2));

z50 = zeros(size(j50,1),size(j50,2));
z51 = zeros(size(j51,1),size(j51,2));
z52 = zeros(size(j52,1),size(j52,2));
z53 = zeros(size(j53,1),size(j53,2));
z54 = zeros(size(j54,1),size(j54,2));
z55 = zeros(size(j55,1),size(j55,2));

NS=[j00 j01 j02 j03 
     j10 j11 j12 j13 
     j20 j21 j22 j23 
     j30 j31 j32 j33];
raysigma = -norm(NS,Inf);


W = (Mx*Mxt + My*Myt)*raysigma;
invW = inv(W);
c22 = j22 + Mxt*invW* Mx;
c23 = j23 + Mxt*invW* My;
c32 = j32 + Myt*invW* Mx;
c33 = j33 + Myt*invW* My;
P = [j00 j01 j02 j03 j04 z05
     j10 j11 j12 j13 j14 z15
     j20 j21 c22 c23 j24 z25
     j30 j31 c32 c33 j34 z35
     j40 j41 j42 j43 j44 z45
     z50 z51 z52 z53 z54 W];
 
lambda = eig(full(J),full(P));
 
 save 'eigenoutput.mat';





