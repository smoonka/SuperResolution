function ratio = alignment(framename, imageb)
camera = textread(framename);
R21=camera(17,1);
R22=camera(18,1);
R23=camera(19,1);
R31=camera(13,1);
R32=camera(14,1);
R33=camera(15,1);
R11 = 0;
R12 = 0;
R13 = 0;
if (R21^2+R31^2<1)
R11=sqrt(1-R21^2-R31^2);
end
if (R22^2+R32^2<1)
R12=sqrt(1-R22^2-R32^2);
end
if (R23^2+R33^2<1)
R13=sqrt(1-R23^2-R33^2);
end
T1=camera(9,1);
T2=camera(10,1);
T3=camera(11,1);
TBB=[R11,R12,R13,T1;R21,R22,R23,T2;R31,R32,R33,T3;0 0 0 1];
TBB = inv(TBB); % extrinsic matrix
fx=camera(33,1)/camera(25,1);
fy=camera(33,1)/camera(27,1);
u0=1080/2;
v0=1920/2;
KBB = [fx,0,u0;0,fy,v0;0,0,1]; % intrinsic matrix
PBB = KBB*[eye(3) zeros(3,1)]*TBB; % projection matrix
camera = textread('clip1-omnicam-calibration.cam');
R21=camera(17,1);
R22=camera(18,1);
R23=camera(19,1);
R31=camera(13,1);
R32=camera(14,1);
R33=camera(15,1);
R11 = 0;
R12 = 0;
R13 = 0;
if (R21^2+R31^2<1)
R11=sqrt(1-R21^2-R31^2);
end
if (R22^2+R32^2<1)
R12=sqrt(1-R22^2-R32^2);
end
if (R23^2+R33^2<1)
R13=sqrt(1-R23^2-R33^2);
end
T1=camera(9,1);
T2=camera(10,1);
T3=camera(11,1);
TO=[R11,R12,R13,T1;R21,R22,R23,T2;R31,R32,R33,T3; 0 0 0 1];
TO = inv(TO); % extrinsic matrix
fx=camera(33,1)/camera(25,1);
fy=camera(33,1)/camera(27,1);
u0=1920/2;
v0=6984/2;
KO = [fx,0,u0;0,fy,v0;0,0,1]; % intrinsic matrix
PO = KO*[eye(3) zeros(3,1)]*TO; % projection matrix
[m,n]=size(PBB);
[U,S,V]=svd(PBB);
r=rank(S);
SR=S(1:r,1:r);
SRc=[SR^-1 zeros(r,m-r);zeros(n-r,r) zeros(n-r,m-r)];
PBBin=V*SRc*U.';
T = PO*PBBin; % Transformation between the frames
imageb=rgb2gray(imageb);
for i = 1:size(imageb,1)
for j = 1:size(imageb,2)
pixb = [i;j;1];
s= T(3,:)*pixb; % depth information or the weight
% transformed pixel coordinates
inew = round(T(1,:)*pixb/s);
jnew = round(T(2,:)*pixb/s);
if (inew > 0 && jnew >0)
imageo(inew,jnew)=imageb(i,j); % aligning the transformed pixels
            fprintf('..\n');
end
end
end
count = 0;
for i = 1:size(imageo,1)
for j = 1:size(imageo,2)
if (imageo(i,j) >0)
count = count+1;
end
            fprintf('..\n');
end
end
fprintf('alignment_done!\n');
ratio = (size(imageb,1)*size(imageb,2)/count);
ratio = sqrt(ratio);
ratio=round(ratio);
