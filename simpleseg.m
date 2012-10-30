function seg=simpleseg(I,initmask,maxits,E,T,alpha)
%--Create a signed distance map(SDF)frommask
phi=bwdist(initmask)-bwdist(1-initmask)-.5;
% mainloop
for its=1:maxits
D=E-abs(I-T);
K=getcurvature(phi);
F=-alpha*D+(1-alpha)*K;
dxplus=shiftR(phi)-phi;
dyplus=shiftU(phi)-phi;
dxminus=phi-shiftL(phi);
dyminus=phi-shiftD(phi);
gradphimaxx=sqrt(max(dxplus,0).^2+max(-dxminus,0).^2);
gradphiminx=sqrt(min(dxplus,0).^2+min(-dxminus,0).^2);
gradphimaxy=sqrt(max(dyplus,0).^2+max(-dyminus,0).^2);
gradphiminy=sqrt(min(dyplus,0).^2+min(-dyminus,0).^2);
gradphimax=sqrt((gradphimaxx.^2)+(gradphimaxy.^2));
gradphimin=sqrt((gradphiminx.^2)+(gradphiminy.^2));
gradphi=(F>0).*(gradphimax)+(F<0).*(gradphimin);
%stabilityCFL
dt=.5/max(max(abs(F.*gradphi)));
%evolvethecurve
phi=phi+dt.*(F).*gradphi;
%reinitialisedistancefuncitonevery50iterations
if(mod(its,50)==0)
phi=bwdist(phi<0)-bwdist(phi>0);
end
%intermediateoutput
if(mod(its,20)==0)
showcontour(I,phi,its);
subplot(2,2,4);surf(phi);shading flat;
end
end
%makemaskfromSDF
seg=phi<=0;% --Getmaskfromlevelset
% --wholematrixderivatives

function shift=shiftD(M)
shift=shiftR(M')';

function shift=shiftL(M)
shift=[M(:,2:size(M,2)) M(:,size(M,2))];

function shift=shiftR(M)
shift=[M(:,1) M(:,1:size(M,2)-1)];

function shift=shiftU(M)
shift=shiftL(M')';

function curvature=getcurvature(phi)
dx=(shiftR(phi)-shiftL(phi))/2;
dy=(shiftU(phi)-shiftD(phi))/2;
dxplus=shiftR(phi)-phi;
dyplus=shiftU(phi)-phi;
dxminus=phi-shiftL(phi);
dyminus=phi-shiftD(phi);
dxplusy=(shiftU(shiftR(phi))-shiftU(shiftL(phi)))/2;
dyplusx=(shiftR(shiftU(phi))-shiftR(shiftD(phi)))/2;
dxminusy=(shiftD(shiftR(phi))-shiftD(shiftL(phi)))/2;
dyminusx=(shiftL(shiftU(phi))-shiftL(shiftD(phi)))/2;
nplusx=dxplus./sqrt(eps+(dxplus.^2)+((dyplusx+dy)/2).^2);
nplusy=dyplus./sqrt(eps+(dyplus.^2)+((dxplusy+dx)/2).^2);
nminusx=dxminus./sqrt(eps+(dxminus.^2)+((dyminusx+dy)/2).^2);
nminusy=dyminus./sqrt(eps+(dyminus.^2)+((dxminusy+dx)/2).^2);

curvature=((nplusx-nminusx)+(nplusy-nminusy)/2);

% --Displaystheimagewithcurvesuperimposed
function showcontour(I,phi,i)
subplot(2,2,3);title('Evolution');
imshow(I,'initialmagnification',200,'displayrange',[0 255]);
hold on;
contour(phi,[0 0],'g','LineWidth',2);
hold off;
title([num2str(i) 'Iterations']);
drawnow;

