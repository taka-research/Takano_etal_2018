close all
clear all;

addpath ~/MATLAB

x=rdmds('XC');
y=rdmds('YC');
z=rdmds('RC');
zf=rdmds('RF');
da=rdmds('RAC');
dz=rdmds('DRF');
dyg=rdmds('DYG');
dxg=rdmds('DXG');
dv=da.*dz;
hc=rdmds('hFacC');
hw=rdmds('hFacW');
hs=rdmds('hFacS');
dv=dv.*hc;
KDOPRemin = 1.6E-8;
beta = 0.858;

[d0,iter]=rdmds('DIC01',NaN);
d1=rdmds('DIC02',NaN);
d2=rdmds('DIC03',NaN);
dop=rdmds('DOP01',NaN,'rec',1);
vel=rdmds('Trans',NaN);

N=size(hc);
Nt=length(iter);

c=squeeze(d0(:,:,:,1,:));
u=squeeze(vel(:,:,:,1,:));
v=squeeze(vel(:,:,:,2,:));
w=squeeze(vel(:,:,:,3,:));
uc=squeeze(d0(:,:,:,2,:));
vc=squeeze(d0(:,:,:,3,:));
wc=squeeze(d0(:,:,:,4,:));
uca=squeeze(d0(:,:,:,5,:));
vca=squeeze(d0(:,:,:,6,:));
wca=squeeze(d0(:,:,:,7,:));
dcx=squeeze(d0(:,:,:,8,:));
dcy=squeeze(d0(:,:,:,9,:));
dcz=squeeze(d0(:,:,:,10,:));
dczi=squeeze(d1(:,:,:,1,:));
cab=squeeze(d1(:,:,:,2,:));
ckpp=squeeze(d1(:,:,:,3,:));
bioa=squeeze(d1(:,:,:,4,:));
carb=squeeze(d1(:,:,:,5,:));
dgas=squeeze(d2(:,:,1,:));
fgco2=squeeze(d2(:,:,2,:));
fwf=squeeze(d2(:,:,3,:));
srelax=squeeze(d2(:,:,4,:));
surForcS=squeeze(d2(:,:,5,:));
sflux=squeeze(d2(:,:,6,:));
ssh=squeeze(d2(:,:,7,:));

dv1=repmat(dv,[1 1 1 Nt]);
for i=1:Nt
   dv1(:,:,1,i)=dv(:,:,1);%+ssh(:,:,i).*da;
end

%advective divergence
u=u.*hw.*repmat(dyg,[1 1 N(3),Nt]).*repmat(dz,[N(1) N(2) 1 Nt]);
v=v.*hs.*repmat(dxg,[1 1 N(3),Nt]).*repmat(dz,[N(1) N(2) 1 Nt]);
w=w.*repmat(da,[1 1 N(3),Nt]);
u(end+1,:,:,:)=u(1,:,:,:);
du = u(2:end,:,:,:)-u(1:end-1,:,:,:);
uca(end+1,:,:,:)=uca(1,:,:,:);
duca = uca(2:end,:,:,:)-uca(1:end-1,:,:,:);
v(:,end+1,:,:)=0;
dvvel=v(:,2:end,:,:)-v(:,1:end-1,:,:);
vca(:,end+1,:,:)=0;
dvca=vca(:,2:end,:,:)-vca(:,1:end-1,:,:);
w(:,:,end+1,:)=0;
dw=-w(:,:,2:end,:)+w(:,:,1:end-1,:);
wca(:,:,end+1,:)=0;
dwca=-wca(:,:,2:end,:)+wca(:,:,1:end-1,:);

udcdx = duca - c.*du;
vdcdy = dvca - c.*dvvel;
wdcdz = dwca - c.*dw;
%dadv = duca + dvca + dwca;
%dadvH=-(duca + dvca);
%dadvV=-dwca;
dadv  = udcdx+vdcdy+wdcdz;
dadvH = -(udcdx+vdcdy);
dadvV = -wdcdz;

%diffusive divergence 
dcx(end+1,:,:,:)=dcx(1,:,:,:);
ddcx = dcx(2:end,:,:,:)-dcx(1:end-1,:,:,:);
dcy(:,end+1,:,:)=0;
ddcy=dcy(:,2:end,:,:)-dcy(:,1:end-1,:,:);
dcz(:,:,end+1,:)=0;
ddcz=-dcz(:,:,2:end,:)+dcz(:,:,1:end-1,:);
dczi(:,:,end+1,:)=0;
ddczi=-dczi(:,:,2:end,:)+dczi(:,:,1:end-1,:);
ckpp(:,:,end+1,:)=0;
dckpp=-ckpp(:,:,2:end,:)+ckpp(:,:,1:end-1,:);

ddifH = -(ddcx + ddcy);
ddifV = -(ddcz + ddczi + dckpp);
ddif = ddcx + ddcy + ddcz + ddczi + dckpp;

%gasflux, upward -> divergence
dgas=dgas.*squeeze(dv1(:,:,1,:));
%dgas2=-fgco2.*da;

%carbon sink due to DOM production 
%production of DOP = 67%, 
dprodDOM=-bioa.*dv1*117*.67;


% 33% is the export production 
wc=zeros(N(1),N(2),N(3)+1,Nt);
dpoc=zeros(N(1),N(2),N(3),Nt);
mask=hc; mask(mask>0)=1;
for i=1:N(1)
 for j=1:N(2)
  for k=1:N(3)-1
   if hc(i,j,k)~=0
    if k==1 & hc(i,j,2)==0
      dpoc(i,j,k,:)=bioa(i,j,k,:)*0.33.*dv1(i,j,k,:)*117;
    else
      km=k+1;
      cexp=bioa(i,j,k,:)*.33*117.*dv1(i,j,k,:);
      zexp=zf(k+1);
      for k2=km:N(3)
         wc(i,j,k2,:)=wc(i,j,k2,:)+mask(i,j,k)*cexp(1,1,1,:).*abs(zf(k2)/zexp).^-beta;
      end
    end
   end
  end
 end
end

% remineralization
dpoc=wc(:,:,1:end-1,:)-wc(:,:,2:end,:);

% remin of CaCO3
dcar=carb.*dv1;

% remin of dom
dremDOM=KDOPRemin*dop*117.*dv1;

% biological pumps
dsoft = dpoc+dprodDOM+dremDOM;
dhard = dcar;

%% freshwater flux
cs=squeeze(c(:,:,1,:));
dfwf = -2.0282./1035.*fwf.*da;
%dsrelax = srelax/35/1000.*2.0.*da;
%dsflux = sflux/35/1000.*2.0.*da;

% carbon tendency (LHS)
cinv = c.*dv1;
dt=31104000*10;
dcdt = zeros(N(1),N(2),N(3),Nt);
dcdt(:,:,:,2:Nt-1) = (cinv(:,:,:,3:Nt)-cinv(:,:,:,1:Nt-2))/(2*dt);

% total flux convergence (RHS)
dc = dadv+ddif+dprodDOM+dpoc+dremDOM+dcar;
for n=1:Nt
   dc(:,:,1,n)=dc(:,:,1,n)+dgas(:,:,n)+dfwf(:,:,n);
end

%%
save dic_tendency_3dim.mat x y z iter zf dv dcdt dfwf dadvH dadvV ddifH ddifV dgas dsoft dhard;


% 185m mass balance
nterms=8;
terms={'horiz adv' 'vert adv' 'horiz mix' 'vert mix' 'gasex' 'fwf' 'soft pump' 'hard pump'};
V=zeros(N(1),N(2));
LHS=zeros(N(1),N(2),Nt);
RHS=zeros(N(1),N(2),Nt,nterms);
for k=1:8
   LHS = LHS + squeeze(dcdt(:,:,k,:));
   RHS(:,:,:,1) = RHS(:,:,:,1) + squeeze(dadvH(:,:,k,:));
   RHS(:,:,:,2) = RHS(:,:,:,2) + squeeze(dadvV(:,:,k,:));
   RHS(:,:,:,3) = RHS(:,:,:,3) + squeeze(ddifH(:,:,k,:));
   RHS(:,:,:,4) = RHS(:,:,:,4) + squeeze(ddifV(:,:,k,:));
   if k==1
      RHS(:,:,:,5) = RHS(:,:,:,5) + squeeze(dgas(:,:,1,:));
      RHS(:,:,:,6) = RHS(:,:,:,6) + squeeze(dfwf(:,:,1,:));
   end
   RHS(:,:,:,7) = RHS(:,:,:,7) + squeeze(dsoft(:,:,k,:));
   RHS(:,:,:,8) = RHS(:,:,:,8) + squeeze(dhard(:,:,k,:));
   V = V + dv(:,:,k);
end
dLHS = LHS./repmat(V,[1 1 Nt]);
dRHS = RHS./repmat(V,[1 1 Nt nterms]);

% zonal mean climatology
dLHS(dLHS==0)=NaN;
dLHSxt=squeeze(nanmean(nanmean(dLHS,3),1));
dRHS(dRHS==0)=NaN;
dRHSxt=squeeze(nanmean(nanmean(dRHS,3),1));

save dic_tendency_185m_8terms.mat LHS RHS dLHS dRHS dLHSxt dRHSxt nterms nterms V;

% 185m mass balance simpler version
nterms=7;
terms={'adv' 'horiz mix' 'vert mix' 'gasex' 'fwf' 'soft pump' 'hard pump'};
LHS=zeros(N(1),N(2),Nt);
RHS=zeros(N(1),N(2),Nt,nterms);
dadv = dadvH + dadvV;
for k=1:7
   LHS = LHS + squeeze(dcdt(:,:,k,:));
   RHS(:,:,:,1) = RHS(:,:,:,1) + squeeze(dadv(:,:,k,:));
   RHS(:,:,:,2) = RHS(:,:,:,2) + squeeze(ddifH(:,:,k,:));
   RHS(:,:,:,3) = RHS(:,:,:,3) + squeeze(ddifV(:,:,k,:));
   if k==1
      RHS(:,:,:,4) = RHS(:,:,:,4) + squeeze(dgas(:,:,1,:));
      RHS(:,:,:,5) = RHS(:,:,:,5) + squeeze(dfwf(:,:,1,:));
   end
   RHS(:,:,:,6) = RHS(:,:,:,6) + squeeze(dsoft(:,:,k,:));
   RHS(:,:,:,7) = RHS(:,:,:,7) + squeeze(dhard(:,:,k,:));
   V = V + dv(:,:,k);
end
dLHS = LHS./repmat(V,[1 1 Nt]);
dRHS = RHS./repmat(V,[1 1 Nt nterms]);

% zonal mean climatology
dLHS(dLHS==0)=NaN;
dLHSxt=squeeze(nanmean(nanmean(dLHS,3),1));
dRHS(dRHS==0)=NaN;
dRHSxt=squeeze(nanmean(nanmean(dRHS,3),1));

save dic_tendency_185m_7terms.mat LHS RHS dLHS dRHS dLHSxt dRHSxt nterms nterms V;

 
