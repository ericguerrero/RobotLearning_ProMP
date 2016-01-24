function [Om_new,mwnew,WK,Swnew,Sy,MK,Sk,GTpos]=ProMP_WMLE(demoY,r,Om_old,WK,mw,Sw,Sy,pk,GGT)
%same as ProMP_MLE, but pk is a column vector with the weights/rewards to
%apply
pk=pk/(sum(pk)+0.00000000000000001); %in case weights are not normalized.
Sy0=Sy;
Sw0=Sw;
mw0=mw;
Om0=Om_old;
 
 
 
%% PREPROCESS DATA: WE DO NOT CONSIDER ROLLOUTS WITH WEIGHT<0.001
regularization=1e-4;

% Build kernels
d=size(demoY,2);
Nt=size(GGT,1)/r;
Ndemos=size(demoY,3);
%Sw=zeros(Nf*r,Nf*r); % weights variance
%WK=zeros(r*Nf,Ndemos); % all demonstration weights
%centers


 
% getGGT;
GTpos=GGT;
 
%% Fit the data with Identity Om (get mw, Sw, and wk)


YYcol=[];
for k=1:Ndemos
    for it=1:Nt
      YYcol=[YYcol;squeeze(demoY(it,:,k))'];  
    end
end

 
% %Get Sy mean over timesteps and demos:
 pmean=mean(squeeze(mean(demoY))')';
 Sy=zeros(d,d);

 for it=1:Nt
     for k=1:Ndemos
     Sy=Sy+(squeeze(demoY(it,:,k))'-pmean)*(squeeze(demoY(it,:,k))'-pmean)';
     end
 end
 Sy=Sy/Nt/Ndemos;
 

%Get big Sy:
SY=kron(eye(Nt),Sy);
 
%% E-STEP\
% get auxiliar matrices
OmN_old=kron(eye(Nt),Om_old);
AA=SY+OmN_old*GTpos*Sw*GTpos'*OmN_old';

%AA=BoundedConditioningCovMatrix(AA,1e9);
AAi=pinv(AA);%AA\eye(Nt*d);
%AAi=pseudoinverse(AA);%AA\eye(Nt*d);
 
BB=Sw;
CC=OmN_old*GTpos*Sw;
ZZ=OmN_old*GTpos;
Smeank=Sw-CC'*AAi*CC;%+1e-6*eye(Nf*r);

% regularization by imposing condition number constraint if necessary
%Smeank2=BoundedConditioningCovMatrix(Smeank,1e9);
Smeank2=(Smeank+Smeank')/2;
dk=0;
MK=[];
YALL=[];
for k=1:Ndemos
% get prior probability p(w|Yk;Om,Sw,mw) for every demo Yk
    YK=YYcol(1+d*Nt*(k-1):d*Nt*k,1);
    wmeank=mw+CC'*AAi*(YK-ZZ*mw);
    MK=[MK wmeank];
    YALL=[YALL YK];
end


%% get new average
mwnew=mw*0;
for k=1:Ndemos
   mwnew= mwnew+MK(:,k)*pk(k,1);
end


 
%% get new Sw and SKM (auxiliar computation)
Swnew=Smeank2;
SKM=zeros(size(Sw));
for k=1:Ndemos

   mk=MK(:,k);
    Swnew=Swnew+pk(k,1)*(mk-mwnew)*(mk-mwnew)';    
    SKM=SKM+(Smeank2+mk*mk')*pk(k,1);
end
 
%% Get new Om

% S1=zeros(r,r);
% S2=zeros(d,r);
% for it=1:Nt
%     phit=GTpos((it-1)*r+1:it*r,:);
%     S1=S1+phit*SKM*phit';
%    for k=1:Ndemos 
%          S2=S2+pk(k,1)*squeeze(demoY(it,:,k))'*(MK(:,k)'*phit');
%    end
% end
% Om_new=S2*pinv(S1);
Om_new=eye(d,r);
%Om_new=Om_old;
 
mw=mwnew;
Sw=Swnew;
Sk=Smeank2;
%Om_new=Om_old;
S4=zeros(d,d);
S5=zeros(d,d);
S6=zeros(d,d);
Om_use=Om_new;
for it=1:Nt
    phit=GTpos((it-1)*r+1:it*r,:);
    for k=1:Ndemos
        ytk=squeeze(demoY(it,:,k))';
        Om_phi_mu=Om_use*phit*MK(:,k);
        S4=S4+pk(k,1)*ytk*(ytk-Om_phi_mu)';
        S5=S5+pk(k,1)*Om_phi_mu*(Om_phi_mu-ytk)';
        
    end
    S6=S6+sum(pk)*Om_use*phit*Smeank2*phit'*Om_use';
end
Sy_new=(S4+S5+S6)/sum(pk)/Nt;
%[vs,ds]=eig(Sy_new);
Sy=real((Sy_new+Sy_new')/2)+eye(d)*0.000001;
    
    
    