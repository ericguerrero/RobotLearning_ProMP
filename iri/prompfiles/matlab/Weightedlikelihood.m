function L=Weightedlikelihood(Sw,Sk,Sy,demoY,Om,MK,GTpos,mw,pk)
Nt=size(demoY,1);
r=size(Om,2);
d=size(demoY,2);
Ndemos=size(demoY,3);
Swi=pinv(Sw);
Syi=pinv(Sy);
L=-1/2*(log(2*pi)+sum(log(eig(Sw)))+Nt*log(2*pi*det(Sy)))*sum(pk);
for k=1:Ndemos
    mk=MK(:,k);
    L=L-0.5*(trace(Swi*Sk)+(mk-mw)'*Swi*(mk-mw))*pk(k,1);
    for i=1:Nt
       phit=GTpos((i-1)*r+1:i*r,:);
       %L=L+squeeze(demoY(i,:,k))*Syi*(0.5*Om*phit*mk-0.5*squeeze(demoY(i,:,k))')*pk(k,1);
       L=L+(squeeze(demoY(i,:,k))*Syi*(0.5*Om*phit*mk-0.5*squeeze(demoY(i,:,k))')+0.5*mk'*phit'*Om'*Syi*squeeze(demoY(i,:,k))')*pk(k,1);
       AUX=phit'*Om'*Syi*Om*phit;
       L=L-0.5*(mk'*AUX*mk+trace(AUX*Sk))*pk(k,1);
   end
    
end
