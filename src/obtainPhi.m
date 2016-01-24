function [PHI] = obtainPhi(Nf,Nt,dof,C,D)
phi = zeros(Nf,Nt);
for i=1:Nt      % trajectory points
    for s=1:Nf  % basis functions
        phi(s,i)=evalexp(i/Nt,C(s),D);
    end
end
PHI = kron(eye(dof),phi);