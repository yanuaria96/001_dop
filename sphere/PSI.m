function psi = PSI(thetai,phii,thetaj,phij)
myu=[0;0;1];
k=2;

ri=[sin(thetai).*cos(phii);sin(thetai).*sin(phii);cos(thetai)];
rj=[sin(thetaj).*cos(phij);sin(thetaj).*sin(phij);cos(thetaj)];
rij=norm(ri-rj);

%psi=1/(rij^2);
%psi=1/(rij^2)/sqrt(Mises_Fisher(k,myu,ri))/sqrt(Mises_Fisher(k,myu,rj));
%psi=1/(rij^2)/sqrt(Valee_Poussin(k,myu,ri))/sqrt(Valee_Poussin(k,myu,rj));
psi=1/(rij^2)/sqrt(PDF(thetai,phii))/sqrt(PDF(thetaj,phij));
end