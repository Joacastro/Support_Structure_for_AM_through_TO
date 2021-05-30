function [R1,R2,e1,e2] = MiGgap_Domains(Rmax,Rmin,Gap)
    Dmax = 2*Rmax; dmin = 2*Rmin;
    % second domain DII_2
    e1 = 0.05;
    [V]    = fsolve(@(X)dom2eq1(Dmax,dmin,e1,X),[pi,Dmax]);
    alp_b  = V(1);
    Db     = V(2);
    if Dmax < (Db+dmin)/2
        [V]    = fzero(@(X)dom2eq2(Dmax,dmin,e1,X),[pi,pi,Dmax]);
        Db     = V(3);
    end
    % Third Domain : DII_3
    Dc     = Dmax + Gap + dmin;
    alp_c  = 2 * acos(1-2*Dmax/Dc);
    bet_c  = 2 * acos(1-2*(Dmax+Gap)/Dc);
    e2 = (bet_c - sin(bet_c) - alp_c + sin(alp_c)) / (2*pi);
    % FINAL RADIUS 
    R1 = Db/2;
    R2 = Dc/2;
    %fprintf(1,'\n R1= %2.4f  R2= %2.4f' ,R1,R2);
    %fprintf(1,'\n e1= %1.4f  e2= %1.4f ',e1,e2);
 %------------------------------------
 function [F] = dom2eq1(Dmax,db,epsi,X)
    alp=X(1);
    Db =X(2);
    F(1) = 2*Dmax -Db*(1-cos(alp/2)); 
    F(2) = 2*pi*(Db^2-db^2)*(1-epsi)-Db^2*(alp-sin(alp))+2*pi*db^2;
function [F] = dom2eq2(Dmax,db,epsi,X)
    alp=X(1);
    bet=X(2);
    Db =X(3);
    F(1) = 2*Dmax -Db*(1-cos(alp/2));
    F(2) = 2*Dmax - (Db-db) - db*(1-cos(bet/2)); 
    F(3) = 2*pi*(Db^2-db^2)*(1-epsi)-Db^2*(alp-sin(alp))+db^2*(bet-sin(bet));