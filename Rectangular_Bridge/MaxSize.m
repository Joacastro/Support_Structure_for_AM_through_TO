function [g,dgdz,MS] = MaxSize(h,dhdy,g,dgdz,DI,MS,penal,nx,ny,Iv)
%% ================= PARAMETERS & COMPUTATION OF DII MATRIX ====================
if ~isfield(MS,'DII')
    Rmax = MS.Rmax;   Rmin = MS.Rmin;           % Max and Min size radius
    [ECtrd]  = Elem_Geometry (nx,ny);           % Areas and Centroids
    [MS.DII,MS.cv]  = MSregions(ECtrd,Rmin,Rmax,nx*ny,ny); %  DII matrix
end 
%% ===================== Maximum Size Constraint using p-mean ==================
p  = MS.p; epsi = MS.epsi;                      % p-exponent & epsilon
if p >= 100
  q = min(penal+0.5,MS.qMax);                   % penalization of voids 
else
  q = min(2-penal,MS.qMax);                     % penalization of voids   
end
N  = nx*ny;                                     % Number of data        
v  = (1-h).^q;       dvdh = -q*(1-h).^(q-1);    % voids vector with penalization
dvdh(h>=0.999) = 0;
gms= epsi-(MS.DII*v+MS.cv); dgmsdv = -MS.DII';  % Amount of voids
s  = gms+1-epsi;                                % Quantity of Interest 
Pm   = (sum(s.^p)/N)^(1/p);                     % Aggregation using p-mean
dPmds= Pm^(1-p)*(s.^(p-1))/N;                   % Sensitivities
Gms  = Pm-1+epsi;                               % Final constraint
dGdz = DI'*(dhdy.*dvdh.*(dgmsdv*dPmds));        % Derivatives of Gms
g    = [ g  ; Gms  ];   dgdz = [dgdz,dGdz(Iv)]; % Final arrays for MMA

%% ======================== FUNCTION: DII Matrix  ==============================
function [DII,cv] = MSregions(Ec,RIN,RMAX,NElem,nely)
d1 = cell(NElem,1); 
xj = Ec(:,1); 
yj = Ec(:,2);
for el = 1:NElem
    xi     = Ec(el,1); yi =Ec(el,2); 
    % linear change of rfil-----------
    ri   = RMAX(1); 
    rf   = RMAX(2);
    Rmax = (ri-rf)/(nely-1)*(yi-1)+rf;
    rim  = RIN(1); 
    rfm  = RIN(2);
    Rin  = (rim-rfm)/(nely-1)*(yi-1)+rfm;
    %---------------------------------
    dist   = sqrt((xi-xj).^2 + (yi-yj).^2);
    [I,~]  = find(dist<=Rmax);                   
    [Iring,~]  = find(dist(I)>Rin);
    Im     = I(Iring);
    I      = Im;
    %Mirroring for the supports -----------------------------
    xs = 0.5;
    if (xi-Rmax) < xs
         yjs         =   yj(I,1) ;
         xjs         =   xj(I,1) + 2*(xs - xj(I,1));
         dists       =   sqrt((xi-xjs).^2 + (yi-yjs).^2);
         [Is,~]      =   find(dists<=Rmax);
         [Isring,~]  =   find(dists(Is)>Rin);
         Is          =   Is(Isring);
         I           =   [Im;I(Is)]; 
    else
         I           =   Im;
    end
    %---------------------------------------------------------------
    d1{el} = [I,zeros(size(I))+el,I*0+1];
end
dII = cell2mat(d1); 
DII = sparse(dII(:,2),dII(:,1),dII(:,3));      
SumDII = (sum(DII,2)); 				       % Extended domain
%SumDII = ones(NElem,1)*SumDII;
DII = spdiags(1./SumDII,0,NElem,NElem)*DII;
cv = ones(NElem,1) - DII*ones(NElem,1);        % Correction of voids
%% ======================== Element coordinates ================================
function [ElemCtrd] = Elem_Geometry(nelx,nely)
ElemCtrd = zeros(nelx*nely,2);
for i1=1:nelx
	for j1=1:nely
        e1=(i1-1)*nely+j1;
        ElemCtrd(e1,1) = i1;
        ElemCtrd(e1,2) = j1;               
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Fernandez E, Koppen S, Kumar P, Alarcon P,   %
% Bauduin S, and Duysinx P.                                                    %
% Please send your comments to: efsanchez@uliege.be                            %
% Disclaimer: The authors reserve all rights but does not guarantee that the   %
% code is free from errors. Furthermore, the author shall not be liable in any %
% event caused by the use of the program.                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%