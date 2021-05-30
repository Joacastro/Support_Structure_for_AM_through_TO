function [Int,Dil,Ero,c] = Rec_Bridge_MG(nelx,nely,volfrac,gap)
%% Continuation Method ========================================================
penal = 1;
% penal_i = 1 ;           % **UD** Initial penal
% penal_f = 3 ;           % **UD** final penal
% beta_i = 1.5 ;          % **UD** Initial beta (check continuation method)
% beta_end = 38 ;         % **UD** Final beta
%% ROBUST DESIGN APPROACH (ERODED-INTERMEDIATE-DILATED)----------------------
rmin            = [5.1,2.6]         ;					% Imposed MinSize (Solid)
ri              = rmin(1)*2/0.749 ;                     % Initial radius
rf              = rmin(2)*2/0.749 ;                     % Final radius
mue             = 0.639           ;                     % Erosion threshold
mui             = 0.5             ;                     % Intermediate threshold
mud             = 0.20            ;                     % Dilated threshold 
beta            = 1.5             ;                     % Initial beta
vRef            = volfrac         ;                     % For scaling volfrac
loopbetaMax     = 50              ;                     % Iter. to increase beta
tdil            = rmin*1.0494     ;                     % Offset distance between the intermediate and dilated design
%% Maximum size constraint
MSi.p    = 125;        MSd.p    = 125;
MSi.Rmin = rmin;       MSi.Rmax = MSi.Rmin+[2.6 1.1];  MSi.qMax = 2;  MSi.epsi = 0.05;
MSd.Rmin = rmin+tdil;  MSd.Rmax = MSi.Rmax + tdil;     MSd.qMax = 2;  MSd.epsi = 0.05;
%% Minimum Gap constraint
Gap = [5,gap];          % Gap = [Min. Gap in the base, Min. Gap in the contact]
%% Min. Gap in Intermediate design
[~,R2,~,e2] = MiGgap_Domains(MSi.Rmax(1),rmin(1),Gap(1));
MGi.Rmin = 0; MGi.Rmax(1,1) = R2;  MGi.qMax=1; MGi.p=60; MGi.epsi(1,1) = e2;
[~,R2,~,e2] = MiGgap_Domains(MSi.Rmax(2),rmin(2),Gap(2));
              MGi.Rmax(1,2) = R2;                        MGi.epsi(1,2) = e2;
%% Min. Gap in Dilated design              
% [~,R2,~,e2] = MiGgap_Domains(MSi.Rmax(1)+tdil(1),rmin(1)+tdil(1),Gap(1)-tdil(1));
% MGd.Rmin = 0;   MGd.Rmax(1,1)=R2;   MGd.epsi(1,1) = e2;    MGd.qMax = 1;  MGd.p = 60;
% 
% [~,R2,~,e2] = MiGgap_Domains(MSi.Rmax(2)+tdil(2),rmin(2)+tdil(2),Gap(2)-tdil(2));
%               MGd.Rmax(1,2) = R2;    MGd.epsi(1,2) = e2;
%% MATERIAL PROPERTIES
E0   = 1;
Emin = 1e-9;
nu   = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS 
% Gap between the branches and the component (GB)
%----------------------------------------------------------
GB = 10;    % Number of elements gap with the component
%----------------------------------------------------------
F = sparse([2:2*(nely+1):2*(nelx+1-GB)*(nely+1)],1,1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
% Gap between the trunks and the component (GT)
% ----------------------------------------------------
GT = 5;    % Number of elements gap with the component
% ----------------------------------------------------
fixeddofs = [1:2:2*(nely+1),2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1-GT), ... 
             2*(nely+1)-1:2*(nely+1):2*(nely+1)*(nelx+1-GT)];
alldofs   = [1:2*(nely+1)*(nelx+1)];
freedofs  = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
% H matrix indeces
% ------------------------------------------------------------
iH = ones(nelx*nely*(2*(ceil(ri)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
% ------------------------------------------------------------
% Filtering regions (for the compute of denominator of H)
% ------------------------------------------------------------
Hsi= zeros(nely,1);           % Filtering regions
Hs = zeros(nelx*nely,nely);   % Index for Filt. Regions
k = 0;
% ------------------------------------------------------------
                                    
for i1 = 1:nelx
  for j1 = 1:nely
      
    e1 = (i1-1)*nely+j1;                 % Number of the element
    
    % linear change of rfil
    % ----------------------------------
    rfil = (ri-rf)/(nely-1)*(j1-1)+rf;
    %-----------------------------------
    
    % Filter in the boundaries
    % -----------------------------------------------------------------
    hsij = j1;                           % y coordinate
    Hs(e1,hsij) = 1;                     % Index w.r.t the bottom edge 
    % -----------------------------------------------------------------
    
    for i2 = max(i1-(ceil(rfil)-1),1):min(i1+(ceil(rfil)-1),nelx)
      for j2 = max(j1-(ceil(rfil)-1),1):min(j1+(ceil(rfil)-1),nely)
        e2    = (i2-1)*nely+j2;
        k     = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,1-sqrt((i1-i2)^2+(j1-j2)^2)/rfil);
        
        % Symmetry condition
        % -----------------------------------------------------
        xs = 0.5;
        if (i2 - (ceil(rfil)-1)) < xs
            i2s = i2 + 2*(xs - i2);
            j2s = j2;
            sHs = max(0,1-sqrt((i1-i2s)^2+(j1 - j2s)^2)/rfil);
            sH(k) = sH(k) + sHs;
        end
        % -----------------------------------------------------
               
        % Filter in the boundaries (take the x midle elemtents to compute Hs)
        %---------------------------------------------------------------------
        if i1 == round(nelx/2) 
            Hsi(hsij,1) = Hsi(hsij,1) + sH(k);            
        end
        %---------------------------------------------------------------------        
      end
    end
  end
end
H = sparse(iH,jH,sH);                      
Hs= Hs*Hsi;                                
H = spdiags(1./Hs,0,nelx*nely,nelx*nely)*H;
%% Passive elements ---------------------------------------------------------
Pas_1 = zeros(nely,nelx);                          % Passive elements matrix
Pas_1(1:1,1:nelx-GB) = 1;                          % Passive solid elements
Pas_0 = ones(nely,nelx);                           % Passive elements matrix
Pas_0(1:1,nelx-GB:nelx) = 0;                       % Passive void elements
Ip  = find(Pas_1(:)==1);                           % Solid Passive Indeces
Ipv = find(Pas_0(:)==0);                           % Void passive indeces
Iv  = find(Pas_1(:)==0);                           % Optimizable design dom.
Nv  = length(Iv);                                  % Num. of variables
%% FOR MMA ------------------------------------------------------------------
xold1   = zeros(Nv,1); xold2 = zeros(Nv,1);         % Parameters for mmasub
low     = zeros(Nv,1); upp   = zeros(Nv,1);         %
movemax = 0.5        ; movemin = 0.1;               % Move limits
cRef    = 1.0        ;                              % To scale obj. 
%% INITIALIZE ITERATION                             %
x       = repmat(volfrac,nely,nelx); 				% size of x matrix
loopbeta= 0; loop = 0; change = 1; 					% initialization of loops
%% ROBUST DESIGN APPROACH (ERODED-INTERMEDIATE-DILATED)----------------------
Tilde  = x                         ;                % Filtered field          
Ero    = Heaviside(Tilde,mue,beta) ;                % eroded
Int    = Heaviside(Tilde,mui,beta) ;                % intermediate
Dil    = Heaviside(Tilde,mud,beta) ;                % dilated   
%% START ITERATION                                  %
while change > 0.001 && loop<(loopbetaMax*10) && penal<=2
  loopbeta  = loopbeta+1;
  loop      = loop+1;  
  mL        = (movemin-movemax)/(2-1)*(penal-1)+movemax;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+Ero(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c  = sum(sum((Emin+Ero.^penal*(E0-Emin)).*ce)); 
  dc = -penal*(E0-Emin)*Ero.^(penal-1).*ce;
  if loopbeta<=3; cRef=abs(c)/30; end
  c  = c/cRef; dc = dc./cRef;
  %% VOLUME CONSTRAINT   
  if rem(loop,10)==0;volfrac=vRef*sum(Dil(:))/sum(Int(:)); end          % Update of volume restricction in the dilated design
  v     = sum(sum(Dil,1),2)/(nelx*nely)/volfrac-1;
  dv    = ones(nelx*nely,1)/(nelx*nely)/volfrac;
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  dEro  = Deriv_Heaviside(Tilde,mue,beta);
  dInt  = Deriv_Heaviside(Tilde,mui,beta);
  dDil  = Deriv_Heaviside(Tilde,mud,beta);
  dc    = H'*(dc(:).*dEro(:)); dc = dc(Iv);   
  dv    = H'*(dv.*dDil(:));    dv = dv(Iv);
  [v,dv,MSi] = MaxSize(Int(:),dInt(:),v,dv,H,MSi,penal,nelx,nely,Iv);       % Maximum size in the intermediate design
  [v,dv,MSd] = MaxSize(Dil(:),dDil(:),v,dv,H,MSd,penal,nelx,nely,Iv);       % Maximum size in the dilated design
  if penal <=1.25+0.25
    [v,dv,MGi] = MaxSize_MinGap(Int(:),dInt(:),v,dv,H,MGi,penal,nelx,nely,Iv);   % Min. Gap in the intermediate design
%    [v,dv,MGd] = MaxSize_MinGap(Dil(:),dDil(:),v,dv,H,MGd,penal,nelx,nely,Iv);   % Min. Gap in the dilated design
  end
  %% UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES  
  xmin = max(zeros(Nv,1),x(Iv)-mL);  xmax=min(ones(Nv,1),x(Iv)+mL);  
  m    = size(v,1); 
  [xnew,~,~,~,~,~,~,~,~,low,upp] = ... 
       mmasub(m,Nv,loopbeta,x(Iv),xmin,xmax,xold1,xold2,c,dc(:),v,dv', ...
       low,upp,1,zeros(m,1),1000*ones(m,1),zeros(m,1),1);
  change = max(abs(xnew(:)-x(Iv)));  xold1=xnew;  xold2=x(Iv); x(Iv)=xnew;
  %% ROBUST DESIGN APPROACH (ERODED-INTERMEDIATE-DILATED)--------------------
  x(Ip)   = 1                                              ;  % Passive Elem
  x(Ipv)  = 0                                              ;
  Tilde(:)= H*x(:); Tilde(Ip)=1; Tilde(Ipv) = 0            ;  % Filtering
  Ero= Heaviside(Tilde,mue,beta); Ero(Ip)=1                ;  % eroded
  Int= Heaviside(Tilde,mui,beta)                           ;  % intermediate
  Dil= Heaviside(Tilde,mud,beta)                           ;  % dilated
  %% PRINT RESULTS & PLOT DENSITIES
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f Mnd.: %2.1f\n',loop,c, ...
    mean(Int(:)),change,sum(4*Int(:).*(1-Int(:)))/(nelx*nely)*100);
  %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
  if loopbeta >= loopbetaMax || change <= 0.01
    penal=penal+0.125;   beta=min(1.5*beta,38.4);           % Continuation
    loopbeta = 0; change=1;                                 % Recount
    fprintf('beta: %g   penal: %g  \n',beta,penal);         % Info. Parameters
  end
end
%% FINAL RESULT : Compliance for the Intermediate Design
sK=reshape(KE(:)*(Emin+Int(:)'.^2.0*(E0-Emin)),64*nelx*nely,1);
K=sparse(iK,jK,sK); K=(K+K')/2; U(freedofs)=K(freedofs,freedofs)\F(freedofs);
ce=sum((U(edofMat)*KE).*U(edofMat),2);
c =sum((Emin+Int(:).^2*(E0-Emin)).*ce); fprintf('Obj(Int): %6.3f.\n',c);
function H = Heaviside(x,mu,beta)
H = (tanh(beta*mu) + tanh(beta*(x-mu)))/(tanh(beta*mu) + tanh(beta*(1-mu)));
function dH = Deriv_Heaviside(x,mu,beta)
dH = (beta*(sech(beta*(x-mu))).^2) /(tanh(beta*mu)+tanh(beta*(1-mu)));