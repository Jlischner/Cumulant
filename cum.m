load spec;
ws = spec(:,7);

elda = -5.956967;
vxc  = -10.455963;
eF = 6.5; %# fermi energy
ReSig = spec(:,8)-vxc;
ImSig = spec(:,10);

#load output;
#vxc = dmu*27.21;
#ws *= 27.21;
#elda = xi*27.21 + 0.00001;
#eF = ef*27.21;
#ReSig = real(Sig)*27.21 - vxc;
#ImSig = imag(Sig)*27.21;

wsN = [-45 :0.1: 2]'; %# freq. grid for cumulant A
Sigma = ReSig + I*ImSig;

A = 1/pi * abs(ImSig) ./ ( (ws-elda-ReSig).^2 + ImSig.^2);

Gamma = abs( ImSig )/pi;
dw = ws(2)-ws(1);

%# compute:
%# alphak = Im dSigma/dw(elda)
%# gammak = Re dSigma/dw(elda)
%# etak   = |Im Sigma(elda)|
%# deltaE = Re Sigma(elda) - Vxc

#elda_save = elda;
#elda = -0.3179;
if( min(ws) < elda && elda < max(ws))
  IndElda = (elda - min(ws))/dw + 1;
  distA   = IndElda - floor(IndElda);
  distB   = ceil(IndElda) - IndElda;
  GammaElda = distB*Gamma(floor(IndElda)) + distA*Gamma(ceil(IndElda));
  DGammaA = Gamma( floor(IndElda)+1) - Gamma( floor(IndElda) -1);
  DGammaA /= 2*dw;
  DGammaB = Gamma( ceil(IndElda)+1) - Gamma( ceil(IndElda) -1);
  DGammaB /= 2*dw;
  DGammaElda = distB*DGammaA + distA*DGammaB;

  deltaE  = distB*ReSig(floor(IndElda)) + distA*ReSig(ceil(IndElda));
  etak = pi*GammaElda;

  DSigmaA = Sigma( floor(IndElda)+1) - Sigma( floor(IndElda) -1);
  DSigmaA /= 2*dw;
  DSigmaB = Sigma( ceil(IndElda)+1) - Sigma( ceil(IndElda) -1);
  DSigmaB /= 2*dw;
  DSigmaElda = distB*DSigmaA + distA*DSigmaB;

  alphak  =  imag(DSigmaElda);
  gammak  = -real(DSigmaElda);
else
  printf("error - bad freq. range!");
  break;
endif;

alphak = 0.;
#deltaE = 0.;
#elda = elda_save;

ts   = [-60: 0.1 : 60]';
dt   = ts(2)-ts(1);
wt   = ws*ts';

%# calculate qp spectral function following Hedin/Almbladh
Cqp   = -I*ts*deltaE - etak*abs(ts) + I*alphak*sign(ts) - gammak;
expCqp= exp(-I*elda*ts) .* exp(Cqp);
exptw = exp(I*ts*ws');

%# compute full spectral function following Gunnarsson/Aryasetiawan
%# have to make choice for chemical potential - middle of gap?
CS = Gamma.*( (ws-elda) < eF )-GammaElda-(ws-elda)*DGammaElda;

%# get rid of divergence at ws=elda
fitrange = 10;
minInd = floor(IndElda) - fitrange;
maxInd = floor(IndElda) + fitrange;
P= polyfit((ws-elda)(minInd:maxInd),CS(minInd:maxInd),2);
PN  = [P(1) 0 0];
CSZ = polyval(PN,ws-elda);

#wCut= 0.01;
wCut = 1.5;
fCut= exp(-(ws-elda).^2/wCut^2);
CSN = CS  .* (1-fCut) + CSZ.*fCut;
CS2 = CSN ./ (ws-elda).^2;

%# Fourier transform to real time
CSt = dw*sum( dmult(CS2, exp(-I*wt)), 1);
CSt = transpose(CSt);
CSt = exp(I*elda*ts) .* CSt;
expC= exp(CSt) .* expCqp;

%# Fourier transform back to freq. space
exptwN = exp(I*ts*wsN');

%# quasiparticle spectral function
Aqp  = sum( dmult(expCqp, exptwN ),1);
Aqp  = transpose(Aqp);
Aqp *= dt/2/pi;

%# full spectral function
Ac = sum( dmult(expC, exptwN ),1);
Ac = transpose(Ac);
Ac *= dt/2/pi;

%# first order expansion
dA = dt* sum( dmult( expCqp .* CSt, exptwN), 1);
dA = transpose(dA);
dA /= 2*pi;

%# first order expansion without working in time
Eqp = elda + deltaE;
Aqpw = (etak*cos(alphak)-(ws-Eqp)*sin(alphak) )./( (ws-Eqp).^2 + etak.^2);
Aqpw *= exp(-gammak)/pi;
dA2 = dw*conv(Aqpw,CS2);

wsC  = dw*[0:length(Aqpw)+length(CS2)-2]';
wsC += 2*ws(1)-elda;
AqpC = (etak*cos(alphak)-(wsC-Eqp)*sin(alphak) )./( (wsC-Eqp).^2 + etak.^2);
AqpC *= exp(-gammak)/pi;
A1   = AqpC + dA2;

%# second order expansion without working in time
CS2conv = dw*conv(CS2,CS2)/2; 
wsconv  = dw*[0:2*length(CS2)-2]';
wsconv += 2*ws(1)-2*elda;

Aqpconv = (etak*cos(alphak)-(wsconv-Eqp)*sin(alphak) )./( (wsconv-Eqp).^2 + etak.^2);
Aqpconv *= exp(-gammak)/pi;
ddA     = dw*conv(Aqpconv, CS2conv);
wsdd    = dw*[0:2*length(wsconv)-2]';
wsdd   += 2*wsconv(1); 
ddAI    = interp1(wsdd,ddA,wsC);
A2      = A1 + ddAI;

%# plasmon pole model
[a,b] = max(Gamma);
wmax  = ws(b);
lambda= 0.002;
dA3   = lambda/pi/(wmax-elda)^2 *Aqpw;
wsPP  = ws-(wmax-elda);
AqpPP = (etak*cos(alphak)-(wsPP-Eqp)*sin(alphak) )./( (wsPP-Eqp).^2 + etak.^2);
AqpPP *= exp(-gammak)/pi;
APP   = AqpPP + dA3;