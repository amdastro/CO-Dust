; read in opacity file

file = 'freedman.dat'
;file = 'dustopacities.dat'

precc = 360000
tempOP = dblarr(precc)
rhoOP = dblarr(precc)
kappaOP = dblarr(precc)
OPENR, 2, file
temp1 = 0.0
temp2 = 0.0
temp3 = 0.0
FOR i=0, precc-1 DO BEGIN
  READF, 2, temp1, temp2, temp3
  rhoOP(i) = temp1
  tempOP(i) = temp2
  kappaOP(i) = temp3
ENDFOR
tempOP = 10^(tempOP)
rhoOP = 10^(rhoOP)
nOP = rhoOP/1.67e-24
kappaOP = 10^(kappaOP) 

CLOSE, 2

file = 'H_He_RHO-T_ideal_eq.d'

precc = 360000
tempEOS = dblarr(precc)
rhoEOS = dblarr(precc)
gammaEOS = dblarr(precc)
eintEOS = dblarr(precc)
prEOS = dblarr(precc)
OPENR, 1, file
temp1 = 0.0
temp2 = 0.0
temp3 = 0.0
temp4 = 0.0
temp5 = 0.0
temp6 = 0.0
temp7 = 0.0
temp8 = 0.0
temp9 = 0.0
FOR i=0, precc-1 DO BEGIN
  READF, 1, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
  rhoEOS(i) = temp1
  tempEOS(i) = temp2
  prEOS(i) = temp3
  eintEOS(i) = temp4
  gammaEOS(i) = temp7
ENDFOR
rhoEOS(where(rhoEOS gt 0.0)) = 0.0
tempEOS(where(tempEOS gt 8.0)) = 0.0
tempEOS = 10^(tempEOS)
rhoEOS = 10^(rhoEOS)
prEOS = 10^(prEOS)
eintEOS = 10^(eintEOS)
nEOS = rhoEOS/1.67e-24
enthEOS = (eintEOS + prEOS)/(rhoEOS)

CLOSE, 1

Msun = 2.0e33
sigSB = 5.67e-5
kb = 1.38e-16
mu = 0.8
mp = 1.67e-24
Lsun = 4.0e33
c = 3.0e10

; mass of shell
M = 1.0e-4*Msun
; time of collision in seconds
t0 = 1.0e6
; velocity of shell
v = 5.0e7
v8 = v/1.0e8
; initial density of shell
n0 = 5.0e13*(t0/(2.0*7.*24.*3600.))^(-3)*(v8^(-1.0))*(M/(1.0e-5*Msun))
; initial radius of shell
R0 = t0*v

Lsh = 1.0e38
vsh = 500.e5
P0 = Lsh*(2./(3.0*!pi*R0^(2.0)*vsh))
tprec = 1000
tempmin = alog(1.0e3)
tempmax = alog(1.0e5)
temprange = dindgen(tprec)*(tempmax-tempmin)/(tprec-1.0) + tempmin
temprange = exp(temprange)
nrange = dblarr(tprec)
eintrange = dblarr(tprec)
kapparange = dblarr(tprec)
FOR j=0, tprec-1 DO BEGIN
  hog = where((tempEOS gt temprange(j)*0.97 and (tempEOS lt temprange(j)*1.03)))
  pig = where((prEOS gt P0*0.90) and (prEOS lt P0*1.1))
  listmatch, hog, pig, dog, cat
  nrange(j) = nEOS(hog(dog[0]))
  eintrange(j) = eintEOS(hog(dog[0]))
  hog = where((tempOP gt temprange(j)*0.97 and (tempOP lt temprange(j)*1.03)))
  pig = where((nOP gt nrange(j)*0.75) and (nOP lt nrange(j)*1.25))
  listmatch, hog, pig, dog, cat
  kapparange(j) = kappaOP(hog(dog[0]))
ENDFOR
taurange = M*kapparange/(4.0*!pi*(R0^(2.0)))
Lcoolrange = (16.*!pi*sigSB*R0^(2.0)*temprange^(4.0))/(9.*(taurange + 2./(3.0*taurange) + 1.33))
tcoolrange = (M/(mp*nrange))*eintrange/(Lcoolrange)

;set_plot, 'ps'
;device, filename = 'tcool.eps', /color

;PLOT, [0.0], [0.0], /nodata, /xlog, /ylog, ystyle=1, yrange = [1.0, 1.0e7], xtickname = ['10!u3!n', '10!u4!n'], color = FSC_COLOR("Black"), xmargin = [12,2], ymargin = [5,2], charthick=4., xcharsize = 1.5, ycharsize = 1.5, xthick=5., ythick=5., xstyle=1, xrange = [1000.,30000.], thick=5, ytitle = 't!dcool!n (s)', xtitle = 'T (K)'
;OPLOT, temprange, tcoolrange2, linestyle=0, thick=5, color = FSC_COLOR("Black")
;XYOUTS, [0.55], [0.89], 'Freedman+08 (No Dust)', charthick=5, /normal, charsize = 1.5
;OPLOT, temprange, tcoolrange, linestyle=0, thick=5, color = FSC_COLOR("Blue")
;XYOUTS, [0.55], [0.84], 'Ferguson+05 (Dust)', charthick=5, /normal, charsize = 1.5
;OPLOT, [1000., 50000], [1.0e6, 1.0e6], linestyle=1, thick=3., color = FSC_COLOR("Black")

;device, /close


; create time array
tmin = alog(1.0e-4*t0)
tmax = alog(50.0*t0)
tprec = 15000
t = dindgen(tprec)*(tmax-tmin)/(tprec-1.0) + tmin
t = exp(t)

nPR = 1.0e15*((t+t0)/(3.0*3600.*24.))^(-3.0)
tempPR = 6000.*((t+t0)/(3.0*3600.*24.))^(-0.5)
R = t*v + R0
; mean density
nbar = M/(4.0*!pi*R^(3.0)*mp)


; define shock luminosity and pressure
Lsh = 1.0e38
vsh = 500.e5
P0 = Lsh*(2./(3.0*!pi*R0^(2.0)*vsh))
tau0 = 3.0*t0
Lsh = Lsh/(1.0 + (t/tau0)^(4.0))
Psh = Lsh*(2./(3.0*!pi*R^(2.0)*vsh))
tempsh = (Lsh/(4.0*!pi*R^(2.0)*sigSB))^(0.25)
tempPR2 = (1.0e38/(4.0*!pi*R^(2.0)*sigSB))^(0.25)

LWD = 1.0e38
tempWD = (3.0*LWD/(16.0*!pi*R^(2.0)*sigSB))^(0.25)

n = dblarr(tprec)
gam = dblarr(tprec)
cs = dblarr(tprec)
pr = dblarr(tprec)
kappa = dblarr(tprec)
eint = dblarr(tprec)
temp = dblarr(tprec)
tau = dblarr(tprec)
Vol = dblarr(tprec)
delta = dblarr(tprec)
Lum = dblarr(tprec)
LumPdV = dblarr(tprec)
dndt = dblarr(tprec)
enth = dblarr(tprec)

; define initial temperature
n(0) = n0
temp(0) = 3.0e4

; calculate initial internal energy and density
hog = where((tempEOS gt temp(0)*0.98 and (tempEOS lt temp(0)*1.02)))
pig = where((prEOS gt Psh(0)*0.95) and (prEOS lt Psh(0)*1.05))
listmatch, hog, pig, dog, cat
eint(0) = eintEOS(hog(dog[0]))
n(0) = nEOS(hog(dog[0]))
enth(0) = enthEOS(hog(dog[0]))

; calculate initial opacity and optical depth
hog = where((tempOP gt temp(0)*0.995 and (tempOP lt temp(0)*1.005)))
pig = where((nOP gt n(0)*0.75) and (nOP lt n(0)*1.25))
listmatch, hog, pig, dog, cat
kappa(0) = kappaOP(hog(dog[0]))
tau(0) = M*kappa(0)/(4.0*!pi*R(0)^(2.0))

radtag = 1
skip = 0
FOR i=1, tprec-1 DO BEGIN
  IF(skip eq 0) THEN BEGIN
  Lrad =  (1./M)*(16.*!pi*sigSB*R(i-1)^(2.0)*temp(i-1)^(4.0))/(9.*(tau(i-1) + 2./(3.0*tau(i-1)) + 1.33))
  Lum(i) = Lrad
  PdV = (Psh(i)-Psh(i-1))/(t(i)-t(i-1))
  PdV = PdV/(mp*n(i-1)) 
  LumPdV(i) = PdV
  enth(i) = enth(i-1) + (PdV-Lrad)*(t(i)-t(i-1))
  
  hog = where((enthEOS gt enth(i)*0.98 and (enthEOS lt enth(i)*1.02)))
  pig = where((prEOS gt Psh(i)*0.95) and (prEOS lt Psh(i)*1.05))
  listmatch, hog, pig, dog, cat
  temp(i) = tempEOS(hog(dog[0]))
  n(i) = nEOS(hog(dog[0]))  
  
  IF(dog[0] eq -1) THEN BEGIN
    hog = where((enthEOS gt enth(i)*0.90 and (enthEOS lt enth(i)*1.1)))
    pig = where((prEOS gt Psh(i)*0.90) and (prEOS lt Psh(i)*1.1))
    listmatch, hog, pig, dog, cat
    temp(i) = tempEOS(hog(dog[0]))
    n(i) = nEOS(hog(dog[0]))
  ENDIF
  IF(dog[0] eq -1) THEN BEGIN
    Print, 'hey'
    STOP
  ENDIF
  ENDIF
  
  IF(radtag eq 1) THEN BEGIN
    IF(temp(i) lt tempsh(i)) THEN BEGIN
      skip = 1
      temp(i) = tempPR(i)
;      temp(i) = tempsh(i)
      hog = where((tempEOS gt temp(i)*0.98 and (tempEOS lt temp(i)*1.02)))
      pig = where((prEOS gt Psh(i)*0.95) and (prEOS lt Psh(i)*1.05))
      listmatch, hog, pig, dog, cat
      n(i) = nEOS(hog(dog[0]))

      IF(dog[0] eq -1) THEN BEGIN
        hog = where((tempEOS gt temp(i)*0.90 and (tempEOS lt temp(i)*1.1)))
        pig = where((prEOS gt Psh(i)*0.80) and (prEOS lt Psh(i)*1.2))
        listmatch, hog, pig, dog, cat
        n(i) = nEOS(hog(dog[0]))
      ENDIF
      IF(dog[0] eq -1) THEN BEGIN
        Print, 'hey'
        STOP
      ENDIF
    ENDIF
  ENDIF
  
  hog = where((tempOP gt temp(i)*0.995 and (tempOP lt temp(i)*1.005)))
  pig = where((nOP gt n(i)*0.75) and (nOP lt n(i)*1.25))
  listmatch, hog, pig, dog, cat
  kappa(i) = kappaOP(hog(dog[0]))
  tau(i) = M*kappa(i)/(4.0*!pi*R(i)^(2.0))
ENDFOR


set_plot, 'ps'
device, filename = 'rho.eps', /color

PLOT, [0.0], [0.0], /nodata, /xlog, /ylog, ystyle=1, yrange = [1.0e9, 3.0e15], color = FSC_COLOR("Black"), xmargin = [12,2], ymargin = [5,2], charthick=4., xcharsize = 1.5, ycharsize = 1.5, xthick=5., ythick=5., xstyle=1, xrange = [100.,max(t)], thick=5, ytitle = 'n (cm!u-3!n)', xtitle = 't (s)'
OPLOT, t, nbar, linestyle=1, thick=10, color = FSC_COLOR("Red")
XYOUTS, [0.5], [0.22], 'n!davg!n', /normal, charthick=5, charsize = 1.5
OPLOT, t, n, linestyle=0, thick=10, color = FSC_COLOR("Blue")
XYOUTS, [0.7], [0.82], 'n', /normal, charthick=5, charsize = 1.5
;OPLOT, t, n2, linestyle=2, thick=5, color = FSC_COLOR("Blue")
OPLOT, t, nPR, linestyle=3, thick=10, color = FSC_COLOR("Brown")
XYOUTS, [0.6], [0.62], 'n!dPR!n', /normal, charthick=5, charsize = 1.5

device, /close


set_plot, 'ps'
device, filename = 'temp.eps', /color

PLOT, [0.0], [0.0], /nodata, /xlog, /ylog, ystyle=1, ytickname = ['10!u3!n', '10!u4!n'], yrange = [300., 1.2e4], color = FSC_COLOR("Black"), xmargin = [14,2], ymargin = [5,2], charthick=4., xcharsize = 1.5, ycharsize = 1.5, xthick=5., ythick=5., xstyle=1, xrange = [100.,max(t)], thick=5, ytitle = 'T (K)', xtitle = 't (s)'
;OPLOT, t, Trad, linestyle=1, thick=2, color = FSC_COLOR("Red")
OPLOT, t, temp, linestyle=0, thick=5, color = FSC_COLOR("Blue")
XYOUTS, [0.3], [0.83], 'T', /normal, charthick=5, charsize = 1.5
;OPLOT, t, temp2, linestyle=2, thick=5, color = FSC_COLOR("Blue")
OPLOT, t, tempPR, linestyle=3, thick=5, color = FSC_COLOR("Brown")
XYOUTS, [0.3], [0.68], 'T!dPR!n', /normal, charthick=5, charsize = 1.5
OPLOT,t, tempsh, linestyle=1, thick=5, color = FSC_COLOR("Red")
XYOUTS, [0.3], [0.52], 'T!dmin!n', /normal, charthick=5, charsize = 1.5

device, /close


set_plot, 'ps'
device, filename = 'temp-rho.eps', /color


PLOT, [0.0], [0.0], /nodata, /xlog, /ylog, ystyle=1, ytickname = ['10!u3!n', '10!u4!n'], yrange = [300., 1.2e4], color = FSC_COLOR("Black"), xmargin = [14,2], ymargin = [5,2], charthick=4., xcharsize = 1.5, ycharsize = 1.5, xthick=5., ythick=5., xstyle=1, xrange = [1.0e9,2.0e15], thick=5, ytitle = 'T (K)', xtitle = 'n (cm!u-3!n)'
;OPLOT, nmean, Trad, linestyle=1, thick=2, color = FSC_COLOR("Red")
OPLOT, n, temp, linestyle=0, thick=5, color = FSC_COLOR("Blue")
XYOUTS, [0.55], [0.4], 'T-n', /normal, charthick=5, charsize = 1.5
;OPLOT, n2, temp2, linestyle=2, thick=5, color = FSC_COLOR("Blue")
OPLOT, nPR, tempPR, linestyle=3, thick=5, color = FSC_COLOR("Brown")
XYOUTS, [0.45], [0.6], 'T!dPR!n-n!dPR!n', /normal, charthick=5, charsize = 1.5

device, /close


n = SMOOTH(n, 1000.)
temp = SMOOTH(temp, 1000.)

n = SMOOTH(n, 1000.)
temp = SMOOTH(temp, 1000.)

n = SMOOTH(n, 100.)
temp = SMOOTH(temp, 100.)
;temp(where(temp lt tempPR)) = tempPR(wherE(temp lt tempPR))


;
;tC1 = 370.*(1.0)*(1.0e10/n)
;tB1 = 4.6e5*(1.0)^(0.2)*(1.0e10/n)
;tcool1 = tC1+tB1
;tC10 = 370.*(10.0)*(1.0e10/n)
;tB10 = 4.6e5*(10.)^(0.2)*(1.0e10/n)
;tcool10 = tC10+tB10
;tC100 = 370.*(100.0)*(1.0e10/n)
;tB100 = 4.6e5*(100.0)^(0.2)*(1.0e10/n)
;tcool100 = tC100+tB100
;tC1000 = 370.*(1000.0)*(1.0e10/n)
;tB1000 = 4.6e5*(1000.0)^(0.2)*(1.0e10/n)
;tcool1000 = tC1000+tB1000

;set_plot, 'x'
;plot, t, tcool1/t, /xlog, /ylog
;oplot, t, tcool10/t, linestyle=1
;oplot, t, tcool100/t, linestyle=2
;oplot, t, tcool1000/t, linestyle=3



;SAVE, filename = 'andrea', n, temp, t

;EcouleV = 5.0e5*((tn/240.)*(n/1.0e10))^(2./3.)

OPENW, 1, '4Andrea_floorSH.dat'
FOR i=0, tprec-1 DO BEGIN
  PRINTF, 1, t(i), n(i), temp(i), Lsh(i)
ENDFOR
CLOSE, 1


OPENW, 1, '4Andrea_PR.dat'
FOR i=0, tprec-1 DO BEGIN
  PRINTF, 1, t(i), nPR(i), tempPR(i)
ENDFOR
CLOSE, 1

STOP

delta(0) = R0/(9.2e3*(v8^(2.0)))
Vol(0) = 4.0*!pi*R0^(2.0)*delta(0)

hog = where((tempEOS gt temp(0)*0.97 and (tempEOS lt temp(0)*1.03)))
pig = where((nEOS gt n(0)*0.90) and (nEOS lt n(0)*1.1))
listmatch, hog, pig, dog, cat

gam(0) = gammaEOS(hog(dog[0]))
pr(0) = prEOS(hog(dog[0]))
eint(0) = eintEOS(hog(dog[0]))
Etot(0) = eint(0)*Vol(0)

hog = where((tempOP gt temp(0)*0.99 and (tempOP lt temp(0)*1.01)))
pig = where((nOP gt n(0)*0.95) and (nOP lt n(0)*1.05))
listmatch, hog, pig, dog, cat
kappa(0) = kappaOP(hog(dog[0]))
tau(0) = kappa(0)*delta(0)*n(0)*mp


STOP


n = dblarr(tprec)
gam = dblarr(tprec)
cs = dblarr(tprec)
pr = dblarr(tprec)
kappa = dblarr(tprec)
eint = dblarr(tprec)
Etot = dblarr(tprec)
temp = dblarr(tprec)
tau = dblarr(tprec)
Vol = dblarr(tprec)
delta = dblarr(tprec)
Lum = dblarr(tprec)
LumPdV = dblarr(tprec)

n(0) = n0
temp(0) = 1.0e4
delta(0) = R0/(9.2e3*(v8^(2.0)))
Vol(0) = 4.0*!pi*R0^(2.0)*delta(0)

hog = where((tempEOS gt temp(0)*0.97 and (tempEOS lt temp(0)*1.03)))
pig = where((nEOS gt n(0)*0.90) and (nEOS lt n(0)*1.1))
listmatch, hog, pig, dog, cat

gam(0) = gammaEOS(hog(dog[0]))
pr(0) = prEOS(hog(dog[0]))
eint(0) = eintEOS(hog(dog[0]))
Etot(0) = eint(0)*Vol(0)

hog = where((tempOP gt temp(0)*0.99 and (tempOP lt temp(0)*1.01)))
pig = where((nOP gt n(0)*0.95) and (nOP lt n(0)*1.05))
listmatch, hog, pig, dog, cat
kappa(0) = kappaOP(hog(dog[0]))
tau(0) = kappa(0)*delta(0)*n(0)*mp

STOP

FOR i=1, tprec-1 DO BEGIN

  hog = where((eintEOS gt eint(i-1)*0.9 and (eintEOS lt eint(i-1)*1.1)))
  pig = where((nEOS gt n(i-1)*0.90) and (nEOS lt n(i-1)*1.1))
  listmatch, hog, pig, dog, cat
  gam(i) = gammaEOS(hog(dog[0]))
  pr(i) = prEOS(hog(dog[0]))
  temp(i) = tempEOS(hog(dog[0]))
  cs(i) = (gam(i)*pr(i)/(n(i-1)*mp))^(0.5)

  
  hog = where((tempOP gt temp(i-1)*0.99 and (tempOP lt temp(i-1)*1.01)))
  pig = where((nOP gt n(i-1)*0.95) and (nOP lt n(i-1)*1.05))
  listmatch, hog, pig, dog, cat
  kappa(i) = kappaOP(hog(dog[0]))
  
  tau(i) = kappa(i)*delta(i-1)*n(i-1)*mp
  Lum(i) = 4.0*!pi*(R(i)^(2.0))*sigSB*(temp(i)^(4.0))/(tau(i) + 1.0/tau(i))  
  dVdt = 4.0*!pi*(2.0*R(i)*v*delta(i-1) + R(i)^(2.0)*cs(i)) 
  LumPdV(i) = 1.0*pr(i)*dVdt
  dEdt = -1.0*LumPdV(i) - 1.0*Lum(i)
  
  Etot(i) = Etot(i-1) + dEdt*(t(i)-t(i-1))
  delta(i) = delta(i-1) + cs(i)*(t(i)-t(i-1))
  Vol(i) = 4.0*!pi*R(i)^(2.0)*delta(i)
  eint(i) = Etot(i)/Vol(i)
  n(i) = n0*((R0/R(i))^(2.0))*(delta(0)/delta(i))
ENDFOR

STOP

set_plot, 'ps'
device, filename = 'densities.eps', /color

PLOT, [0.0], [0.0], /nodata, /xlog, /ylog, ystyle=1, yrange = [1.0e6, 1.0e15], color = FSC_COLOR("Black"), xmargin = [12,2], ymargin = [5,2], charthick=4., xcharsize = 1.5, ycharsize = 1.5, xthick=5., ythick=5., xstyle=1, xrange = [min(t),max(t)], thick=5, ytitle = 'n (/cc)', xtitle = 't (s)'
OPLOT, t, nmean, linestyle=1, thick=2, color = FSC_COLOR("Red")
OPLOT, t, n, linestyle=0, thick=5, color = FSC_COLOR("Blue")
;OPLOT, t, n2, linestyle=2, thick=5, color = FSC_COLOR("Blue")
OPLOT, t, nPR, linestyle=3, thick=5, color = FSC_COLOR("Brown")

device, /close


set_plot, 'ps'
device, filename = 'temperatures.eps', /color

PLOT, [0.0], [0.0], /nodata, /xlog, /ylog, ystyle=1, yrange = [1.0e2, 1.2e4], color = FSC_COLOR("Black"), xmargin = [12,2], ymargin = [5,2], charthick=4., xcharsize = 1.5, ycharsize = 1.5, xthick=5., ythick=5., xstyle=1, xrange = [min(t),max(t)], thick=5, ytitle = 'T (K)', xtitle = 't (s)'
;OPLOT, t, Trad, linestyle=1, thick=2, color = FSC_COLOR("Red")
OPLOT, t, temp, linestyle=0, thick=5, color = FSC_COLOR("Blue")
;OPLOT, t, temp2, linestyle=2, thick=5, color = FSC_COLOR("Blue")
OPLOT, t, TPR, linestyle=3, thick=5, color = FSC_COLOR("Brown")

device, /close


set_plot, 'ps'
device, filename = 'Trho.eps', /color

PLOT, [0.0], [0.0], /nodata, /xlog, /ylog, ystyle=1, yrange = [1.0e2, 1.2e4], color = FSC_COLOR("Black"), xmargin = [12,2], ymargin = [5,2], charthick=4., xcharsize = 1.5, ycharsize = 1.5, xthick=5., ythick=5., xstyle=1, xrange = [1.0e10,1.0e15], thick=5, ytitle = 'T (K)', xtitle = 'n (/cc)'
;OPLOT, nmean, Trad, linestyle=1, thick=2, color = FSC_COLOR("Red")
OPLOT, n, temp, linestyle=0, thick=5, color = FSC_COLOR("Blue")
;OPLOT, n2, temp2, linestyle=2, thick=5, color = FSC_COLOR("Blue")
OPLOT, nPR, TPR, linestyle=3, thick=5, color = FSC_COLOR("Brown")

device, /close


STOP
END