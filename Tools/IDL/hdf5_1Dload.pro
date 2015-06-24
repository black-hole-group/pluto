
;+
;
; NAME:      HDF5_1DLOAD
;
; AUTHOR:    A. Mignone, O. Tesileanu, C. Zanni
;
; REQUIRES:  HDF5 support.
;
; PURPOSE:   read a one-dimensional HDF5 file and store its content on a
;            dedicated structure called 'Data1D'. For each level variables 
;            will be stored sequentially using array concatenation, starting
;            from the finest level to the base one. 
;            One can then plot the data using the psym keyword, e.g.,
;
;            IDL> plot,Data1D.xx, Data1D.rho,psym=1
;            
;             
;
; SYNTAX:    HDF5_1DLOAD, nout, level, Data1D
;
; KEYWORDS:  None
;
; LAST MODIFIED
;
;   May 7, 2011 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO HDF5_1DLOAD, nout, maxlev, Data1D

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR

 qrho = [0.0] & qpr  = [0.0]
 qvx  = [0.0] & qvy  = [0.0] & qvz  = [0.0]
 qbx  = [0.0] & qby  = [0.0] & qbz  = [0.0]
 xx = [0.0]

 pload,/hdf5,dim=1,NOUT,lev=MAXLEV
 indx = WHERE(AMRBoxes EQ MAXLEV)

 qrho = [rho(indx)] & qpr = [pr(indx)]
 qvx  = [v1(indx)]  & qvy = [v2(indx)] & qvz  = [v3(indx)]
 qbx  = [b1(indx)]  & qby = [b2(indx)] & qbz  = [b3(indx)]
 xx = [x1(indx)]

 AMRBoxesMax = AMRBoxes
 xmax        = x1
 FOR n=MAXLEV-1,0,-1 DO BEGIN
   print,"Building lev ",n
   pload,/hdf5,dim=1,NOUT,lev=n
   indx = WHERE((AMRBoxes EQ n) AND $
                REBIN(AMRBoxesmax,n1) EQ n)
   qrho = [qrho,rho(indx)] & qpr = [qpr,pr(indx)]
   qvx  = [qvx, v1(indx)]  & qvy  = [qvy, v2(indx)] & qvz  = [qvz, v3(indx)]
   qbx  = [qbx, b1(indx)]  & qby  = [qby, b2(indx)] & qbz  = [qbz, b3(indx)]
   xx = [xx,x1(indx)]
 ENDFOR

 Data1D = {xx:xx, rho:qrho, pr: qpr, $
                  vx:qvx, vy:qvy, vz:qvz, $
                  bx:qbx, by:qby, bz:qbz}
END

