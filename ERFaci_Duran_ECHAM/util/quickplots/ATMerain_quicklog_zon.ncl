load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("grow_up_2d")
function grow_up_2d(x)
begin

  msg = 1e20    ; missing value

; Get dimension sizes of input array
  dims    = dimsizes(x)
  nlev    = dims(0)
  nlev1   = nlev + 1

; Create new variable to hold larger array
  newdims    = dims
  newdims(0) = nlev1
  xnew = new(newdims,typeof(x))

; Copy attributes from old to new
  copy_VarAtts(x,xnew)

; Copy dimension names from old to new
  xnew!0 = x!0
  xnew!1 = x!1

; Get the level array and increase by one
  lev = x&$x!0$
  lev1 = new(nlev1,typeof(lev))

; Copy the old levels to new variable
  lev1(0:nlev-1) = lev
  lev1(nlev)     = 1.

; Copy coordinate variable
  xnew&$x!0$ = lev1
  xnew&$x!1$ = x&$x!1$

; Copy old x to new x
  xnew(0:nlev-1,:) = x
  xnew(nlev,:)     = msg

  xnew@_FillValue = msg
  return(xnew)

end

undef("BWdefine_zon_draw_strings")
 procedure BWdefine_zon_draw_strings (nlati,Code,pltdate,comment,wks)
;************************************************
;  Define min max mean for black and white plot
;  Draw text on plot using NDC coordinates.
;************************************************
local txres

begin
;***compute weights for mean
  nlats=nlati/2
  gs=doubletofloat(gaus(nlats))
  gweight=gs(:,1)
;***compute mean 
  meanV=dim_avg_wgt(Code,gweight,1)

  minV=min(Code)
  maxV=max(Code)
  meanV=dim_avg(meanV)

  statV="min="+minV+"   mean="+meanV(0)+"   max="+maxV

  txres               = True                    ; text mods desired
  txres@txFontHeightF = 0.014                   ; font smaller. default big
  gsn_text_ndc(wks,pltdate,0.9,0.1,txres)
  gsn_text_ndc(wks,comment,0.17,0.1,txres)
  txres@txFontHeightF = 0.016 
  gsn_text_ndc(wks,statV,0.5,0.22,txres)
end
;***end procedure 
undef("define_zon_draw_strings")
 procedure define_zon_draw_strings (nlati,Code,pltdate,comment,wks)
;************************************************
;  Define min max mean
;  Draw text on plot using NDC coordinates.
;************************************************
local txres

begin
;***compute weights for mean
  nlats=nlati/2
  gs=doubletofloat(gaus(nlats))
  gweight=gs(:,1)
;***compute mean 
  meanV=dim_avg_wgt(Code,gweight,1)

  minV=min(Code)
  maxV=max(Code)
  meanV=dim_avg(meanV)

  statV="min="+minV+"   mean="+meanV(0)+"   max="+maxV

  txres               = True                    ; text mods desired
  txres@txFontHeightF = 0.014                   ; font smaller. default big
  gsn_text_ndc(wks,pltdate,0.93,0.08,txres)
  gsn_text_ndc(wks,comment,0.17,0.08,txres)
  txres@txFontHeightF = 0.016 
  gsn_text_ndc(wks,statV,0.6,0.22,txres)
end
;***end procedure 

;************************************************
begin

  values    = asciiread("var.txt",9,"string")
  name = values(0)
  typ = values(1)
  run = values(2)
  meantime = values(3)
  comment = values(4)
  workdir = values(5)        ; plot & input directory
  nlati = stringtointeger(values(6))
  eratime = values(7)
  pltdate = values(8)       ; plot date

  print(workdir)
  print(eratime)

;***define plot
   pltType = "png"                        ; x11, eps, ps, pdf, png
   pltPath = workdir+"/atmlog_"+typ+"_"+run+"_"

;************************************************
;   Temperature code130 zonal
;************************************************

if (isfilepresent("Ubusy_c130.nc")) then

;***read code 130
    fili    = workdir+"/Ubusy_c130.nc"
    f       = addfile (fili , "r") 
    c130     = f->var130(0,:,:,0)
    c130     = c130 
    delete(c130@long_name)                  ; don't draw long_name by default
    delete(c130@units)                      ; don't draw unit by default
    printVarSummary(c130)

;***open plot
    pltName = pltPath+"c130"                      ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@gsnFrame             = False            ; don't advance frame yet
    BWres@tiMainString         = "Temperature [C] "+typ
    BWres@tiYAxisString       = " "               ; Y-Axis no label
    BWres@gsnCenterString      = run+meantime
    BWres@gsnContourZeroLineThicknessF = 2. 	  ; thickness of zero contour
    BWres@gsnContourNegLineDashPattern = 1        ; negative contours to dash 
    BWres@trYLog   = True                         ;  logarithmic Y-axis
    BWres@tmXMajorGrid      = True                ; X Axis grid lines
    BWres@tmXMajorGridLineColor = "grey"          ; color X Axis grid lines
    BWres@tmYMajorGrid      = True                ; Y Axis grid lines
    BWres@tmYMajorGridLineColor = "grey"          ; color Y Axis grid lines
    BWres@tmYLMode = "Explicit"
;    BWres@tmYLValues = (/0.01,0.04,0.1,0.2,0.5,1,10,50,100,500,1000/)
;    BWres@tmYLLabels = (/0.01,0.04,0.1,0.2,0.5,1,10,50,100,500,1000/) 
    BWres@tmYLValues = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    BWres@tmYLLabels = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    BWres@tmYRMode = "Automatic"

;    BWres@tmYRMode = "Explicit"
;    BWres@tmYUseLeft = False
;    BWres@tmYRValues = (/15578,13239,11066,9102,7406,5964,4752,3743,2914,2235,1685,1245,901,637,440,296,193,122,74,43,23,11,4,1/)/100.
;    BWres@tmYRLabels = "" + (/15578,13239,11066,9102,7406,5964,4752,3743,2914,2235,1685,1245,901,637,440,296,193,122,74,43,23,11,4,1/)
;    BWres@tmYRValues = (/10000,1000,100/)
;    BWres@tmYRLabels = "" + (/10000,1000,100/)

    BWres@cnInfoLabelOn        = False            ; turn off contour label box
    BWres@cnLineLabelFormat    = "#f"             ; no unnecessary 0
    BWres@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    BWres@cnMinLevelValF       =   180.          ; set min contour level
    BWres@cnMaxLevelValF       =   300.          ; set max contour level
    BWres@cnLevelSpacingF      =      5.          ; set contour spacing

    BWres@vpWidthF  = 0.8		; set width of plot 
    BWres@vpHeightF = 0.5		; set height of plot


    plot = gsn_csm_pres_hgt(wks,c130, BWres)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    BWdefine_zon_draw_strings (nlati,c130,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
    delete (wks)
;;    delete(BWres@cnLevels)
  end if

;************************************************
;   Temperature code130 zonal model-ERAin
;************************************************

if (isfilepresent("Ubusy_d130.nc")) then

;***read code 130
    fili    = workdir+"/Ubusy_d130.nc"
    f       = addfile (fili , "r") 
    d130    = f->var130(0,:,:,0)
    d130new = grow_up_2d(d130)
    
    delete(d130new@long_name)                  ; don't draw long_name by default
    delete(d130new@units)                      ; don't draw unit by default
    printVarSummary(d130new)

;***open plot
    pltName = pltPath+"c130_era"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@gsnFrame             = False            ; don't advance frame yet
    res@tiMainString         = "Temperature [C] "+typ
    res@tiYAxisString        = " "              ; Y-Axis no label
    res@gsnCenterString      = run+meantime+" - Erain"+eratime

    res@trYLog                = True            ; logarithmic Y-axis
    res@tmXMajorGrid          = True            ; X Axis grid lines
    res@tmXMajorGridLineColor = "grey"          ; color X Axis grid lines
    res@tmYMajorGrid          = True            ; Y Axis grid lines
    res@tmYMajorGridLineColor = "grey"          ; color Y Axis grid lines
    res@tmYLMode = "Explicit"
;;    res@tmYLValues = (/100,200,300,400,500,600,700,800,900,1000/)
;;    res@tmYLLabels = (/100,200,300,400,500,600,700,800,900,1000/) 
    res@tmYLValues = (/0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    res@tmYLLabels = (/0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    res@tmYRMode = "Automatic"

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.07           ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-16,-8,-4,-2,-1,1,2,4,8,16/)
    res@cnFillColors         = (/2,4,6,7,9,0,11,12,14,15,17/)

    res@vpWidthF  = 0.8		; set width of plot 
    res@vpHeightF = 0.5		; set height of plot

    plot = gsn_csm_pres_hgt(wks,d130new, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_zon_draw_strings (nlati,d130new,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if


;************************************************
;   U code131 zonal
;************************************************

if (isfilepresent("Ubusy_c131.nc")) then

;***read code 131
    fili    = workdir+"/Ubusy_c131.nc"
    f       = addfile (fili , "r") 
    c131     = f->var131(0,:,:,0)
    
    delete(c131@long_name)                  ; don't draw long_name by default
    delete(c131@units)                      ; don't draw unit by default
    printVarSummary(c131)

;***open plot
    pltName = pltPath+"c131"                    ; plot name
    wks     = gsn_open_wks(pltType, pltName) 


;***create plot
    BWres                      = True             ; plot mods desired
    BWres@gsnFrame             = False            ; don't advance frame yet
    BWres@tiMainString         = "U [m/s] "+typ
    BWres@tiYAxisString       = " "               ; Y-Axis no label
    BWres@gsnCenterString      = run+meantime
    BWres@gsnContourZeroLineThicknessF = 2. 	  ; thickness of zero contour
    BWres@gsnContourNegLineDashPattern = 1        ; negative contours to dash 
    BWres@trYLog   = True                         ; logarithmic Y-axis
    BWres@tmXMajorGrid      = True                ; X Axis grid lines
    BWres@tmXMajorGridLineColor = "grey"          ; color X Axis grid lines
    BWres@tmYMajorGrid      = True                ; Y Axis grid lines
    BWres@tmYMajorGridLineColor = "grey"          ; color Y Axis grid lines
    BWres@tmYLMode = "Explicit"
;;    BWres@tmYLValues = (/100,200,300,400,500,600,700,800,900,1000/)
;;    BWres@tmYLLabels = (/100,200,300,400,500,600,700,800,900,1000/) 
;    BWres@tmYLValues = (/0.01,0.04,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
;    BWres@tmYLLabels = (/0.01,0.04,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    BWres@tmYLValues = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    BWres@tmYLLabels = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    BWres@tmYRMode = "Automatic"

    BWres@cnInfoLabelOn        = False            ; turn off contour label box
    BWres@cnLineLabelFormat    = "#f"             ; no unnecessary 0
    BWres@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    BWres@cnMinLevelValF       =   -50.           ; set min contour level
    BWres@cnMaxLevelValF       =    50.           ; set max contour level
    BWres@cnLevelSpacingF      =     5.           ; set contour spacing

    BWres@vpWidthF  = 0.8		; set width of plot 
    BWres@vpHeightF = 0.5		; set height of plot


    plot = gsn_csm_pres_hgt(wks,c131, BWres)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    BWdefine_zon_draw_strings (nlati,c131,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
    delete(BWres@cnLevels)
  end if

;************************************************
;   U code131 zonal model-ERA40
;************************************************

if (isfilepresent("Ubusy_d131.nc")) then

;***read code 131
    fili    = workdir+"/Ubusy_d131.nc"
    f       = addfile (fili , "r") 
    d131     = f->var131(0,:,:,0)
    d131new = grow_up_2d(d131)

    delete(d131new@long_name)                  ; don't draw long_name by default
    delete(d131new@units)                      ; don't draw unit by default
    printVarSummary(d131new)

;***open plot
    pltName = pltPath+"c131_era"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@gsnFrame             = False            ; don't advance frame yet
    res@tiMainString         = "U [m/s] "+typ
    res@tiYAxisString       = " "               ; Y-Axis no label
    res@gsnCenterString      = run+meantime+" - Erain"+eratime

    res@trYLog   = True           ; no logarithmic Y-axis
    res@tmXMajorGrid      = True                 ; X Axis grid lines
    res@tmXMajorGridLineColor = "grey"           ; color X Axis grid lines
    res@tmYMajorGrid      = True                 ; Y Axis grid lines
    res@tmYMajorGridLineColor = "grey"           ; color Y Axis grid lines
    res@tmYLMode = "Explicit"
;    res@tmYLValues = (/100,200,300,400,500,600,700,800,900,1000/)
;    res@tmYLLabels = (/100,200,300,400,500,600,700,800,900,1000/) 
    res@tmYLValues = (/0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    res@tmYLLabels = (/0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    res@tmYRMode = "Automatic"

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.07           ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-16,-8,-4,-2,-1,1,2,4,8,16/)
    res@cnFillColors         = (/2,4,6,7,9,0,11,12,14,15,17/)


    res@vpWidthF  = 0.8		; set width of plot 
    res@vpHeightF = 0.5		; set height of plot

    plot = gsn_csm_pres_hgt(wks,d131new, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_zon_draw_strings (nlati,d131new,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;************************************************
;   V code132 zonal
;************************************************

if (isfilepresent("Ubusy_c132.nc")) then

;***read code 132
    fili    = workdir+"/Ubusy_c132.nc"
    f       = addfile (fili , "r") 
    c132     = f->var132(0,:,:,0)
    
    delete(c132@long_name)                  ; don't draw long_name by default
    delete(c132@units)                      ; don't draw unit by default
    printVarSummary(c132)

;***open plot
    pltName = pltPath+"c132"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 


;***create plot
    BWres                      = True             ; plot mods desired
    BWres@gsnFrame             = False            ; don't advance frame yet
    BWres@tiMainString         = "V [m/s] "+typ
    BWres@tiYAxisString       = " "               ; Y-Axis no label
    BWres@gsnCenterString      = run+meantime
    BWres@gsnContourZeroLineThicknessF = 2. 	  ; thickness of zero contour
    BWres@gsnContourNegLineDashPattern = 1        ; negative contours to dash 
    BWres@trYLog   = True                         ; logarithmic Y-axis
    BWres@tmXMajorGrid      = True                ; X Axis grid lines
    BWres@tmXMajorGridLineColor = "grey"          ; color X Axis grid lines
    BWres@tmYMajorGrid      = True                ; Y Axis grid lines
    BWres@tmYMajorGridLineColor = "grey"          ; color Y Axis grid lines
    BWres@tmYLMode = "Explicit"
;    BWres@tmYLValues = (/100,200,300,400,500,600,700,800,900,1000/)
;    BWres@tmYLLabels = (/100,200,300,400,500,600,700,800,900,1000/) 
;;    BWres@tmYLValues = (/0.01,0.04,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
;;    BWres@tmYLLabels = (/0.01,0.04,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    BWres@tmYLValues = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    BWres@tmYLLabels = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    BWres@tmYRMode = "Automatic"

    BWres@cnInfoLabelOn        = False            ; turn off contour label box
    BWres@cnLineLabelFormat    = "#f"
    BWres@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    BWres@cnMinLevelValF       =   -10.           ; set min contour level
    BWres@cnMaxLevelValF       =    10.           ; set max contour level
    BWres@cnLevelSpacingF      =     0.5          ; set contour spacing

    BWres@vpWidthF  = 0.8		; set width of plot 
    BWres@vpHeightF = 0.5		; set height of plot


    plot = gsn_csm_pres_hgt(wks,c132, BWres)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    BWdefine_zon_draw_strings (nlati,c132,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
;;    delete(BWres@cnLevels)
  end if

;************************************************
;   V code132 zonal model-ERAin
;************************************************

if (isfilepresent("Ubusy_d132.nc")) then

;***read code 132
    fili    = workdir+"/Ubusy_d132.nc"
    f       = addfile (fili , "r") 
    d132     = f->var132(0,:,:,0)
    d132new = grow_up_2d(d132)

    delete(d132new@long_name)                  ; don't draw long_name by default
    delete(d132new@units)                      ; don't draw unit by default
    printVarSummary(d132)

;***open plot
    pltName = pltPath+"c132_era"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@gsnFrame             = False            ; don't advance frame yet
    res@tiMainString         = "V [m/s] "+typ
    res@tiYAxisString       = " "               ; Y-Axis no label
    res@gsnCenterString      = run+meantime+" - Erain"+eratime

    res@trYLog   = True           ; no logarithmic Y-axis
    res@tmXMajorGrid      = True                 ; X Axis grid lines
    res@tmXMajorGridLineColor = "grey"           ; color X Axis grid lines
    res@tmYMajorGrid      = True                 ; Y Axis grid lines
    res@tmYMajorGridLineColor = "grey"           ; color Y Axis grid lines
    res@tmYLMode = "Explicit"
;    res@tmYLValues = (/100,200,300,400,500,600,700,800,900,1000/)
;    res@tmYLLabels = (/100,200,300,400,500,600,700,800,900,1000/) 
    res@tmYLValues = (/0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    res@tmYLLabels = (/0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    res@tmYRMode = "Automatic"

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.07           ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1./)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)


    res@vpWidthF  = 0.8		; set width of plot 
    res@vpHeightF = 0.5		; set height of plot

    plot = gsn_csm_pres_hgt(wks,d132new, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_zon_draw_strings (nlati,d132new,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

end
