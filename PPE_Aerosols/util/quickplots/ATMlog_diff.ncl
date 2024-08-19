load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

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

  values    = asciiread("var.txt",8,"string")
  typ = values(0)
  runA = values(1)
  runB = values(2)
  meantime = values(3)
  comment = values(4)
  workdir = values(5)        ; plot & input directory
  nlati = stringtointeger(values(6))
  pltdate = values(7)       ; plot date

  print(workdir)

;***define plot
   pltType = "png"                        ; x11, eps, ps, pdf, png
   pltPath = workdir+"/atmlog_"+typ+"_"+runA+"_"+runB+"_"

;************************************************
;   Temperature code130 zonal 
;************************************************

if (isfilepresent("Ubusy_d130.nc")) then

;***read code 130
    fili    = workdir+"/Ubusy_d130.nc"
    f       = addfile (fili , "r") 
    d130     = f->var130(0,:,:,0)
    d130     = d130 
    delete(d130@long_name)                  ; don't draw long_name by default
    delete(d130@units)                      ; don't draw unit by default
    printVarSummary(d130)

;***open plot
    pltName = pltPath+"c130"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@gsnFrame             = False            ; don't advance frame yet
    res@tiMainString         = "Temperature [C] "+typ
    res@tiYAxisString        = " "               ; Y-Axis no label
    res@gsnCenterString      = meantime

    res@trYLog   = True           ; no logarithmic Y-axis
    res@tmXMajorGrid      = True                 ; X Axis grid lines
    res@tmXMajorGridLineColor = "grey"           ; color X Axis grid lines
    res@tmYMajorGrid      = True                 ; Y Axis grid lines
    res@tmYMajorGridLineColor = "grey"           ; color Y Axis grid lines
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    res@tmYLLabels = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    res@tmYRMode = "Automatic"

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.07           ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-8,-4,-2,-1,-.5,.5,1,2,4,8/)
    res@cnFillColors         = (/2,4,6,7,9,0,11,12,14,15,17/)

    res@vpWidthF  = 0.8		; set width of plot 
    res@vpHeightF = 0.5		; set height of plot

    plot = gsn_csm_pres_hgt(wks,d130, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_zon_draw_strings (nlati,d130,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if


;************************************************
;   U code131 zonal model-ERA40
;************************************************

if (isfilepresent("Ubusy_d131.nc")) then

;***read code 131
    fili    = workdir+"/Ubusy_d131.nc"
    f       = addfile (fili , "r") 
    d131     = f->var131(0,:,:,0)
    d131     = d131 
    delete(d131@long_name)                  ; don't draw long_name by default
    delete(d131@units)                      ; don't draw unit by default
    printVarSummary(d131)

;***open plot
    pltName = pltPath+"c131"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@gsnFrame             = False            ; don't advance frame yet
    res@tiMainString         = "U [m/s] "+typ
    res@tiYAxisString       = " "               ; Y-Axis no label
    res@gsnCenterString      = meantime

    res@trYLog   = True           ; no logarithmic Y-axis
    res@tmXMajorGrid      = True                 ; X Axis grid lines
    res@tmXMajorGridLineColor = "grey"           ; color X Axis grid lines
    res@tmYMajorGrid      = True                 ; Y Axis grid lines
    res@tmYMajorGridLineColor = "grey"           ; color Y Axis grid lines
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    res@tmYLLabels = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    res@tmYRMode = "Automatic"

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.07           ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-8,-4,-2,-1,-.5,.5,1,2,4,8/)
    res@cnFillColors         = (/2,4,6,7,9,0,11,12,14,15,17/)


    res@vpWidthF  = 0.8		; set width of plot 
    res@vpHeightF = 0.5		; set height of plot

    plot = gsn_csm_pres_hgt(wks,d131, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_zon_draw_strings (nlati,d131,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;************************************************
;   V code132 zonal 
;************************************************

if (isfilepresent("Ubusy_d132.nc")) then

;***read code 132
    fili    = workdir+"/Ubusy_d132.nc"
    f       = addfile (fili , "r") 
    d132     = f->var132(0,:,:,0)
    d132     = d132 
    delete(d132@long_name)                  ; don't draw long_name by default
    delete(d132@units)                      ; don't draw unit by default
    printVarSummary(d132)

;***open plot
    pltName = pltPath+"c132"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@gsnFrame             = False            ; don't advance frame yet
    res@tiMainString         = "V [m/s] "+typ
    res@tiYAxisString       = " "               ; Y-Axis no label
    res@gsnCenterString      = meantime

    res@trYLog   = True           ; no logarithmic Y-axis
    res@tmXMajorGrid      = True                 ; X Axis grid lines
    res@tmXMajorGridLineColor = "grey"           ; color X Axis grid lines
    res@tmYMajorGrid      = True                 ; Y Axis grid lines
    res@tmYMajorGridLineColor = "grey"           ; color Y Axis grid lines
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    res@tmYLLabels = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    res@tmYRMode = "Automatic"

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.07           ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)


    res@vpWidthF  = 0.8		; set width of plot 
    res@vpHeightF = 0.5		; set height of plot

    plot = gsn_csm_pres_hgt(wks,d132, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_zon_draw_strings (nlati,d132,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame
    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if


end
