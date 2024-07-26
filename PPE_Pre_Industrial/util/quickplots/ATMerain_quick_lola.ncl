load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("define_draw_strings")
 procedure define_draw_strings (nlati,Code,pltdate,comment,wks)
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
  zonal_mean=dim_avg(Code)
  meanV=dim_avg_wgt(zonal_mean,gweight,1)

  minV=min(Code)
  maxV=max(Code)

  statV="min="+minV+"   mean="+meanV+"   max="+maxV

  txres               = True                    ; text mods desired
  txres@txFontHeightF = 0.014                   ; font smaller. default big
  gsn_text_ndc(wks,pltdate,0.85,0.1,txres)
  gsn_text_ndc(wks,comment,0.1,0.1,txres)
  txres@txFontHeightF = 0.016 
  gsn_text_ndc(wks,statV,0.5,0.22,txres)
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

;***define plot
   pltType = "png"                        ; x11, eps, ps, pdf, png
   pltPath = workdir+"/atmlola_"+typ+"_"+run+"_"

;************************************************
;   Temperature 850 hPa code130
;************************************************

if (isfilepresent("Ubusy_c130lola.nc")) then

;***read code 130
    fili    = workdir+"/Ubusy_c130lola.nc"
    f       = addfile (fili , "r") 
    c130     = f->var130(0,0,:,:)
    c130     = c130 
    delete(c130@long_name)                  ; don't draw long_name by default
    delete(c130@units)                      ; don't draw unit by default
    printVarSummary(c130)

;***open plot
    pltName = pltPath+"c130"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
   i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 


;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Temperature [C] 850 hPa "+typ
    res@gsnCenterString      = run+meantime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30/)
    res@cnFillColors         = (/2,4,5,6,7,8,9,10,11,12,14,15,16,17/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c130, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c130,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if
;************************************************
;   Temperature 850 hPa code130 (model-ERAin)
;************************************************

if (isfilepresent("Ubusy_d130lola.nc")) then

;***read code 130
    fili    = workdir+"/Ubusy_d130lola.nc"
    f       = addfile (fili , "r")
    c130     = f->var130(0,0,:,:)
    c130     = c130
    delete(c130@long_name)                  ; don't draw long_name by default
    delete(c130@units)                      ; don't draw unit by default
    printVarSummary(c130)

;***open plot
    pltName = pltPath+"c130_era"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName)

   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map
   i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map


;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Temperature [C] 850 hPa "+typ
    res@gsnCenterString      = run+meantime+" - Erain"+eratime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-16,-8,-4,-2,-1,1,2,4,8,16/)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0

    plot = gsn_csm_contour_map_ce(wks,c130, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c130,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;************************************************
;  Geopotential height  500 hPa code156 (model-ERAin)
;************************************************

if (isfilepresent("Ubusy_d156lola.nc")) then

;***read code 156
    fili    = workdir+"/Ubusy_d156lola.nc"
    f       = addfile (fili , "r") 
    c156     = f->var156(0,0,:,:)
    c156     = c156 
    delete(c156@long_name)                  ; don't draw long_name by default
    delete(c156@units)                      ; don't draw unit by default
    printVarSummary(c156)

;***open plot
    pltName = pltPath+"c156_era"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
   i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 


;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Geopotential height [dam] 500 hPa "+typ
    res@gsnCenterString      = run+meantime+" - Erain"+eratime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-10,-8,-4,-2,-1,1,2,4,8,10/)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c156, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c156,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;************************************************
;   Geopotential height [dam] 500 hPa code156
;************************************************

if (isfilepresent("Ubusy_c156lola.nc")) then

;***read code 156
    fili    = workdir+"/Ubusy_c156lola.nc"
    f       = addfile (fili , "r") 
    c156     = f->var156(0,0,:,:)
    c156     = c156 
    delete(c156@long_name)                  ; don't draw long_name by default
    delete(c156@units)                      ; don't draw unit by default
    printVarSummary(c156)

;***open plot
    pltName = pltPath+"c156"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 


;***create plot
    BWres                      = True              ; plot mods desired
    BWres@gsnFrame             = False             ; don't advance frame yet
    BWres@tiMainString         = "Geopotential height [dam] 500 hPa "+typ
    BWres@gsnCenterString      = run+meantime

    BWres@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
    BWres@cnMinLevelValF       =  480.             ; set min contour level
    BWres@cnMaxLevelValF       =  600.             ; set max contour level
    BWres@cnLevelSpacingF      =  16.              ; set contour spacing
    BWres@cnInfoLabelOn        = False             ; turn off contour label box

    BWres@mpCenterLonF         = 0                 ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c156, BWres) ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c156,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame


  end if
;************************************************
;   code149---Velocity potential 200 hPa code149
;************************************************

if (isfilepresent("Ubusy_c149lola.nc")) then

;***read code 149
    fili    = workdir+"/Ubusy_c149lola.nc"
    f       = addfile (fili , "r") 
    c149     = f->var149(0,0,:,:)
    c149     = c149 
    delete(c149@long_name)                  ; don't draw long_name by default
    delete(c149@units)                      ; don't draw unit by default
    printVarSummary(c149)

;***open plot
    pltName = pltPath+"c149"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 


;***create plot
    BWres                      = True              ; plot mods desired
    BWres@gsnFrame             = False             ; don't advance frame yet
    BWres@tiMainString         = "Velocity potential [10~S~6~N~m~S~2~N~/s] 200 hPa "+typ
    BWres@gsnContourZeroLineThicknessF = 2. 	  ; thickness of zero contour
    BWres@gsnContourNegLineDashPattern = 1        ; negative contours to dash 
    BWres@gsnCenterString      = run+meantime

    BWres@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
    BWres@cnMinLevelValF       =  -20.             ; set min contour level
    BWres@cnMaxLevelValF       =  20.             ; set max contour level
    BWres@cnLevelSpacingF      =  2.              ; set contour spacing
    BWres@cnInfoLabelOn        = False             ; turn off contour label box

    BWres@mpCenterLonF         = 0                 ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c149, BWres) ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c149,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame


  end if
end
