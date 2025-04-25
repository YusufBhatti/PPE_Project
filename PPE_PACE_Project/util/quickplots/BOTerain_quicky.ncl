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
  gsn_text_ndc(wks,comment,0.18,0.1,txres)
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

;  print(eratime)
;  print(workdir)
;  print(pltdate)

;***define plot
   pltType = "png"                        ; x11, eps, ps, pdf, png
   pltPath = workdir+"/bot_"+typ+"_"+run+"_"

;************************************************
;   Sea Level Pressure code151
;************************************************

if (isfilepresent("Ubusy_c151.nc")) then

;***read code 151
    fili    = workdir+"/Ubusy_c151.nc"
    f       = addfile (fili , "r") 
    c151     = f->var151(0,:,:)
    c151     = c151 - 1000.
    delete(c151@long_name)                  ; don't draw long_name by default
    delete(c151@units)                      ; don't draw unit by default
    printVarSummary(c151)

;***open plot
    pltName = pltPath+"c151"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
   i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 


;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Sea Level Pressure (hPa-1000) "+typ
    res@gsnCenterString      = run+meantime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels            = (/-25,-20,-15,-10,-5,0,5,10,15,20,25/)
    res@cnFillColors    = (/2,4,6,7,8,9,10,11,12,14,15,16/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c151, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c151,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;************************************************
;   Sea Level Pressure code151 (model-ERAin)
;************************************************

if (isfilepresent("Ubusy_d151.nc")) then

;***read code 151 (model -ERAin)
    fili     =  workdir+"/Ubusy_d151.nc"
    print(fili)

    f       = addfile (fili , "r")  
;   printVarSummary(f)
 
    d151     = f->var151(0,:,:)       
    delete(d151@long_name)            ; don't draw long_name by default
    delete(d151@units)                ; don't draw unit by default

;***open plot
    pltName = pltPath+"c151_era"              ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Sea Level Pressure (hPa-1000)"+typ
    res@gsnCenterString      = run+meantime+" - Erain"+eratime

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-10,-5,-2,-1,1,2,5,10/)
;    res@cnFillColors         = (/4,6,7,9,0,11,12,14,15,16/)
    res@cnFillColors         = (/4,6,8,9,0,11,12,14,15,17/)
    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,d151, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,d151,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if

;************************************************
;   Total cloud cover code164 
;************************************************

if (isfilepresent("Ubusy_c164.nc")) then

;***read code 164
    fili     =  workdir+"/Ubusy_c164.nc"
    print(fili)

    f       = addfile (fili , "r")  
   printVarSummary(f)
 
    c164     = f->var164(0,:,:)       
    delete(c164@long_name)            ; don't draw long_name by default

;***open plot
    pltName = pltPath+"c164"              ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 
    i = NhlNewColor(wks,0.6,0.6,0.6)               
    i = NhlNewColor(wks,0.4,0.4,0.4)  
             
;***create plot
;***resource list for color contour leves
    res1                      = True
    res1@gsnFrame             = False             ; don't advance frame
    res1@tiMainString         = "Total cloud cover [%] "+typ
    res1@gsnCenterString      = run+meantime

    res1@cnLevelSelectionMode = "ExplicitLevels"  ; set explicit contour levels
    res1@cnLevels             = (/10,30,50,70,90/)
    res1@cnFillColors         = (/3,6,9,18,19,20/)


    res1@pmLabelBarOrthogonalPosF = 0.2           ; move farther to plot
    res1@lbOrientation        = "horizontal"      ; vertical label bar
    res1@gsnSpreadColors      = False             ; use full range of colors

    res1@mpCenterLonF         = 0.                ; center plot at lon 0

;***resource list for contour lines
    res2                      = True
    res2@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
    res2@cnMinLevelValF       =  10.              ; set min contour level
    res2@cnMaxLevelValF       =  90.              ; set max contour level
    res2@cnLevelSpacingF      =  10.              ; set contour spacing
    res2@cnLineLabelsOn       = False             ; no contour labels
    res2@cnInfoLabelOn        = False             ; turn off contour label box

    plot = gsn_csm_contour_map_overlay(wks,c164,c164,res1,res2)

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c164,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

end if
;************************************************
;   2 m Temperature code167 (model-ERAin)
;************************************************

if (isfilepresent("Ubusy_d167.nc")) then

;***read code 167 (model -ERAin)
    fili     =  workdir+"/Ubusy_d167.nc"
    print(fili)

    f       = addfile (fili , "r")  
;   printVarSummary(f)
 
    d167     = f->var167(0,:,:)       
    delete(d167@long_name)            ; don't draw long_name by default
    delete(d167@units)                ; don't draw unit by default

;***open plot
    pltName = pltPath+"c167_era"              ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "2 m Temperature [C] "+typ
    res@gsnCenterString      = run+meantime+" - Erain"+eratime

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-5,-2,-1,-0.5,0.5,1,2,5/)
;    res@cnFillColors         = (/4,6,7,9,0,11,12,14,15,16/)
    res@cnFillColors         = (/4,6,8,9,0,11,12,14,15,17/)
    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,d167, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,d167,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if
;************************************************
;   Surface temperature code169
;************************************************

if (isfilepresent("Ubusy_c169.nc")) then

;***read code 169
    fili    = workdir+"/Ubusy_c169.nc"
    f       = addfile (fili , "r") 
    c169     = f->var169(0,:,:)
    c169     = c169 
    delete(c169@long_name)                  ; don't draw long_name by default
    delete(c169@units)                      ; don't draw unit by default
    printVarSummary(c169)

;***open plot
    pltName = pltPath+"c169"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
   i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 


;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Surface temperature [C] "+typ
    res@gsnCenterString      = run+meantime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels            = (/-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30/)
    res@cnFillColors    = (/2,4,5,6,7,8,9,10,11,12,14,15,16,17/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c169, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c169,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if
;************************************************
;   Zonal wind stress code180
;************************************************

if (isfilepresent("Ubusy_c180.nc")) then

;***read code 180
    fili    = workdir+"/Ubusy_c180.nc"
    f       = addfile (fili , "r") 
    c180     = f->var180(0,:,:)
    c180     = c180 
    delete(c180@long_name)                  ; don't draw long_name by default
    delete(c180@units)                      ; don't draw unit by default
    printVarSummary(c180)

;***open plot
    pltName = pltPath+"c180"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
   i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 


;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Zonal wind stress [mN/m~S~2~N~] "+typ
    res@gsnCenterString      = run+meantime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels            = (/-200,-100,-50,-20,0,20,50,100,200/)
    res@cnFillColors    = (/2,4,6,8,9,10,11,12,14,15,16/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c180, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c180,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if
;************************************************
;   Zonal wind stress code180  ERAin
;************************************************

if (isfilepresent("Ubusy_c180e.nc")) then

;***read code 180
    fili    = workdir+"/Ubusy_c180e.nc"
    f       = addfile (fili , "r") 
    c180     = f->var180(0,:,:)
    c180     = c180 
    delete(c180@long_name)                  ; don't draw long_name by default
    delete(c180@units)                      ; don't draw unit by default
    printVarSummary(c180)

;***open plot
    pltName = pltPath+"c180_era"                     ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
   i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 


;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Zonal wind stress [mN/m~S~2~N~] "+typ
    res@gsnCenterString      = "Erain"+eratime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels            = (/-200,-100,-50,-20,0,20,50,100,200/)
    res@cnFillColors    = (/2,4,6,8,9,10,11,12,14,15,16/)
    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c180, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c180,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if
;************************************************
;   Total precipitation code4 
;************************************************

if (isfilepresent("Ubusy_c4.nc")) then

;***read code 4
    fili     =  workdir+"/Ubusy_c4.nc"
    print(fili)

    f       = addfile (fili , "r")  
   printVarSummary(f)
 
    c4     = f->var4(0,:,:)       
    delete(c4@long_name)            ; don't draw long_name by default
    delete(c4@units)                ; don't draw unit by default

;***open plot
    pltName = pltPath+"c004"              ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 
    i = NhlNewColor(wks,0.6,0.6,0.6)               
    i = NhlNewColor(wks,0.4,0.4,0.4)  
             
;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Total precipitation [mm/d] "+typ
    res@gsnCenterString      = run+meantime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels            = (/1,2,4,8,12/)
    res@cnFillColors    = (/0,18,9,6,3,17/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c4, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c4,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)

end if
;************************************************
;   Liquid water + ice content code231_150 
;************************************************

if (isfilepresent("Ubusy_c231_150.nc")) then

;***read code 231_150
    fili     =  workdir+"/Ubusy_c231_150.nc"
    print(fili)

    f       = addfile (fili , "r")  
   printVarSummary(f)
 
    c231     = f->var231(0,:,:)       
    delete(c231@long_name)            ; don't draw long_name by default
    delete(c231@units)                ; don't draw unit by default

;***open plot
    pltName = pltPath+"c231"              ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 
    i = NhlNewColor(wks,0.6,0.6,0.6)               
    i = NhlNewColor(wks,0.4,0.4,0.4)  
             
;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Liquid water + ice content [g/m~S~2~N~] "+typ
    res@gsnCenterString      = run+meantime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels            = (/10,50,100,200,300/)
    res@cnFillColors    = (/3,6,8,18,19,20/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c231, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c231,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)

end if


;************************************************
;   Column water vapour code230 
;************************************************

if (isfilepresent("Ubusy_c230.nc")) then

;***read code 230 
   fili     =  workdir+"/Ubusy_c230.nc"
   print(fili)

   f       = addfile (fili , "r")  
   printVarSummary(f)
 
   c230     = f->var230(0,:,:)       
   delete(c230@long_name)            ; don't draw long_name by default
   delete(c230@units)                ; don't draw unit by default

;***open plot
    pltName = pltPath+"c230"              ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
     res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Column water vapour [kg/m~S~2~N~] "+typ
    res@gsnCenterString      = run+meantime

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/10,20,30,40,50/)
    res@cnFillColors         = (/0,18,9,6,3,17/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c230, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c230,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if
 

;************************************************
;   Column water vapour code230 (ERAin=code137) (Model-ERAin)
;************************************************

if (isfilepresent("Ubusy_d230.nc")) then

;***read code 230 (model -ERAin)
    fili     =  workdir+"/Ubusy_d230.nc"
    print(fili)

    f       = addfile (fili , "r")  
;   printVarSummary(f)
 
    d230     = f->var230(0,:,:)       
    delete(d230@long_name)            ; don't draw long_name by default
    delete(d230@units)                ; don't draw unit by default

;***open plot
    pltName = pltPath+"c230_era"              ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 

;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Column water vapour [kg/m~S~2~N~] "+typ
    res@gsnCenterString      = run+meantime+" - Erain"+eratime

    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-10,-5,-2,-1,1,2,5,10/)
    res@cnFillColors         = (/4,6,7,9,0,11,12,14,15,16/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,d230, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,d230,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if

;************************************************
;   code140/code229 --only land----Ws/Wsmax  
;************************************************

if (isfilepresent("Ubusy_c140.nc")) then

;***read code 140
    fili     =  workdir+"/Ubusy_c140.nc"
    print(fili)

    f       = addfile (fili , "r")  
   printVarSummary(f)
 
    c140     = f->var140(0,:,:)       
    delete(c140@long_name)            ; don't draw long_name by default
    delete(c140@units)                ; don't draw unit by default

;***open plot
    pltName = pltPath+"c140"              ; plot name
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 
    i = NhlNewColor(wks,0.6,0.6,0.6)               
    i = NhlNewColor(wks,0.4,0.4,0.4)  
             
;***create plot
    res                      = True               ; plot mods desired
    res@gsnFrame             = False              ; don't advance frame yet
    res@tiMainString         = "Ws/Wsmax [%] "+typ
    res@gsnCenterString      = run+meantime
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.2            ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels            = (/10,30,50,70,90/)
    res@cnFillColors    = (/18,9,8,7,4,2/)

    res@mpFillOn             = False              ; turn off gray continents
    res@mpCenterLonF         = 0                  ; Centers the plot at lon 0 

    plot = gsn_csm_contour_map_ce(wks,c140, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (nlati,c140,pltdate,comment,wks)

    frame(wks)                                    ; now advance frame

    delete(res@cnLevels)
    delete (res@cnFillColors)

end if
end
