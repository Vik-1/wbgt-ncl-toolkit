;===============================================================
; Enhanced WBGT Calculator with Improved Globe Temperature Physics
; Includes wind speed effects and radiation balance corrections
;===============================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  ; for wetbulb_stull
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/heat_stress.ncl"  ; for heat_wbgt_inout

;---------------------------------------------------------------
; Improved Globe Temperature Calculation (Liljegren-inspired)
;---------------------------------------------------------------
undef("calculate_globe_temperature")
function calculate_globe_temperature(Ta:numeric, SW:numeric, WS:numeric)
local D, k_air, nu, alpha_sp, epsilon, epsilon_a, sigma, \
      TaK, Re, h_c, left, Tg_prev, max_iter, tol, iter, \
      f_val, f_prime, Tg_next, nlat, nlon, ilat, ilon
begin
    ;=== Physical Constants ===
    D        = 0.05     ; Globe diameter (m)
    k_air    = 0.025    ; Thermal conductivity of air (W/m/K)
    nu       = 1.5e-5   ; Kinematic viscosity of air (m²/s)
    alpha_sp = 0.95     ; Solar absorptivity
    epsilon  = 0.95     ; Longwave emissivity
    epsilon_a= 0.85     ; Atmospheric emissivity
    sigma    = 5.67e-8  ; Stefan-Boltzmann constant

    ;=== Initialize Output ===
    Tg = new(dimsizes(Ta), typeof(Ta), Ta@_FillValue)
    nlat = dimsizes(Ta&lat)
    nlon = dimsizes(Ta&lon)

    ;=== Main Calculation Loop ===
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            ; Skip missing values
            if(ismissing(Ta(ilat,ilon)) .or. \
               ismissing(SW(ilat,ilon)) .or. \
               ismissing(WS(ilat,ilon))) then
                continue
            end if

            ; Convert to Kelvin
            TaK = Ta(ilat,ilon) + 273.15

            ;=== Calculate Convective Coefficient ===
            WS_val = where(WS(ilat,ilon).lt.0, 0.0, WS(ilat,ilon))  ; Ensure non-negative
            Re = WS_val * D / nu  ; Reynolds number
            h_c = 0.0014 * (Re^0.6) * (k_air / D)  ; Convective heat transfer coefficient

            ;=== Radiation Balance ===
            ; Absorbed solar + atmospheric longwave
            left = alpha_sp * SW(ilat,ilon) + epsilon * epsilon_a * sigma * TaK^4

            ;=== Newton-Raphson Solver ===
            Tg_prev = TaK  ; Initial guess
            max_iter = 100
            tol = 0.01     ; Convergence tolerance (K)
            
            do iter = 0, max_iter-1
                ; Energy balance function: 
                ; f(Tg) = εσTg⁴ + h_c(Tg - TaK) - left
                f_val = epsilon*sigma*Tg_prev^4 + h_c*(Tg_prev - TaK) - left
                
                ; Derivative: f'(Tg) = 4εσTg³ + h_c
                f_prime = 4.0 * epsilon * sigma * Tg_prev^3 + h_c
                
                ; Avoid division by zero
                if(abs(f_prime).lt.1e-7) then
                    break
                end if
                
                ; Newton-Raphson update
                Tg_next = Tg_prev - f_val / f_prime
                
                ; Check convergence
                if(abs(Tg_next - Tg_prev).lt.tol) then
                    Tg(ilat,ilon) = Tg_next - 273.15  ; Convert to °C
                    break
                end if
                
                Tg_prev = Tg_next
                
                ; Use last iteration if max reached
                if(iter.eq.max_iter-1) then
                    Tg(ilat,ilon) = Tg_next - 273.15
                end if
            end do
        end do
    end do

    copy_VarCoords(Ta, Tg)
    return(Tg)
end

;---------------------------------------------------------------
; Main Program
;---------------------------------------------------------------
begin
    ;=== Input Files ===
    T_file  = addfile("air_temperature.nc", "r")       ; °C
    RH_file = addfile("relative_humidity.nc", "r")     ; %
    SW_file = addfile("solar_radiation.nc", "r")       ; W/m²
    WS_file = addfile("wind_speed.nc", "r")            ; m/s

    ;=== Read Variables ===
    T  = T_file->T
    RH = RH_file->RH
    SW = SW_file->SW
    WS = WS_file->WS

    ;=== Ensure Consistent Missing Values ===
    if(isatt(T, "_FillValue")) then
        fillval = T@_FillValue
        RH = where(ismissing(RH), fillval, RH)
        SW = where(ismissing(SW), fillval, SW)
        WS = where(ismissing(WS), fillval, WS)
    end if

    time = T&time
    lat  = T&lat
    lon  = T&lon

    ;=== Output Array ===
    wbgt = new(dimsizes(T), typeof(T), T@_FillValue)

    ;=== Compute WBGT at Each Time Step ===
    do t = 0, dimsizes(time)-1
        print("Processing time step: " + t + "/" + (dimsizes(time)-1))
        
        ; Natural wet-bulb temperature (°C)
        tw = wetbulb_stull(T(t,:,:), RH(t,:,:), (/0,0/), False)
        
        ; Globe temperature (°C)
        tg = calculate_globe_temperature(T(t,:,:), SW(t,:,:), WS(t,:,:))
        
        ; Air temperature (°C)
        ta = T(t,:,:)
        
        ; Outdoor WBGT calculation
        wbgt(t,:,:) = heat_wbgt_inout(tw, tg, ta, (/0,0/), 2)
    end do

    ;=== Output to NetCDF ===
    output_file = "wbgt_improved.nc"
    system("rm -f " + output_file)  ; Remove if exists
    ncout = addfile(output_file, "c")
    
    ; Global attributes
    ncout@title = "WBGT Calculation with Improved Physics"
    ncout@source = "NCL script with Liljegren-inspired globe temperature"
    ncout@history = "Created " + systemfunc("date")
    
    ; Define dimensions
    dimNames = (/"time", "lat", "lon"/)
    dimSizes = dimsizes(wbgt)
    dimUnlim = (/True, False, False/)
    filedimdef(ncout, dimNames, dimSizes, dimUnlim)
    
    ; Define variables
    filevardef(ncout, "time", typeof(time), "time")
    filevardef(ncout, "lat",  typeof(lat),  "lat")
    filevardef(ncout, "lon",  typeof(lon),  "lon")
    filevardef(ncout, "WBGT", typeof(wbgt), dimNames)
    
    ; Variable attributes
    ncout->time@long_name = "Time"
    ncout->time@units     = time@units
    
    ncout->lat@long_name = "Latitude"
    ncout->lat@units     = lat@units
    
    ncout->lon@long_name = "Longitude"
    ncout->lon@units     = lon@units
    
    ncout->WBGT@long_name = "Wet Bulb Globe Temperature"
    ncout->WBGT@units = "degC"
    ncout->WBGT@description = "Outdoor WBGT calculated with improved physics"
    
    ; Write data
    ncout->time = time
    ncout->lat  = lat
    ncout->lon  = lon
    ncout->WBGT = wbgt
    
    print(" ")
    print("✅ Successfully created WBGT file: " + output_file)
    print("   Dimensions: time=" + dimsizes(time) + \
          ", lat=" + dimsizes(lat) + ", lon=" + dimsizes(lon))
end
