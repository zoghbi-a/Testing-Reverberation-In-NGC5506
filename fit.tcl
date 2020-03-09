source ~/codes/xspec/az.tcl
lmod relxill

# ---------------------------------------------- #
# fit po+zga to 5-7 keV to get enegy of the line #
proc fit_1 {sfile {suff ""}} {
    para err 5; para lev 5
    da $sfile
    ign 0.0-5., 7.-**
    mo po+zga & 1 & 1e-2& 6.4 & 0 -1 & 0.00618& 1e-4
    fit 1000
    az_calc_errors [az_free_params] fits/fit_1$suff 1.0
}

# ---------------------------------------------------- #
# fit zedge*po+zga to 4-8 keV to get enegy of the line #
proc fit_2 {sfile {suff ""}} {
    para err 5; para lev 5
    da $sfile
    ign 0.0-4., 8.-**
    mo zedge*po+zga& 7.1 .1 6.9 6.9 7.6 7.6 & 0.1 & 0.00618 & 1 & 1e-2& 6.4 & 0 -1 & 0.00618& 1e-4
    fit 1000
    az_calc_errors [az_free_params] fits/fit_2$suff 1.0
}

# ----------------------------------------- #
# scan 0.3-10 residuals after adding xillver  #
proc fit_3 {sfile {suff ""}} {
    para err 5; para lev 1; para st 40; query yes
    da $sfile
    ign 0.0-0.3, 10.-**
    mo tbabs*ztbabs*(ztbabs*po + xillver + brems) & 0.041 -1& .1 .1 .001 .001 & 0.00618 & \
        2 & 0.00618 & 1.7 & 1e-2 & \
        =6 & 1 -1 & 1000 -1 & 0 -1 & 0.00618 & 30 -1 & -1 -1 & 1e-4 & 0.2 & 1e-3
    fit 1000
    set xcm fits/fit_3$suff.xcm; rm $xcm > /dev/null 2>&1; save all $xcm
    az_scan_en_norm $xcm 6 {16 2 10 200} {log 19 1e-6 1e-4 100} 1 0.00618
}

# ---------------------------------------------------------------- #
# scan 2-10 residuals; similar to fit_3 but fits to 2-10 keV only  #
proc fit_3a {sfile {suff ""}} {
    para err 5; para lev 1; para st 40; query yes
    da $sfile
    ign 0.0-2.0, 10.-**
    mo tbabs*(ztbabs*po + xillver) & 0.041 -1& \
        2 & 0.00618 & 1.7 & 1e-2 & \
        =6 & 1 -1 & 1000 -1 & 0 -1 & 0.00618 & 30 -1 & -1 -1 & 1e-4
    fit 1000
    set xcm fits/fit_3a$suff.xcm; rm $xcm > /dev/null 2>&1; save all $xcm
    az_scan_en_norm $xcm 4 {6 2 10 200} {log 9 1e-6 1e-4 100} 1 0.00618
}

# ---------------------------------------------------------------- #
# simulate stat improvements in fit_3a when adding a gaussian line  #
proc fit_3a_sim {sfile {suff ""}} {
    query yes
    az_sim_dchi2 10000 fits/fit_3a${suff}.xcm 4
}

# ------------------------------------------- #
# xillver + broad gaussian; start from fit_3  #
proc fit_4 {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_3$suff
    add 6 zga & 6.7 & 0.2 & 0.00618& 1e-5
    fit 1000
    az_calc_errors [az_free_params] fits/fit_4$suff 1.0
}

# --------------------------------- #
# similar to fit_4, but use fluxes  #
proc fit_4a {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_4$suff
    add 4 cflux & 7 & 10 & -10
    new 10 1 -1
    add 6 cflux & 6.1 & 6.7 & -12
    new 21 1 -1
    add 8 cflux & 4 & 9 & -11
    new 28 1 -1
    set ifree [az_free]; freez $ifree
    thaw 8 13 24; fit 1000
    thaw $ifree
    az_calc_errors [az_free_params] fits/fit_4a$suff 1.0
}

# ------------------------------------------------ #
# xillver + two narrow gaussians; start from fit_3 #
proc fit_5 {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_3$suff
    add 6 zga & 6.7 & 0 -1 & 0.00618& 1e-5
    add 7 zga & 6.9 & 0 -1 & 0.00618& 1e-5
    fit 1000
    az_calc_errors [az_free_params] fits/fit_5$suff 1.0
}

# ----------------- #
# fit_5 with fluxes #
proc fit_5a {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_5$suff
    add 4 cflux & 7 & 10 & -10
    new 10 1 -1
    add 6 cflux & 6.1 & 6.7 & -12
    new 21 1 -1
    add 8 cflux & 4 & 9 & -11
    new 28 1 -1
    add 10 cflux & 4 & 9 & -11
    new 35 1 -1
    set ifree [az_free]; freez $ifree
    thaw 8 13 24 31; fit 1000
    thaw $ifree
    az_calc_errors [az_free_params] fits/fit_5a$suff 1.0
}

# ----------------------------------- #
# xillver + relxill; start from fit_3 #
proc fit_6 {sfile {suff ""}} {
    para err 5; para lev 5; para st 5; query yes
    @fits/fit_3$suff
    add 4 relxill & 3. .1 2 2 10 10 & =6& 999 -1 1 1 1000 1000 & ,-1& 30 .1 & \
        2. .1 1.3 1.3 800 800 & =8+1.0& 0.00618& 1.8 & 3& 1 -1& 1000 -1& -1 -1& 1e-4
    new 20=14; new 22=14
    new 23=16; new 24=17; new 27=10
    fit 1000
    stepp best 15 0 3 40
    
    az_calc_errors [az_free_params] fits/fit_6$suff 1.0
}

# ----------------------------------------------------------------- #
# xillver + relxill, similar to fit_6, but fit all spectra together #
proc fit_7 {sfile {suff ""}} {
    para err 35; para lev 1; para st 35; query yes
    @fits/fit_6__1
    add 3 tbpcf & =7 & 1 -1 & 0.00618 & /*
    del 4
    
    foreach i {2 3 4 5 6 7 8 9} { da $i:$i spec_$i.grp }
    ign 1-9:0.0-0.3,10.-**
    new 7 3 -1
    new 16 3
    new 31 0.5
    new 12 50
    
    az_untie_spec {2 3 4 5 6 7 8 9} {4 12 15 20 22 }
    foreach i [az_find 15] {new [expr $i+6]=$i; new [expr $i+8]=$i}
    foreach i [az_find 12] {new $i 100 .1 1.3 1.3 800 800}
    renorm; fit 1000

    # 3195/1519
    set xcm fits/fit_7a.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    
    # add apec #
    add 7 apec & .8 & 1 -1 & 0.00618 & 3e-5 & /*
    fit
    # 2829/1517
    set xcm fits/fit_7b.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    
    # thaw cf
    thaw 5; fit
    # 2640/1516
    set xcm fits/fit_7c.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    
    
    # add diskbb at the lowest energies
    add 8 bb & 0.67 & 1e-3; fit 
    # 2149/1514
    set xcm fits/fit_7d.xcm; rm $xcm >/dev/null 2>&1; save all $xcm    
    az_calc_errors [az_free_params] fits/fit_7d 1.0
    
    # add wa;
    add 3 zxipc & .3 & -1 & 1 -1 & 0.00618
    untie [az_find 4] [az_find 5]; fit
    foreach i [az_find 5] {stepp best $i -1 1 30}
    # 1957/1496
    set xcm fits/fit_7e.xcm; rm $xcm >/dev/null 2>&1; save all $xcm  
    
    # add two lines around 2 keV
    add 10 zga & 1.65 & 0.1 & 0.00618 & 5e-5; fit
    add 11 zga & 2.2 & 0.07 & 0.00618& 3e-5; fit
    foreach i [az_find 5] {stepp best $i -1 1 30}
    # 1765/1490
    set xcm fits/fit_7f.xcm; rm $xcm >/dev/null 2>&1; save all $xcm  
    
    # add a line at 6.9 keV
    add 12 zga & 6.9 & 0 -1 & 0.00618 & 1e-5;fit
    # 1685/1488
    set xcm fits/fit_7g.xcm; rm $xcm >/dev/null 2>&1; save all $xcm  
    az_calc_errors [az_free_params] fits/fit_7g 1.0
    
}

# ----------------------------------------------------- #
# xillver + relxill, similar to fit_7, but fit 2-10 keV #
proc fit_8 {sfile {suff ""}} {
    para err 35; para lev 1; para st 40; query yes
    @fits/fit_6__1
    add 3 tbfeo & =8 & & & 0.00618 & /*
    del 8; del 4; del 2
    
    foreach i {2 3 4 5 6 7 8 9} { da $i:$i spec_$i.grp }
    ign 1-9:0.0-2.,10.-**
    new 6 3 -1
    new 15 3
    new 11 150
    
    az_untie_spec {2 3 4 5 6 7 8 9} {14 19 21}
    foreach i [az_find 14] {new [expr $i+6]=$i; new [expr $i+8]=$i}
    new 16=4; new 23=4
    renorm; fit 1000

    
    stepp best 10 10 80 40
    # 1432/1148
    set xcm fits/fit_8a.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    
    
    # untie nh #
    untie [az_find 2]
    # 1396/1140
    set xcm fits/fit_8b.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    
    # thaw A_Fe
    thaw 4; fit
    new 10 ,,,3 3 89 89
    # 1328/1139
    set xcm fits/fit_8c.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    az_calc_errors [az_free_params] fits/fit_8c 1.0
    # fit ok, but incl is 89; try the following ...
    
    ####
    #add nustar data
    da 10:10 spec_a.grp
    ign 10:0.0-3.,79.-**
    az_untie_spec {10} {2 14 19 21}
    new 281=275; new 283=275
    fit
    stepp best 10 10 89 20
    # 2154/1914
    set xcm fits/fit_8d.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    
    # add also b spectra
    da 11:11 spec_b.grp; ign 11:0.0-3.,79.-**
    foreach i {292 304 309 310 311 312} {new $i=[expr $i-29]}
    add 1 cons & 1 -1 & /*
    untie 301; thaw 301
    fit
    # 2919/2678
    set xcm fits/fit_8e.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    az_calc_errors [az_free_params] fits/fit_8e 1.0
    chain len 100000; chain burn 1000000; chain walk 460; para walk 40
    chain run fits/fit_8e.fits
    
    # nh of nustar is below the rest; set it to the average
    new 273=(3+33+63+93+123+153+183+213+243)/9.
    # 2938/2679
    set xcm fits/fit_8e1.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    az_calc_errors [az_free_params] fits/fit_8e1 1.0
    
    # untie rin:
    foreach i [az_find 12] {new $i [tcloutr par 12]}; new 312=282; fit 100
    # 2898/2669
    set xcm fits/fit_8f.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2898 2669 2919 2678: p=0.0228; barely significant
    
    # from fit_8e; do separate low and high flux rin
    new 42 [tcloutr par 12]
    foreach i {162 192 222} {new $i=42}
    fit 1000
    # 2913/2677
    set xcm fits/fit_8g.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2913 2677 2919 2678: p=0.0189; again, barely significant
    az_calc_errors {12 42} fits/fit_8g 1.0
    
    
    # from fit_8e, test separate theta for low/high flux
    foreach i [az_find 28] {new $i=[expr $i-17]}
    new 41 [tcloutr par 11]
    foreach i {161 191 221} {new $i=41}
    fit 1000
    # 2916/2677
    set xcm fits/fit_8h.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2916 2677 2919 2678: p=0.097; not significant
    
    
    # from fit_8e; replace relxill with gaussian; useful later for cov spectra
    add 4 zga & 6.7 & 0.1 & 0.00618& 1e-5 
    del 5
    foreach i [az_find 13] {new $i=[expr $i-2]}
    freez 18
    untie [az_find 10]; new 210=190
    fit 1000
    set xcm fits/fit_8z.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
}

# -------------------------------------- #
# xillver + relxilllp, start from fit_8e #
proc fit_9 {sfile {suff ""}} {
    para err 35; para lev 1; para st 40; query yes
    @fits/fit_8e
    
    add 4 relxilllp & 50 .1 3 3 500 500 & ,-1 & 61& 190 .1 6 6 998 998 & 1000 &\
        0.00618& 1.8 & 3 & 1.56 & 1000& ,-1& 1 & 1e-4 & /*
    foreach i [az_find 13] {new $i=[expr $i+15]}
    untie [az_find 13]
    del 5; del 5
    foreach i [az_find 20] {new $i=[expr $i-7]}
    untie [az_find 19]
    new 283=256; new 289=262; new 21=15; new 25=9; new 15=5
    thaw 9
    renor; fit 1000
    # 2933/2687
    set xcm fits/fit_9a.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
        az_calc_errors [az_free_params] fits/fit_9a 1.0
    #chain len 100000; chain burn 1000000; chain walk 460; para walk 40
    #chain run fits/fit_9a.fits
    
    # 9a and untie rin:
    foreach i [az_find 10] {new $i [tcloutr par 10]}; new 280=253; fit 100
    # 2920/2678
    set xcm fits/fit_9b.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2920 2678 2933 2687: p=0.22; not significant
    
    # 9a and untie h:
    foreach i [az_find 7] {new $i [tcloutr par 7]}; new 277=250; fit 100
    # 2915/2678
    set xcm fits/fit_9c.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2915 2678 2933 2687: p=0.057; barely significant
    
    # from fit_9a; do separate low and high flux rin
    new 37 [tcloutr par 10]
    foreach i {145 172 199} {new $i=37}
    fit 1000
    # 2933/2686
    set xcm fits/fit_9d.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # no imporvement
    
    # from fit_9a; do separate low and high flux height
    new 34 [tcloutr par 7]
    foreach i {142 169 196} {new $i=34}
    fit 1000
    # 2933/2686
    set xcm fits/fit_9e.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    az_calc_errors {7 34} fits/fit_9e 1.0
    # not significant
    # h_low = 93-15+70; h_hi = 96-13+80
}

# ------------------------------------ #
# compare spec fit fit_9a to lag model #
prov fit_9_vs_lag {} {
    setpl ener
    mkdir -p plots
    
    # show the 6.7 keV feature
    @fits/fit_9a
    az_plot u plots/fit_9a 9a 1 1
    
    new 17 0; new 18 0; freez 1-297
    thaw 19 46 73 100 127 154 181 208 235 262; fit
    az_plot u plots/fit_9a_r0 9a_r0 1 1
    
    # lag model
    @fits/fit_9a
    freez 1-297; thaw 19 46 73 100 127 154 181 208 235 262
    new 9 15.05; new 10 4.4 -1 2 2; fit
    save all fits/fit_9a__vs_lag
    
    az_plot u plots/fit_9a__lag 9a__lag 1 1
    
    new 17 0; new 18 0; fit
    az_plot u plots/fit_9a_r0__lag 9a_r0__lag 1 1
}


# ---------------------------------------- #
# xillver + grid25_aout, start from fit_8e #
proc fit_10 {sfile {suff ""}} {
    para err 35; para lev 1; para st 40; query yes
    @fits/fit_8e
    
    add 6 atable{grid25_aout.fits} & 1e23 & 3.4& 0.00618& 1e-4 & /*
    untie [az_find 26]; new 366=332
    del 4
    foreach i [az_find 13] {new $i=[expr $i-6]}
    new 5 1 -1; new 18 62 -1
    fit 1000
    # 2978/2680
    set xcm fits/fit_10a.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    az_calc_errors [az_free_params] fits/fit_10a 1.0
    
    # from fit_10a, untie nh of grid25_aout
    untie [az_find 9]; new 209=189; fit
    # 2977/2671
    set xcm fits/fit_10b.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    
    # from fit_10a, untie xi of grid25
    untie [az_find 10]; new 210=190; fit
    # 2954/2671
    set xcm fits/fit_10c.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2954 2671 2978 2680 -> 0.01; somewhat significant
    
    # from fit_10a; fix xi for low/hi flux intervals
    new 30 [tcloutr par 10] .1
    foreach i {110 130 150} {new $i=30}
    fit 1000
    # 2969/2679
    set xcm fits/fit_10d.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2969 2679 2978 2680 -> 0.0044; *significant*
    # comparing to the all-free fit in fit_10c; ftest 2954 2671 2969 2679 ->0.09
    # i.e. two values for the low/high flux is perferred over a signal value for all
    # and over each having its own.
    
    # 10a; but put grid25_aout inside the absorber.
    add 5 tbfeo & & 1 -1& 1 -1& 0.00618 & /*
    foreach i [az_find 9] {new $i=[expr $i-6]}
    fit 1000
    # 2952/2680
    set xcm fits/fit_10e.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
}

# ---------------------------------- #
# xillver + apec, start from fit_10a #
proc fit_20 {sfile {suff ""}} {
    para err 35; para lev 1; para st 40; query yes
    @fits/fit_10a
    
    add 5 apec & 7 & 1 -1& 0.00618& 1e-3& /*
    del 6
    untie [az_find 12]; new 212=192
    fit 1000
    # 2849/2681
    set xcm fits/fit_20a.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    az_calc_errors [az_free_params] fits/fit_20a 1.0
    
    # from fit_20a, untie kt
    untie [az_find 9]; new 209=189; fit
    # 2827/2672
    set xcm fits/fit_20b.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2827 2672 2849 2681 -> 0.014
    
    # from fit_20a; fix kt for low/hi flux intervals
    new 29 [tcloutr par 9] .1
    foreach i {109 129 149} {new $i=29}
    fit 1000
    # 2842/2680
    set xcm fits/fit_20c.xcm; rm $xcm >/dev/null 2>&1; save all $xcm
    # ftest 2842 2680 2849 2681 -> 0.010; somwhat significant
    # comparing to the all-free fit in fit_20c; ftest 2827 2672 2842 2680  ->0.078
    # i.e. two values for the low/high flux is perferred over a signal value for all
    # and over each having its own.
}

# -------------------------------- #
# base model for the added spectra #
# based on fit_8z in the spectral modeling
proc fit_add {{suff ""}} {
    # suff examples: 16lch, 16lch_hi, 16lch_s7
    da spec_add_${suff}.grp.b
    ign 1, [tcloutr nchan]
    mo TBabs(TBfeo(zgauss + powerlaw) + xillver) & 0.041 -1& 2.8 & 1 -1 & 1.09 -1& \
        0.00618& 6.755 -1& 0.2 -.1 .1 .1 .6 .6 & 0.00618 & 1e-5 .1 1e-9 1e-9 & 1.85 & \
        1e-2 &=10&=4& 1000 -1& 0 -1& 0.00618& 61.57 -1& -1 -1 & 4.821e-4 -1
    renorm; fit 1000; thaw 6 7; fit
    az_calc_errors [az_free_params] fits/spec_add_${suff} 1.0
}

# -------------------------------------------- #
# base relativistc model for the added spectra #
# based on fit_8e in the spectral modeling
proc fit_add_rel {{suff ""}} {
    # suff examples: 16lch, 16lch_hi, 16lch_s7
    da spec_add_${suff}.grp.b
    ign 1, [tcloutr nchan]
    mo TBabs(TBfeo(relxill + powerlaw) + xillver) & 0.041 -1 & 2.8 & 1 -1& 1.56 -.1 & \
        0.00618 & 3 -1 & =6 & 999& .998 -1& 61.6 -1& 189.9 -1 1.3 1.3 400 400& =p8+1& \
        0.00618 & 1.8 & 3 -1& =4& 1000 -1& -1 -1 & 6.4e-5 & =14& 1.6e-2& =14 & =4& =17& \
        0 -1& 6.18e-3& =10& -1 -1& 3.09e-4 -1
    renorm; fit 1000
    az_calc_errors [az_free_params] fits/spec_add_rel_${suff} 1.0
}

# --------------------------------------------- #
# base Compton-thin model for the added spectra #
# based on fit_10d in the spectral modeling
proc fit_add_thin {{suff ""}} {
    # suff examples: 16lch, 16lch_hi, 16lch_s7
    da spec_add_${suff}.grp.b
    ign 1, [tcloutr nchan]
    mo TBabs(ztbabs*powerlaw + atable{grid25_aout.fits} + xillver) & 0.041 -1 & 3 & \
        0.00618 & 1.8 & 2e-2 & 1.1e21 -1 & 3.6 & 0.00618& 0.15 & =4 & 1 -1& 1000 -1& 0 -1 & \
        0.00618 & 62 -1 & -1 -1& 5.47162E-04 -1
    fit 1000
    az_calc_errors [az_free_params] fits/spec_add_thin_${suff} 1.0
}

# ----------------------------------------------- #
# base Compton-thin model that is inside absorber #
# based on fit_10e in the spectral modeling
proc fit_add_thin2 {{suff ""}} {
    # suff examples: 16lch, 16lch_hi, 16lch_s7
    da spec_add_${suff}.grp.b
    ign 1, [tcloutr nchan]
    mo TBabs(ztbabs*(powerlaw + atable{grid25_aout.fits}) + xillver) & 0.041 -1 & 3 & \
        0.00618 & 1.8 & 2e-2 & 1.1e21 -1 & 3.6 & 0.00618& 0.15 & =4 & 1 -1& 1000 -1& 0 -1 & \
        0.00618 & 62 -1 & -1 -1& 5.46816E-04 -1
    fit 1000
    az_calc_errors [az_free_params] fits/spec_add_thin2_${suff} 1.0
}

# ----------------------------------------------- #
# base thermal plasma model for the added spectra #
# based on fit_20c in the spectral modeling
proc fit_add_therm {{suff ""}} {
    # suff examples: 16lch, 16lch_hi, 16lch_s7
    da spec_add_${suff}.grp.b
    ign 1, [tcloutr nchan]
    mo TBabs(ztbabs*powerlaw + apec + xillver) & 0.041 -1 & 3 & \
        0.00618 & 1.8 & 2e-2 & 9 & 1 -1& 0.0061 & 1e-3& =4 & 1 -1& 1000 -1& 0 -1 & \
        0.00618 & 62 -1 & -1 -1& 5.38018E-04 -1
    fit 1000
    az_calc_errors [az_free_params] fits/spec_add_therm_${suff} 1.0
}

# ------------------------------------------------- #
# use the model from fit_add to fit the cov spectra #
proc fit_cov {{suff ""} {bsuff 16lch} {FQ {1 2}}} {
    # suff: "", "_lo", "_hi", "_s7" etc
    # bsuff: binning suffix; e.g. 16lch
    # FQ is the frequencies to fit
    para err 20; query yes; setpl ene
    # loop over frequency bins #
    foreach ii $FQ {
        @fits/spec_add_${bsuff}${suff}
        freez 2 6 7 10; del 5; new 9 1e-9 -1
        
        da cov${suff}_f${ii}.pha
        ign 1, [tcloutr nchan]
        az_calc_errors [az_free_params] fits/cov_add${suff}_f${ii} 1.0
        cpd /cps;az_plot u plots/cov_add${suff}_f${ii} ${suff}_f${ii} 1 1
        
        # thaw gamma #
        thaw 10; fit
        az_calc_errors [az_free_params] fits/cov_add${suff}_f${ii}_m1 1.0
        cpd /cps;az_plot u plots/cov_add${suff}_f${ii}_m1 ${suff}_f${ii}_m1 1 1
        
        # thaw gaussian norm #
        thaw 9; fit
        az_calc_errors [az_free_params] fits/cov_add${suff}_f${ii}_m2 1.0
        cpd /cps;az_plot u plots/cov_add${suff}_f${ii}_m2 ${suff}_f${ii}_m2 1 1
        
        # thaw gaussian energy #
        thaw 6; fit
        az_calc_errors [az_free_params] fits/cov_add${suff}_f${ii}_m3 1.0
        cpd /cps;az_plot u plots/cov_add${suff}_f${ii}_m2 ${suff}_f${ii}_m3 1 1
    }
}


# ----------------------------------------------------- #
# use the model from fit_add_rel to fit the cov spectra #
proc fit_cov_rel {{suff ""} {bsuff "16lch"} {FQ {1 2}}} {
    # suff: "", "_lo", "_hi", "_s7" etc
    # bsuff: binning suffix; e.g. 16lch
    # FQ is the frequencies to fit
    para err 20; query yes; setpl ene
    # loop over frequency bins #
    foreach ii $FQ {
        @fits/spec_add_rel_${bsuff}${suff}
        freez 2 14; del 5; new 19 1e-9 -1 1e-9 1e-9
        
        da cov${suff}_f${ii}.pha
        ign 1, [tcloutr nchan]
        az_calc_errors [az_free_params] fits/cov_add_rel${suff}_f${ii} 1.0
        cpd /cps;az_plot u plots/cov_add_rel${suff}_f${ii} ${suff}_f${ii} 1 1
        
        # thaw gamma #
        thaw 14; fit
        az_calc_errors [az_free_params] fits/cov_add_rel${suff}_f${ii}_m1 1.0
        cpd /cps;az_plot u plots/cov_add_rel${suff}_f${ii}_m1 ${suff}_f${ii}_m1 1 1
        
        # thaw refl. norm #
        new 19 1e-6 .1; fit
        az_calc_errors [az_free_params] fits/cov_add_rel${suff}_f${ii}_m2 1.0
        cpd /cps;az_plot u plots/cov_add_rel${suff}_f${ii}_m2 ${suff}_f${ii}_m2 1 1
        
    }
}

# ----------------------------------------------------- #
# use the model from fit_add_thin to fit the cov spectra #
proc fit_cov_thin {{suff ""} {bsuff "16lch"} {FQ {1 2}}} {
    # suff: "", "_lo", "_hi", "_s7" etc
    # bsuff: binning suffix; e.g. 16lch
    # FQ is the frequencies to fit
    para err 20; query yes; setpl ene
    # loop over frequency bins #
    foreach ii $FQ {
        @fits/spec_add_thin_${bsuff}${suff}
        freez 2 4 7; del 5; new 9 1e-9 -1 1e-9 1e-9
        
        da cov${suff}_f${ii}.pha
        ign 1, [tcloutr nchan]
        az_calc_errors [az_free_params] fits/cov_add_thin${suff}_f${ii} 1.0
        cpd /cps;az_plot u plots/cov_add_thin${suff}_f${ii} ${suff}_f${ii} 1 1
        
        # thaw gamma #
        thaw 4; fit
        az_calc_errors [az_free_params] fits/cov_add_thin${suff}_f${ii}_m1 1.0
        cpd /cps;az_plot u plots/cov_add_thin${suff}_f${ii}_m1 ${suff}_f${ii}_m1 1 1
        
        # thaw refl. norm #
        new 9 1e-6 .1; fit
        az_calc_errors [az_free_params] fits/cov_add_thin${suff}_f${ii}_m2 1.0
        cpd /cps;az_plot u plots/cov_add_thin${suff}_f${ii}_m2 ${suff}_f${ii}_m2 1 1
        
    }
}

# ------------------------------------------------------- #
# use the model from fit_add_thin2 to fit the cov spectra #
proc fit_cov_thin2 {{suff ""} {bsuff "16lch"} {FQ {1 2}}} {
    # suff: "", "_lo", "_hi", "_s7" etc
    # bsuff: binning suffix; e.g. 16lch
    # FQ is the frequencies to fit
    para err 20; query yes; setpl ene
    # loop over frequency bins #
    foreach ii $FQ {
        @fits/spec_add_thin2_${bsuff}${suff}
        freez 2 4 7; del 5; new 9 1e-9 -1 1e-9 1e-9
        
        da cov${suff}_f${ii}.pha
        ign 1, [tcloutr nchan]
        az_calc_errors [az_free_params] fits/cov_add_thin2${suff}_f${ii} 1.0
        cpd /cps;az_plot u plots/cov_add_thin2${suff}_f${ii} ${suff}_f${ii} 1 1
        
        # thaw gamma #
        thaw 4; fit
        az_calc_errors [az_free_params] fits/cov_add_thin2${suff}_f${ii}_m1 1.0
        cpd /cps;az_plot u plots/cov_add_thin2${suff}_f${ii}_m1 ${suff}_f${ii}_m1 1 1
        
        # thaw refl. norm #
        new 9 1e-6 .1; fit
        az_calc_errors [az_free_params] fits/cov_add_thin2${suff}_f${ii}_m2 1.0
        cpd /cps;az_plot u plots/cov_add_thin2${suff}_f${ii}_m2 ${suff}_f${ii}_m2 1 1
        
    }
}

# ----------------------------------------------------- #
# use the model from fit_add_therm to fit the cov spectra #
proc fit_cov_therm {{suff ""} {bsuff "16lch"} {FQ {1 2}}} {
    # suff: "", "_lo", "_hi", "_s7" etc
    # bsuff: binning suffix; e.g. 16lch
    # FQ is the frequencies to fit
    para err 20; query yes; setpl ene
    # loop over frequency bins #
    foreach ii $FQ {
        @fits/spec_add_therm_${bsuff}${suff}
        freez 2 4 6; del 5; new 9 1e-9 -1 1e-9 1e-9
        
        da cov${suff}_f${ii}.pha
        ign 1, [tcloutr nchan]
        az_calc_errors [az_free_params] fits/cov_add_therm${suff}_f${ii} 1.0
        cpd /cps;az_plot u plots/cov_add_therm${suff}_f${ii} ${suff}_f${ii} 1 1
        
        # thaw gamma #
        thaw 4; fit
        az_calc_errors [az_free_params] fits/cov_add_therm${suff}_f${ii}_m1 1.0
        cpd /cps;az_plot u plots/cov_add_therm${suff}_f${ii}_m1 ${suff}_f${ii}_m1 1 1
        
        # thaw refl. norm #
        new 9 1e-6 .1; fit
        az_calc_errors [az_free_params] fits/cov_add_therm${suff}_f${ii}_m2 1.0
        cpd /cps;az_plot u plots/cov_add_therm${suff}_f${ii}_m2 ${suff}_f${ii}_m2 1 1
        
    }
}


# ------------------------------------------ #
# plot the components with error using mcmc  #
# used to get the reflection fraction
proc plot_components {xcm {nsim 400}} {
    # xcm: file in fits/ containing the fit; e.g. spec_add_17lch, cov_add_hi_f1_m1
    @fits/$xcm
    fit
    chain len 10000; chain burn 50000; chain walk 50; para walk 35; chain unload 1
    rm fits/${xcm}.fits  >/dev/null 2>&1
    chain run fits/${xcm}.fits
    setpl ener; setpl add
    rm tmp_*qdp
    for {set i 1} {$i<=$nsim} {incr i} {
        az_rand_pars; setpl co wd tmp_$i; pl d; setpl del 1
    }
    exec cat {*}[glob tmp*qdp] | grep -Ev "^R|\!|\@" > fits/${xcm}.plot.dat
    rm tmp*qdp
}


# ---------------------------------------------------- #
# scan 2-10 cov residuals; start from model _m1 above  #
proc cov_scan {suff} {
    # suff e.g.: f2, hi_f1, hi_f2, lo_f2
    para lev 1; para st 40; query yes
    az_scan_en_norm fits/cov_add_${suff}.xcm 1 {log 1 2 10 200} {log 4 1e-5 5e-5 200} 1 0.00618
}

# ---------------------------------------------------------------- #
# simulate stat improvements in fit_3a when adding a gaussian line  #
proc cov_scan_sim {suff} { 
    query yes
    az_sim_dchi2 10000 fits/cov_add_${suff}.xcm 1
}
