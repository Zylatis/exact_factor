#!/usr/bin/gnuplot


set terminal cairolatex standalone pdf size  4.5in,3.5in  dashed dl 1 transparent


set style line 1 pointtype 1
set style line 2 pointtype 4
set style line 3 pointtype 2
set style line 4 pointtype 19
set style line 5 pointtype 14
#~ set style line 6 pointtype 2

set style line 99 lt 0 lc rgb "black"  lw 1

set style line 7 lc  "web-blue" lw 2 #pointtype 1
set style line 8 lc rgb  "black" lw 2 #pointtype 19
set style line 9 lc rgb  "red" lw 2 #pointtype 1
set style line 10 lc rgb  "brown" lw 2 #pointtype 1
set style line 11 lc rgb  "blue" lw 2 #pointtype 1
set style line 13 lc rgb "black" lw 2 dashtype 2

set style line 14 lc rgb  "brown" lw 2 dashtype 4
set style line 15 lc rgb  "blue" lw 2 dashtype 4
set style line 16 lc rgb  "black" lw 2 dashtype 4


set termoption dashlength 5
dR = 0.0090045
set output "SM_setup.tex"
set key font "Helvetica, 25"
set xlabel "R"
set ylabel ""
set key at 1.9,-0.05
set grid ls 99	
set ytics format "\\large %g"
set xtics format "\\large %g"
set key spacing 0.6
plot [-8:8][-0.4:0.5] "../NACV/NACV1-12_2k.txt" u ($0*dR-9):($1/10) w l lw 2 title "\\large $d_{12}/10$", \
   "../NACV/NACV2-12_2k.txt" u ($0*dR-9):($1/20) w l lc "red" lw 2 title "\\large $d^{2}_{12}/20$" ,\
 "../BOPES/1_bopes_2k.dat" u ($0*dR-9):($2) w l lc "black" lw 2 title "\\large$\\epsilon_{BO}$", \
  "../BOPES/2_bopes_2k.dat"u ($0*dR-9):($2) w l lc "black" lw 2 title ""

##################################################################

set output "mask_example.tex"

plot [0:8][0:1.5] "../output_git_run_maskA/1global_mask.txt" u ($0*dR-9):($1) w l ls 9 title "" ,\
"../output_git_run_maskA/1nd.txt" u ($0*dR-9):($1) w l title "" ls 8

##################################################################
#~ set lmargin 0
#~ set output "no_mask.tex"
#~ set key at 8, 4.5
#~ set ytics format "\\large %g"
#~ set xtics format "\\large %g"
#~ plot [-9:9][-1:5] "../output_git_run_nomask/57TDPES.txt" u ($0*dR-9):($1) w l ls 11 title "\\large $\\epsilon(R,t)$" ,\
#~ "../output_git_run_nomask/57GCOC.txt" u ($0*dR-9):($1/10) w l title "\\large $\\nabla \\chi/\\chi$" ls 10 ,\
"../output_git_run_nomask/57nd.txt" u ($0*dR-9):($1) w l title "\\large $\|\\chi|^2$" ls 8 
##################################################################

#~ set output "maskuentalk.tex"
#~ unset key
#~ Rmin=730
#~ Rmax=1500
#~ set key at 900, -0.8
#~ plot [-3:5][-1:1.1]"../output_git_run_maskuen/61TDPES.txt" u ($0*dR-9):($1) w l ls 11 title "TDPES" ,\
#~ "../output_git_run_maskuen/61GCOC.txt" u ($0*dR-9):($1/10) w l ls 10 title "Re(GCOC)" ,\
#~ "../output_git_run_maskuen/61global_mask.txt" u ($0*dR-9):($1) w l lc "dark-green" lw 3 title "mask" ,\
"../output_git_run_maskuen/61nd.txt" u ($0*dR-9):($1) w l ls 8 title "$\|\\chi|^2$"

##################################################################

#~ set terminal pngcairo enhanced size 1500,1080 dashed font "Helvetica,25" 
set terminal cairolatex standalone pdf size 4.5in,3.5in  dashed dl 1 transparent
set xtics
set ytics
set format y
set title ""	
set output "mask_dt.tex"
set key at 3, -0.25
t = 340
set title "Mask on $\\partial_t C(R,t)$ at $t = 340$au"
plot [0:4][-1:1.5	] "../output_git_run_maskdt/".t."TDPES.txt" u  ($0*dR-9):($1) w l ls 11 title "$\\epsilon(R,t)$" ,\
"../output_git_run_maskA/".(t/10)."nd.txt" u ($0*dR-9):($1) w l title "" ls 13,\
"../output_git_run_maskdt/".t."nd.txt" u ($0*dR-9):($1) w l title "$\|\\chi(R,t)|^2$" ls 9,\
"../output_git_run_maskA/".(t/10	)."TDPES.txt" u ($0*dR-9):($1) w l title "best case"  ls 13



##################################################################

# PLOTTING MASK, TDPES, AND BOPES (NO GCOC FOR NOW?)

set terminal cairolatex standalone pdf size 9.5in,5in  dashed transparent
ali_folder = "../ext_exact_tdnd/"
unset format y
set xlabel ""
set output "mask_gcoc_grid.tex" 
set multiplot layout 2,3 
i = "1"
set xtics -8, 4,8
set title "t = ".(i*10)."a.u."
N=0.742399
dR = 0.0090045
plot [:][0:0.3]"../output_git_run_maskA/".i."TDPES.txt" u ($0*dR-9):($1) w l title "" lw 3 lc "blue" ,\
"../output_git_run_maskA/".i."global_mask.txt" u ($0*dR-9):($1/4) w l title "" lc black lw 2,\
"../BOPES/2_bopes_2k.dat" u ($0*dR-9):($2) w l title "" lc "gray"  ,\
"../BOPES/1_bopes_2k.dat" u ($0*dR-9):($2) w l title ""  lc "gray" ,\
"../output_git_run_maskA/".i."nd.txt" u ($0*dR-9):($1/10) w l title ""   lc "red" lw 3,\
ali_folder.i."_tdnd.dat"  u ($0*dR-9):($1/10)  w l title "" dashtype 4 lc "black" 



set format y ""
i = "10"
set title "t = ".(i*10)."a.u."
plot [:][0:0.3]"../output_git_run_maskA/".i."TDPES.txt" u ($0*dR-9):($1) w l title "" lw 3 lc "blue" ,\
"../output_git_run_maskA/".i."global_mask.txt" u ($0*dR-9):($1/4) w l title "" lc black lw 2,\
"../BOPES/2_bopes_2k.dat" u ($0*dR-9):($2) w l title "" lc "gray"  ,\
"../BOPES/1_bopes_2k.dat" u ($0*dR-9):($2) w l title ""  lc "gray" ,\
"../output_git_run_maskA/".i."nd.txt" u ($0*dR-9):($1/10) w l title ""   lc "red" lw 3 ,\
ali_folder.i."_tdnd.dat"  u ($0*dR-9):($1/10)  w l title "" dashtype 4 lc "black" 


i = "40"
set title "t = ".(i*10)."a.u."
plot [:][0:0.3]"../output_git_run_maskA/".i."TDPES.txt" u ($0*dR-9):($1) w l title "" lw 3 lc "blue" ,\
"../output_git_run_maskA/".i."global_mask.txt" u ($0*dR-9):($1/4) w l title "" lc black lw 2,\
"../BOPES/2_bopes_2k.dat" u ($0*dR-9):($2) w l title "" lc "gray"  ,\
"../BOPES/1_bopes_2k.dat" u ($0*dR-9):($2) w l title ""  lc "gray" ,\
"../output_git_run_maskA/".i."nd.txt" u ($0*dR-9):($1/10) w l title ""   lc "red" lw 3 ,\
ali_folder.i."_tdnd.dat"  u ($0*dR-9):($1/10)  w l title "" dashtype 4 lc "black" 



################# bottom row ################# 
i = "100"

set title "t = ".(i*10)."a.u."
unset format y
plot [:][0:0.3]"../output_git_run_maskA/".i."TDPES.txt" u ($0*dR-9):($1) w l title "" lw 3 lc "blue" ,\
"../output_git_run_maskA/".i."global_mask.txt" u ($0*dR-9):($1/4) w l title "" lc black lw 2,\
"../BOPES/2_bopes_2k.dat" u ($0*dR-9):($2) w l title "" lc "gray"  ,\
"../BOPES/1_bopes_2k.dat" u ($0*dR-9):($2) w l title ""  lc "gray" ,\
"../output_git_run_maskA/".i."nd.txt" u ($0*dR-9):($1/10) w l title ""   lc "red" lw 3,\
ali_folder.i."_tdnd.dat"  u ($0*dR-9):($1/10)  w l title "" dashtype 4 lc "black" 



set format y ""
i = "120"
set title "t = ".(i*10)."a.u."
plot [:][0:0.3]"../output_git_run_maskA/".i."TDPES.txt" u ($0*dR-9):($1) w l title "" lw 3 lc "blue" ,\
"../output_git_run_maskA/".i."global_mask.txt" u ($0*dR-9):($1/4) w l title "" lc black lw 2,\
"../BOPES/2_bopes_2k.dat" u ($0*dR-9):($2) w l title "" lc "gray"  ,\
"../BOPES/1_bopes_2k.dat" u ($0*dR-9):($2) w l title ""  lc "gray" ,\
"../output_git_run_maskA/".i."nd.txt" u ($0*dR-9):($1/10) w l title ""   lc "red" lw 3,\
ali_folder.i."_tdnd.dat"  u ($0*dR-9):($1/10)  w l title "" dashtype 4 lc "black" 



i = "149"
set title "t = ".(i*10)."a.u."
plot [:][0:0.3]"../output_git_run_maskA/".i."TDPES.txt" u ($0*dR-9):($1) w l title "" lw 3 lc "blue" ,\
"../output_git_run_maskA/".i."global_mask.txt" u ($0*dR-9):($1/4) w l title "" lc black lw 2,\
"../BOPES/2_bopes_2k.dat" u ($0*dR-9):($2) w l title "" lc "gray"  ,\
"../BOPES/1_bopes_2k.dat" u ($0*dR-9):($2) w l title ""  lc "gray" ,\
"../output_git_run_maskA/".i."nd.txt" u ($0*dR-9):($1/10) w l title ""   lc "red" lw 3 ,\
ali_folder.i."_tdnd.dat"  u ($0*dR-9):($1/10)  w l title "" dashtype 4 lc "black" 


unset multiplot


##################################################################

#~ ##################################################################

set terminal cairolatex standalone pdf size 6.5in,4in  dashed transparent

set output "mask_gcoc_smalldt.tex" 
set key at 0.5, 1.1
set key samplen 2
set xlabel "R"

set xtics
set ytics
set format y
unset lmargin
set multiplot layout 1,2 title ""
set xtics -0.4,0.4, 5
set ytics 0, 0.2, 1.3
Rmin = -0.5
Rmax = 1.
t = 117
set title "\\LARGE $\\epsilon(R,t)$"
plot  [Rmin:Rmax]"../output_git_run_maskA/".t."TDPES.txt" u ($0*dR-9):($1) w l title "\\large $dt=0.1$ au" lc "blue" lw 2,\
 "../output_git_run_maskA_smalldt/".t."TDPES.txt" u ($0*dR-9):($1) w l  title "\\large $ dt=0.01$ au"  lc "red"   lw 2

#~ unset ytics
set ytics -30, 10, 30
#~ set ytics auto
set title "\\LARGE Re($\\nabla \\chi/\\chi$)"  
  plot  [Rmin:Rmax][-31:31]  "../output_git_run_maskA_smalldt/".t."GCOC.txt" u ($0*dR-9):($1) w l dashtype 1  lc "red" lw 2 title "",\
  "../output_git_run_maskA/".t."GCOC.txt" u ($0*dR-9):($1)  w l title "" lc "blue"  lw 2

#~ set title "\\LARGE Re($|\\chi|^2$)"  
  #~ plot  [Rmin:Rmax]  "../output_git_run_smalldt/117nd.txt" u ($0*dR-9):($1/30) w l dashtype 1  lc "red" lw 2 title "",\
  #~ "../output_git_run_maskA/117nd.txt" u ($0*dR-9):($1/30)  w l title "" lc "blue"  lw 2
  
unset multiplot
##################################################################


set terminal cairolatex standalone pdf size 4.5in,3.5in  dashed transparent
set output "mask_uen_grid2.tex" 
set title ""
set xtics -5, 1,5
set ytics auto
set key at -0.8,0.965
set title "$t = 790$a.u."
t = 790
plot  [-4.5:2][-0.05:1]	"../output_git_run_maskuen_tightermask/".t."nd.txt" u ($0*dR-9):($1) w l ls 10 title "$\\kappa = 10^{-5}$ \\tiny $(t_{max} = 798a.u.)$" ,\
"../output_git_run_maskuen_tightestmask/".t."nd.txt" u ($0*dR-9):($1) w l ls 11 title "$\\kappa = 10^{-3}$ \\tiny $(t_{max} = 878a.u.)$"  ,\
"../output_git_run_maskA/".(t/10)."TDPES.txt" u ($0*dR-9):($1) w l ls 16 title "" ,\
"../output_git_run_maskuen_tightermask/".t."TDPES.txt" u ($0*dR-9):($1) w l ls 14 title "" ,\
"../output_git_run_maskuen_tightestmask/".t."TDPES.txt" u ($0*dR-9):($1) w l ls 15 title "" ,\
"../output_git_run_maskA/".(t/10)."nd.txt" u ($0*dR-9):($1) w l ls 8 title "best case" 

	
##################################################################


##################################################################

set terminal cairolatex standalone pdf size 8in,3.5in  dashed dl 1 transparent
set output "linear_gcoc_grid.tex" 
set xtics
set ytics
set format y

set multiplot layout 1,3 #title " Artificial analytic GCOC at t = 10au"
set xtics -9, 3, 9
set key at 8.5, 73
set title "$t=0$ value"
plot  "../output_git_run_freechi_linear/1TDPES.txt" u ($0*dR-9):($1) w l ls 11 title "$\\epsilon(R,t)$" ,\
 "../output_git_run_freechi_linear/1GCOC.txt" u ($0*dR-9):($1) w l ls 10 title "Re$(\\nabla\\chi/\\chi)$" 

unset key
set title "RHS max"
plot  "../output_git_run_freechi_linear_RHS/1TDPES.txt" u ($0*dR-9):($1) w l ls 11  ,\
 "../output_git_run_freechi_linear_RHS/1GCOC.txt" u ($0*dR-9):($1) w l   ls 10

set title "centre max"
plot  [:][-5:]"../output_git_run_freechi_linear_gauss/1TDPES.txt" u ($0*dR-9):($1) w l ls 11 ,\
"../output_git_run_freechi_linear_gauss/1GCOC.txt" u ($0*dR-9):($1) w l  ls 10


unset multiplot

