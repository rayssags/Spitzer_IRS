; PROGRAM TO FIT THE CONTINUUM AND THE IONIC LINES OF THE SELECTED SOURCES 
;
; The sources were extracted from A. Hernán-Caballero and E. Hatziminaoglou (2011).
; Site: http://www.denebola.org/atlas/?p=data
;
; IMPORTANTE: The sources must have passed through redshift's correction before using
; this program, what can be done with check_and_correct_drd5.py
;
; The code used to fit is the PAHFIT. Reference: Smith, J.D.T., Draine B.T., et al.,
; 2007, ApJ, 656, 770
; One of its programs was modified by Anelise Audibert - pahfit_plot
;
; PARAMETER: directory's path of the sources (placed at line 20) 
;
; Code written by Carla Martinez Canelo- 2015, July
; (based on runall from Anelise Audibert) 

pro fit_pahfit

CD, '/home/rayssa/Win/Programs/Spitzer/'
files = FILE_SEARCH('*.tbl')

foreach element, files, idx do begin
 name = STRMID(element, 0, STRLEN(element)-4)
 readcol, element, lambda, flux, error, delimiter="|", skipline=1
 fit=pahfit(lambda,flux,error,/PLOT_PROGRESS,REPORT=name+'_pahfit.txt')
; plot, lambda, flux/lambda
	
 readcol,'components_model.dat', lam_mod , cont_mod, feat_mod, lines_mod, stars_mod, dust_mod, bestfit_mod, ext_mod
 openw,1,name+'_mod_components.dat',WIDTH=250
 printf,1,'#',name
 printf,1,'#','lam_mod',STRING(9B),'cont_mod',STRING(9B),'feat_mod',STRING(9B),'lines_mod', $
              STRING(9B),'stars_mod',STRING(9B),'dust_mod',STRING(9B),'bestfit_mod',STRING(9B),'ext_mod'

 for k=0,n_elements(lam_mod)-1 do begin
   printf,1, lam_mod[k], cont_mod[k], feat_mod[k], lines_mod[k], stars_mod[k], dust_mod[k], bestfit_mod[k], ext_mod[k]
  endfor
 
 close,1
 
endforeach
 
end 
