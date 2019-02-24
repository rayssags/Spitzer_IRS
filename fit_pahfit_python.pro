pro fit_pahfit_python

CD, '/home/rayssa/Win/Programs/Spitzer/3_PAHfit_input/'
files = FILE_SEARCH('*.tbl')

foreach element, files, idx do begin
 name = STRMID(element, 0, STRLEN(element)-17)
 readcol, element, lambda, flux, error, skipline=1, delimiter="|"
 fit=pahfit(lambda,flux,error,/PLOT_PROGRESS,REPORT=name+'_pahfit_output.txt')
	
 readcol,'components_model.dat', lam_mod , cont_mod, feat_mod, lines_mod, stars_mod, dust_mod, bestfit_mod, ext_mod
 
 openw,1,'../4_PAHfit_output/'+name+'_mod_components.dat',WIDTH=250
 
 printf,1,'#',name
; printf,1,'#',decoded_params.reduced_chi_sq ;testar se não vai dar problema
 printf,1,'lam_mod','|','cont_mod','|','feat_mod','|','lines_mod', '|','stars_mod','|','dust_mod','|','bestfit_mod','|','ext_mod'

 for k=0,n_elements(lam_mod)-1 do begin
   printf,1, lam_mod[k],'|', cont_mod[k],'|', feat_mod[k],'|', lines_mod[k],'|', stars_mod[k],'|', dust_mod[k],'|', bestfit_mod[k],'|', ext_mod[k]
  endfor
 
 close,1
 
endforeach
 
end
