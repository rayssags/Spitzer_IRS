;+
; NAME:
;
;   PAHFIT
;
; VERSION:
;
;   1.2 (Oct, 2006)
;
; PURPOSE:
;
;   Decompose low-resolution Spitzer IRS Spectra into dust continuua,
;   PAH and dust features, starlight, atomic and molecular lines, and
;   silicate absorption, returning the fit parameters for each.
;
; REFERENCE:
; 
;   Smith, J.D., Draine B.T., et al., ApJ, 2006 (astro-ph/0610913)
;
; CATEGORY:
;
;   IRS Spectral fitting, PAH Emission, Silicate
;
; CALLING SEQUENCE:
;
;   fit=pahfit(lambda, flux, flux_uncertainty,[REDSHIFT=,
;              STARLIGHT_TEMPERATURE=,
;              CONTINUUM_TEMPERATURES=, LINES=, DUST_FEATURES=,
;              /PLOT_PROGRESS, /NO_EXTINCTION, /NO_FIT,PARINFO=,/NO_FIX,
;              /VALUES_ONLY, /NO_MEGAJANSKY_SR,REPORT=,_EXTRA=e])
;
; INPUTS:
;
;    lambda: The wavelength in microns, IN THE OBSERVED FRAME (see
;       REDSHIFT).  This should ideally span the range ~5-~35um,
;       i.e. IRS SL+LL.
;
;    flux: The flux intensity in MJy/sr (or another suitably scaled
;       f_nu unit).  The spectrum must be continguous, pre-stitched
;       from the 4 IRS low-resolution individual order spectra.
;
; OPTIONAL INPUTS:
;
;    flux_uncertainty: The flux uncertainty in the same units as the
;       flux.  If no uncertainties are passed, equal weighting is
;       assigned each point, and no formal parameter uncertainties are
;       returned.
;
; KEYWORD PARAMETERS:
;
;    REDSHIFT: The redshift cz in km/s.  If not passed, the redshift
;       is assumed to be 0.0, which may adversely affect the fit.
;
;    STARLIGHT_TEMPERATURE: Temperature to use for the effective
;       stellar continuum blackbody (default=5000K).
;
;    CONTINUUM_TEMPERATURES: Array of continuum temperatures for
;       modified blackbodies to fit.
;       (default=300.D,200.D,135.D,90.D,65.D,50.D,40.D,35.D)
;
;    LINES: Array of structures for atomic lines to fit.  Each
;       structure should be of the form {WAVELENGTH:0.0D, NAME:''}
;       where name is the line name (e.g. [NeII]).  Defaults to a
;       range of lines of H_2, Ar, S, Ne, O, Si, and Fe.
;
;    DUST_FEATURES: Array of structures for resolved dust features to
;       fit.  Each structure should be of the form {WAVELENGTH: 0.0D,
;       FRAC_FWHM: 0.0D}, where wavelength is the central wavelength,
;       and FRAC_FWHM is the fractional FWHM of the Drude profile.
;       Defaults to a collection of features fitting the most common
;       PAH bands (see Smith, Draine, et al., 2006), tuned to a high
;       S/N sample of SINGS galaxies.
;
;    PLOT_PROGRESS: If set, update the progress of the fitter function
;       by plotting the spectrum with fit intead of printing the
;       Chi-Sq. Pass XSIZE= and YSIZE= to override the default
;       progress window size.
;
;    NO_EXTINCTION: Fix the extinction optical depth to 0.0.
;
;    SCREEN_EXTINCTION: Rather than the default fully-mixed extinction
;       geometry, assume a foreground screen.
;
;    CHIAR_TIELENS: Use the galactic center extinction curve of Chiar
;       & Tielens, 2006, which extends to ~30um.
;
;    KEMPER_VRIEND_TIELENS: Use the silicate profile of Kemper,
;       Vriend, & Tielens, integrated into an extinction curve with an
;       exponent 1.7 power-law and artificial 18um silicate feature
;       Drude profile (default).
;
;    NO_FIT: Don't actually run a fit, just populate PARINFO.  Useful
;       to return the default PARINFO array.
;
;    PARINFO: The parameter info structure array, as input/output.  If
;       passed in, this structure will be used to seed the starting
;       conditions, and (optionally) fix or limit individual
;       parameters, presumably differently than the default.  Keywords
;       which control how PARINFO is interpreted are VALUES_ONLY, and
;       NO_FIX.  On output, PARINFO will contain the updated parameter
;       structure.  This provides a powerful means to alter the
;       fitting conditions, zero out or fix individual parameters,
;       place different limits on parameter values, etc.  Note that
;       parameters are assigned by name, and only parameters which
;       match by name will be updated.  Other parameters which do not
;       match (either the default set, or a modified set, passed in
;       with LINES/DUST_FEATURE/CONTINUUM_TEMPERATURES) will be
;       included unchanged in the fit.  Since dust feature names
;       include the wavelength, updates to the default dust feature
;       set will necessarily result in parameters which do not get
;       set.  The solution to this is to first run PAHFIT with NO_FIT
;       to update the parinfo structure, and then modify the structure
;       to suit, calling PAHFIT again for the real fit.
;
;     NO_FIX: When passing a parinfo structure, disregard any
;       differences in whether a parameter's value is fixed.
;
;     VALUES_ONLY: When passing a parinfo structure in, only set the
;       values of parameters, not whether they are limited, whether
;       they are fixed, and the limits themselves.  Starting values
;       are clipped to the limiting range, if set, for each parameter.
;
;     NO_MEGAJANSKY_SR: By default, input flux intensities are assumed
;       to have units of MJy/sr, and output line and feature strengths
;       are computed in W/m^2/sr, by integrating over the line or
;       feature (f_nu dnu).  If another f_nu unit is used instead, set
;       this keyword to avoid incorrectly converting the units.  See
;       Restrictions, below.
;
;     REPORT: If set to a filename, generate a text report of the
;        fitting outcome, in the given filename.  If set otherwise
;        (e.g. /REPORT), print the report in the IDL console.
;
;     _EXTRA: Other keyword parameters to MPFITFUN can be passed.  In
;        particular, passing the YFIT keyword a named variable will
;        return the best fitting model.
;      
; OUTPUT:
;
;    The decoded fit parameters in the rest frame of the target, as a
;    nested structure with the following members:
;
;      STARLIGHT: A structure with two (four, with uncertainties) tags:
;         TEMPERATURE(_UNC): The stellar continuum temperature used (K).
;         TAU(_UNC): The effective optical depth, roughly the fraction
;            of the solid angle covered by stellar surfaces (unitless).
;
;      DUST_CONTINUUM: An array of structures with two (four, with
;         uncertainties) tags:
;         TEMPERATURE(_UNC): The (fixed) input temperatures for each
;            dust continuum component (K).
;         TAU(_UNC): Dust component scaling parameters (unitless).
;
;      DUST_FEATURES: An array of structures for each dust feature,
;         with five (nine, with uncertainties) tags:
;         WAVELENGTH(_UNC): The (fixed) central wavelength of the
;            feature (microns).
;         CENTRAL_INTEN(_UNC): The central intensity of the Drude
;            profile feature (MJy/sr, or other input units of FLUX
;            above).
;         FWHM(_UNC): The (fixed) fractional full width half-maximum
;            of the feature (no units), relative to the central
;            wavelength.
;         INT_STRENGTH(_UNC): The integrated feature strength, in
;            W/m^2/sr for input units of MJy/sr.
;         EQW: The feature equivalent width, in microns.  Note that,
;            for broad features which can cover regions with no
;            continuum, the continuum is replaced by its
;            profile-averaged value in this calculation
;
;      LINES: An array of structures for each line, with five (nine,
;         with uncertainties) tags:
;         WAVELENGTH(_UNC): The (fixed) central wavelength of the
;            feature (microns).
;         CENTRAL_INTEN(_UNC): The central intensity of the Drude
;            profile (MJy/sr, or other input units of FLUX above).
;         FWHM(_UNC): The (fixed) fractional full width half-maximum
;            of the line (no units), relative to the central
;            wavelength.
;         INT_STRENGTH(_UNC): The integrated line strength, in
;            W/m^2/sr for input units of MJy/sr.
;         EQW: The line equivalent width, in microns.
;
;      EXTINCTION: A structure with two (four, with uncertainties)
;         tags:
;         TAU_9_7(_UNC): The optical depth of the mixed (or screen)
;            dust extinction, at 9.7um (the local extinction maximum).
;         TYPE: Either 0, for the Chiar's and Tielens based curve, or
;            1, for the Kemper et al.-based curve.
;         BETA(_UNC): The (fixed) mixing parameter between the
;            silicate profile and power-law extinction curve components.
; 
;      FINAL_FIT: The final fitted function, evaluated at each of the
;         *rest-wavelength* (not observed frame) points.
;         
;      USE_UNCERTAINTY: A boolean, indicating that, if set, each
;         parameter PARAM in the set above has a matching field named
;         PARAM_UNC which gives the formal uncertainty in the same
;         units.  Fixed parameters will have an uncertainy of 0.0.
;         Requires statistical uncertaintiest to be passed in as the
;         argument FLUX_UNCERTAINTY.
;
;      REDUCED_CHI_SQUARE: The final chi-square divided by degrees of
;         freedom.
;
;      COVARIANCE: The NxN covariance matrix (where N is the number of
;         parameters in the fit).  Necessary for the
;         pahfit_combine_feature_strength() function, to estimate
;         uncertainties of coadded feature strengths correctly for
;         correlated features.  Most off-diagonal entries will be
;         blank (e.g., for fixed parameters).
;
; COMMON BLOCKS:
; 
;    PAHFIT_PARAMS: Passes additional information for iteration
;       update/etc.  See pahfit_params.pro.
;
;
; RESTRICTIONS:
;
;   The passed spectrum must be in the observed frame, and the
;   returned fit is in the rest frame.
;   
;   Units: The native flux intensity units assumed are MJy/sr.  Any
;   other f_nu unit can be used, but the Starlight TAU parameter is
;   physically meaningful only for these units.  If other units are
;   used, pass NO_MEGAJANSKY_SR to avoid line strength units being
;   converted into W/m^2/sr (by multiplication by a factor of 1.e-20).
;   Instead, they will have units [UNIT]*Hz, where [UNIT] is the f_nu
;   unit passed.
;
; PROCEDURE:
;
;   Requires Craig Markwardt's MPFIT library
;   (http://cow.physics.wisc.edu/~craigm/idl/fitting.html)
;
; EXAMPLE:
;
;   fit=pahfit(rest_lam,flux,flux_unc,REDSHIFT=cz,/PLOT_PROGRESS,
;              XSIZE=1000,YSIZE=600, /SCREEN)
;
; MODIFICATION HISTORY:
;
;  2006-02-01 (JD Smith): Added Chiars & Tielens and screen dust geom.
;  2005-07-01 (JD Smith & Bruce Draine): Written.
;-
;##############################################################################
;
; LICENSE
;
;  Copyright (C) 2005-2006 J.D. Smith
;
;  This file is part of PAHFIT.
;
;  PAHFIT is free software; you can redistribute it and/or modify it
;  under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2, or (at your option)
;  any later version.
;
;  PAHFIT is distributed in the hope that it will be useful, but
;  WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;  General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with PAHFIT; see the file COPYING.  If not, write to the Free
;  Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
;  MA 02110-1301, USA.
;
;##############################################################################

;=========================================================================
;  pahfit_gaussian -- Return a Gaussian profile characterized by
;                     central wavelength, central intensity (in
;                     MJy/sr), and fractional FWHM.  
;=========================================================================
function pahfit_gaussian,lambda,lam_0,central_inten,frac_fwhm
  return,central_inten*exp(-(lambda-lam_0)^2*2.7725887222397811D/ $
                           (frac_fwhm*lam_0)^2)
end


;=========================================================================
;  pahfit_lorentzian -- Return a Lorentzian profile characterized by
;                       central wavelength, central intensity (in
;                       MJy/sr), and fractional FWHM.
;=========================================================================
function pahfit_lorentzian,lambda,lam_0,central_inten,frac_fwhm
  g=(frac_fwhm*lam_0/2.D)^2
  return,central_inten*g/((lambda-lam_0)^2+g)
end


;=========================================================================
;  pahfit_drude -- Return a Drude profile characterized by central
;                  wavelength, central intensity (in MJy/sr), and
;                  fractional FWHM.
;=========================================================================
function pahfit_drude,lambda,lam_0,central_inten,frac_fwhm
  return,central_inten*frac_fwhm^2/((lambda/lam_0-lam_0/lambda)^2+frac_fwhm^2)
end


;=========================================================================
;  pahfit_feature_strength -- Compute strength of one or more Drude or
;                             Gaussian profiles, passed as a structure
;                             or separarate central wavelength (um),
;                             intensity and fwhm arguments.  If
;                             MEGAJANSKY_SR is passed, then the return
;                             units are W/m^2/sr.
;=========================================================================
function pahfit_feature_strength, lam_0_in,central_inten,frac_fwhm, $
                                  MEGAJANSKY_SR=mj,LAM_O_UNC=lam_0_unc, $
                                  CENTRAL_INTEN_UNC=central_inten_unc, $
                                  FRAC_FWHM_UNC=frac_fwhm_unc,UNCERTAINTY=unc,$
                                  DRUDE=drude,GAUSSIAN=gaussian
  if n_params() eq 1 then begin 
     ;; Assume a structure with all the appropriate keywords.
     st=lam_0_in
     if size(st,/TYPE) ne 8 then $
        message,'Must pass decoded structure or individual Drude parameters.'
     lam_0=st.wavelength
     central_inten=st.central_inten
     frac_fwhm=st.fwhm
     if arg_present(unc) then begin 
        lam_0_unc=st.wavelength_unc
        central_inten_unc=st.central_inten_unc
        frac_fwhm_unc=st.fwhm_unc
     endif 
  endif else lam_0=lam_0_in
  
  ;; Factor of ~0.67 between the Gaussian and Drude
  strength=(2.9979246D14/lam_0)/2*central_inten*frac_fwhm * $
           (keyword_set(drude)?!DPI:sqrt(!DPI/ALOG(2.D)))

  if arg_present(unc) && n_elements(lam_0_unc) ne 0 then $
        unc=sqrt(lam_0_unc^2 * (strength/lam_0)^2 + $
                 central_inten_unc^2 * (strength/central_inten)^2 + $
                 frac_fwhm_unc^2 * (strength/frac_fwhm)^2)
     
  if keyword_set(mj) then begin
     strength*=1D-20            ;return units are now W/m^2(/sr)
     if n_elements(unc) ne 0 then unc*=1D-20
  endif 
  return,strength
end


;=========================================================================
;  pahfit_blackbody_nu -- Return the blackbody flux intensity B_nu at
;                         the given temperature in MJy/sr
;=========================================================================
function pahfit_blackbody_nu,lambda,T
  return,3.97289D13/lambda^3/(exp(1.4387752D4/lambda/T)-1.D)
end


;=========================================================================
;  pahfit_exinction -- Return the overall extinction factor in a
;                      uniformly mixed model, or, optionally a
;                      foreground SCREEN model.
;=========================================================================
function pahfit_extinction,lambda,tau_9_7,ext_curve,SCREEN=screen
  if tau_9_7 eq 0.0 then return,replicate(1.0D, n_elements(lambda))
  tau_lam=tau_9_7*ext_curve
  if keyword_set(screen) then return,exp(-tau_lam) else $
     return,(1.D - exp(-tau_lam))/tau_lam
end


;=========================================================================
;  pahfit_tau_ct -- Extinction curve from the galactic center
;                   measurements of Chiar & Tielens, 2006, ApJ, 637,
;                   774.  Call via PAHFIT_TAU.
;=========================================================================
function pahfit_tau_ct,lambda
  dir=file_dirname((routine_info('PAHFIT',/SOURCE,/FUNCTIONS)).PATH)
  readfmt,filepath(ROOT=dir,'chiar+tielens_2005.dat'),$
          'F6.3,1X,F5.3',lam_ct,a_gcen,SKIPLINE=14
  
  ext=dblarr(n_elements(lambda),/NOZERO)
  lam_mx=max(lam_ct,mpos)
  wh=where(lambda gt lam_mx,cnt,COMPLEMENT=wh2,NCOMPLEMENT=cnt2)
  if cnt gt 0 then $
     ext[wh]=a_gcen[mpos]*(lam_mx/lambda[wh])^1.7D
  
  if cnt2 gt 0 then $
     ext[wh2]=interpol(a_gcen,lam_ct,lambda[wh2])
  
  ext/=interpol(ext,lambda,9.7) ;normalize to 9.7um
  return,ext
end


;=========================================================================
;  pahfit_tau_kvt -- Custom extinction curve built from silicate
;                    profile of Kemper, Vriend, & Tielens, 2004, ApJ ,
;                    609, 826.  Constructed as a weighted sum of two
;                    components: silicate profile, and an exponent 1.7
;                    power-law.
;=========================================================================
function pahfit_tau_kvt,lambda,beta
  ;; Relative silicate profile from 
  kvt_prof=[[8.0 , .06   ], $
            [8.2 , .09   ], $
            [8.4 , .16   ], $
            [8.6 , .275  ], $
            [8.8 , .415  ], $
            [9.0 , .575  ], $
            [9.2 , .755  ], $
            [9.4 , .895  ], $
            [9.6 , .98   ], $
            [9.7 , .99   ], $
            [9.75, 1.00  ], $   ;this profile has fwhm = 11.03 - 8.71 = 2.32um
            [9.8 , .99   ], $
            [10.0, .94   ], $
            [10.2, .83   ], $
            [10.4, .745  ], $
            [10.6, .655  ], $
            [10.8, .58   ], $
            [11.0, .525  ], $
            [11.2, .43   ], $
            [11.4, .35   ], $
            [11.6, .27   ], $
            [11.8, .20   ], $
            [12.0, .13   ], $ 
            [12.2, .09   ], $ 
            [12.4, .06   ], $ 
            [12.6, .045  ], $ 
            [12.7, .04314]]
  
  ext=dblarr(n_elements(lambda))
  
  kvt_mn=min(kvt_prof[0,*],MAX=kvt_mx,SUBSCRIPT_MAX=kvt_mx_loc,kvt_mn_loc)
  lowwh=where(lambda lt kvt_mn,lowcnt)
  if lowcnt gt 0 then $
     ext[lowwh]=kvt_prof[1,kvt_mn_loc]*exp(2.03D*(lambda[lowwh]-kvt_mn))
  
  on_prof=where(lambda gt kvt_mn AND lambda lt kvt_mx,cnt)
  if cnt ne 0 then $
     ext[on_prof]=interpol(kvt_prof[1,*],kvt_prof[0,*],lambda[on_prof])
  
  ;; Extend the KVT profile by fitting a quadratic
  fit=where(lambda ge kvt_mx-2.5 AND lambda le kvt_mx)
  anchor=where(lambda ge kvt_mx+2. AND lambda le kvt_mx+3.)
  fit=[fit,anchor]
  
  extend=where(lambda gt kvt_mx AND lambda lt kvt_mx+2.,cnt)
  if cnt gt 0 then ext[extend]=interpol(ext[fit],lambda[fit],lambda[extend],$
                                        /LSQUADRATIC)
  
  ;if cnt gt 0 then $
  ;   ext[between]=spline(lambda[[o1,o2]],ext[[o1,o2]],lambda[between])
  ;interpol(ext[[o1,o2]],lambda[[o1,o2]],lambda[between], $
  ;                         /LSQUADRATIC)

  ;; Add a Drude profile for the 18um silicate features
  ext=ext>0.0D + pahfit_drude(lambda,18.,.4,.247)
  
  ;; Form linear combination of modified silicate and powerlaw.
  ext=(1.D - beta)*ext + beta*(9.7D/lambda)^1.7D
  return,ext
end


;=========================================================================
;  pahfit_tau -- Return the extinction curve interpolated from the
;                saved fixed curve, which is recomputed if necessary.
;=========================================================================
function pahfit_tau,lambda,beta, $
                    CHIAR_TIELENS=ct, $
                    KEMPER_VRIEND_TIELENS=kvt
  common pahfit_tau_save, pahfit_tau_type,pahfit_kvt_beta, $
     pahfit_tau_lambda, pahfit_tau
  
  type=keyword_set(ct)
  if n_elements(pahfit_tau_type) eq 0 || $
     type ne pahfit_tau_type || $
     (type eq 0 && n_elements(beta) gt 0 && $
      n_elements(pahfit_kvt_beta) gt 0 && pahfit_kvt_beta ne beta) $
  then begin 
     pahfit_tau_type=type
     ;; Create tau curve on a well-sampled grid from 2--40um
     if n_elements(pahfit_tau_lambda) eq 0 then $
        pahfit_tau_lambda=2.+findgen(500)/499.*38.
     if type eq 1 then begin 
        pahfit_tau=pahfit_tau_ct(pahfit_tau_lambda)
     endif else begin 
        if n_elements(beta) ne 0 then pahfit_kvt_beta=beta
        pahfit_tau=pahfit_tau_kvt(pahfit_tau_lambda,pahfit_kvt_beta)
     endelse 
  endif 
  return,interpol(pahfit_tau,pahfit_tau_lambda,lambda)
end


;=========================================================================
;  pahfit_components -- Compute and return via keywords the component
;                       contributions to a given decoded model
;                       parameter set, and return the full model.
;                       EXTINCTION_FAC is the extinction pre-factor
;                       multiplicative term.  Input EXTINCTION_CURVE,
;                       or it will be recomputed based on the
;                       preferences.
;=========================================================================
function pahfit_components,lambda,decoded_params, $
                           DUST_CONTINUUM=dust_continuum, $
                           TOTAL_DUST_CONTINUUM=dc_tot,STARLIGHT=stars, $
                           TOTAL_CONTINUUM=tot_cont, $
                           DUST_FEATURES=dust_features, $
                           TOTAL_DUST_FEATURES=df_tot, $
                           LINES=lines,TOTAL_LINES=lines_tot, $
                           EXTINCTION_FAC=ext, $
                           EXTINCTION_CURVE=ext_curve, $
                           DISABLE_EXTINCTION=de
  ;; Extinction
  ext_pars=decoded_params.extinction
  
  if ~keyword_set(ext_curve) then $
     ext_curve=pahfit_tau(lambda,ext_pars.beta, $
                          CHIAR_TIELENS=ext_pars.type eq 1, $
                          KEMPER_VRIEND_TIELENS=ext_pars.type eq 0)
  
  ext=pahfit_extinction(lambda,ext_pars.tau_9_7,ext_curve, $
                        SCREEN=ext_pars.screen)
  
  do_ext=~keyword_set(de) 
  ;; Starlight
  star_par=decoded_params.starlight
  stars=star_par.tau*pahfit_blackbody_nu(lambda,star_par.temperature)
  if do_ext then stars*=ext
  
  ;; Dust continuum
  dc_pars=decoded_params.dust_continuum
  dust_continuum=make_array(n_elements(lambda),n_elements(dc_pars), $
                           /DOUBLE,/NOZERO)
  for i=0,n_elements(dc_pars)-1 do begin 
     this_dc=dc_pars[i].tau*(9.7D/lambda)^2* $
             pahfit_blackbody_nu(lambda,dc_pars[i].temperature)
     if do_ext then this_dc*=ext
     dust_continuum[0,i]=this_dc
     if i eq 0 then dc_tot=this_dc else dc_tot+=this_dc
  endfor 
  
  tot_cont=stars+dc_tot
  
  ;; Lines
  line_pars=decoded_params.lines
  lines=make_array(n_elements(lambda),n_elements(line_pars), $
                   /DOUBLE,/NOZERO)
  for i=0,n_elements(line_pars)-1 do begin 
     this_line=pahfit_gaussian(lambda,line_pars[i].wavelength, $
                               line_pars[i].central_inten, $
                               line_pars[i].fwhm)
     if do_ext then this_line*=ext
     lines[0,i]=this_line
     if i eq 0 then lines_tot=this_line else lines_tot+=this_line
  endfor 
  
  ;; Dust Features
  df_pars=decoded_params.dust_features
  dust_features=make_array(n_elements(lambda),n_elements(df_pars), $
                            /DOUBLE,/NOZERO)
  for i=0,n_elements(df_pars)-1 do begin 
     this_df=pahfit_drude(lambda,df_pars[i].wavelength, $
                          df_pars[i].central_inten,df_pars[i].fwhm)
     if do_ext then this_df*=ext
     dust_features[0,i]=this_df
     if i eq 0 then df_tot=this_df else df_tot+=this_df
  endfor 
  
  return,stars+dc_tot+lines_tot+df_tot
end


;=========================================================================
;  pahfit_feature_eqw -- Return the equivalent width
;                        (integral[(I_nu-I_cont)/I_cont d_lam]) in um
;                        of the passed profile (Drude or Lorentzian).
;                        Set BROAD_LIMIT to replace the continuum in
;                        the EQW calculation by its profile-averaged
;                        value for fractional FWHM greater than that
;                        limit.  Useful when the EQW requires
;                        extrapolation to regions where the continuum
;                        can vanish, such that f_line/f_continuum
;                        diverges.  Default value BROAD_LIMIT=.05
;=========================================================================
function pahfit_feature_eqw,decoded_params,lam_0_in,central_inten,frac_fwhm, $
                          BROAD_LIMIT=bl,DRUDE=drude
  
  if n_params() le 2 then begin 
     ;; Assume a structure with all the appropriate keywords.
     st=lam_0_in
     if size(st,/TYPE) ne 8 then $
        message,'Must pass decoded structure or individual Drude parameters.'
     lam_0=st.wavelength
     central_inten=st.central_inten
     frac_fwhm=st.fwhm
     if arg_present(unc) then begin 
        lam_0_unc=st.wavelength_unc
        central_inten_unc=st.central_inten_unc
        frac_fwhm_unc=st.fwhm_unc
     endif 
  endif else lam_0=lam_0_in
  
  if n_elements(bl) eq 0 then bl=.05D
  
  nlines=n_elements(lam_0) 
  eqw=fltarr(nlines)
  
  drude=keyword_set(drude) 
  lmin=(lam_0-frac_fwhm*lam_0*6)>0.
  lmax=lam_0+frac_fwhm*lam_0*6
  
  for i=0,nlines-1 do begin 
     lam=findgen(100)/99*(lmax[i]-lmin[i])+lmin[i]
     yfit=pahfit_components(lam,decoded_params,TOTAL_CONTINUUM=continuum)
     if drude then $
        lnu=pahfit_drude(lam,lam_0[i],central_inten[i],frac_fwhm[i]) $
     else  $
        lnu=pahfit_gaussian(lam,lam_0[i],central_inten[i],frac_fwhm[i])
     
     if n_elements(bl) gt 0 && frac_fwhm[i] gt bl then begin 
        ilam=int_tabulated(lam,lnu)
        weighted_cont=int_tabulated(lam,lnu*continuum)/ilam
        eqw[i]=ilam/weighted_cont
     endif else eqw[i]=int_tabulated(lam, lnu/continuum)
  endfor 
  return,eqw
end 


;=========================================================================
;  pahfit_function -- The full model fitting function; return the
;                     model given the wavelength vector and parameter
;                     set, parameter counts and adopted silicate
;                     profile.
;=========================================================================
function pahfit_function, lambda, P
  @pahfit_params
  
  ;; parameter organization: [stars, dust continua, lines,
  ;;                          dust features, extinction]
  
  ;; starlight (B_nu)
  model=P[0]*pahfit_blackbody_nu(lambda,P[1])
  off=2
  
  ;; add dust continua (nu^2 B_nu)
  for i=0,n_cont_temps-1 do begin 
     model+=P[off]*(9.7D/lambda)^2*pahfit_blackbody_nu(lambda,P[off+1])
     off+=2
  endfor 
  
  ;; Add lines and dust features (Drude: lam_central, central_inten, fwhm)
  for i=0,n_lines-1 do begin 
     model+=pahfit_gaussian(lambda,P[off],P[off+1],P[off+2])
     off+=3
  endfor 
  for i=0,n_dust_features-1 do begin 
     model+=pahfit_drude(lambda,P[off],P[off+1],P[off+2])
     off+=3
  endfor 

  ;; Cut by extinction (mixed, or screen)
  if P[off] ne 0.0D then $      ;tau_9.7, with pre-computed ext_curve
     model*=pahfit_extinction(lambda,P[off],ext_curve, $
                              SCREEN=P[off+3])
  
  return,model
end


;=========================================================================
;  pahfit_add_extra_info -- Supplement the decoded structure with
;                           additional information.
;=========================================================================
pro pahfit_add_extra_info,decoded,RED_CHI_SQ=rcsq, $
                          COVAR=covar,STRENGTH_EQW=st_eqw,FINAL_FIT=ff, $
                          _EXTRA=e
  ;; Add Chi-Square, covariance, strengths, and equivalent widths, if
  ;; asked for
  if n_elements(ff) ne 0 then decoded=create_struct(decoded,'FINAL_FIT',ff)
  if n_elements(rcsq) ne 0 then $
     decoded=create_struct(decoded,'REDUCED_CHI_SQ',rcsq)

  if keyword_set(st_eqw) then pahfit_add_strength_eqw,decoded,_EXTRA=e
  
  if n_elements(covar) gt 0 then $
     decoded=create_struct(decoded,'COVARIANCE',covar)
end


;=========================================================================
;  pahfit_decode -- Decode the raw parameter vector
;=========================================================================
function pahfit_decode,p,parinfo,line_names,PERROR=perr,COVAR=covar

  err=n_elements(perr) gt 0
  cov=n_elements(covar) gt 0
  
  ;; Starlight
  wh=where(parinfo.PAHFIT_TYPE eq 0b)
  stars={temperature:p[wh[1]],tau:p[wh[0]]}
  if err then stars=create_struct(stars,'temperature_unc',perr[wh[1]], $
                                  'tau_unc',perr[wh[0]])
  ;; Dust continuum
  wh=where(parinfo.PAHFIT_TYPE eq 1b,nc)
  nc/=2
  dust={temperature:0.0D,tau:0.0D}
  if err then dust=create_struct(dust,'temperature_unc',0.0D, $
                                 'tau_unc',0.0D)
  dust=replicate(dust,nc)
  inds=indgen(nc)*2
  dust.tau=p[wh[inds]]
  dust.temperature=p[wh[inds+1]]
  if err then begin 
     dust.tau_unc=perr[wh[inds]]
     dust.temperature_unc=perr[wh[inds+1]]
  endif 
  
  ;; Lines
  wh=where(parinfo.PAHFIT_TYPE eq 2b,nl)
  nl/=3
  lines={name:'',wavelength:0.0D,central_inten:0.0D,fwhm:0.0D, $
         int_strength:0.0D,eqw:0.0D}
  if err then lines=create_struct(lines,'wavelength_unc',0.0D, $
                                  'central_inten_unc',0.0D,'fwhm_unc',0.0D, $
                                  'int_strength_unc',0.0D)
  
  if cov then lines=create_struct(lines,'parinfo_covar_index',0L)
  
  lines=replicate(lines,nl)
  inds=indgen(nl)*3
  lines.wavelength=p[wh[inds]]
  lines.name=line_names
  lines.central_inten=p[wh[inds+1]]
  if cov then lines.parinfo_covar_index=wh[inds+1] 
    
  lines.fwhm=p[wh[inds+2]]
  if err then begin 
     lines.wavelength_unc=perr[wh[inds]]
     lines.central_inten_unc=perr[wh[inds+1]]
     lines_fwhm_unc=perr[wh[inds+2]]
  endif 
  
  ;; Dust features
  wh=where(parinfo.PAHFIT_TYPE eq 3b,nf)
  nf/=3
  features={wavelength:0.0D,central_inten:0.0D,fwhm:0.0D,int_strength:0.0D, $
            eqw:0.0D}
  if err then features=create_struct(features,'wavelength_unc',0.0D, $
                                     'central_inten_unc',0.0D,'fwhm_unc',0.0D,$
                                     'int_strength_unc',0.0D)
  if cov then features=create_struct(features,'parinfo_covar_index',0L)
     
  features=replicate(features,nf)
  inds=indgen(nf)*3
  features.wavelength=p[wh[inds]]
  features.central_inten=p[wh[inds+1]]
  if cov then features.parinfo_covar_index=wh[inds+1]
  features.fwhm=p[wh[inds+2]]
  if err then begin
     features.wavelength_unc=perr[wh[inds]]
     features.central_inten_unc=perr[wh[inds+1]]
     features.fwhm_unc=perr[wh[inds+2]]
  endif
  
  ;; Extinction
  wh=where(parinfo.PAHFIT_TYPE eq 4b)
  extinction={tau_9_7:p[wh[0]],beta:p[wh[1]],type:p[wh[2]],screen:p[wh[3]]}
  if err then $
     extinction=create_struct(extinction,'tau_9_7_unc', perr[wh[0]], $
                              'beta_unc',perr[wh[1]])

  return,{STARLIGHT:stars, DUST_CONTINUUM:dust, LINES:lines, $
          DUST_FEATURES:features, EXTINCTION: extinction, $
          USE_UNCERTAINTY:err}
end


;=========================================================================
;  pahfit_add_strength_eqw -- Enhance a decoded structure by adding
;                             strength (+strength_unc), and equivalent
;                             width to the individual line and dust
;                             features.
;=========================================================================
pro pahfit_add_strength_eqw,decoded,_EXTRA=e
  for j=2,3 do begin 
     dust=j eq 3                ;drude's for dust features
     ;; Strength for all the features
     strength=pahfit_feature_strength(decoded.(j),UNCERTAINTY=strength_unc, $
                                      DRUDE=dust,_EXTRA=e)
     
     ;; Widths for all the features
     eqw=pahfit_feature_eqw(decoded,decoded.(j),DRUDE=drude)
     
     decoded.(j).INT_STRENGTH=strength
     if n_elements(strength_unc) ne 0 then $
        decoded.(j).INT_STRENGTH_UNC=strength_unc
     decoded.(j).EQW=eqw
  endfor 
end


;=========================================================================
;  pahfit_iterproc -- Report on the fitter's progress
;=========================================================================
pro pahfit_iterproc,funct, pars, iteration, chi_sq, QUIET=quiet, $
                    PARINFO=parinfo,_EXTRA=e
  @pahfit_params
  if plot_progress then begin 
     wset,plot_dbwin
     decoded=pahfit_decode(pars,parinfo,line_names)
     title=string(FORMAT='(%"Iteration %3d, Chi-Sq=%9.3g, Tau_9_7=%9.3g")', $
                  iteration, chi_sq, decoded.extinction.tau_9_7)
     pahfit_plot,decoded,pf_lambda, pf_intensity, pf_errors, TITLE=title, $
                 /FAST_PLOT,_EXTRA=e
     wset,plot_win
     device,COPY=[0,0,!D.X_SIZE,!D.Y_SIZE,0,0,plot_dbwin] 
  endif else begin 
     if keyword_set(quiet) then return
     print,FORMAT='(%"Iteration %3d: Chi-Sq=%8.3g")',iteration,chi_sq
  endelse 
end


;=========================================================================
;  pahfit_report -- Print a report summarizing the fit parameters.
;=========================================================================
pro pahfit_report,decoded_params,filename
  info=['===============================================================',$
        ' PAHFIT v1.2 JD Smith & Bruce Draine', $
        '  '+systime(), $
        '===============================================================']
  info=[info, 'Reduced Chi-Square (statistical errors only): '+ $
        string(FORMAT='(F8.3)',decoded_params.reduced_chi_sq)]
  
  info=[info,string(10b),'Starlight: ']
  stars=decoded_params.starlight
  
  unc=decoded_params.USE_UNCERTAINTY
  if unc then f='%8.3g ( %8.3g )' else f='%8.3g'
  if unc then lam_f='%8.3f ( %8.3f )' else lam_f='%6.3g'

  form='(%"  T_star: %8.3g  tau_star: '+f+'")'
  if unc then $
     info=[info,string(FORMAT=form,stars.temperature,stars.tau,stars.tau_unc)] $
  else info=[info,string(FORMAT=form,stars.temperature,stars.tau)]
  
  dust=decoded_params.dust_continuum
  form='(%"  T_dust: %8.3g  tau_dust: '+f+'")'
  info=[info,string(10b),'Dust Continuum: ']
  for i=0,n_elements(dust)-1 do begin 
     if unc then $
        info=[info,string(FORMAT=form,dust[i].temperature, $
                          dust[i].tau,dust[i].tau_unc)] $
     else info=[info,string(FORMAT=form,dust[i].temperature,dust[i].tau)]
  endfor
  
  lines=decoded_params.lines
  form='(%"  %11s lam: ' + lam_f + ' cen_inten: ' + f + $
       '\n             fwhm: '+f+'     power: '+f+ $
       '\n              eqw: %8.3g")'
  
  info=[info,string(10b),'Lines: ']
  for i=0,n_elements(lines)-1 do begin 
     if unc then begin 
        info=[info, $
              string(FORMAT=form,lines[i].name, $
                     lines[i].wavelength,lines[i].wavelength_unc, $
                     lines[i].central_inten, lines[i].central_inten_unc, $
                     lines[i].fwhm*lines[i].wavelength, $
                     lines[i].fwhm_unc*lines[i].wavelength, $
                     lines[i].int_strength,lines[i].int_strength_unc, $
                     lines[i].eqw)] 
     endif else begin 
        info=[info, $
              string(FORMAT=form,lines[i].name, $
                     lines[i].wavelength, $
                     lines[i].central_inten, $
                     lines[i].fwhm*lines[i].wavelength, $
                     lines[i].int_strength, $
                     lines[i].eqw)]
     endelse 
  endfor
  
  df=decoded_params.dust_features
  info=[info,string(10b),'Dust Features: ']
  for i=0,n_elements(df)-1 do begin 
     name="DF_"+string(FORMAT='(F0.1)',df[i].wavelength)
     if unc then $
        info=[info,$
              string(FORMAT=form,name, $
                     df[i].wavelength,df[i].wavelength_unc, $
                     df[i].central_inten, df[i].central_inten_unc, $
                     df[i].fwhm*df[i].wavelength, $
                     df[i].fwhm_unc*df[i].wavelength, $
                     df[i].int_strength,df[i].int_strength_unc, $
                     df[i].eqw)] $
     else $
        info=[info, $
              string(FORMAT=form,name, $
                     df[i].wavelength, $
                     df[i].central_inten, $
                     df[i].fwhm*df[i].wavelength, $
                     df[i].int_strength, $
                     df[i].eqw)]
  endfor
  
  extinction=decoded_params.extinction
  info=[info,string(10b),'Extinction (mixed): ']
  form='(%"       Tau(9.7um): '+f+' Beta: %6.2g")'
  if unc then $
     info=[info,string(FORMAT=form,extinction.tau_9_7,extinction.tau_9_7_unc, $
                       extinction.beta)] $
  else info=[info,string(FORMAT=form,extinction.tau_9_7, extinction.beta)]
  
  if n_elements(filename) ne 0 then begin
     openw,un,filename,/GET_LUN
     printf,un,transpose(info)
     free_lun,un
  endif else print,transpose(info)
end


;=========================================================================
;  pahfit -- Fit a mixed dust and line model to a rest-frame IRS
;            low-resolution spectrum.  Pass an observed frame spectrum
;            along with redshift to ensure correct line widths.
;=========================================================================
function pahfit,lambda_0,intensity,errors, REPORT=report, $
                STARLIGHT_TEMPERATURE=star_temp, $
                CONTINUUM_TEMPERATURES=cont_temps, $
                LINES=lines, DUST_FEATURES=dust_features, $
                REDSHIFT=cz, PLOT_PROGRESS=plot_prog, $
                XSIZE=xs, YSIZE=ys, NO_EXTINCTION=no_extinct, $
                SCREEN_EXTINCTION=screen_extinction, $
                CHIAR_TIELENS=ct, KEMPER_VRIEND_TIELENS=kvt,$
                PARINFO=parinfo_in, NO_FIX=no_fix, VALUES_ONLY=vo, $
                NO_FIT=nf, NO_MEGAJANSKY_SR=no_mjy,_REF_EXTRA=e
  
  @pahfit_params
  
  if n_params() eq 0 then $
     message,'Usage: fit=pahfit(observed_wavelength, intensity, uncertainty)'
     
  
  if n_params() lt 3 then begin
     if keyword_set(nf) then begin 
        if n_elements(intensity) eq 0 then $ ; just for parameter estimation
           intensity=replicate(1.,n_elements(lambda_0))
     endif else $
        message,'Usage: fit=pahfit(observed_wavelength, intensity, uncertainty)'
  endif
  
  plot_progress=keyword_set(plot_prog)
  
  if plot_progress then begin 
     if n_elements(xs) eq 0 then xs=500
     if n_elements(ys) eq 0 then ys=400
     if !D.X_SIZE ne xs || !D.Y_SIZE ne ys || !D.WINDOW lt 0 then  $
        window,TITLE='PAHFIT: Fitter Progress',XSIZE=xs,YSIZE=ys,_EXTRA=e
     plot_win=!D.WINDOW
     window,/FREE,/PIXMAP,XSIZE=xs,YSIZE=ys,_EXTRA=e
     plot_dbwin=!D.WINDOW
     wset,plot_win
  endif 
  
  if n_elements(cz) eq 0 then begin 
     if ~keyword_set(nf) then $
        message,'No redshift specified, assuming 0',/CONTINUE
     cz=0.0
  endif 
  
  if cz ne 0.0 then $
     lambda=lambda_0/(1.+cz/299792.46) $
  else lambda=lambda_0
  
  no_errors=n_elements(errors) eq 0
  if no_errors then errors=replicate(1.D,n_elements(lambda))
  
  pf_lambda=lambda & pf_intensity=intensity & pf_errors=errors
    
  ;;--------- Parameter set defaults
  ;; Fixed single temperature for stellar spectrum
  if n_elements(star_temp) eq 0 then star_temp=5000.D
  
  ;; Fixed set of temperatures for modified dust blackbodies
  if n_elements(cont_temps) eq 0 then $
     cont_temps=[300.D,200.D,135.D,90.D,65.D,50.D,40.D,35.D]
  
  ;; Fixed set of lines (specified by central wavelengths only)
  if n_elements(lines) eq 0 then begin
     lines=[{WAVELENGTH: 5.5115D,  NAME: "H2 S(7)"}, $
;; Severely blended with PAH 6.2, and usually weak, enable by-hand?            
            {WAVELENGTH: 6.1088D,  NAME: "H2 S(6)"}, $
            {WAVELENGTH: 6.9091D,  NAME: "H2 S(5)"}, $
            {WAVELENGTH: 6.985274D,NAME: "[ArII]"},  $
            {WAVELENGTH: 8.0258D,  NAME: "H2 S(4)"}, $
            {WAVELENGTH: 8.99138D, NAME: "[ArIII]"}, $
            {WAVELENGTH: 9.6649D,  NAME: "H2 S(3)"}, $
            {WAVELENGTH: 10.5105D, NAME: "[SIV]"}, $
            {WAVELENGTH: 12.2785D, NAME: "H2 S(2)"}, $
            {WAVELENGTH: 12.813D,  NAME: "[NeII]"}, $
            {WAVELENGTH: 15.555D,  NAME: "[NeIII]"}, $
            {WAVELENGTH: 17.0346D, NAME: "H2 S(1)"}, $
            {WAVELENGTH: 18.713D,  NAME: "[SIII] 18"}, $
            {WAVELENGTH: 25.91D,   NAME: "[OIV]"}, $
            {WAVELENGTH: 25.989D,  NAME: "[FeII]"}, $
            {WAVELENGTH: 28.2207D, NAME: "H2 S(0)"}, $
            {WAVELENGTH: 33.480D,  NAME: "[SIII] 33"}, $
            {WAVELENGTH: 34.8152D, NAME: "[SiII]"}]   
;;            {WAVELENGTH: 35.349D,  NAME: "[FeII]"}]
  endif
  line_names=lines.NAME
  
  ;; Fixed set of dust features, as central wavelengths and fractional FWHM
  if n_elements(dust_features) eq 0 then begin
     dust_features=[{WAVELENGTH: 5.27D,  FRAC_FWHM: 0.034D}, $
                    {WAVELENGTH: 5.70D,  FRAC_FWHM: 0.035D}, $
                    {WAVELENGTH: 6.22D,  FRAC_FWHM: 0.030D}, $
                    {WAVELENGTH: 6.69D,  FRAC_FWHM: 0.07D},  $
                    {WAVELENGTH: 7.42D,  FRAC_FWHM: 0.126D}, $
                    {WAVELENGTH: 7.60D,  FRAC_FWHM: 0.044D}, $
                    {WAVELENGTH: 7.85D,  FRAC_FWHM: 0.053D}, $
                    {WAVELENGTH: 8.33D,  FRAC_FWHM: 0.05D},  $
                    {WAVELENGTH: 8.61D,  FRAC_FWHM: 0.039D}, $
                    {WAVELENGTH: 10.68D, FRAC_FWHM: 0.02D},  $
                    {WAVELENGTH: 11.23D, FRAC_FWHM: 0.012D}, $
                    {WAVELENGTH: 11.33D, FRAC_FWHM: 0.032D}, $
                    {WAVELENGTH: 11.99D, FRAC_FWHM: 0.045D}, $
                    {WAVELENGTH: 12.62D, FRAC_FWHM: 0.042D}, $
                    {WAVELENGTH: 12.69D, FRAC_FWHM: 0.013D}, $
                    {WAVELENGTH: 13.48D, FRAC_FWHM: 0.04D},  $
                    {WAVELENGTH: 14.04D, FRAC_FWHM: 0.016D}, $
                    {WAVELENGTH: 14.19D, FRAC_FWHM: 0.025D}, $
                    {WAVELENGTH: 15.9D,  FRAC_FWHM: 0.02D},  $
                    {WAVELENGTH: 16.45D, FRAC_FWHM: 0.014D}, $
                    {WAVELENGTH: 17.04D, FRAC_FWHM: 0.065D}, $
                    {WAVELENGTH: 17.375D,FRAC_FWHM: 0.012D}, $
                    {WAVELENGTH: 17.87D, FRAC_FWHM: 0.016D}, $
                    {WAVELENGTH: 18.92D, FRAC_FWHM: 0.019D},  $
                    {WAVELENGTH: 33.1D,  FRAC_FWHM: 0.05D}]
     ;;{WAVELENGTH: 18.88D, FRAC_FWHM: 0.04D}, $
     ;;{WAVELENGTH: 35.2D,  FRAC_FWHM: 0.04D} $
  endif 
  
  n_cont_temps=n_elements(cont_temps)
  n_lines=n_elements(lines)
  n_dust_features=n_elements(dust_features)
  
  
  ;;------------- Create PARINFO parameter array with defaults and limits
  med_inten=double(median(intensity))>0.0D
  
  par={VALUE:0.0D,FIXED:0b,LIMITED:[0b,0b],LIMITS:[0.0D,0.0D],PARNAME:'', $
       PAHFIT_TYPE:0b}
  
  ;;------------  Stellar continuum
  pars=replicate(par,2)
  pars.PAHFIT_TYPE=0b
  tau_guess=interpol(intensity,lambda,5.5,/LSQUADRATIC)/ $
            pahfit_blackbody_nu(5.5,star_temp)>0.0D
  pars[0].PARNAME='tau_star'
  pars[0].VALUE=tau_guess
  pars[0].LIMITED[0]=1b
  pars[0].LIMITS[0]=0.0D        ; depth greater than zero
  pars[1].PARNAME='T_star'
  pars[1].VALUE=star_temp    
  pars[1].FIXED=1b
  
  parinfo=pars
     
  ;;------------  Thermal dust continuum components
  pars=replicate(par,2*n_cont_temps)
  pars.PAHFIT_TYPE=1b
  names=string(FORMAT='(F0.1)',cont_temps)+'K'
  ind=lindgen(n_cont_temps)*2
  
  ;; Guess the tau from the flux at wavelength maximum
  mx=max(lambda,min=mn)
  max_lam=mn>2898./cont_temps<mx
  tau_guess=fltarr(n_cont_temps)
  for i=0,n_cont_temps-1 do begin 
     tau_guess[i]=interpol(intensity,lambda,max_lam[i],/LSQUADRATIC)/ $
                  pahfit_blackbody_nu(max_lam[i],cont_temps[i])* $
                  (max_lam[i]/9.7D)^2/5.>0.0D
  endfor 
  
  pars[ind].PARNAME='tau_dust['+names+']'
  pars[ind].VALUE=tau_guess     ;tau_dust
  pars[ind].LIMITED[0,*]=1b 
  pars[ind].LIMITS[0,*]=0.0D    ;non-negative taus
  ind++
  pars[ind].PARNAME='T_dust['+names+']'
  pars[ind].VALUE=cont_temps    ;fixed dust temperatures
  pars[ind].FIXED=1b
  
  parinfo=[parinfo,pars]
  
  ;;------------ Lines
  pars=replicate(par,3*n_lines)
  pars.PAHFIT_TYPE=2b
  ind=lindgen(n_lines)*3
  pars[ind].PARNAME='line_lambda['+lines.NAME+']'
  pars[ind].VALUE=lines.WAVELENGTH 
  pars[ind].LIMITED=1b
  pars[ind].LIMITS=transpose([[lines.WAVELENGTH-0.05D], $
                              [lines.WAVELENGTH+0.05D]])
  ind++
  pars[ind].PARNAME='line_central_inten['+lines.NAME+']'
  pars[ind].VALUE=med_inten/2.
  pars[ind].LIMITED[0,*]=1b
  pars[ind].LIMITS[0,*]=0.0D    ;non-negative central_inten
  ind++
  pars[ind].PARNAME='line_frac_fwhm['+lines.NAME+']'
  pars[ind].LIMITED=1b
  ;; Relative Line width limits, vary based on redshift/module
  obs_wav=lines.WAVELENGTH*(1.+cz/299792.46) ; observed frame wavelengths
  ;; Approximate observed frame breaks between module orders
  breaks=[0.,7.55,14.6,20.7,100000.]
  width=[.053,.1,.14,.34]       ;FWHM, in um
  for i=0,n_elements(breaks)-2 do begin 
     wh=where(obs_wav gt breaks[i] AND obs_wav le breaks[i+1],cnt)
     if cnt eq 0 then continue
     pars[ind[wh]].VALUE=width[i]/obs_wav[wh]
     pars[ind[wh]].LIMITS= $
        rebin(transpose(pars[ind[wh]].VALUE),2,cnt) * $
        rebin([.90,1.10],2,cnt,/SAMPLE)
  endfor
  parinfo=[parinfo,pars]
  
  ;;------------  Dust features
  pars=replicate(par,3*n_dust_features)
  pars.PAHFIT_TYPE=3b
  ind=lindgen(n_dust_features)*3
  names=string(FORMAT='(F0.1)',dust_features.WAVELENGTH)
  pars[ind].PARNAME='dust_feature_lambda['+names+']'
  pars[ind].VALUE=dust_features.WAVELENGTH 
  ;;pars[ind].LIMITED=1b
  ;;pars[ind].LIMITS[0,*]=transpose(pars[ind].VALUE)-.1
  ;;pars[ind].LIMITS[1,*]=transpose(pars[ind].VALUE)+.1
  pars[ind].FIXED=1b
  ind++
  pars[ind].PARNAME='dust_feature_central_inten['+names+']'
  pars[ind].VALUE=med_inten/2.
  pars[ind].LIMITED[0,*]=1b
  pars[ind].LIMITS[0,*]=0.0D    ;>0
  ind++
  pars[ind].PARNAME='dust_feature_frac_fwhm['+names+']'
  pars[ind].VALUE=dust_features.FRAC_FWHM
  pars[ind].FIXED=1b
  ;;pars[ind].LIMITED=1b
  ;;pars[ind].LIMITS[0,*]=transpose(pars[ind].VALUE)*.9
  ;;pars[ind].LIMITS[1,*]=transpose(pars[ind].VALUE)*1.1
  
  parinfo=[parinfo,pars]
  
  ;;------------  Mixed Dust Attenuation
  pars=replicate(par,4)
  pars.PAHFIT_TYPE=4b
  pars[0].PARNAME='tau_9.7'
  pars[0].LIMITED[0]=1b
  pars[0].LIMITS[0]=0.0D        ;tau>0
  if keyword_set(no_extinct) then begin 
     pars[0].FIXED=1b
     pars[0].VALUE=0.0D
  endif else pars[0].VALUE=keyword_set(screen_extinction)?0.1D:0.5D
  pars[1].PARNAME='Beta_Extinct'
  pars[1].VALUE=0.1D
  pars[1].FIXED=1b
  
  pars[2].PARNAME='Curve_Type'
  pars[2].VALUE=keyword_set(ct)?1:0
  pars[2].FIXED=1
  
  pars[3].PARNAME='Screen'
  pars[3].VALUE=keyword_set(screen_extinction) 
  pars[3].FIXED=1
  
  parinfo=[parinfo,pars]
  
  ;; Assign globals: saved tau(lambda) for the working wavelength
  ext_curve=pahfit_tau(lambda,pars[1].VALUE,CHIAR_TIELENS=keyword_set(ct), $
                       KEMPER_VRIEND_TIELENS=keyword_set(kvt))
  
  ;; For an input PARINFO, assign the parameters by matching PARNAME
  for i=0,n_elements(parinfo_in)-1 do begin 
     wh=where(parinfo.parname eq parinfo_in[i].parname,cnt)
     if cnt eq 0 then continue
     if keyword_set(vo) then begin 
        parinfo[wh[0]].value=parinfo_in[i].value
        for k=0,1 do begin 
           if parinfo[wh[0]].limited[k] then begin 
              old=parinfo[wh[0]].value
              if k eq 0 then parinfo[wh[0]].value>=parinfo[wh[0]].limits[k] $
              else parinfo[wh[0]].value<=parinfo[wh[0]].limits[k]
              if parinfo[wh[0]].value ne old then begin 
                 message,'Warning, parameter pinned at '+ $
                         (['lower','upper'])[k]+' limit: ' + $
                         parinfo[wh[0]].parname,/CONTINUE
              endif 
           endif 
        endfor 
     endif else begin 
        if keyword_set(no_fix) then fixed=parinfo[wh[0]].fixed
        st=parinfo[wh[0]]
        struct_assign,parinfo_in[i],st,/NOZERO
        parinfo[wh[0]]=st
        if keyword_set(no_fix) then parinfo[wh[0]].fixed=fixed
     endelse 
  endfor 
  
  ;; If requested, don't fit, just update the parameters list
  if keyword_set(nf) then begin 
     if arg_present(parinfo_in) then parinfo_in=parinfo ;for output
     return,-1
  endif 
  
  ;; Call the LM fitter
  fit=mpfitfun('pahfit_function',lambda,intensity,errors, $
               PARINFO=parinfo,STATUS=status, ERRMSG=errmsg,PERROR=perror, $
               ITERPROC='pahfit_iterproc',/QUIET,BESTNORM=bn,DOF=df, $
               COVAR=covar,YFIT=yfit,_EXTRA=e)
  
  if plot_progress then wdelete,plot_dbwin
  
  if status LE 0 then message, errmsg
  
  parinfo.VALUE=fit             ;returned fit parameters
  if arg_present(parinfo_in) then parinfo_in=parinfo ;for parameter output

  ;; Decode the parameter vector, computing strengths and EQW's, and
  ;; adding error estimates, reduced chi-square, and covariance.
  decoded=pahfit_decode(fit,parinfo,lines.NAME,PERROR=perror, COVAR=covar)  
  pahfit_add_extra_info,decoded,COVAR=covar, $
                        MEGAJANSKY_SR=~keyword_set(no_mjy), $
                        /STRENGTH_EQW,RED_CHI_SQ=bn/df, $
                        FINAL_FIT=yfit
  
  ;; Save decoded info to text file, if requested
  if keyword_set(report) then begin 
     if size(report,/TYPE) eq 7 then pahfit_report,decoded,report $
     else pahfit_report,decoded
  endif 
  
  return, decoded
end
