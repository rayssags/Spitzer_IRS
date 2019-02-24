def pah_plot_ptbr(word):
	fit_pahfit_ptbr = {"lam_mod":"comprimento de onda",
                   "cont_mod":"contínuo",
                   "feat_mod":"PAHs",
                   "lines_mod":"linhas atômicas e moleculares",
                   "stars_mod":"contínuo estelar",
                   "dust_mod":"contínuo da poeira",
                   "bestfit_mod":"modelo",
                   "ext_mod":"extinção"}
	return fit_pahfit_ptbr[word]

def clean_and_PAHfit(path='/home/rayssa/Win/Programs/Spitzer/1_Raw/*.tbl', z=0.00502):
    import glob
    import pandas as pd
    files = glob.glob(path)
    redshift_cf = (1+z)
    a = path.find('*')
    b = len(path)
    c = a - b
    e = path[:c]
    for filename in files:
        c = pd.read_csv(filename, comment='\\', quotechar="|", sep='\s+', skip_blank_lines=True,
                        float_precision='high')
        c = c.drop(c.index[[0,1]])
        col = c.keys()
        col_new = []
           
        for i in range(len(col)): col_new.append(list(col)[i].replace(" ", ""))
        c.columns = col_new
                
        filename_clean = "/home/rayssa/Win/Programs/Spitzer/2_Clean/"+ filename.replace(e,"")[:-8] +"_clean.tbl"
        c.to_csv(filename_clean,mode = 'w', header=1, index=0, sep='|')

        c = c.apply(pd.to_numeric)     
        c['WAVELENGTH']=c['WAVELENGTH']/redshift_cf
        c['FLUX']=c['FLUX']*redshift_cf
        
        filename_pahfit = "/home/rayssa/Win/Programs/Spitzer/3_PAHfit_input/"+ filename.replace(e,"")[:-8]+"_pahfit_input.tbl"
        c.to_csv(filename_pahfit,mode = 'w', header=1, index=0, sep='|')

        
    return print('As tabelas em {} foram limpas e corrigidas'.format(path))

def PAHdb_to_pandas(path='/home/rayssa/Win/Programs/Spitzer/5_PAHdb_output/*.tbl'):
    import glob
    import pandas as pd
    files = glob.glob(path)
    for pahdb_output in files:
        f= open(pahdb_output,"r")
        lines = f.readlines()
        for i in [0,2,3]: lines[i] = '#'+lines[i]
        f= open(pahdb_output,"w")
        f.writelines(lines)

def PAHdb_to_plot(path='/home/rayssa/Win/Programs/Spitzer/5_PAHdb_output/*.tbl'):
    import glob
    import pandas as pd
    files = glob.glob(path)

    a = path.find('*')
    b = len(path)
    c = a - b
    e = path[:c]
    
    for pahdb_output in files:
        df = pd.read_csv(pahdb_output, delim_whitespace=True, skip_blank_lines=True, comment='#', float_precision='high')
        df.rename(columns = {df.keys()[0]:'Wavenumber'}, inplace = True)
        df['Wavenumber'] = df['Wavenumber'].apply(pd.to_numeric)

        filename_pahdb = "/home/rayssa/Win/Programs/Spitzer/6_PAHdb_final/"+ pahdb_output.replace(e,"")[:-16] +"_pahdb.tbl"

        pahdb = pd.DataFrame()
        pahdb['Wavelength'] = 1e4/df['Wavenumber']
        pahdb['Flux'] = df.iloc[:, 1:].sum(axis=1)
        pahdb.to_csv(filename_pahdb, mode = 'w', header=1, index=0, sep='|')



