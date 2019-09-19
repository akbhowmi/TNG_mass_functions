import numpy
import arepo_package

Nbins = 15
log_mass_min = 6
log_mass_max = 12

basePath='/n/ghernquist/Illustris/Runs/Illustris-1/'

run='L205n2500TNG'

basePath='/n/hernquistfs3/IllustrisTNG/Runs/'+run+'/output/'

for desired_redshift,col in zip([0.,3., 5.],['blue','red','green']):
    object_type='subhalo'
#    basePath = '/ufrc/lblecha/aklantbhowmick/arepo_runs_aklant/L25_n128/output/'
    category='bh'
    centers,HMF,dHMF,output_redshift=arepo_package.get_mass_function(category,object_type,desired_redshift,basePath,Nbins,log_mass_min,log_mass_max,list_all=False)
    numpy.save('./'+run+'_'+object_type+'_'+category+'_z%.2f.npy'%desired_redshift,[centers,HMF,dHMF,output_redshift])
