#! /usr/bin/env python
# -*- coding: utf-8 -*-
import shutil, os,  subprocess, fileinput, string, sys, stat, json

server = os.uname()[1]  # hostname

# path to ICES config
#pathToICES = '../../../config/ices.json'  # this should be the ONLY hard-coded path
#pathToICES = '../../../../../../www-v4.1/htdocs/ICES/config/ices.json'  # this should be the ONLY hard-coded path
pathToICES ='../../config/ices.json'

iceyfp = open(pathToICES)
assert (iceyfp), "Couldn't open "+pathToICES
cfg = json.load(iceyfp)  # read in the configuration object from the JSON file
iceyfp.close()

modelFP = open(cfg['SERVERS'][server]['DOCROOT']+"/"+cfg['ICESROOT']+"/config/"+cfg['TYPES']['Coma']['CONFIG'])
Models = json.load(modelFP)
modelFP.close()

# Gobal parameters

# path to not installed cases
#pathToNotInstalledCases = '/Volumes/HD-2TB/NetBeansProjects/ices-2012/www/htdocs/Models/Hybrid2/CasestoInstall'
pathToNotInstalledCases = cfg['SERVERS'][server]['DOCROOT']+"/"+cfg['MODELS']+"/"+Models['Hybrid2']['CasesToInstall']

# path to installed cases
#pathToInstalledCases = '/Volumes/HD-2TB/NetBeansProjects/ices-2012/www/htdocs/Data/Coma/Hybrid2'
pathToInstalledCases = cfg['SERVERS'][server]['DOCROOT']+"/"+cfg['DATA']+"/"+cfg['TYPES']['Coma']['DATA']+"/Hybrid2"

# path to folder with file 'aikef.py'
#pathToProgram = '/Volumes/HD-2TB/NetBeansProjects/ices-2012/www/htdocs/Models/Hybrid2'
pathToProgram = cfg['SERVERS'][server]['DOCROOT']+"/"+cfg['MODELS']+"/Hybrid2"

# name of orbit file
#orbitFile = 'orbit-input.txt'
orbitFile = Models['Hybrid2']['Trajectory']

# name of output file
#orbitFileOut = 'orbit-output.txt'
orbitFileOut = Models['Hybrid2']['Output']

# path to OpenMPI
# MacPorts
#pathToOpenMPI = '/opt/local/lib/openmpi'
pathToOpenMPI = Models['Hybrid2']['pathToOpenMPI']
# Apple MacOSX (default): 
# pathToOpenMPI = '/usr/local/openmpi-1.2.3_64'

# internal variables
doInstall = False
doRun = False

print 'Running ... '

# check if arguments are available
if ( len(sys.argv) > 2 and len(sys.argv)<6 ):   # 3 to 5 is acceptable
  # check for second argument
  if(sys.argv[2] == 'install' ):
    doInstall = True 
    print ' - set install'
  elif(sys.argv[2] == 'run' ):
    doRun = True
    print ' - set run'
  else :
    print 'Error: Unkown option'
    print 'first argument: case'
    print 'second argument: install or run'
    print 'third argument (on run only): path to runtime directory'  # *JK* - i.e. where the output orbit file will be created and written
    print 'fourth argument (on run only): orbit trajectory input file'  # *JK* added 6/14/2012
    sys.exit()
  # check for first argument
  case = sys.argv[1]
  
  #############################################
  # in case of install
  if(doInstall==True):
    # check for Makefile
    if os.path.exists('Makefile') :
      print ' - found Makefile'
    else:
      print 'Error: Could not find Makefile'
      sys.exit()
    if os.path.exists(os.path.join(pathToNotInstalledCases,case+'.tar')) :
      print ' - found case file'
    else :
      print 'Error: Could not find Case File'
      print os.path.join(pathToNotInstalledCases,case+'.tar')
      sys.exit()
  #############################################
  # in case of run
  if(doRun==True):
    if os.path.exists(os.path.join(pathToInstalledCases,case)) :
      print ' - found case folder'
      # For use on the web site, we need to make each runtime independant of the others.  argv[3] will contain that path
      RuntimePath = os.path.join(pathToInstalledCases,case,'bin')  ## *JK* the default place where the runtime path will be
      if( len(sys.argv)>=4 ):
          if( os.path.exists(sys.argv[3]) ):
               RuntimePath = sys.argv[3]
               print ' - Using alternate runtime path: '+RuntimePath
          else:
               print ' - WARNING: Alternate runtime path ('+sys.argv[3]+') not found; using standalone path: '+RuntimePath
      else:
          print " - Using runtime path: "+RuntimePath
      if( len(sys.argv)>=5 ):
      	  if( os.path.exists(os.path.join(sys.argv[3],sys.argv[4])) ):
      	  	  TrajectoryFile = os.path.join(sys.argv[3],sys.argv[4])
      	  else:
      	  	  if( os.path.exists(os.path.join(sys.argv[3],'orbit-input.txt')) ):
	      	  	  TrajectoryFile = os.path.join(sys.argv[3],'orbit-input.txt')
	      	  else:
	      	  	  print 'Error: Could not find input Trajectory file'
	      	  	  sys.exit()
    else :
      print 'Error: Could not find case folder'
      print 'please install the case first'
      sys.exit()
else:
  print 'Error: Number of arguments is wrong'
  print 'first argument: case'
  print 'second argument: install or run'
  sys.exit()


#######################################################
# Install the Case
if (doInstall==True):
  print ' - install case'
  if os.path.exists(os.path.join(pathToInstalledCases,case)) :
    print ' - remove old files'
    shutil.rmtree(os.path.join(pathToInstalledCases,case))
  print ' - create case directory'
  os.mkdir(os.path.join(pathToInstalledCases,case))
  print ' - unpack case file '+case+'.tar'
  subprocess.call(['tar','xf',os.path.join(pathToNotInstalledCases,case+'.tar'),'-C',os.path.join(pathToInstalledCases,case)])
  print ' - copy makefile'
  shutil.copy('Makefile',os.path.join(os.path.join(pathToInstalledCases,case,'Makefile')))
  print ' - create folders'
  os.mkdir(os.path.join(pathToInstalledCases,case,'bin'))
  os.mkdir(os.path.join(pathToInstalledCases,case,'data'))
  os.mkdir(os.path.join(pathToInstalledCases,case,'obj'))
  print ' - move state files'
  shutil.move(os.path.join(pathToInstalledCases,case,'State'),os.path.join(pathToInstalledCases,case,'bin','State'))
  print ' - run make: please wait ....'
  os.chdir(os.path.join(pathToInstalledCases,case))
  subprocess.call('make')
  if os.path.exists(os.path.join(pathToInstalledCases,case,'bin','aikef_mpi')) :
    print ' - installation of '+ case +' was successful'
  else:
    print 'Error: Can not find aikef_mpi'
    sys.exit()
    
#######################################################
# Run the Case
if (doRun==True):
  print ' - run case'
  print ' - link to State file directory'
  if os.path.exists(os.path.join(pathToInstalledCases,case,'bin','State')) :
     if( RuntimePath != os.path.join(pathToInstalledCases,case,'bin') ) :
          os.symlink(os.path.join(pathToInstalledCases,case,'bin','State'),os.path.join(RuntimePath,'State'))
     # Else, leave things as they are, and be sure NOT to unlink
  else :
     print 'Error: State directory does not exist: '+os.path.join(pathToInstalledCases,case,bin,'State')
     sys.exit()
     
  if os.path.exists(os.path.join(RuntimePath,'orbit.txt')) :
    print ' - remove old orbit file'
    os.remove(os.path.join(RuntimePath,'orbit.txt'))  # *JK* Using RuntimePath for all instead of hard-coded path
# *JK* We no longer have a fixed orbit file
#  print ' - search for new orbit file'
#   if os.path.exists(os.path.join(pathToProgram,orbitFile)) :
#     print ' - found orbit file'
#   else : 
#     print 'Error: Could not found orbit file'
#     sys.exit()
  print ' - prepare orbit file: '+orbitFile+' Trajectory file: '+TrajectoryFile
  try:
    # Open new orbit file
    file = open(os.path.join(RuntimePath,'orbit.txt'), 'w')   # *JK* Using RuntimePath for all instead of hard-coded path
    searchForStart = False
    # Read orbit-raw file line by line
    for line in fileinput.input(TrajectoryFile):  ## copy from the trajectory file specfied in arg4
          # Check if '#START' was found
          if searchForStart==False :
	          if string.find(line, '#START') != -1:
	               searchForStart = True
          else :
	          # '#Start' was already found
	          # put new lines into orbit file for A.I.K.E.F.
	     file.writelines(line)
    
    if searchForStart==False : #No '#START' in file copy file
          shutil.copy(os.path.join(pathToProgram,orbitFile),os.path.join(RuntimePath,'orbit.txt'))  ## RuntimePath will reflect the proper directory, whether default or overridden
          
    # Close File
    file.close()
  except IOError:
    print 'Error: preparing orbit file failed!'
    sys.exit()
  print ' - run A.I.K.E.F.: Please wait ...' 
  # If we use -wdir in mpirun, we don't have to chdir(); RuntimePath is where we'll be running from; current directory is preserved after the run
  # print 'Running: '+os.path.join(pathToOpenMPI,'bin','mpirun')+' '+os.path.join(pathToInstalledCases,case,'bin','aikef_mpi')+' '+'-np 4'+' '+'-wdir '+RuntimePath
  subprocess.call([os.path.join(pathToOpenMPI,'bin','mpirun'),os.path.join(pathToInstalledCases,case,'bin','aikef_mpi'),'-np 4','-wdir '+RuntimePath])
  if( len(sys.argv)<4 ):   # *JK* if we've not specified an alternate runtime directory, then move the results to the program directory
     print ' - move orbit file'
     shutil.move(os.path.join(pathToInstalledCases,case,'bin','orbit.txt'),os.path.join(pathToProgram,orbitFileOut))
     print ' - Results are found in: '+os.path.join(pathToProgram,orbitFileOut)
  else :
     os.rename(os.path.join(RuntimePath,'orbit.txt'),os.path.join(RuntimePath,orbitFileOut))
     print ' - Results are found in: '+os.path.join(RuntimePath,orbitFileOut)
     ## Only zap the link if its actually a link...NO DIRS!!!
     if( stat.S_ISLNK(os.lstat(os.path.join(RuntimePath,'State')).st_mode) != 0 ) :  ## must use lstat() to detect symlinks; stat() follows symlinks
          os.remove(os.path.join(RuntimePath,'State'))

print ' - finish.'
sys.exit()
