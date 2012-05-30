#!/usr/bin/env python
# ======================================================================

# Globally useful modules, imported here and then accessible by all
# functions in this file:

import matplotlib
# Force matplotlib to not use any Xwindows backend:
matplotlib.use('Agg')

# Fonts, latex:
matplotlib.rc('font',**{'family':'serif', 'serif':['TimesNewRoman']})
matplotlib.rc('text', usetex=True)

import numpy,atpy,pylab,sys,getopt,pyfits

# ======================================================================

def make_simulated_data(argv):

  USAGE = """
  NAME
    make_simulated_data.py

  PURPOSE
    Given a FITS table of 10 year AGN lightcurves (from Greg) and some lens 
    IDs, make simple simulated LSST data. 

  COMMENTS
    Use atpy to keep control of tables. Output is in BJB-compliant format.
    Sampling patterns are generated internally, later on may want to import 
    one from LSST opsim output.

  USAGE
    make_simulated_data.py [flags] [options]

  FLAGS
          -h --help    Print this message [0]
          -v --verbose Be verbose [0]

  INPUTS
          lctable      FITS table of lightcurves, in GD format
          ID1 ID2 etc  List of lens IDs to use.
          
  OPTIONAL INPUTS
          --catalog    Lens catalog [def=oguri_qso_LSST_fstar.fits]
          -m ms        Make all sources have the same magnitude
          
  OUTPUTS
          stdout       Useful information

  EXAMPLES

  BUGS

  HISTORY
    2011-09-09 started Marshall (UCSB)    
  """

  # --------------------------------------------------------------------

  try:
      opts, args = getopt.getopt(argv[1:], "hvm:c:", ["help","verbose","source-magnitude","catalog="])
  except getopt.GetoptError, err:
      # print help information and exit:
      print str(err) # will print something like "option -a not recognized"
      print USAGE
      sys.exit(2)

  catalog = 'oguri_qso_LSST_fstar.fits'
  ms = 0
  verbose = False
  for o,a in opts:
      if o == "-v":
          vb = True
      elif o in ("-h", "--help"):
          print USAGE
          sys.exit()
      elif o in ("-c", "--catalog"):
          catalog = a
      elif o in ("-m", "--source-magnitude"):
          ms = float(a)
      else:
          assert False, "unhandled option"

  print "args =",args
  if len(args) > 1:
    lcfile = args[0]
    lensID = numpy.array(args[1:],dtype=numpy.int)
  else:
    print USAGE
    sys.exit()

  print "Making simulated data for",len(lensID),"lenses",lensID
  print "from lightcurves in",lcfile

  # --------------------------------------------------------------------

  # Universal variables:
  campaign = 4 # years
  imnames = ['A','B','C','D']
  startdate = 100

  # Read in lens and image data in expected format:

  lc = atpy.Table(lcfile)

  if vb: print "Read in",len(lc),"lightcurve structures from ",lcfile

  # lc[i][0] = array of mags, mean = 20
  # lc[i][1] = array of epochs (in days), from 0 to 3650
  # lc[i][2] = tau
  # lc[i][3] = sigma
  
  # Loop over desired quads, pulling out time delays and image mags:

  hdulist = pyfits.open(catalog)
  d = hdulist[1].data
  hdulist.close()

  # >>> cols = hdulist[1].columns
  # >>> cols.names
  # ['LENSID', 'FLAGTYPE', 'NIMG', 'ZLENS', 'VELDISP', 'ELLIP', 'PHIE', 
  # 'GAMMA', 'PHIG', 'ZSRC', 'XSRC', 'YSRC', 'MAGI_IN', 'MAGI', 'IMSEP', 
  # 'XIMG', 'YIMG', 'MAG', 'DELAY', 'KAPPA', 'FSTAR', 'DD', 'DDLUM', 
  # 'ABMAG_I', 'APMAG_I', 'KCORR', 'DS', 'DDS', 'SIGCRIT', 'DSLUM', 'L_I', 
  # 'REFF', 'REFF_T']

  if vb: print "Read in",len(d),"mock lenses from ",catalog
  
  # Pull out lenses we want:
  
  for j in range(len(lensID)):

    if vb: 
      print "***************************** lens ID ",lensID[j],"*****************************"

    # Search through lens data to find lensID, and pull out row of data:
    
    index = numpy.where(d.field('lensid') == lensID[j])
    if len(index) == 0:
      print "ERROR: lens",lensID[j]," not found in",catalog
      return
    index = int(index[0])
    dd = d[index]

    # Now assign variables from this row:

    nim = dd.field('nimg')
    ID = dd.field('lensid')
    if vb: print "Lens number",ID," has",nim," images"

    if ms == 0:
      ms  = dd.field('magi_in')
    # Otherwise ms was given on command line...
    mui = dd.field('mag')
    ti  = dd.field('delay')

    if vb: print "Re-ordering by arrival time..."
    index = range(nim)
    index = numpy.argsort(ti)
    ti = ti[index]
    mui = mui[index]
    for i in range(nim):
      # NB. TTake integer part of time delay, to avoid interpolation
      # when we shift the light curvess later:
      ti[i] = int(ti[i])
      if vb: print "Image",i,": mu,t =",mui[i],ti[i]

    mi = numpy.zeros(nim)     
    for i in range(nim):
      mi[i] = ms - 2.5*numpy.log10(numpy.abs(mui[i]))
    print "Mean image magnitudes:",mi

    # Now compute 4-image lightcurves for each (tau,sigma) pair:
    for kt in range(3):
      for ks in range(3):
        k = ks*3 + kt
        tau = lc[k][2]
        taustring = "%s" % int(tau)
        sigma = lc[k][3]
        sigmastring = "%s" % int(sigma*1000)
        
        print "AGN properties: tau = ",tau,", sigma = ",sigma
        
        # Loop over observing strategies - note truth file is (1,12):
        for cadence in (1,4,8):
          cadencestring = str(cadence)
          for season in (12,4,8):
            seasonstring = str(season)

            # Write filename and open file:
            outfile = 'oguri_qso_ID='+str(ID)+'_tau='+taustring+'_sigma='+sigmastring
            outfile += '_cadence='+cadencestring+'_season='+seasonstring
            
            GO = False
            if (cadence > 1 and season < 12):
              # Noisy mock data with realistic sampling:
              outfile += '_lightcurve.txt'
              nsigma = 5
              sigma_f = 0.2
              GO = True
            elif (cadence == 1 and season == 12):  
              # Effectively noiseless reference curves:
              outfile += '_truthcurve.txt'
              nsigma = 10000.0
              sigma_f = 0.0
              GO = True
            
            if GO:
               print "Output filename: "+outfile

               # Read in lightcurve:
               t0 = lc[k][1]
               m0 = lc[k][0]

               # Define sampling pattern an array of times:
               t_obs = sampling_pattern(startdate,cadence,season,campaign)

               # Now make shifted lightcurves, one for each image.
               # Subtract 20, and add mean mags for each image - and offset 
               # in time as well. Note the integer time delays instituted above.
               # Then sample according to the pattern:

               x = numpy.zeros([nim,4,len(t_obs)])

               for i in range(nim):

                 t = t0 + ti[i]
                 m = m0 - 20 + mi[i]

                 m_true = sample_from_lightcurve(t_obs,t,m)

                 m_obs,m_err = addnoise(m_true,nsigma=nsigma)

                 x[i,0,:] = t_obs
                 x[i,1,:] = m_obs
                 x[i,2,:] = m_err
                 x[i,3,:] = i

               print "Calculated",len(t_obs),"-point lightcurve"
               print "Writing to file..."

               write_bjb_lcfile(outfile,x,mi,ti,imnames,sigma,tau,sigma_f)

               print "...done"
            
           # End of loop over seasons.
        # End of loop over cadence.
      # End of loop over tau values.
    # End of loop over sigma values.       
  
    print "On to the next lens!"
    
  # End of loop over lens IDs   
    
  print "All done."
        
  return
    
# ======================================================================
# Generate a vector of noisy magnitudes, aand estimates of the uncertainties
# on them...

def addnoise(truem,mlimit=24,nsigma=5):

  # Convert magnitude to flux, in units such that a mlimit mag source
  # has flux f0 = 1:
  truef = 10.0**(0.4*(mlimit - truem))

  # Assume mlimit is n-sigma, such that the noise rms is f0/nsigma.
  # Generate Gaussian noise vector:
  sigma_f = 1.0/float(nsigma)
  noise = numpy.random.normal(0.0, sigma_f, len(truef))

  # Add noise to fluxes:
  f = truef + noise

  # Back to mags:
  m = mlimit - 2.5*numpy.log10(f)

  # Now: compute the rms uncertainty on m: bright mags are less uncertain!
  # Approximation from quadrature formula:
  m_err = 0.2 * (2.5/numpy.log(10)) / f

  # Note that we use f, not truef here: because we don't know truef!
  # If f fluctuates low, it gets a large m_err: good. But if f fluctuates 
  # high, we upweight that magnitude point. This is the residual skewness,
  # that is not captured by the quadrature approximation.

  return m,m_err

# ======================================================================
# Return a sampling pattern - a numpy array of epochs.

def sampling_pattern(startdate,cadence,season,campaign): 
    
   # Start the pattern:
   pattern = numpy.array([])
   
   # How many observations per season?
   nobs = int(season*30/cadence)
   
   # Populate with observations:
   for year in range(campaign):
     for obs in range(nobs):
       epoch = startdate + year*365 + obs*cadence
       pattern = numpy.append(pattern,epoch)
       
   return pattern

# ======================================================================
# Sample brightnesses from a given lightcurve (t,m) lightcurve.

def sample_from_lightcurve(t_obs,t,m):

  #  Simple linear interpolation:
  m_obs = numpy.interp(t_obs,t,m)
    
  return m_obs
  
# ======================================================================
# Write a BJB-compliant lightcurve file.

def write_bjb_lcfile(filename,x,mag,dt,imnames,sigma,tau,noise):
    
  FILE = open(filename,"w")

  # First write header:
  FILE.write("# Simulated lensed AGN lightcurve\n")
  FILE.write("#   (Dobler & Marshall 2011)\n")
  FILE.write("# \n")
  FILE.write("# True image parameters:\n")
  FILE.write("#   name   dt/days   mean mag\n")
  for i in range(len(imnames)):
    FILE.write("#     %s     %.1f      %f\n" % (imnames[i],dt[i],mag[i]))
  FILE.write("# True source parameters:\n")
  FILE.write("#     sigma = %f\n" % sigma)
  FILE.write("#       tau = %.3f\n" % tau)
  FILE.write("# True noise parameters:\n")
  FILE.write("#   sigma_f = %.1f\n" % noise)
  FILE.write("# \n")
  FILE.write("# Data columns:\n")
  FILE.write("#   t/days        mag   mag_err   image\n")
  
  # Now write data, line by line:
  for i in range(len(imnames)):
    for j in range(len(x[0,0,:])):
      t_obs = x[i,0,j]
      m_obs = x[i,1,j]
      m_err = x[i,2,j]
      im    = x[i,3,j]
      FILE.write("     %.1f  %f  %f     %d\n" % (t_obs,m_obs,m_err,im))
        
  FILE.close()

  return
  
# ======================================================================
# If called as a script, the python variable __name__ will be set to 
# "__main__" - so test for this, and execute the main program if so.
# Writing it like this allows the function plot_oguri_lens to be called
# from the python command line as well as from the unix prompt.

if __name__ == '__main__':
  make_simulated_data(sys.argv)

# ======================================================================


