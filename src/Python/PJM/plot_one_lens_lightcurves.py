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

import numpy,pylab,sys,getopt,string

# ======================================================================

def plot_one_lens_lightcurves(argv):

  USAGE = """
  NAME
    plot_one_lens_lightcurves.py

  PURPOSE
    Given a BJB lighcurve file, plot all 4 image lightcurves with sensible
    offsets.

  COMMENTS
    Optionally include predicted curves? True curves? Can plot many files...

  USAGE
    plot_one_lens_lightcurves.py [flags] [options] lcfile

  FLAGS
          -h --help    Print this message [0]
          -v --verbose Be verbose [0]
          --eps        Postscript output [0, png is default]

  INPUTS
          lcfile       BJB format lightcurve file (or list of files)

  OPTIONAL INPUTS
          -t --truth   Overlay true lightcurve in grey (on each plot)
          
  OUTPUTS
          stdout       Useful information

  EXAMPLES


  BUGS

  HISTORY
    2011-09-12 started Marshall (KIPAC)
  """

  # --------------------------------------------------------------------

  try:
      opts, args = getopt.getopt(argv[1:], "hvt:", ["help","eps","truth="])
  except getopt.GetoptError, err:
      # print help information and exit:
      print str(err) # will print something like "option -a not recognized"
      print USAGE
      sys.exit(2)

  vb = False
  eps = False
  truefile = 0
  lcfiles = []
  for o,a in opts:
      if o == "-v":
          vb = True
      elif o in ("-h", "--help"):
          print USAGE
          sys.exit()
      elif o in ("--eps"):
          eps = True
      elif o in ("-t","--truth"):
          truefile = a
          lcfiles.append(truefile)
      else:
          assert False, "unhandled option"

  if len(args) > 0:
    nlens = len(args)
    lcfiles += args
  else:
    print USAGE
    sys.exit()

  print "Making plots of",nlens,"lens systems"
  print lcfiles

  imnames = ['A','B','C','D']
  colors = ['blue','green','orange','red']
  # --------------------------------------------------------------------

  # Loop over lightcurve files:

  for k in range(len(lcfiles)):

    lcfile = lcfiles[k]

    if vb: 
      print "***************************** "+lcfile+" *****************************"
    if truefile != 0 and k ==0: 
      truth = True
      if vb: print "This is the true lightcurve"
    else:
      truth = False

    # Read in lightcurve data in expected format:

    x = numpy.loadtxt(lcfile)

    t = x.T[0]
    m = x.T[1]
    merr = x.T[2]
    im = x.T[3]
    im = im.astype(int)
    
    nim = numpy.max(im)+1
    nepochs = len(t)/nim
    print "Read in",nepochs,"epochs of data, for",nim,"images"

    tt = numpy.zeros([nim,nepochs])
    mm = numpy.zeros([nim,nepochs])
    dm = numpy.zeros([nim,nepochs])
    
    for i in range(nim):
      index = numpy.where(im == i)
      tt[i,:] = t[index]
      mm[i,:] = m[index]
      dm[i,:] = merr[index]
      
    # Figure out offsets for display - measure mean mag and stdev for each 
    # image, and offset accordingly.
    
    mmean = numpy.zeros(nim)    
    mstdev = numpy.zeros(nim)
    for i in range(nim):
      mmean[i] = numpy.average(mm[i,:])
      mstdev[i] = numpy.std(mm[i,:])
    originalmmean = mmean.copy()
       
    # Compute offsets for the first file, and apply to subsequent files:
    if k == 0:
      index = numpy.flipud(numpy.argsort(mmean))
      offset = numpy.zeros([nim])
      for i in range(nim):
        ii = index[i]
        if i == 0: 
          offset[ii] = 0.0
        else:
          jj = index[i-1]
          diff = (mmean[jj]-mstdev[jj]) - (mmean[ii]+mstdev[ii])
          offset[ii] = int(10*diff)/10.0 - 0.5
        mmean[ii] += offset[ii]
        mm[ii,:] += offset[ii]
        # print "updated",ii,"th mean to ",mmean[ii]
        # print "mmean is now",mmean
    else:
      for i in range(nim):
        mmean[i] += offset[i]
        mm[i,:] += offset[i]
    # Reset index to get new means!
    index = numpy.flipud(numpy.argsort(mmean))

    # Fix labels:
    if k==0:
      for i in range(nim):
        if offset[i] < 0.0:
          plus = '$-$'
        else:
          plus = '$+$'
        imnames[i] = imnames[i]+plus+str(numpy.abs(offset[i]))+' mag'
        if vb: print "Mag offset applied: ",imnames[i]
    
    # Report on this file's means:
    if vb: 
      for i in range(nim):
        print originalmmean[i],'updated to',mmean[i]
    
    if truth:
      tt0 = tt.copy()  
      mm0 = mm.copy()  
      dm0 = dm.copy()  

    # Now set overall plot limits - once:
    if k==0:
      tmin = 0.0
      tmax = 4*365 + 200
      mmax = numpy.max(mm[index[0],:]) + 0.5
      mmin = numpy.min(mm[index[nim-1],:]) - 0.5

    # Skip plotting if we are looking at the truth file:
    if not truth:

      # Start figure:
      fig = pylab.figure(figsize=(8,6))

      # Plot lightcurves, and label each one:
      for i in range(nim):
      
        pylab.plot(tt[i,:], mm[i,:], color=colors[i], \
                 marker='o', markersize=3, markeredgewidth=0, \
                 linestyle='')
      
        if truefile != 0: 
          pylab.plot(tt0[i,:], mm0[i,:], color='gray', alpha=0.5, \
                 marker='', linestyle='-')
        
        index = numpy.where(tt[i,:] < 400) 
        mpos = numpy.average(mm[i,index]) - numpy.std(mm[i,index])
        pylab.annotate('\\bf '+imnames[i],(tmin+40,mpos-0.1), \
                 color=colors[i], fontsize=12, fontweight='bold')

      # Set axes labels:
      pylab.xlabel("t / days",fontsize=20)
      pylab.ylabel("AB magnitude",fontsize=20)
      # Set axis limits:
      pylab.axis([tmin,tmax,mmax,mmin])
      # Add a grid:
      pylab.grid(color='grey', linestyle='--', linewidth=0.5)

      # Plot graph to file:
      pieces = string.split(lcfile,'.')
      stem = string.join(pieces[:-1],'.')
      if eps:
        epsfile = stem+".eps"
        pylab.savefig(epsfile)
        print "Figure saved to file:",epsfile
      else:
        pngfile = stem+".png"
        pylab.savefig(pngfile)
        print "Figure saved to file:",pngfile

  return

# ======================================================================

# If called as a script, the python variable __name__ will be set to 
# "__main__" - so test for this, and execute the main program if so.
# Writing it like this allows the function plot_one_lens_lightcurves to be called
# from the python command line as well as from the unix prompt.

if __name__ == '__main__':
  plot_one_lens_lightcurves(sys.argv)

# ======================================================================


