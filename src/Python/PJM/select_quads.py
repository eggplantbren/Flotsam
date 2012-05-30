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

import numpy,pylab,sys,getopt,pyfits

# ======================================================================

def select_quads(argv):

  USAGE = """
  NAME
    select_quads.py

  PURPOSE
    Select 3 quads from LSST catalog, and plot.

  COMMENTS
    FITS table from Dobler has lens galaxy sizes in it. Best to use pyfits 
    utilities.

  USAGE
    select_quads.py [flags] [options]

  FLAGS
          -h --help    Print this message [0]
          -v --verbose Be verbose [0]

  INPUTS

  OPTIONAL INPUTS
          --catalog    Lens catalog [def=oguri_qso_LSST_fstar.fits]
          
  OUTPUTS
          stdout       Useful information

  EXAMPLES

  BUGS

  HISTORY
    2011-03-05 started Marshall (UCSB)    
  """

  # --------------------------------------------------------------------

  try:
      opts, args = getopt.getopt(argv[1:], "hvl:c:", ["help","catalog=","logfile="])
  except getopt.GetoptError, err:
      # print help information and exit:
      print str(err) # will print something like "option -a not recognized"
      print USAGE
      sys.exit(2)

  catalog = 'oguri_qso_LSST_fstar.fits'
  verbose = False
  for o,a in opts:
      if o == "-v":
          vb = True
      elif o in ("-h", "--help"):
          print USAGE
          sys.exit()
      elif o in ("-c", "--catalog"):
          catalog = a
      else:
          assert False, "unhandled option"

  # --------------------------------------------------------------------

  # Read in lens and image data in expected format:

  hdulist = pyfits.open(catalog)
  d = hdulist[1].data
  hdulist.close()

  if vb: print "Read in",len(d),"rows of lens and image data from ",catalog

  # >>> cols = hdulist[1].columns
  # >>> cols.names
  # ['LENSID', 'FLAGTYPE', 'NIMG', 'ZLENS', 'VELDISP', 'ELLIP', 'PHIE', 
  # 'GAMMA', 'PHIG', 'ZSRC', 'XSRC', 'YSRC', 'MAGI_IN', 'MAGI', 'IMSEP', 
  # 'XIMG', 'YIMG', 'MAG', 'DELAY', 'KAPPA', 'FSTAR', 'DD', 'DDLUM', 
  # 'ABMAG_I', 'APMAG_I', 'KCORR', 'DS', 'DDS', 'SIGCRIT', 'DSLUM', 'L_I', 
  # 'REFF', 'REFF_T']

  # Pull out quads:
  
  nimg  = d.field('nimg')
  index = numpy.where(nimg == 4)

  lensid = d.field('lensid')[index]
  zd     = d.field('zlens')[index]
  sigma  = d.field('veldisp')[index]
  zs     = d.field('zsrc')[index]
  fstar  = d.field('fstar')[index]
  if vb: print "Found",len(lensid),"quads"

  zdlow,zdhigh = 0.4,0.6
  zslow,zshigh = 1.6,2.4
  sigmalow,sigmahigh = 150,300
  
  # Now select within these ranges:
  index = numpy.where((zd > zdlow) * (zd < zdhigh) * \
                      (zs > zslow) * (zs < zshigh) * \
                      (sigma > sigmalow) * (sigma < sigmahigh))
  if vb: print len(index[0]),"are in 1-sigma volume:"
  
  print lensid[index]
                      
  return
                      
# ======================================================================

# If called as a script, the python variable __name__ will be set to 
# "__main__" - so test for this, and execute the main program if so.
# Writing it like this allows the function plot_oguri_lens to be called
# from the python command line as well as from the unix prompt.

if __name__ == '__main__':
  select_quads(sys.argv)

# ======================================================================


