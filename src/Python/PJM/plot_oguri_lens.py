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

def plot_oguri_lens(argv):

  USAGE = """
  NAME
    plot_oguri_lens.py

  PURPOSE
    Given a lensID, pull out the information on that lens from the 
    Oguri catalog and compute some basic quantities, 
    then plot them on the sky.

  COMMENTS
    FITS table from Dobler has lens galaxy sizes in it. Best to use pyfits 
    utilities.

  USAGE
    plot_oguri_lens.py [flags] [options] lensID

  FLAGS
          -h --help    Print this message [0]
          -v --verbose Be verbose [0]
          --eps        Postscript output [0, png is default]
          --showdisk   Plot a seeing disk for comparison [0]

  INPUTS
          lensID       Lens ID number

  OPTIONAL INPUTS
          --catalog    Lens catalog [def=oguri_qso_LSST_fstar.fits]
          
  OUTPUTS
          stdout       Useful information

  EXAMPLES

    First double in LSST catalog:
          plot_oguri_lens.py 18060 

    First quad in LSST catalog:
          plot_oguri_lens.py 4376

  BUGS
    Image opening angles greater than 180deg are not possible - need to 
    replace the largest with 360-sum(phi1:phi3)

  HISTORY
    2010-06-13 started Marshall (KIPAC)
    2011-03-05 adaptd for FITS catalogs Marshall (UCSB)    
  """

  # --------------------------------------------------------------------

  try:
      opts, args = getopt.getopt(argv[1:], "hvl:c:", ["help","eps","catalog=","logfile="])
  except getopt.GetoptError, err:
      # print help information and exit:
      print str(err) # will print something like "option -a not recognized"
      print USAGE
      sys.exit(2)

  catalog = 'oguri_qso_LSST_fstar.fits'
  vb = False
  showdisk = False
  eps = False
  for o,a in opts:
      if o == "-v":
          vb = True
      elif o in ("-h", "--help"):
          print USAGE
          sys.exit()
      elif o in ("-c", "--catalog"):
          catalog = a
      elif o in ("--eps"):
          eps = True
      else:
          assert False, "unhandled option"

  if len(args) > 0:
    lensID = numpy.array(args,dtype=numpy.int)
  else:
    print USAGE
    sys.exit()

  print "Making plots of",len(lensID),"lens systems"

  imnames = ['A','B','C','D']
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

  # --------------------------------------------------------------------

  # Loop over lenses:

  for k in range(len(lensID)):

    if vb: 
      print "***************************** lens ID ",int(lensID[k]),"*****************************"

    # --------------------------------------------------------------------

    # Search through lens data to find lensID, and pull out row of data:
    
    index = numpy.where(d.field('lensid') == lensID[k])
    if len(index) == 0:
      print "ERROR: lens",lensID[k]," not found in",catalog
      return
    index = int(index[0])
    dd = d[index]

    # Now assign variables from this row:

    xd = 0.0
    yd = 0.0
    nim = dd.field('nimg')
    id = dd.field('lensid')
    if vb: print "Lens number",id," has",nim," images"


    zd     = dd.field('zlens')
    sigma  = dd.field('veldisp')
    zs     = dd.field('zsrc')
    if vb: print "Geometry: zd =",zd,", zs =",zs,", sigma =",sigma,"km/s"

    ms     = dd.field('magi_in')
    m3     = dd.field('magi')
    dtheta = dd.field('imsep')
    epsd   = dd.field('ellip')
    phid   = dd.field('phie')
    gammax = dd.field('gamma')
    phix   = dd.field('phig')
    xs     = dd.field('xsrc')
    ys     = dd.field('ysrc')
    type   = dd.field('flagtype')
    Dd     = dd.field('dd')
    DLd    = dd.field('ddlum')
    absmd  = dd.field('abmag_i')
    md     = dd.field('apmag_i')
    Kcorrd = dd.field('kcorr')
    Ds     = dd.field('ds')
    Dds    = dd.field('dds')
    Sigmac = dd.field('sigcrit')
    DLs    = dd.field('dslum')
    # Effective radius in arcsec:
    size   = dd.field('reff_t')
   
    if vb: print "Lens light: md =",md,", size = ",size,"arcsec, epsd =",epsd,", phid =",phid,"deg"
    
    if vb: print "Lens external shear: gammax =",gammax,", phix =",phix,"deg"

    xi  = dd.field('ximg')
    yi  = dd.field('yimg')
    mui = dd.field('mag')
    ti  = dd.field('delay')
    fstari = dd.field('fstar')
    ri = numpy.zeros(nim)
    phi = numpy.zeros(nim)
    for i in range(nim):
      ri[i]  = numpy.sqrt(xi[i]*xi[i]+yi[i]*yi[i])
    for i in range(nim):
      if i == (nim-1):
        j = 0
      else:
        j = i + 1
      phi[i] = numpy.arccos((xi[j]*xi[i] + yi[j]*yi[i])/(ri[j]*ri[i]))*180.0/3.141592654
#     for i in range(nim):
#       if vb: print "Image",i,": x,y,r,phi,mu,t,fstar =",xi[i],yi[i],ri[i],phi[i],mui[i],ti[i],fstari[i]

    # Good - got lens!

    # Sort images by opening angle:

#     print "Re-ordering by image opening angle..."
#     index = range(nim)
#     index = numpy.flipud(numpy.argsort(numpy.abs(phi)))
#     ti = ti[index]
#     xi = xi[index]
#     yi = yi[index]
#     ri = ri[index]
#     phi = phi[index]
#     mui = mui[index]
#     fstari = fstari[index]
#     for i in range(nim):
#       if vb: print "Image",i,": x,y,r,phi,mu,t,fstar =",xi[i],yi[i],ri[i],phi[i],mui[i],ti[i],fstari[i]

    # Sort images by arrival time:

    if vb: print "Re-ordering by arrival time..."
    index = range(nim)
    index = numpy.argsort(ti)
    ti = ti[index]
    xi = xi[index]
    yi = yi[index]
    ri = ri[index]
    phi = phi[index]
    mui = mui[index]
    fstari = fstari[index]
    for i in range(nim):
      if vb: print "Image",i,": x,y,r,phi,mu,t,fstar =",xi[i],yi[i],ri[i],phi[i],mui[i],ti[i],fstari[i]

    # Compute image magnitudes:
    mi = numpy.zeros(nim)
    lfi = numpy.zeros(nim)
    for i in range(nim):
      mi[i] = ms - 2.5*numpy.log10(numpy.abs(mui[i]))
      lfi[i] = 0.4*(24-mi[i])
    if vb: print "Lens, image magnitudes:",md,mi
    lfd = 0.4*(24-md)
    if vb: print "Lens, image log fluxes:",lfd,lfi

    # Compute time delays relative to first image:
    dt = numpy.zeros(nim-1)
    for i in range(nim-1):
      dt[i] = ti[i+1] - ti[0]
    if vb: print "Time delays:",dt

#     # Work out cuspiness/opening angle:
# 
#     if nim == 2:
#     #  r0 = numpy.sqrt(xi[0]*xi[0] + yi[0]*yi[0])
#     #  r1 = numpy.sqrt(xi[1]*xi[1] + yi[1]*yi[1])
#     #  theta = numpy.arccos((xi[0]*xi[1] + yi[0]*yi[1])/(r0*r1))*180.0/3.141592654
#       print "Double opening angle =",phi[1],"deg"
#     elif nim == 4:
#       crossness = (phi[3]+phi[2])/(phi[0]+phi[1])
#       cuspiness = 1.0 - crossness
#       foldiness = 1.0 - phi[3]/phi[2]
#       if vb: 
#         print "Quad crossness =",crossness
#         print "Quad cuspiness =",cuspiness
#         print "Quad foldiness =",foldiness

    # ------------------------------------------------------------------
    # Compute caustics and critical curves:




    # ------------------------------------------------------------------

    # Start figure:
    fig = pylab.figure(figsize=(8,8))
    # BUG: apsect needs to be equal.

    # Axes limits, useful sizes:
    xmax = 1.99
    dm = 1.0/10

    # Plot command sets its own axes. 'bp' = blue pentagons
    # pylab.plot(xi, yi, 'bp')
    pylab.plot(xi, yi, color='blue', \
               marker='+', markersize=10, markeredgewidth=2, \
               linestyle='')
    pylab.plot(xs, ys, color='lightblue', \
               marker='+', markersize=10, markeredgewidth=2, \
               linestyle='')
    pylab.plot(xd, yd, color='orange', \
               marker='+', markersize=10, markeredgewidth=2, \
               linestyle='')

    # Ellipse to represent lens size:
    bd = size*numpy.sqrt(1.0 - epsd)
    ad = size/numpy.sqrt(1.0 - epsd)
    ell = Ellipse((xd,yd),(ad,bd),orientation=phid, alpha=0.2, fc='orange')    
    pylab.gca().add_patch(ell)
#     cir = pylab.Circle((xd,yd), radius=dm*lfd, alpha=0.2, fc='orange')
#     pylab.gca().add_patch(cir)
    
    # Circles to represent image brightness:
    for i in range(nim):
      cir = pylab.Circle((xi[i],yi[i]), radius=dm*lfi[i], alpha=0.2, fc='blue')
      pylab.gca().add_patch(cir)
    # Label images with image name:
      text = imnames[i]
      pylab.annotate(text, (xi[i]-0.3,yi[i]-0.02), fontsize=12)
    # Label images with f* value:
      text = "%4.2f" % fstari[i]
      pylab.annotate(text, (xi[i]-0.3,yi[i]-0.12), fontsize=12)

    # Circle to represent seeinf:
    if showdisk:
      cir = pylab.Circle((1.5,-1.5), radius=0.7/2.0, alpha=0.1, fc='grey')
      pylab.gca().add_patch(cir)

#     # Legend giving opening angle/cuspiness:
#     if nim == 2:
#       text = "$\Delta\phi$ = %5.1f deg" % phi[0]
#     elif nim == 4:
#       text = "C = %4.2f, F = %4.2f" % (cuspiness, foldiness)
#     else:
#       text = "Naked cusp"
#     pylab.annotate(text, (10,20), xytext=None, fontsize=14, \
#                      xycoords='axes points',textcoords='axes points')

    # Legend giving time delays:
    if nim == 2:
      text = "$\Delta t_{\\rm AB}$ = %5.1f days" % dt[0]
      pylab.annotate(text, (15,20), xytext=None, fontsize=18, \
                     xycoords='axes points',textcoords='axes points')
    elif nim == 4:
      text = "$\Delta t_{\\rm AB}$ = %5.1f days" % dt[0]
      pylab.annotate(text, (15,60), xytext=None, fontsize=18, \
                     xycoords='axes points',textcoords='axes points')
      text = "$\Delta t_{\\rm AC}$ = %5.1f days" % dt[1]
      pylab.annotate(text, (15,40), xytext=None, fontsize=18, \
                     xycoords='axes points',textcoords='axes points')
      text = "$\Delta t_{\\rm AD}$ = %5.1f days" % dt[2]
      pylab.annotate(text, (15,20), xytext=None, fontsize=18, \
                     xycoords='axes points',textcoords='axes points')

    # Legend giving lens, source redshift:
    text1 = "$z_d$ = %5.2f" % zd
    text2 = "$z_s$ = %5.2f" % zs
    pylab.annotate(text1, (15,430), xytext=None, fontsize=18, \
                     xycoords='axes points',textcoords='axes points')
    pylab.annotate(text2, (15,410), xytext=None, fontsize=18, \
                     xycoords='axes points',textcoords='axes points')

    # Plot title:
    title = "Lensed QSO, ID="+str(id)
    pylab.title(title,fontsize=20)
    # Set axes labels:
    pylab.xlabel("x / arcsec",fontsize=20)
    pylab.ylabel("y / arcsec",fontsize=20)
    # Set axis limits:
    pylab.axis([-xmax,xmax,-xmax,xmax])
    # Add a grid:
    pylab.grid(color='grey', linestyle='--', linewidth=0.5)

    # Plot graph to file:
    if eps:
      epsfile = "oguri_qso_ID="+str(id)+".eps"
      pylab.savefig(epsfile)
      print "Figure saved to file:",epsfile
    else:
      pngfile = "oguri_qso_ID="+str(id)+".png"
      pylab.savefig(pngfile)
      print "Figure saved to file:",pngfile

    # Interactive graphics:
    # pylab.show()

  sys.exit()

# ======================================================================

# Fake an ellipse using an N-sided polygon
def Ellipse((x,y), (rx, ry), orientation=0, resolution=100, **kwargs):
    theta = 2*numpy.pi*pylab.arange(resolution)/resolution
    xs = x + rx * numpy.cos(theta)
    ys = y + ry * numpy.sin(theta)
    xr =  xs * numpy.cos(orientation) - ys * numpy.sin(orientation)
    yr =  xs * numpy.sin(orientation) + ys * numpy.cos(orientation)
    return pylab.Polygon(zip(xr, yr), **kwargs)

# ======================================================================

# If called as a script, the python variable __name__ will be set to 
# "__main__" - so test for this, and execute the main program if so.
# Writing it like this allows the function plot_oguri_lens to be called
# from the python command line as well as from the unix prompt.

if __name__ == '__main__':
  plot_oguri_lens(sys.argv)

# ======================================================================


