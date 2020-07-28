import os

from pylab import *

class spectra:
  def __init__(self, filename):
      specfile = open(filename, 'rb')
      
      # header
      self.z = fromfile(specfile, dtype=float64, count=1)[0]
      self.om = fromfile(specfile, dtype=float64, count=1)[0]     # cosmo parameters
      self.ol = fromfile(specfile, dtype=float64, count=1)[0]
      self.ob = fromfile(specfile, dtype=float64, count=1)[0]
      self.h = fromfile(specfile, dtype=float64, count=1)[0]
      self.box = fromfile(specfile, dtype=float64, count=1)[0]    # boxsize kpc/h
      self.xh = fromfile(specfile, dtype=float64, count=1)[0]     # hydrogen mass fraxtion
      self.nbins = fromfile(specfile, dtype=int32, count=1)[0]    # pixel per LOS
      self.nlos = fromfile(specfile, dtype=int32, count=1)[0]     # number of LOSs
      
      self.a = 1.0/(1.0+self.z)
            
      # postions
      self.dirlos = fromfile(specfile, dtype=int32, count=self.nlos)   # direction 0=x, 1=y, 2=z of the LOSs
      self.xlos = fromfile(specfile, dtype=float64, count=self.nlos)   # coordinates of a point through which the line passes
      self.ylos = fromfile(specfile, dtype=float64, count=self.nlos)
      self.zlos = fromfile(specfile, dtype=float64, count=self.nlos)
      
      # coordinate and velocity of pixel center
      self.pixpos = fromfile(specfile, dtype=float64, count=self.nbins) # first pixel center is at simulation coordinate 0 and has pixpos = 0
      self.pixvel = fromfile(specfile, dtype=float64, count=self.nbins) # pixpos * a * H(a)
      
      self.dxpix = self.box/self.nbins
      self.dxpix_phys = self.dxpix*self.a
      
      assert self.box > 1000.0 # otherwise likely not kpc/h units
      
      self.vmax = (self.pixvel[-1]-self.pixvel[0])/(self.nbins-1)*self.nbins
      self.dvpix = self.pixvel[1]-self.pixvel[0]
      
      # gas overdensity Delta
      self.rhoH_over_rhoHmean = 10.0**fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
      
      # HI Lyman-alpha
      self.nHI_frac = fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
      self.temp_HI = fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
      self.vel_HI = fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
      
      self.tau_HI = fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
            
      curpos = specfile.tell()
      specfile.seek(0,os.SEEK_END)
      endpos = specfile.tell()
      
      self.have_HeII = False
      
      if curpos != endpos:
        specfile.seek(curpos,os.SEEK_SET)
        
        ## HeII Lyman-alpha
        self.nHeII_frac = fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
        self.temp_HeII = fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
        self.vel_HeII = fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
        self.tau_HeII = fromfile(specfile, dtype=float64, count=self.nbins*self.nlos).reshape(self.nlos,self.nbins)
      
        self.have_HeII = True
      
      curpos = specfile.tell()
      if curpos != endpos:
        print(("WARNING: not reading the whole file", curpos, "of", endpos))
      
      specfile.close()

# ---
