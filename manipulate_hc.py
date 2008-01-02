#
# manipulate a HC type data profile
#
import math
import pylab as p
import matplotlib
from matplotlib.backends.backend_gtk import FigureCanvasGTK, NavigationToolbar   
matplotlib.use('GTK')


class ManipulateXYData:
    """

    handles x y data that specifies HC-type profiles
    
    use like:

    >>>

    mp = ManipulateXYData(filename,mode)
    p.connect('button_press_event', mp.on_click)
    p.connect('button_release_event', mp.on_release)

    <<<

    INPUT

    filename: data name
    
    mode: 1: data are viscosity
          2: data are density scaling factors


    xtol amd ytol are relative tolerances

    inspired by http://www.scipy.org/Cookbook/Matplotlib/Interactive_Plotting
    
    """
    
    # initialize class
    def __init__(self, filename, mode, xtol = None, ytol = None):
        if xtol == None:
            xtol = 0.1
        if ytol == None:
            ytol = 0.1
        self.xtol = xtol
        self.ytol = ytol
        

        self.visc_norm = 1.e21  # some constants
        self.radius_km = 6371.
        self.cmb_km = 2891.


        self.figure = p.figure()
        self.axis = self.figure.add_subplot(111)

        #
        # 1: viscosity 
        # 2: density 
        self.plot_mode = mode
        #
        # read data 
        self.read_data(filename,mode)
        #
        # convert to plotting format
        self.convert_data(self.plot_mode,False)

        self.datax0 = self.datax # copy for restore
        self.datay0 = self.datay

        self.moving = -1        # if point are being move


        self.zlabels = [300,660,1750] # for plot
        #
        self.verbose = True     # progress output

        self.xl = []
        self.yl = []

        # start a plot
        self.redraw_plot()        


    def distance(self, x1, x2, y1, y2):
        """
        return the distance between two points
        """
        return(math.sqrt( (x1 - x2)**2 + (y1 - y2)**2 ))


    def __call__(self, event):
        #
        # get the x and y data coords
        #
        x, y = event.xdata, event.ydata

        if event.inaxes:
            print 'generic call?: data coords', x,y

    def on_click(self, event):
        # 
        # mouse click event
        # 
        if event.inaxes:        # within region
            # data coordinates
            x, y = event.xdata, event.ydata
            if self.axis == event.inaxes:
                #
                # look for close points
                #
                cps = []
                i=0
                data_cont = zip(self.datax,self.datay) # reformat
                for xd,yd in data_cont: # compute distances for
                    # those points that are
                    # within range compute
                    # tolerance
                    if xd == 0:         # compute tolerances
                        xts = 0.5
                    else:
                        xts = abs(xd)*self.xtol
                    if yd == 0:
                        yts = 0.5
                    else:
                        yts = abs(yd)*self.ytol
                    #
                    # if close, compute distance
                    if  (abs(x-xd) < xts) and  (abs(y-yd) < yts) :
                        cps.append( (self.distance(xd,x,yd,y),xd,yd,i) )
                    i=i+1
                if cps:             # if we found some point close enough, sort them and use the closest
                    cps.sort()
                    dist, xd, yd, i = cps[0] # closest

                if event.button == 2: # center mouse click: add point to list
                    if not cps or dist > 1:
                        if self.verbose:
                            print 'adding data point %7.2f, %7.2f ' % (x, y)
                        self.datax.append(x)
                        self.datay.append(y)
                        self.redraw_plot()
                    else:
                        if self.verbose:
                            print 'there is already a point at %7.2f, %7.2f ' % (x, y)
                else:
                    # 
                    # left or right
                    # 
                    if cps:
                        if event.button == 1: 
                            # left mouse button, move this data point
                            self.moving = i
                            if self.verbose:
                                print 'moving data point %5i ' % i, 'from %7.2f, %7.2f ' % (xd, yd)
                        elif event.button == 3: 
                    # right click: removing this data point
                            if self.verbose:
                                print 'removing data point %7.2f, %7.2f ' % (self.datax[i],self.datay[i])
                            del self.datax[i]
                            del self.datay[i]
                            self.redraw_plot()
                    else:
                        if self.verbose:
                            print 'did not find data close to click  %7.2f, %7.2f' % (x,y)


    def on_release(self, event):
        # mouse has been released
        if event.inaxes:
            xd, yd = event.xdata, event.ydata
            if self.axis == event.inaxes:
                if self.moving > -1: # are we actually moving a point?
                    if self.verbose:
                        print 'assigning %7.2f, %7.2f to data point %5i' % (xd, yd, self.moving)
                    i=0;xn,yn=[],[]
                    data_cont=zip(self.datax,self.datay)
                    # this could be dealt with smarter
                    self.datax, self.datay = [], []
                    for x,y in data_cont: # replace the self.moving-th point with the current location
                        if i==self.moving:
                            self.datax.append(xd);self.datay.append(yd)
                        else:
                            self.datax.append(x);self.datay.append(y)
                        i+=1
                    self.redraw_plot()
                    self.moving = -1 # reset

    def redraw_plot(self):      # refresh the plot
        """

        redraw a plot

        """
        self.sortlevels() # sort data and get layer entries
        
        # get the figure handle
        p.axes(self.axis)
        self.axis.clear()
        if self.plot_mode == 1:  # viscosity 
            xmin,xmax = 1e-3,1e3
            if min(self.datax) < xmin:
                xmin /= 10.
            if max(self.datax) > xmax:
                xmax *= 10.
            self.axis.semilogx(self.datax,self.datay,'o') # plot actual profile
            self.axis.semilogx(self.xl,self.yl,linewidth=3,color='red') # plot layers
            self.axis.set_xlabel('viscosity [1e21 Pas]')
            self.add_pmantle_ornaments(xmin,xmax,True)

        elif self.plot_mode == 2: # density scaling factor
            xmin,xmax = -0.1,0.4
            if min(self.datax) < xmin:
                xmin = self.datax *0.8
            if max(self.datax) > xmax:
                xmax = self.datay *1.2
            self.axis.plot(self.datax,self.datay,'o') # plot actual profile
            self.axis.plot(self.xl,self.yl,linewidth=3,color='blue')
            self.axis.set_title('left mouse: move center: add right: remove point')
            self.axis.set_xlabel('scale factor')
            self.add_pmantle_ornaments(xmin,xmax,False)
# what is the renderer?
#        self.axis.draw('GTKAgg')
        p.draw()
  

    def add_pmantle_ornaments(self,xmin,xmax,uselogx):
        """
        add ornaments typical for the earth's mantle to the plot
        """
        self.axis.grid(True)
        self.axis.set_title('left mouse: move center: add right: remove point')
        self.axis.set_ylabel('depth [km]')
        x = [xmin,xmax];

        if self.plot_mode == 1:
            xoff = xmin*2.5
        else:
            xoff = 0.025*(xmax-xmin)
  
        for z in self.zlabels: # add a few labels
            y=[-z,-z];
            self.axis.text(xmin+xoff,-z+10.,str(z)+' km',fontstyle='italic')
            if uselogx:
                self.axis.semilogx(x,y,linewidth=2,color='black',linestyle='--')
            else:
                self.axis.plot(x,y,linewidth=2,color='black',linestyle='--')
        self.axis.set_ylim([-self.cmb_km, 0])
        self.axis.set_xlim([xmin,xmax])

  

    def reset_data(self,event):
        if self.verbose:
            print 'resetting to original data'
        self.datax = self.datax0
        self.datay = self.datay0
        self.redraw_plot()

    def sortlevels(self):  
        """
        
        sort through a list of weirdly formatted viscosity file
        values and add data point to make a plot look nice

        also, assign layer plot data 

        """
        # sort the z and eta vectors
        data = []
        for zl, el in zip(self.datay, self.datax):
            data.append((zl, el))
        data.sort();
        z,eta = [], []
        for zl,el in data:
            z.append(zl); eta.append(el)

        zn, en =[], []
        n = len(z)
        if n:
            if z[0] > -self.cmb_km:
                zn.append(-self.cmb_km)
                en.append(eta[0])
            i=0
            while i < n:
                zn.append(z[i])
                en.append(eta[i])
                i += 1
            if n and z[n-1] < 0:
                zn.append(0)
                en.append(eta[n-1])
        self.datax = en
        self.datay = zn


        # convert the point-based data to one that can be plotted as 
        # layers
        data = []
        i=0; n = len(self.datax);
        while i < n:
            if i > 0 and self.datay[i] != self.datay[i-1]:
               data.append((self.datay[i],self.datax[i-1]))
            data.append((self.datay[i],self.datax[i]))
            i += 1
        self.xl,self.yl = [],[]
        for yl,xl in data:
            self.xl.append(xl)
            self.yl.append(yl)

    def read_data(self,filename,mode):
        """
        
        read HC profile data from file and return datax, datay
        
        """
        self.datax,self.datay = [],[]
        f=open(filename,'r')
        for line in f:
            val = line.split()
            if len(val) != 2:
                print 'error file ', filename, ' appears to be in wrong format'
                print 'expecting'
                if self.mode == 1:
                    print 'radius[non_dim] viscosity[Pas]'
                elif self.mode == 2:
                    print 'radius[non_dim] density_scale'
                else:
                    print 'unknown'
                exit()
            self.datax.append(val[0])
            self.datay.append(val[1])
        f.close()



    def convert_data(self,mode,reverse):
        """ convert input data to plotting format and vice versa """
        tmpx, tmpy = [],[]
        i=0
        for i in range(0,len(self.datax)):
            if mode == 1:       # viscosity
                if not reverse:
                    tmpx.append(float(self.datay[i])/self.visc_norm)
                    tmpy.append(-(1.-float(self.datax[i]))*self.radius_km)
                else:
                    tmpy.append(float(self.datax[i])*self.visc_norm)
                    tmpx.append((self.radius_km+float(self.datay[i]))/self.radius_km)
            elif mode == 2:     # density 
                if not reverse:
                    tmpx.append(float(self.datay[i]))
                    tmpy.append(-(1.-float(self.datax[i]))*self.radius_km)
                else:
                    tmpy.append(float(self.datax[i]))
                    tmpx.append((self.radius_km+float(self.datay[i]))/self.radius_km)
        self.datax,self.datay = tmpx,tmpy
 

    def save_and_exit(self,event):
        if self.verbose:
            print 'saving modified data'
        filename = 'tmp.dat'
        if self.plot_mode == 1:
            print 'saving modified viscosity profile data to ', filename
        elif self.plot_mode == 2:
            print 'saving modified density profile data to ', filename
            
        #
        # convert data back
        self.convert_data(self.plot_mode,True)
        f=open(filename,'w')
        i=0
        for i in range(0,len(self.datax)):
            ostring = "%8.5f\t%12.7e\n" % (self.datax[i], self.datay[i])
            f.writelines(ostring)
        f.close()
        p.close(self.figure)

    def exit(self,event):
        if self.verbose:
            print 'exiting without saving'
        p.close(self.figure)
