#! /usr/bin/env python
#
# plot and mofiy a HC viscosity file
#
from optparse import OptionParser
import math
import pylab as p
import matplotlib
from matplotlib.backends.backend_gtk import FigureCanvasGTK, NavigationToolbar   
matplotlib.use('GTK')
from manipulate_hc import *


parser = OptionParser()
parser.add_option("--vf", "--viscosity-file", type="string", dest="filename", default="visc.C", \
                      help="HC type viscosity file to display and modify")
(options, args) = parser.parse_args()
filename = options.filename

print 'using ', filename, ' as HC viscosity profile '


mp = ManipulateXYData(filename,1)
#
# connect main data window plot handling routines with mouse click procedures
p.connect('button_press_event', mp.on_click)
p.connect('button_release_event', mp.on_release)
# add some GUI neutral buttons
reset_button =     matplotlib.widgets.Button(p.axes([0.85, 0.025, 0.1, 0.04]), 'Reset',  hovercolor='gray')
cancel_button =    matplotlib.widgets.Button(p.axes([0.75, 0.025, 0.1, 0.04]), 'Cancel',  hovercolor='gray')
done_button =      matplotlib.widgets.Button(p.axes([0.65, 0.025, 0.1, 0.04]), 'Save',   hovercolor='gray')
# button handling routines
reset_button.on_clicked(mp.reset_data)
done_button.on_clicked(mp.save_and_exit)
cancel_button.on_clicked(mp.exit)

p.show()
    



