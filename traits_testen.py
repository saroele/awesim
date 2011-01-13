
print 'hello world'
#from enthought.traits.api import HasTraits, Str, Int
#import enthought.traits.ui
#
#class SimpleEmployee(HasTraits):
#    first_name = Str
#    last_name = Str
#    department = Str
#
#    employee_number = Str
#    salary = Int
#
#sam = SimpleEmployee()
#sam.configure_traits()


#from enthought.traits.api import *
#from enthought.traits.ui.api import *
#
#class Camera( HasTraits ):
#   """ Camera object """
#   gain = Enum(1, 2, 3,
#      desc="the gain index of the camera",
#      label="gain", )
#   exposure = CInt(10,
#      desc="the exposure time, in ms",
#      label="Exposure", )
#
#   def capture(self):
#      """ Captures an image on the camera and returns it """
#      print "capturing an image at %i ms exposure, gain: %i" % (
#                    self.exposure, self.gain )
#
#if  __name__ == "__main__":
#   camera = Camera()
#   camera.configure_traits()
#   camera.capture()


#from enthought.traits.api import *
#from enthought.traits.ui.api import *
#
#class Check( HasTraits ):
#   """ user check """
#   input = Str
#
#
## def capture(self):
##    """ Captures an image on the camera and returns it """
##    print "capturing an image at %i ms exposure, gain: %i" % (
##                  self.exposure, self.gain )
#
#if  __name__ == "__main__":
#   check = Check(input='y')
#   check.configure_traits()
#   if check.input =='y' or check.input=='Y':
#       print 'True'
#   else:
#       print 'False'


from threading import Thread
from time import sleep
from enthought.traits.api import *
from enthought.traits.ui.api import View, Item, ButtonEditor

class TextDisplay(HasTraits):
    string =  String()

    view= View( Item('string',show_label=False, springy=True, style='custom' ))


class CaptureThread(Thread):
    def run(self):
        self.display.string = 'Camera started\n' + self.display.string
        n_img = 0
        while not self.wants_abort:
            sleep(.5)
            n_img += 1
            self.display.string = '%d image captured\n' % n_img \
                                                    + self.display.string
        self.display.string = 'Camera stopped\n' + self.display.string

class Camera(HasTraits):
    start_stop_capture = Button()
    display = Instance(TextDisplay)
    capture_thread = Instance(CaptureThread)

    view = View( Item('start_stop_capture', show_label=False ))

    def _start_stop_capture_fired(self):
        if self.capture_thread and self.capture_thread.isAlive():
            self.capture_thread.wants_abort = True
        else:
            self.capture_thread = CaptureThread()
            self.capture_thread.wants_abort = False
            self.capture_thread.display = self.display
            self.capture_thread.start()

class MainWindow(HasTraits):
    display = Instance(TextDisplay, ())

    camera = Instance(Camera)

    def _camera_default(self):
        return Camera(display=self.display)

    view = View('display', 'camera', style="custom", resizable=True)


if __name__ == '__main__':
    MainWindow().configure_traits()



import wx

import matplotlib
# We want matplotlib to use a wxPython backend
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

from enthought.traits.api import Any, Instance
from enthought.traits.ui.wx.editor import Editor
from enthought.traits.ui.wx.basic_editor_factory import BasicEditorFactory

class _MPLFigureEditor(Editor):

    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()
        
    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # The panel lets us add additional controls.
        panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        # matplotlib commands to create a canvas
        mpl_control = FigureCanvas(panel, -1, self.value)
        sizer.Add(mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
        toolbar = NavigationToolbar2Wx(mpl_control)
        sizer.Add(toolbar, 0, wx.EXPAND)
        self.value.canvas.SetMinSize((10,10))
        return panel

class MPLFigureEditor(BasicEditorFactory):

    klass = _MPLFigureEditor


if __name__ == "__main__":
    # Create a window to demo the editor
    from enthought.traits.api import HasTraits
    from enthought.traits.ui.api import View, Item
    from numpy import sin, cos, linspace, pi

    class Test(HasTraits):

        figure = Instance(Figure, ())

        view = View(Item('figure', editor=MPLFigureEditor(),
                                show_label=False),
                        width=400,
                        height=300,
                        resizable=True)

        def __init__(self):
            super(Test, self).__init__()
            axes = self.figure.add_subplot(111)
            t = linspace(0, 2*pi, 200)
            axes.plot(sin(t)*(1+0.5*cos(11*t)), cos(t)*(1+0.5*cos(11*t)))

    Test().configure_traits()


from enthought.traits.api import *
from enthought.traits.ui.api import *

class EchoBox(HasTraits):
    input =  Int()
    output = Str(desc='square of input if input<500', label='Square')
    other = Int
    knopke=Button()
    nothing=5.67
    lijstje=List(['een','twee','drie'],desc='mijn lijstje', label= 'a list')
    lijstje2=List(['een','twee','drie'],desc='mijn lijstje2', label= 'list 2')
    
    def __init__(self,value):
        self.input=value
        self.other=value-10

    def _input_changed(self):
        self.other=self.input-10
        if self.input<500:
            self.output = str(self.input**2)
        else:
            self.output='too big for me'
            
    def _knopke_fired(self):
        self.input+=5

view1 = View(HSplit(Group(Item(name='input'),
                         Item(name='output'),
                         Item('knopke', show_label=False),
                         label = 'Interactive input',
                         show_border = True),
                    Group(Item('lijstje'),
                          Item('lijstje2', editor=CheckListEditor(values=['een','twee','drie'])))),
                    buttons = [OKButton, CancelButton],
                    resizable=True,
                    title='Mijn eerste interface!!')

e=EchoBox(600)
print 'Nu kunt ge spelen...'
e.configure_traits(view=view1)

print 'Output is %s \nOther is %d'  % (e.output,e.other)