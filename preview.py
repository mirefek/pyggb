import os
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

from random_constr import Construction

class DisplayWindow(Gtk.Window):

    def __init__(self, width, height, datadir):
        super(DisplayWindow, self).__init__()

        self.datadir = datadir
        self.fnames = sorted([
            fname for fname in os.listdir(datadir)
            if fname.endswith(".txt")
        ])
        self.index = 0
        self.construction = Construction(display_size = (width, height))
        self.load_construction()

        self.darea = Gtk.DrawingArea()
        self.darea.connect("draw", self.on_draw)
        self.darea.set_events(Gdk.EventMask.KEY_PRESS_MASK)
        self.add(self.darea)

        self.connect("key-press-event", self.on_key_press)

        self.set_title("Display")
        self.resize(width, height)
        self.set_position(Gtk.WindowPosition.CENTER)
        self.connect("delete-event", Gtk.main_quit)
        self.show_all()

    def load_construction(self):
        if self.index == len(self.fnames): self.index = 0
        elif self.index < 0: self.index += len(self.fnames)
        fname = self.fnames[self.index]
        print(self.index, fname)
        self.construction.load(os.path.join(self.datadir, fname))
        self.construction.generate(max_attempts = 0)

    def on_draw(self, wid, cr):

        self.construction.render(cr)

    def on_key_press(self,w,e):

        keyval = e.keyval
        keyval_name = Gdk.keyval_name(keyval)
        #print(keyval_name)
        regenerate = False
        if keyval_name in ("Up", "Down", "Left", "Right"):
            regenerate = True
            if keyval_name in ("Up", "Left"): self.index -= 1
            else: self.index += 1
            self.load_construction()
            
        elif keyval_name == "space":
            regenerate = True
        elif keyval_name == "Escape":
            Gtk.main_quit()
        else:
            return False

        if regenerate:
            self.construction.generate(max_attempts = 0)
            self.darea.queue_draw()

if __name__ == "__main__":
    DisplayWindow(400, 300, "ggb-benchmark/true")
    Gtk.main()
