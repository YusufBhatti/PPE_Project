#! /usr/bin/env python
#-----------------------------------------------------------------------------
from twisted.internet import gtk2reactor 

gtk2reactor.install()

from twisted.internet import reactor
from twisted.internet import protocol
from twisted.protocols import basic

from time import ctime

import pygtk

pygtk.require('2.0')

import gtk

PORT = 19876

class LogProtocol(basic.LineReceiver):
    def connectionMade(self):
        clnt = self.clnt = self.transport.getPeer().host
        textbuffer = self.factory.textview.get_buffer()
        enditer = textbuffer.get_end_iter()
        textbuffer.insert(enditer,'[%s] %s %s\n' % (ctime(), ' connected from ', clnt))
#        print '[%s] %s %s' % (ctime(), ' connected from ', clnt)
#        print '... connected from', clnt

    def lineReceived(self,line):
#        print '[%s] %s' % (ctime(), line)
        textbuffer = self.factory.textview.get_buffer()
        enditer = textbuffer.get_end_iter()
        textbuffer.insert(enditer,'[%s] %s\n' % (ctime(), line))

class LogFactory(protocol.ServerFactory):
    def close_application(self, widget):
        gtk.main_quit()

    def __init__(self):
        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        window.set_resizable(True)  
        window.set_default_size(600, 300)
        window.connect("destroy", self.close_application)
        window.set_title("logServer")
        window.set_border_width(0)

        box = gtk.VBox(False, 10)
        box.set_border_width(10)
        window.add(box)
        box.show()

        sw = gtk.ScrolledWindow()
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        textview = gtk.TextView()
        textview.set_editable(False)
        textview.set_cursor_visible(False)
        
        sw.add(textview)
        sw.show()
        textview.show()

        box.pack_start(sw)

        button = gtk.Button("close")
        button.connect("clicked", self.close_application)
        box.pack_start(button, False, False, 0)
        button.set_flags(gtk.CAN_DEFAULT)
        button.grab_default()
        button.show()
        window.show()

        self.textview = textview

    protocol = LogProtocol

#-----------------------------------------------------------------------------

# print 'waiting for connection ...'
reactor.listenTCP(PORT, LogFactory())
reactor.run()

#-----------------------------------------------------------------------------
