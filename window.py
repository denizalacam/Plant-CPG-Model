#! /usr/bin/python

import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

try:
    from flask import Flask
    print '# web server available'

except:
    print '# web server not available'


for k in plt.rcParams.keys():

	if k.split('.')[0] == 'keymap':
		plt.rcParams[k] = ''


class window(object):

	title = "Window"
	figsize = (5, 3)
	debug = False

	def __init__(self, position=None, infoWindow=None):

		self.fig = plt.figure(type(self).title, figsize=type(self).figsize)

                self.key_func_dict = dict(escape=window.quit)
		self.fig.canvas.mpl_connect('key_press_event', self.on_key)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focusIn)
		self.info = infoWindow


		if not position == None:

			try:
				self.fig.canvas.manager.window.wm_geometry(position)

			except:
				pass


	def on_key(self, event):

                if event.key == 'escape':
		        self.key_func_dict[event.key](self)

                else:
			if type(self).debug:
                                self.key_func_dict[event.key](self)

			else:
		        	try:
                                	self.key_func_dict[event.key](self)
		        	except:
                                	self.key_func_dict[event.key] = lambda x: None



	def echo(self, string):
		if self.info == None: 	print string
		else: 			self.info.set(string)


	def focusIn(self, string):
		pass


        def quit(self):
            print "exiting ..."
            exit(0)


	"""
	if websupport:
		@...
		def toHtml(self):
			return m.fig_to_html(self.fig)
		pass
	"""


if __name__ == "__main__":
	win = window()
	plt.show()
