import sys
import os
import PySimpleGUI as sg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import csv

os.chdir('..')
path = os.getcwd()+"/src/"
sys.path.append(os.getcwd())

import correa


def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

AppFont = 'Any 16'
sg.theme('LightGrey')

layout_main_welcome = [[sg.Text("Correa main menu")], [sg.Button("SINGLE CONTOUR MODE")], [sg.Button("BATCH MODE")]]

layout_main_single = [[sg.Text("Correa single contour mode")], [sg.Button("LOAD CONTOUR"), sg.Button("LOAD FOCAL POINT"),sg.Button("PLOT"), sg.Button("INFO")],[sg.Button("CLOSE")]]


layout_main_batch =  [[sg.Text("Correa batch contour mode")], [sg.Button("LOAD CONTOUR"), sg.Button("LOAD FOCAL POINT"),sg.Button("PLOT"), sg.Button("INFO")],[sg.Button("CLOSE")]]

layout_main = layout_main_welcome
window_main = sg.Window("Correa", layout_main, resizable=True)

loaded = False


# Create an event loop
while True:
	event_main, values_main = window_main.read()
    # End program if user closes window or
	if event_main == "SINGLE CONTOUR MODE":
		window_single = wg.Window("Correa single contour mode", layout_main_single, resizable=True)
		window_main.close()# = sg.Window("Correa", layout_main_single, resizable = True)
  
	if event_main == "BATCH MODE":
		window_main = sg.Window("Correa", layout_main_batch, resizable = True)
    
	if event_main == "LOAD CONTOUR":
		filename = sg.popup_get_file('Enter the file containing the contour you wish to process.')
		if ( values_main['focal'] == True):
			focal_file = sg.popup_get_file('Enter the file contaning the focal point you wish to use.')
		else:
			p = correa.create_polygon(filename)
		loaded = True


		
	#if (event_main == "ANALYSE") & (loaded == True):
	#	correa.analyse_polygon(p)
  
	#if (event_main == "ANALYSE") & (loaded == False):
	#	sg.popup_ok("You have not loaded a contour please load one.", title = "Warning: no contour")
	
	if (event_main == "PLOT") & (loaded == True):
		if values_main['focal'] == False:
			_VARS = {'window': False}
			layout = [[sg.Canvas(key='figCanvas')],
			[sg.Button('Exit', font=AppFont)]]
			_VARS['window'] = sg.Window('Correa: display contour',
								layout,
								finalize=True,
								resizable=True,
								element_justification="right",
								size=(500,500), 
								modal=False)
			x = []
			y = []
			for i in range(len(p.vertices())):
				x.append(p.vertices()[i][0])
				y.append(p.vertices()[i][1])
			# make fig and plot
			fig = plt.figure()
			plt.plot(x, y)
			# Instead of plt.show
			draw_figure(_VARS['window']['figCanvas'].TKCanvas, fig)
			# MAIN LOOP
			while True:
				event_plot, values_plot = _VARS['window'].read(timeout=200)
				if event_plot == sg.WIN_CLOSED or event_plot == 'Exit':
					break
			_VARS['window'].close()
		else:
			window_plot = sg.popup_error("oops not ready yet")#, [[sg.Text("this functioinality is not ready yet.")], [sg.Button("Close")]])
	
	if (event_main == "PLOT") & (loaded == False):
		sg.popup_cancel("You have not loaded a contour, please load one.", title = "Warning: no contour")
  
	if (event_main == "INFO") & (loaded == True):
		polygon_size = p.size()
		layout_info = [[sg.Text(f"information for contour {filename}")], [sg.Text(f"Number of points in contour {polygon_size}.\nWilmore energy is {p.willmore()}")]]
		window_info = sg.Window("Correa information", layout_info, modal = False)
  
		while True:
			event_info, values_info = window_info.read()
			if event_info == "CLOSE" or event_info == sg.WIN_CLOSED:
				break

	if (event_main == "INFO") & (loaded == False):
		sg.popup_ok("You have not loaded a contour, please load one.", title = "Warning: no contour")

	
	if event_main == "CLOSE" or event_main == sg.WIN_CLOSED:
		break

if loaded == True:
    del(p)

window_main.close()




