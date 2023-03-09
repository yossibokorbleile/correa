import sys
import os
import PySimpleGUI as sg



os.chdir('..')
path = os.getcwd()+"/src/"
sys.path.append(os.getcwd())

import correa

layout = [[sg.Text("Hello from PySimpleGUI")], [sg.Button("LOAD")],[sg.Button("CLOSE")]]


# Create the window

window = sg.Window("Demo", layout)
right_click_menu=sg.MENU_RIGHT_CLICK_EDITME_VER_EXIT

MENU_RIGHT_CLICK_EDITME_VER_EXIT = ['', ['Edit Me']]

# Create an event loop
while True:
	event, values = window.read()
    # End program if user closes window or
    # presses the OK button
    
	if event == 'Edit Me':
		sg.execute_editor(__file__)
	if event == "LOAD":
		correa.analyse_polygon('/Users/yossi/AAU/ToMaCo/contours/contours/01kPa/24h_Trypsin24h_2-FITC_001_png_contour.csv')
		
		
	if event == "CLOSE" or event == sg.WIN_CLOSED:
		break

window.close()