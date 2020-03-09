# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import wx
#import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from relpermtable import RelPermTable as RelPermCl

class KrPanel(wx.Panel):
       
    def __init__(self, parent):
        super().__init__(parent)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.list_ctrl = wx.ListCtrl(
                self, size = (-1, 300),
                style = wx.LC_REPORT | wx.BORDER_SUNKEN
        )
        
        self.list_ctrl.InsertColumn(0,  'Sw', width=100)
        self.list_ctrl.InsertColumn(1, 'Krw', width=100)
        self.list_ctrl.InsertColumn(2, 'Krg', width=100)
        self.list_ctrl.InsertColumn(3,  'Pc', width=100)
        main_sizer.Add(self.list_ctrl, 0, wx.ALL | wx.EXPAND, 5)
        
        self.SetSizer(main_sizer)
        
    def load_KrPc_data(self, in_dfr):
        for index, row in in_dfr.iterrows():
            self.list_ctrl.InsertItem(index, "{:.3f}".format(row["Sg"]))
            self.list_ctrl.SetItem(index, 1, "{:.6f}".format(row["Krg"]))
            self.list_ctrl.SetItem(index, 2, "{:.6f}".format(row["Krw"]))
            self.list_ctrl.SetItem(index, 3, "{:.6f}".format(row["Pc"]))
        
class PlotPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        right_sizer = wx.BoxSizer(wx.VERTICAL)
        self.fig, self.ax1 = plt.subplots()
        self.ax1.grid()
        self.canvas = FigureCanvas(self, -1, self.fig)
        right_sizer.Add(self.canvas, 0, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(right_sizer)
        self.unsc_plot = False
        self.scal_plot = False
            
    def PlotKrPcData(self, df_in, df_in2 = None, in_swco = 0, in_sgco = 0):
        if not self.unsc_plot:
            self.ax1.plot(df_in['Sg'], df_in['Krg'], 'r-', df_in['Sg'], df_in['Krw'], 'g-')
            self.unsc_plot = True
        if df_in2 is not None and not self.scal_plot:
            self.line3, self.line4 = self.ax1.plot(df_in2['Sg'], df_in2['Krg'], 'k--', df_in2['Sg'], df_in2['Krw'], 'b--')
            self.scal_plot = True
        elif df_in2 is not None and self.scal_plot:
            self.line3.set_xdata(df_in2['Sg'])
            self.line3.set_ydata(df_in2['Krg'])
            self.line4.set_xdata(df_in2['Sg'])
            self.line4.set_ydata(df_in2['Krw'])

        self.ax1.plot((1-in_swco), 0, 'go')
        self.ax1.annotate('1-SWL', ((1-in_swco), 0.05))
        self.ax1.plot((in_sgco), 0, 'go', label = 'SGL')
        self.ax1.annotate('SGL', (in_sgco, 0.05))
        self.canvas.draw()

#Class for the right bottom panel for entry of EPS data
class EPSPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        lbl_SGCR = wx.StaticText(self, label="SGCR", pos=(20, 20))
        self.ent_SGCR = wx.TextCtrl(self, value="", pos=(100,20), size=(200,-1))
        lbl_SWCR = wx.StaticText(self, label="SWCR", pos=(20, 50))
        self.ent_SWCR = wx.TextCtrl(self, value="", pos=(100,50), size=(200,-1))
        lbl_SWU = wx.StaticText(self, label="SWU", pos=(20, 80))
        self.ent_SWU = wx.TextCtrl(self, value="", pos=(100,80), size=(200,-1))
        lbl_SGU = wx.StaticText(self, label="SGU", pos=(20, 110))
        self.ent_SGU = wx.TextCtrl(self, value="", pos=(100,110), size=(200,-1))
        lbl_KRG = wx.StaticText(self, label="KRG", pos=(20, 140))
        self.ent_KRG = wx.TextCtrl(self, value="", pos=(100,140), size=(200,-1))
        lbl_KRW = wx.StaticText(self, label="KRW", pos=(20, 170))
        self.ent_KRW = wx.TextCtrl(self, value="", pos=(100,170), size=(200,-1))
        lbl_KRGR = wx.StaticText(self, label="KRGR", pos=(20, 200))
        self.ent_KRGR = wx.TextCtrl(self, value="", pos=(100,200), size=(200,-1))
        lbl_KRWR = wx.StaticText(self, label="KRWR", pos=(20, 230))
        self.ent_KRWR = wx.TextCtrl(self, value="", pos=(100,230), size=(200,-1))
        self.calc_button = wx.Button(self, label="Calculate", pos=(110,260))

class KrFrame(wx.Frame):
    def __init__(self):
        super().__init__(parent=None,
             title='Gas Relative Permeability Editor', size=(900, 800))
        self.sp = wx.SplitterWindow(self)
        self.rightSplitter = wx.SplitterWindow(self.sp) #Another splitter to split right panel into two vertical ones
        self.leftSplitter = wx.SplitterWindow(self.sp)
        self.panel01 = KrPanel(self.leftSplitter)
        self.panel02 = PlotPanel(self.rightSplitter)
        self.panel03 = EPSPanel(self.rightSplitter) #Third panel for scaled end point entry
        self.panel03.calc_button.Bind(wx.EVT_BUTTON, self.on_btn)
        self.panel04 = KrPanel(self.leftSplitter)
        self.rightSplitter.SplitHorizontally(self.panel02, self.panel03, 400) #Splitting right panel into two horizontally
        self.leftSplitter.SplitHorizontally(self.panel01, self.panel04, 400)
        self.sp.SplitVertically(self.leftSplitter, self.rightSplitter, 450)
        self.create_menu()
        self.Show()
        
    def on_btn(self, event):
        try:
            self.set_SGCR = float(self.panel03.ent_SGCR.GetValue())
            self.set_SWCR = float(self.panel03.ent_SWCR.GetValue())
            self.set_SGU = float(self.panel03.ent_SGU.GetValue())
            self.set_SWU = float(self.panel03.ent_SWU.GetValue())
            self.set_KRG = float(self.panel03.ent_KRG.GetValue())
            self.set_KRW = float(self.panel03.ent_KRW.GetValue())
            self.set_KRGR = float(self.panel03.ent_KRGR.GetValue())
            self.set_KRWR = float(self.panel03.ent_KRWR.GetValue())
            self.object.three_point_scaling(set_SGCR = self.set_SGCR, set_SWCR = self.set_SWCR, set_SGU = self.set_SGU, set_SWU = self.set_SWU,
                                            set_KRG = self.set_KRG, set_KRW = self.set_KRW, set_KRGR = self.set_KRGR, set_KRWR = self.set_KRWR)
            self.panel04.load_KrPc_data(self.object.sc_df)
            self.panel02.PlotKrPcData(self.object.df, self.object.sc_df, 1 - self.set_SWU, 1 - self.set_SGU)
        except ValueError:
            wx.MessageBox('Check input!', caption='Error', style=wx.OK | wx.ICON_EXCLAMATION)
        
    def create_menu(self):
        menu_bar = wx.MenuBar()
        file_menu = wx.Menu()
        open_file_menu_item = file_menu.Append(
                wx.ID_ANY, 'Open File', 'Open a Relative Permeability File'
        )
        save_file_menu_item = file_menu.Append(
                wx.ID_ANY, 'Save Scaled File', 'Save Scaled Relative Permeability File'
        )
        exit_menu_item = file_menu.Append(
                wx.ID_EXIT, "Exit", "Exit the application"
        )
        menu_bar.Append(file_menu, '&File')
        self.Bind(
                event = wx.EVT_MENU,
                handler = self.on_open_file,
                source = open_file_menu_item
        )
        self.Bind(
                event = wx.EVT_MENU,
                handler = self.on_save_file,
                source = save_file_menu_item
        )
        self.Bind(
                event = wx.EVT_MENU,
                handler = self.on_exit,
                source = exit_menu_item
        )
        self.SetMenuBar(menu_bar)
        
    def on_open_file(self, event):
        title = "Choose a Relative Permeability file:"
        dlg = wx.FileDialog(self, title, "", "",
                           "Eclipse include files (*.INC) | *.INC", style = wx.FD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.object = RelPermCl(dlg.GetPath())
            self.object.three_point_scaling()
            self.panel01.load_KrPc_data(self.object.df)
            self.panel04.load_KrPc_data(self.object.sc_df)
            self.panel02.PlotKrPcData(self.object.df, in_swco = self.object.sgco, in_sgco = self.object.swco)
            self.panel03.ent_SGCR.SetValue(str(self.object.sgcr))
            self.panel03.ent_SWCR.SetValue(str(self.object.swcr))
            self.panel03.ent_SWU.SetValue(str(self.object.swmax))
            self.panel03.ent_SGU.SetValue(str(self.object.sgmax))
            self.panel03.ent_KRG.SetValue(str(self.object.krg_max))
            self.panel03.ent_KRW.SetValue(str(self.object.krw_max))
            self.panel03.ent_KRGR.SetValue(str(self.object.krgr))
            self.panel03.ent_KRWR.SetValue(str(self.object.krwr))
        dlg.Destroy()
        
    def on_save_file(self, event):
        title = "Save Scaled Relative Permeability file:"
        dlg = wx.FileDialog(self, title, "", "",
                           "Eclipse include files (*.INC) | *.INC", style = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            self.object.do_save_data(dlg.GetPath())
        dlg.Destroy()
        
    def on_exit(self, e):
        self.Close()
        
if __name__ == '__main__':
    app = wx.App(False)
    frame = KrFrame()
    app.MainLoop()
    del app