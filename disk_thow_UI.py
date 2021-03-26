from tkinter import *
import tkinter as tk
from tkinter import ttk
import numpy as np
from disc_throw_onlyfuncs import compute
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import os

class UI():
    def __init__(self):
        self.app_wind = Tk()

    def __init__(self,v0,angle0,om0):
        self.app_wind = Tk()
        self.initv0   = v0
        self.initangle = angle0
        self.initom0 = om0
        self.solName = '%g_%g_%g.npy'%(float(self.initv0),float(self.initangle),float(self.initom0))
        self.sol = np.load('./solutions/%s'%self.solName)

    def updateVelGraph(self):
        self.sol = np.load('./solutions/%s'%self.solName)
        self.pltVel.clf()
        self.velPltAx = self.pltVel.add_subplot(111)
        self.velPltAx.plot(self.sol[-1, :], self.sol[3, :], label='v_x')
        self.velPltAx.plot(self.sol[-1, :], self.sol[4, :], label='v_y')
        self.velPltAx.plot(self.sol[-1, :], self.sol[5, :], label='v_z')
        self.velPltAx.set_xlabel('t')
        self.velPltAx.set_ylabel('v_i')
        self.pltVel.legend()
        self.pltVelCan.draw()
        self.pltVelTB.update()

    def updateTrajGraph(self):
        self.sol = np.load('./solutions/%s'%self.solName)
        self.pltTraj.clf()
        self.trajPltAx = self.pltTraj.add_subplot(111, projection='3d')
        self.trajPltAx.plot(self.sol[0, :], self.sol[1, :],self.sol[2, :])
        self.trajPltAx.set_xlabel('x')
        self.trajPltAx.set_ylabel('y')
        self.trajPltAx.set_zlabel('z')
        self.pltTrajCan.draw()
        self.pltTrajTB.update()

    def updateOMVel(self):
        self.lsSols = os.listdir('./solutions/')
        self.solOm.set(self.solName)
        self.opVel['menu'].delete(0, 'end')
        for choice in self.lsSols:
            self.opVel['menu'].add_command(label=choice, command=tk._setit(self.solOm,choice))

    def on_opVel_change(self,*args):
        self.solName = self.solOm.get()
        self.updateVelGraph()

    def on_opTraj_change(self,*args):
        self.solName = self.solOmTr.get()
        self.updateTrajGraph()

    def on_tab_change(self,event):
        tab = event.widget.tab('current')['text']
        if tab == 'Solver':
            pass

        elif tab == 'Trajectory':
            self.updateTrajGraph()
            self.updateOMVel()

        # update velocity graph
        elif tab == 'Velocity':
            self.updateVelGraph()
            self.updateOMVel()


    def solverBtClick(self):
        self.v0 = float(self.v0En.get())
        self.angl = float(self.angleEn.get())
        self.init_rotation = float(self.init_rotEn.get())
        solution = compute(self.v0, self.angl, self.init_rotation)
        self.sol = np.append(solution.y,solution.t.reshape(1,-1),axis=0)
        self.solName = '%g_%g_%g.npy'%(self.v0, self.angl, self.init_rotation)
        np.save('./solutions/%s'%self.solName,self.sol)


    def runUI(self):
        #general stuff
        self.app_wind.title('Solver of the disk-golf disk throw')
        # self.app_wind.geometry('500x500')

        # header
        self.header = Frame(self.app_wind)
        Label(self.header, font=('Verdana',15),text='Solver of the disk-golf disk throw').pack()
        Label(self.header, font=('Verdana',10), text='authors: Jan Hrůza, Tomáš Hlavatý').pack()
        self.header.pack()

        # tabs
        self.tabs = ttk.Notebook(self.app_wind)
        self.tab1 = ttk.Frame(self.tabs)
        self.tab2 = ttk.Frame(self.tabs)
        self.tab3 = ttk.Frame(self.tabs)
        self.tabs.add(self.tab1, text='Solver')
        self.tabs.add(self.tab3, text='Trajectory')
        self.tabs.add(self.tab2, text='Velocity')

        # solver tab
        self.solver_fr = Frame(self.tab1)
        # column 0
        Label(self.solver_fr, font=('Verdana',12,'bold'),text='Parameters of the throw').grid(row=2,columnspan=2)
        Label(self.solver_fr, text='starting velocity').grid(row=3,column=0)
        Label(self.solver_fr, text='angle of the throw').grid(row=4,column=0)
        Label(self.solver_fr, text='init rotation').grid(row=5,column=0)
        self.solver_bt = Button(master=self.solver_fr, width=20,font=('Helvetica', 15), text='Solve',bg='grey', command=self.solverBtClick)
        self.solver_bt.grid(row=6, column=0,columnspan=2)
        # column 1
        self.v0En = Entry(self.solver_fr)
        self.angleEn = Entry(self.solver_fr)
        self.init_rotEn = Entry(self.solver_fr)
        self.v0En.insert(0, self.initv0)
        self.angleEn.insert(0,self.initangle)
        self.init_rotEn.insert(0, self.initom0)
        self.v0En.grid(row=3, column=1)
        self.angleEn.grid(row=4, column=1)
        self.init_rotEn.grid(row=5, column=1)
        self.solver_fr.pack(fill=BOTH,expand=True)

        # trajectory tab
        self.traj_fr = Frame(self.tab3)
        Label(self.traj_fr,text='option menu: v0_theta_omega.npy').grid(row=0,column=0)
        # option list
        self.lsSols = os.listdir('./solutions/')
        self.solOmTr = StringVar(self.traj_fr)
        self.solOmTr.set(self.lsSols[0])
        self.solOmTr.trace('w', self.on_opTraj_change)
        self.opTr = OptionMenu(self.traj_fr, self.solOmTr,*self.lsSols)
        self.opTr.grid(row=0, column=1)
        # graph
        self.pltTraj = matplotlib.pyplot.figure(figsize=(7, 5), dpi=100)
        self.pltTrajCan = FigureCanvasTkAgg(self.pltTraj,master=self.traj_fr)
        self.pltTrajCan.get_tk_widget().grid(row=1, columnspan=2)
        self.pltTrajTB = NavigationToolbar2Tk(self.pltTrajCan, self.tab3)
        self.updateTrajGraph()
        self.traj_fr.pack()

        # velocitty solutions tab
        self.solution_fr = Frame(self.tab2)
        Label(self.solution_fr,text='option menu: v0_theta_omega.npy').grid(row=0,column=0)
        #option list
        self.lsSols = os.listdir('./solutions/')
        self.solOm = StringVar(self.solution_fr)
        self.solOm.set(self.lsSols[0])
        self.solOm.trace('w',self.on_opVel_change)
        self.opVel = OptionMenu(self.solution_fr, self.solOm , *self.lsSols)
        self.opVel.grid(row=0,column=1)
        # graph
        self.pltVel = matplotlib.pyplot.figure(figsize=(7, 5), dpi=100)
        self.pltVelCan = FigureCanvasTkAgg(self.pltVel,master=self.solution_fr)
        self.pltVelCan.get_tk_widget().grid(row=1,columnspan=2)
        self.pltVelTB = NavigationToolbar2Tk(self.pltVelCan, self.tab2)
        self.updateVelGraph()
        self.solution_fr.pack()
        # packing of the frames
        self.tabs.pack(fill=BOTH,expand=True,)
        self.tabs.bind('<<NotebookTabChanged>>', self.on_tab_change)
        self.app_wind.mainloop()


UI_ = UI(10,20,-30)
UI_.runUI()


# window = Tk()

