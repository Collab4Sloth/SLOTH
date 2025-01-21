import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
import os
import xml.etree.ElementTree as ET
from cycler import cycler
import sys
from scipy.interpolate import interp1d

class SLOTHVisu:

    def __init__(self):
        self._field = "c"
        self._pathPVD = "./"
        self._pathPVTU = "./"
        self._xData = {}
        self._yData = {}
        self._plotHeight = 6
        self._plotLength = 8
        self._timestepsToPlot = [0]
        self._figName = "save"
        self._figTitle = ""
        self._figYLabel = self._field
        self._figXLabel = "Position $x$ [m]"
        self._figFontSize = 10
        self._figMarker = "o"
        self._figMarkerSize = 10
        self._figLw = 2
        self._figLs = "-"
        self._figLabelSize = 10
        self._figLegendFontSize = 10
        self._problemName = "Problem1"
        self.fig = ""
        self.ax = ""
        self._xlim = None
        plt.rc(
            "axes",
            prop_cycle=cycler(
                color=[
                    "black",
                    "red",
                    "blue",
                    "orange",
                    "green",
                    "purple",
                    "brown",
                    "pink",
                    "grey",
                ]
            ),
        )
        if "ipykernel" in sys.modules:
            self._isNotebook = True
        else:
            self._isNotebook = False

    def set_problemName(self, val):
        self._problemName = val

    def set_figName(self, val):
        self._figName = val

    def set_xlim(self, val):
        self._xlim = val

    def set_figLegendFontSize(self, val):
        self._figLegendFontSize = val

    def set_figLabelSize(self, val):
        self._figLabelSize = val

    def set_figMarkerSize(self, val):
        self._figMarker = val

    def set_figMarker(self, val):
        self._figMarker = val

    def set_figLw(self, val):
        self._figLw = val

    def set_figLs(self, val):
        self._figLs = val

    def set_figFontSize(self, val):
        self._figFontSize = val

    def set_figYLabel(self, val):
        self._figYLabel = val

    def set_figXLabel(self, val):
        self._figXLabel = val

    def set_figTitle(self, val):
        self._figTitle = val

    def set_timestepsToPlot(self, val):
        self._timestepsToPlot = val

    def get_timestepsToPlot(self):
        return self._timestepsToPlot

    def set_pathPVD(self, path):
        self._pathPVD = path

    def get_pathPVD(self):
        return self._pathPVD

    def set_pathPVTU(self, path):
        self._pathPVTU = path

    def get_pathPVTU(self):
        return self._pathPVTU

    def set_field(self, field):
        self._field = field

    def get_field(self):
        return self._field

    def set_xData(self, data):
        self.xData = data

    def get_xData(self):
        return self.xData

    def set_yData(self, data):
        self.yData = data

    def get_yData(self):
        return self.yData

    def set_plotHeight(self, value):
        self._plotHeight = value

    def get_plotHeight(self):
        return self._plotHeight

    def set_plotLength(self, value):
        self._plotLength = value

    def get_plotLength(self):
        return self._plotLength

    def readPVTU(self, chemin_pvtu, nom_champ):
        try:
            
            dataset = pv.read(chemin_pvtu)

            
            if nom_champ not in dataset.point_data.keys():
                print(f"Champ '{nom_champ}' non trouvé dans le fichier.")
                return

            
            x = dataset.points[:, 0] 
            valeurs = dataset[nom_champ]  

            
            indices = np.argsort(x)
            x = x[indices]
            valeurs = valeurs[indices]
        except Exception as e:
            print(f"Erreur : {e}")
        return x, valeurs

    def loopOnFolder(self):
        chemin = self._pathPVTU
        sous_dossiers = []
        try:
            # Utilisation de os.walk pour parcourir le répertoire
            for root, dirs, _ in os.walk(chemin):
                for dossier in dirs:
                    sous_dossiers.append(os.path.join(root, dossier))
            return sous_dossiers
        except Exception as e:
            print(f"Error : {e}")
            return []

    def getAllData(self):
        folderWithDataToPlot = self.loopOnFolder()
        self.fig, self.ax = self.createPlot()
        for case in folderWithDataToPlot:
            x, val = self.readPVTU(case + "/data.pvtu", self._field)
            self.ax.plot(x, val, "-o", label= self._field)
            self.ax.legend()
        plt.show()

    def readPVD(self):
        # Parse le fichier XML .pvd
        try:
            tree = ET.parse(self._pathPVD + self._problemName + ".pvd") 
            root = tree.getroot()

            timesteps = []
            PVTUFiles = []

            for dataset in root.findall(".//DataSet"):
                timestep = float(
                    dataset.attrib["timestep"]
                )  # Récupère l'attribut timestep
                file_name = dataset.attrib["file"]  # Récupère l'attribut file

                timesteps.append(timestep)
                PVTUFiles.append(file_name)
            return timesteps, PVTUFiles

        except Exception as e:
            print(f"Error, can't read this .pvd file : {e}")
            return [], []

    def createPlot(self):
        self.fig, self.ax = plt.subplots()
        self.fig.set_size_inches(self._plotLength, self._plotHeight)
        self.ax.tick_params(axis="x", labelsize=self._figLabelSize)
        self.ax.tick_params(axis="y", labelsize=self._figLabelSize)
        self.ax.set_xlabel(self._figXLabel, fontsize=self._figFontSize)
        self.ax.set_ylabel(self._figYLabel, fontsize=self._figFontSize)
        if self._figTitle != "":
            self.ax.set_title(self._figTitle, fontsize=self._figFontSize)
        self.ax.grid(True)
        if self._xlim != None:
            self.ax.set_xlim(self._xlim)
        

    def addDataToPlot(self, x, y, lab=""):
        self.ax.plot(
            x,
            y,
            linestyle=self._figLs,
            linewidth=self._figLw,
            marker=self._figMarker,
            label=lab,
        )

    def finishPlot(self):
        self.ax.legend(fontsize=self._figLegendFontSize)
        self.fig.savefig(self._figName + ".png",bbox_inches='tight')
        if not(self._isNotebook):
            self.fig.show()
            plt.close(self.fig)       

    def plotWantedValues(self):
        ts, files = self.readPVD()
        self.createPlot()
        for i in range(0, len(ts), 1):
            if ts[i] in self._timestepsToPlot:
                x, val = self.readPVTU(self._pathPVD + files[i], self._field)
                self.addDataToPlot(x, val, "t=" + str(ts[i]) + "s")
        self.finishPlot()

    def extractWantedValues(self):
        ts, files = self.readPVD()
        for i in range(0, len(ts), 1):
            if ts[i] in self._timestepsToPlot:
                self._xData[str(ts[i])], self._yData[str(ts[i])] = self.readPVTU(
                    self._pathPVD + files[i], self._field
                )
        return (self._xData, self._yData)
    
    def plotWantedValuesWithAnalyticalSolution(self,pos,lambda_function):
        ts, files = self.readPVD()
        self.createPlot()
        for i in range(0, len(ts), 1):
            if ts[i] in self._timestepsToPlot:
                x, val = self.readPVTU(self._pathPVD + files[i], self._field)
                self.addDataToPlot(x, val, "t=" + str(ts[i]) + "s")
                lines = self.ax.get_lines()
                color = lines[-1].get_color()
                self.ax.plot(pos,[lambda_function(x_,ts[i]) for x_ in pos],color = color, lw = self._figLw, ls = "-")
        self.finishPlot()

    
    def plotWithTransformation(self,lambda_function):
        ts, files = self.readPVD()
        self.createPlot()
        for i in range(0, len(ts), 1):
            if ts[i] in self._timestepsToPlot:
                x, val = self.readPVTU(self._pathPVD + files[i],self._field)
                self.ax.plot(x,lambda_function(val), label = "t=" + str(ts[i]) + "s", lw = self._figLw, ls = self._figLs)
        self.finishPlot()

    def plotWithTransformationAndResize(self,lambda_function):
        ts, files = self.readPVD()
        self.createPlot()
        x_plot = []
        for i in range(0, len(ts), 1):
            if ts[i] in self._timestepsToPlot:
                x, val = self.readPVTU(self._pathPVD + files[i], self._field)
                y_plot = lambda_function(val)
                x_plot = x[:len(y_plot)]
                self.ax.plot(x_plot, y_plot, label = "t=" + str(ts[i]) + "s", lw = self._figLw, ls = self._figLs)
        self.finishPlot()

    def plotErrorComparaison(self,xData1,yData1,xData2,yData2):
        self.createPlot()
        keys1 = [float(key) for key in xData1.keys()]
        keys2 = [float(key) for key in xData2.keys()]
        err = []
        for i in range(0,len(keys1),1):
            if keys1[i] in keys2:
                err = abs((yData1[str(keys1[i])][:] - yData2[str(keys1[i])][:]) / yData1[str(keys1[i])][:])
                self.ax.semilogy(xData1[str(keys1[i])], err, label = "t=" + str(keys1[i]) + "s", lw = self._figLw, ls = self._figLs)
        self.finishPlot()
    
    def plotComparaison(self,xData1,yData1,xData2,yData2):
        self.createPlot()
        keys1 = [float(key) for key in xData1.keys()]
        keys2 = [float(key) for key in xData2.keys()]
        for i in range(0,len(keys1),1):
            if keys1[i] in keys2:
                self.ax.plot(xData1[str(keys1[i])], yData1[str(keys1[i])], label = "t=" + str(keys1[i]) + "s", lw = self._figLw, ls = "solid")
                lines = self.ax.get_lines()
                color = lines[-1].get_color()
                self.ax.plot(xData2[str(keys1[i])], yData2[str(keys1[i])], label = "t=" + str(keys1[i]) + "s", lw = self._figLw, ls = "dashdot",color=color)
        self.finishPlot()

    def plotIsoValues(self, targetIso):
        ts, files = self.readPVD()
        self.createPlot()
        ts_ = []
        pos = []
        for i in range(0, len(ts), 1):
            if ts[i] in self._timestepsToPlot:
                x, val = self.readPVTU(self._pathPVD + files[i], self._field)
                f = interp1d(val, x, kind='linear', bounds_error=False, fill_value="extrapolate")
                ts_.append(ts[i])
                pos.append(f(targetIso))
        self.addDataToPlot(ts_, pos)
        self.finishPlot()

    
    def plotIsoValuesForMultipleSim(self,x_,y_,targetIso, label = 'None'):
        self.createPlot()
        if label == 'None':
            label = ['' for i in range(len(x_))]
        for i in range(len(x_)):
            pos = []
            time = []
            for key in x_[i].keys():
                f = interp1d(y_[i][key], x_[i][key], kind='linear', bounds_error=False, fill_value="extrapolate")
                pos.append(float(f(targetIso)))
                time.append(float(key))
            self.addDataToPlot(time, pos, lab=label[i])
        self.finishPlot()

    def plotErrorIsoValues(self,x_ref,y_ref,x_,y_,targetIso, label = 'None'):
        self.createPlot()
        self.ax.semilogy([], [], label = '', lw = 0, ls = self._figLs)
        if label == 'None':
            label = ['' for i in range(len(x_))]
        for i in range(len(x_)):
            error = []
            time = []
            for key in x_[i].keys():
                if (key in x_ref.keys() and float(key) !=0):
                    f = interp1d(y_[i][key], x_[i][key], kind='linear', bounds_error=False, fill_value="extrapolate")
                    f_ref = interp1d(y_ref[key], x_ref[key], kind='linear', bounds_error=False, fill_value="extrapolate")
                    curr_error = abs(float(f(targetIso)) - float(f_ref(targetIso))) / abs(float(f_ref(targetIso)))
                    error.append(curr_error)
                    time.append(float(key))
            self.ax.semilogy(time, error, label = label[i], lw = self._figLw, ls = self._figLs)
            print(label[i] + "    " + str(error[-1]))
        self.finishPlot()


    def plotMultipleData(self,x,y,ls = None):
        if ls == None :
            ls = ['-' for i in range(len(x))] 
        self.createPlot()
        for key in x[0].keys():
            ts = float(key)
            if ts in self._timestepsToPlot:
                self.ax.plot(x[0][key], y[0][key], label = "t=" + key + "s", lw = self._figLw, ls = "solid")
                lines = self.ax.get_lines()
                color = lines[-1].get_color()
                for i in range(len(x)):
                    self.ax.plot(x[i][key], y[i][key], lw = self._figLw, ls = ls[i], color=color)
        self.finishPlot()

