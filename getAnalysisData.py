#!/usr/bin/env python

########################################
########################################
# This file is part of Ph.D. work by Alexandra Strang on Adelie population dynamics in the Ross Sea
# Copyright (C) 2025 Alexandra Strang and Dean Anderson 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ADDITIONAL NOTES ABOUT THIS SCRIPT.
########################################
########################################

import os
import json
import numpy as np
from numba import njit
from osgeo import gdal, ogr, osr
import resource
from pathlib import Path

osr.UseExceptions()
gdal.UseExceptions()


class Params:
    def __init__(self):
        ## DATA FOLDERS
        self.alexandraDataPath = os.path.join(os.getenv('ADELIEPROJDIR', default = '.'), 
            'Alexandra_Data')
        self.deanInputDataPath = os.path.join(os.getenv('ADELIEPROJDIR', default = '.'), 
            'Dean_Data', 'Results_GuanoTerrain', 'Rasters_2m')
        self.deanOutputDataPath = os.path.join(os.getenv('ADELIEPROJDIR', default = '.'), 
            'Dean_Data', 'Results_GuanoTerrain', 'AnalysisData')
        # ## MAKE NEW RESULTS DIRECTORY IF DOESN'T EXIST
        if not os.path.isdir(self.deanOutputDataPath):
            os.makedirs(self.deanOutputDataPath)
        ## PATH TO TERRAIN DATA
        self.terrain_folder = os.path.join(self.alexandraDataPath, 'colony_terrain_rasters')
        self.terrainVariables = ['aspect', 'clipped', 'Roughness', 'slope', 'TRI']
        self.colonies = ['Crozier', 'Royds']
        ## OUTPUT DATA PATHS AND FILENAMES
        self.analysisDataFName = os.path.join(self.deanOutputDataPath, 'analysisData.json')

class ProcessRasters(object):
    def __init__(self, params):
        """
        CLASS OBJ TO PROCESS DATA
        """
        self.params = params

        #############################
        ## RUN FUNCTIONS
        self.sampleRasters()

        ## END RUNNING 
        #############################

    def sampleRasters(self):
        """
        ## LOOP THRU COLONIES AND ASSOCIATED DATA, GET DATA, WRITE JSON
        """
        self.dataDict = {}
        for col in self.params.colonies:
            ## GET COUNT DATA FIRST
            countPath = os.path.join(self.deanInputDataPath, '{}PenguinCounts_2m.tif'.format(col))
            
        


def main():
    params = Params()
    processor = ProcessRasters(params)

    maxMem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print('Max Mem Usage in KB', maxMem)

if __name__ == '__main__':
    main()
