# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es), May 2013
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import sqlite3

from pwem import Domain
from pwem.emlib.image import ImageHandler
from pwem.viewers.viewer_chimera import Chimera
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from ccp4.protocols.protocol_coot import (CootRefine, COOTPDBTEMPLATEFILENAME,
                                          OUTPUTDATABASENAMESWITHLABELS,
                                          DATABASETABLENAME)

# TODO: very likely this should inherit from ProtocolViewer
# not from XmippViewer. But then I get an empty form :-(


class CootRefineViewer(Viewer):
    """ Visualize the output of the coot protocol """
    _label = 'coot viewer'
    _targets = [CootRefine]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **args):
            # TODO if input volume is not mrc this will not work.
        # Construct the coordinate file and visualization
        bildFileName = os.path.abspath(self.protocol._getExtraPath("axis.bild"))
        dims = []
        samplings = []
        if len(self.protocol.inputVolumes) == 0:
            if self.protocol.pdbFileToBeRefined.get().getVolume() is not None:
                dim = self.protocol.pdbFileToBeRefined.get().getVolume().getDim()[0]
                sampling = self.protocol.pdbFileToBeRefined.get().getVolume().\
                    getSamplingRate()
                dims.append(dim)
                samplings.append(sampling)
        else:
            for i in range(len(self.protocol.inputVolumes)):
                dim = self.protocol.inputVolumes[i].get().getDim()[0]
                sampling = self.protocol.inputVolumes[i].get().\
                    getSamplingRate()
                dims.append(dim)
                samplings.append(sampling)
        if len(dims) != 0 and len(samplings) != 0:
            Chimera.createCoordinateAxisFile(max(dims),
                                     bildFileName=bildFileName,
                                     sampling=max(samplings))
        else:
            dim = 150.
            sampling = 1.
            Chimera.createCoordinateAxisFile(dim,
                                             bildFileName=bildFileName,
                                             sampling=sampling)

        fnCmd = self.protocol._getExtraPath("chimera.cxc")
        f = open(fnCmd, 'w')
        f.write("open %s\n" % bildFileName)
        f.write("cofr 0,0,0\n")

        outputsVol = []
        if len(self.protocol.inputVolumes) == 0:
            if self.protocol.pdbFileToBeRefined.get().getVolume() is not None:
                outputVol = self.protocol.pdbFileToBeRefined.get().getVolume()
                outputsVol.append(outputVol)
        else:
            for i in range(len(self.protocol.inputVolumes)):
                outputVol = self.protocol.inputVolumes[i].get()
                outputsVol.append(outputVol)

        count = 2
        if len(outputsVol) != 0:
            for outputVol in outputsVol:
                outputVolFileName = os.path.abspath(
                        ImageHandler.removeFileType(outputVol.getFileName()))
                x, y, z = outputVol.getOrigin(force=True).getShifts()
                f.write("open %s\n" % outputVolFileName)
                f.write("volume #%d  style surface voxelSize %f\n"
                        "volume #%d  origin %0.2f,%0.2f,%0.2f\n"
                        % (count, outputVol.getSamplingRate(), count, x, y, z))
                count += 1

        # counter = 1
        # template = self.protocol._getExtraPath(COOTPDBTEMPLATEFILENAME)
        databasePath = self.protocol._getExtraPath(OUTPUTDATABASENAMESWITHLABELS)
        conn = sqlite3.connect(databasePath)
        c = conn.cursor()
        sql = 'SELECT fileName FROM %s where saved = 1 order by id' % \
              DATABASETABLENAME
        c.execute(sql)
        for row in c:
            pdbFileName = os.path.abspath(row[0])
            if not pdbFileName.endswith(".mrc"):
                f.write("open %s\n" % pdbFileName)

        f.close()
        conn.close()
        # run in the background
        chimeraPlugin = Domain.importFromPlugin('chimera', 'Plugin', doRaise=True)
        chimeraPlugin.runChimeraProgram(chimeraPlugin.getProgram(), fnCmd + "&")
        return []
