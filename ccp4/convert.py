# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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
"""
This module contains converter functions that will serve to:
1. define ccp4 environ
TODO:
2. Read/Write CCP4 specific files
"""

import os
import pyworkflow.utils as pwutils
from ccp4 import Plugin


def runCCP4Program(program, args="", extraEnvDict=None, cwd=None):
    """ Internal shortcut function to launch a CCP4 program. """
    env = Plugin.getEnviron()
    if extraEnvDict is not None:
        env.update(extraEnvDict)
    pwutils.runJob(None, program, args, env=env, cwd=cwd)


def validVersion(major=7, minor=0.056, greater=True):
    """ Return ccp4 version as string. Example: 7.0.056"""

    versionFileName = os.path.join(Plugin.getHome(), 'lib',
                                   'ccp4','MAJOR_MINOR')

    if not os.path.exists(versionFileName):
        return False

    with open(versionFileName,"r") as f:
        _major, _minor = f.readline().split(".",1)
        _major = int(_major)
        _minor = float(_minor)
        if greater:
            if _major > major or \
                (_major == major and _minor >= minor):
                return True
        else:
            if _major == major and _minor == minor:
                return True
    return False
