# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *
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

import pyworkflow.utils as pwutils
from pyworkflow import VERSION_1_2
from pwem.objects import Volume, EMObject
from pwem.objects import AtomStruct
from pwem.emlib.image import ImageHandler
from pwem.convert import Ccp4Header
from ccp4 import Plugin
from ccp4.convert import (runCCP4Program, validVersion)
from pwem.protocols import EMProtocol
from pyworkflow.protocol.constants import STATUS_FINISHED
from pyworkflow.protocol.params import (MultiPointerParam, PointerParam,
                                        BooleanParam, StringParam)
from pyworkflow.utils.properties import Message
from ccp4.constants import CCP4_BINARIES
import sqlite3


# template for new atomic models, first ID is the coot model id
# the second id increases for each time the model is saved
COOTPDBTEMPLATEFILENAME = "coot_%06d_Imol_%04d_version_%04d.pdb" # protId, modelID, counter
# filename for coot script file
COOTSCRIPTFILENAME = "cootScript.py"
# database that stores filenames and corresponding models
OUTPUTDATABASENAMESWITHLABELS = "outpuDataBaseNameWithLabels.sqlite"
#table with the information
DATABASETABLENAME = 'pdb'

TYPE_3DMAP = 0
TYPE_ATOMSTRUCT = 1


class CootRefine(EMProtocol):
    """
    Interactive refinement and correction environment for macromolecular
    atomic models guided by cryo-EM density maps. The protocol enables
    users to inspect, modify, validate, and iteratively improve structural
    models within a graphical molecular interpretation workflow.

    AI Generated:

    Coot Refinement (CootRefine) - User Manual
        Overview

        The Coot Refinement protocol provides an interactive environment for
        manual and semi-automatic refinement of macromolecular atomic models
        against cryo-EM density maps. Its main purpose is to help users improve
        structural interpretation by visually inspecting the agreement between
        the experimental density and the current atomic model while allowing
        direct editing and local rebuilding of problematic regions.

        In cryo-EM workflows, automated refinement procedures often produce
        models that still require human inspection. Flexible loops, side-chain
        conformations, poorly resolved segments, and local sequence assignment
        errors frequently need expert intervention. This protocol integrates
        the Coot molecular graphics system into Scipion-oriented workflows so
        that model correction, rebuilding, and validation can be performed in
        an iterative and biologically meaningful manner.

        Biological Context and Typical Applications

        Structural refinement is not simply a geometric optimization process.
        The biological interpretation of a cryo-EM structure strongly depends
        on the quality of the fitted atomic model. Small tracing errors may
        alter conclusions regarding ligand binding, catalytic mechanisms,
        conformational states, or intermolecular interfaces.

        This protocol is especially useful after automated fitting or
        refinement stages, when the user needs to manually inspect uncertain
        regions. Typical applications include rebuilding flexible loops,
        correcting backbone tracing, improving side-chain placement, refining
        domain interfaces, adjusting chain assignments, or validating local
        geometry in medium- or high-resolution reconstructions.

        The protocol also supports comparative interpretation workflows in
        which multiple reference structures are loaded together with one or
        more cryo-EM maps. This allows users to compare conformations, inspect
        structural variability, or transfer biological interpretation from
        related structures.

        Input Maps and Atomic Structures

        The protocol accepts one or more cryo-EM density maps together with a
        principal atomic structure that will be refined interactively. Optional
        auxiliary atomic models can also be loaded as references. These
        reference structures are useful for comparative modeling, domain
        placement, conformational analysis, or rebuilding poorly resolved
        regions using homologous structures as guides.

        When several maps are provided, users can visually compare different
        reconstructions or focus refinement on maps generated under different
        processing strategies. This becomes particularly valuable when working
        with flexible systems, multi-state reconstructions, or local refinement
        approaches.

        Density Map Normalization

        The protocol optionally normalizes input maps before visualization.
        This step is important because molecular graphics tools often expect
        density values within stable intensity ranges. Proper normalization
        improves contour visualization, map interpretation, and local fitting
        behavior.

        In biological practice, normalization is especially useful when maps
        originate from heterogeneous processing pipelines or when density
        amplitudes vary significantly between reconstructions. Consistent map
        scaling facilitates visual comparison and reduces ambiguity during
        model rebuilding.

        Interactive Refinement Workflow

        The refinement process is intentionally interactive and user driven.
        Users inspect the density map, evaluate the local agreement between
        the map and the model, and manually apply corrections where necessary.
        This workflow is particularly important in cryo-EM because flexible or
        poorly resolved regions frequently require expert interpretation beyond
        what automated refinement algorithms can reliably provide.

        The protocol supports iterative correction cycles in which the model
        can be repeatedly edited, refined, and exported back into the Scipion
        workflow. This allows refinement to become part of a larger structural
        analysis strategy that may include validation, classification,
        comparison of conformational states, or additional automated refinement
        stages.

        Chain Navigation and Segment Refinement

        The protocol incorporates tools that facilitate rapid navigation across
        chains and residue ranges during interactive sessions. This is
        especially useful when refining large macromolecular assemblies,
        ribosomes, membrane complexes, or multi-chain protein systems.

        Segment-oriented refinement workflows allow users to focus on local
        regions of interest instead of refining the entire structure at once.
        From a biological perspective, local rebuilding is often preferable
        because many cryo-EM maps contain regions with variable local
        resolution. Stable cores may already be correct while flexible loops,
        interfaces, or peripheral domains require additional attention.

        The protocol also supports iterative traversal through amino acid
        segments, enabling users to systematically inspect and refine extended
        regions of the structure in a controlled and reproducible manner.

        Model Tracking and Data Management

        An important aspect of the protocol is the management of intermediate
        structural states generated during refinement. Multiple versions of the
        model can be produced throughout the interactive session, allowing
        users to preserve refinement history and recover previous structural
        states if necessary.

        This traceability is biologically important because refinement is often
        exploratory. Users may test alternative interpretations for ambiguous
        regions, compare different backbone traces, or evaluate multiple
        conformational hypotheses before selecting the final biologically
        meaningful model.

        Integration With Validation Workflows

        The protocol is designed to work together with downstream validation
        and refinement tools commonly used in structural biology workflows.
        Interactive rebuilding is frequently followed by stereochemical
        validation, real-space refinement, or comparison against independent
        reconstructions.

        In practical cryo-EM studies, users often alternate between automated
        refinement and manual correction several times before obtaining a final
        publication-quality model. This protocol acts as the human-guided
        component of that iterative refinement cycle.

        Outputs and Their Interpretation

        The protocol produces refined atomic structures that preserve the
        biological interpretation introduced during interactive editing.
        Normalized maps generated during preprocessing may also be exported for
        further visualization or downstream processing.

        Because the refinement process is interactive, the final structural
        quality depends strongly on the expertise of the user and on the local
        interpretability of the density map. Users should always visually
        inspect the agreement between density and model and complement manual
        refinement with independent validation tools.

        Practical Recommendations

        In routine cryo-EM practice, it is advisable to begin refinement by
        inspecting the global agreement between the model and the map before
        focusing on local details. Large tracing errors or domain misplacements
        should be corrected first because they can bias subsequent local
        refinement.

        Flexible loops and poorly resolved peripheral regions should be treated
        conservatively. Over-interpretation of weak density may introduce
        biologically misleading structural features. When uncertainty remains,
        preserving simpler or partially unresolved models is often preferable
        to introducing speculative coordinates.

        Comparative reference structures are particularly valuable when working
        at intermediate resolution. However, users should avoid forcing the map
        to match an expected conformation if the experimental density suggests
        genuine structural differences.

        Final Perspective

        Interactive refinement remains one of the most biologically important
        stages of cryo-EM structural interpretation. Although automated methods
        continue to improve, expert-guided rebuilding is still essential for
        resolving ambiguous regions, validating biological plausibility, and
        ensuring that the final structural model accurately reflects the
        experimental data. This protocol provides a flexible environment for
        integrating human expertise into iterative cryo-EM model refinement
        workflows.
    """
    _label = 'coot refinement'
    _program = ""
    _version = VERSION_1_2
    COOT = CCP4_BINARIES['COOT']
    COOTINI='coot.txt'
    EDITOR='editor.py'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolumes', MultiPointerParam, pointerClass="Volume",
                      label='Input Volume/s', allowsNull=True,
                      help="Set of volumes to process")
        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize', important=True,
                      help='If set to True, particles will be normalized in '
                           'the way COOT prefers it. It is recommended to '
                           '*always normalize your particles* if the maximum '
                           'value is higher than 1.')
        form.addParam('pdbFileToBeRefined', PointerParam,
                      pointerClass="AtomStruct",
                      label='Atomic structure to be refined',
                      help="PDBx/mmCIF file to be refined. This PDBx/mmCIF object "
                           "will be saved after refinement, will be saved")
        form.addParam('inputPdbFiles', MultiPointerParam,
                      pointerClass="AtomStruct",
                      allowsNull=True,
                      label='Other reference atomic structures',
                      help="Other PDBx/mmCIF files used as reference. These PDBx/mmCIF "
                           "objects will not be saved")
        form.addParam('extraCommands', StringParam,
                      default='',
                      condition='False',
                      label='Extra commands for chimera viewer',
                      help="""Add extra commands in cmd file. Use for testing
                      """)
        form.addParam('doInteractive', BooleanParam, default=True,
                      label='Interactive', condition='False',
                      help="""It makes coot an interactive protocol""")
        # TODO: when is this used, just for testing?
        form.addParam('phythonscript', StringParam, default="",
                      label='pythonScript', condition='False',
                      help="""calls coot with '--python string'""")
        form.addParam('inputProtocol', PointerParam, allowsNull=True,
                      default=None,
                      condition='False',
                      label="Input protocols", important=True,
                      pointerClass='PhenixProtRunMolprobity, '
                               'PhenixProtRunRSRefine',
                      help="Father protocol. This is used for trazability "
                       "when coot is launched by a viewer ")

        form.addSection(label='Help')
        form.addLine('Press "w" in coot to transfer the pdb file from coot '
                     'to scipion.\nYou may also execute (Calculate -> '
                     'Scripting -> Python) the command scipion_write(imol, '
                     '[pdblabel]).\nExample: scipion_write(0,"new_name")\n'
                     'imol is the PDB id.\npdblabel is an optional parameter '
                     'that assign that label to the produced pdb. By default '
                     'the label is outcoot0001\n'
                     'Press "x" in coot to change from one chain to '
                     'the previous one.\nPress "X" in coot to change from one '
                     'chain to the next one.\nPress "U" in coot to initiate '
                     'global variables.\nYou have to set in advance the '
                     'protocolDirectory/extra/coot.ini text file:\n['
                     'myvars]\nimol: '
                     '0\naa_main_chain: '
                     'X\naa_auxiliary_chain: XX\naaNumber: 160\nstep: 15\nIn '
                     'this case global variables will initiate in '
                     'aminoacid number 160\nand each shift over the '
                     'sequence will include a segment of 15 aminoacids.\n'
                     'Press "z" in coot to refine those upstream 15 '
                     'aminoacids included in each step.\nPress "Z" in coot ' \
                     'to refine those downstream 15 aminoacids included in ' \
                     'each step.\nPress "E" in coot to print the ' \
                     'environment.\nPress "e" in coot to finish your ' \
                     'project. Then your project will not be interactive ' \
                     'anymore.')

        # --------------------------- INSERT steps functions ---------------

    def _insertAllSteps(self):
        convertId = self._insertFunctionStep('convertInputAndSaveToDBStep')
        self.step = self._insertFunctionStep('runCootStep',
                                             prerequisites=[convertId],
                                             interactive=self.doInteractive)


    # --------------------------- STEPS functions --------------------------

    def convertInputAndSaveToDBStep(self):
        """ init database to store the last name of the files used and
            convert 3D maps to MRC '.mrc' format.
            This step  is run once even if the protocol is relunched
        """
        databasePath = self._getExtraPath(OUTPUTDATABASENAMESWITHLABELS)

        # create database and table
        # this table will be used to record the last version of any file
        # save in coot
        conn = sqlite3.connect(databasePath)

        # create table
        # saved = 0, means this file need to be converted to a scipion object
        # type  = 0-> Map, 1 -> atom struct
        # predefined macros for type:
        # TYPE_3DMAP = 0
        # TYPE_ATOMSTRUCT = 1

        sqlCommand = """create table if not exists %s 
                               (id integer primary key AUTOINCREMENT,
                                modelId integer,
                                fileName text,
                                labelName text,
                                type int,
                                saved integer default 1
                                )""" % (DATABASETABLENAME)

        conn.execute(sqlCommand)

        # create view to retrieve the id of the last copy
        # for each model id  (coot call imol to this id)
        sqlCommand = """CREATE VIEW lastid AS 
                        SELECT modelId, max(id) as id
                        FROM %s
                        GROUP BY modelId
        """ % DATABASETABLENAME

        conn.execute(sqlCommand)

        inVolumes, norVolumesNames = self._getVolumesList()

        sqlCommand = """INSERT INTO %s (modelId, fileName, labelName, type, saved) values
                        (%d, '%s', '%s', %d, %d)"""

        #process main atomic Structure
        counter=0
        pdbFileToBeRefined = self.pdbFileToBeRefined.get().getFileName()
        base = os.path.basename(pdbFileToBeRefined)
        conn.execute(sqlCommand % (DATABASETABLENAME, counter, pdbFileToBeRefined,
                                   os.path.splitext(base)[0], TYPE_ATOMSTRUCT, 1 # saved
                                   )
                     )
        counter += 1

        # Process another atom structures
        for pdb in self.inputPdbFiles:
            fileName = pdb.get().getFileName()
            base = os.path.basename(fileName)
            conn.execute(sqlCommand % (DATABASETABLENAME, counter,
                                       fileName,
                                       os.path.splitext(base)[0], TYPE_ATOMSTRUCT,
                                       1 # saved
                                       )
                         )
            counter += 1

        # Process 3D maps
        # normalize them if needed
        ih = ImageHandler()
        for inVol, norVolName in zip(inVolumes, norVolumesNames):
            inVolName = inVol.getFileName()

            if inVolName.endswith(".mrc"):
                inVolName += ":mrc"
            if norVolName.endswith(".mrc"):
                norVolName += ":mrc"
            if not ih.existsLocation(norVolName):
                if True:  # self.doNormalize:
                    img = ImageHandler()._img
                    img.read(inVolName)
                    mean, dev, min, max = img.computeStats()
                    img.inplaceMultiply(1./max)
                    img.write(norVolName)
                else:
                    ImageHandler().convert(inVolName, norVolName)
                Ccp4Header(norVolName, readHeader=True).copyCCP4Header(
                    inVol.getOrigin(force=True).getShifts(),
                    inVol.getSamplingRate(), originField=Ccp4Header.START)

            conn.execute(sqlCommand%(DATABASETABLENAME, counter, norVolName[:-4],
                                     os.path.basename(norVolName)[:-8], TYPE_3DMAP, 0 # saved
                                     )
                         )
            counter += 1

        conn.commit()

    def runCootStep(self):

        databasePath = self._getExtraPath(OUTPUTDATABASENAMESWITHLABELS)
        createScriptFile(0,  # imol
                         self._getExtraPath(COOTSCRIPTFILENAME),  # save script in extra otherwise is lost
                         # when continue
                         self._getExtraPath(COOTPDBTEMPLATEFILENAME), # default template name for
                                                                      # fro newPDBs
                         self.extraCommands.get(),
                         self._getExtraPath(self.COOTINI),  # coot.ini
                         self._getExtraPath(self.EDITOR),   # editor.py
                         databasePath,
                         table_name=DATABASETABLENAME,
                         protId=self.getObjId()
                         )

        args = ""

        #  extraCommands option is only used for tests
        if self.extraCommands.get() != '':
            args += " --no-graphics "
        args += " --script " + self._getExtraPath(COOTSCRIPTFILENAME)
        if len(self.phythonscript.get()) > 1:
            args += " --script {phythonscript}".format(
                phythonscript=self.phythonscript.get())
        # script with auxiliary files
        self._log.info('Launching: ' + Plugin.getProgram(self.COOT) + ' ' + args)

        # run in the background
        try:
            runCCP4Program(Plugin.getProgram(self.COOT), args)
        except Exception as e:
            print(pwutils.redStr("=============================================================="))
            print(pwutils.redStr("ERROR: Coot execution failed"))
            print(pwutils.redStr("Last SAVED atomic model will be shown if protocol is continued"))
            print(pwutils.redStr("=============================================================="))
 
        self.createOutput()

    def createOutput(self):
        """ Copy the PDB structure and register the output object.
        """
        databasePath = self._getExtraPath(OUTPUTDATABASENAMESWITHLABELS)
        getModels(databasePath, DATABASETABLENAME)

        # open database
        conn = sqlite3.connect(databasePath)
        if not _checkTableExists(conn, DATABASETABLENAME):
            conn.close()
            return

        c = conn.cursor()

        # read atom struct filename and label in a loop
        c.execute('SELECT fileName, labelName '
                  'FROM %s '
                  'WHERE saved = 0 AND type=%d' % (DATABASETABLENAME, TYPE_ATOMSTRUCT))
        for row in c:
            pdbFileName = row[0]
            pdbLabelName = row[1]
            pdb = AtomStruct()
            pdb.setFileName(pdbFileName)

            outputs = {str(pdbLabelName) : pdb}
            self._defineOutputs(**outputs)
            self._defineSourceRelation(self.inputPdbFiles, pdb)

        # files has been saved
        sql = 'UPDATE %s SET saved = 1 WHERE saved=0 ' \
              'AND type=%d' % (DATABASETABLENAME, TYPE_ATOMSTRUCT)
        c.execute(sql)

        # check if normalized files are saved
        sql = "SELECT count(*) " \
              "FROM %s " \
              "WHERE saved = 0 " \
              "  AND type = %d LIMIT 1" % (DATABASETABLENAME, TYPE_3DMAP)
        c.execute(sql)
        result = c.fetchone()[0]

        # update  saved parameter
        if result>0:
            sql = 'UPDATE %s SET saved = 1 WHERE saved=0 ' \
                  'AND type=%d' % (DATABASETABLENAME, TYPE_3DMAP)
            c.execute(sql)
        conn.commit()
        conn.close()

        # save normalized  vols...
        if result>0:
            inVolumes, norVolumesNames = self._getVolumesList()
            counter = 1
            for inVol, norVolName in zip(inVolumes, norVolumesNames):
                outVol = Volume()
                sampling = inVol.getSamplingRate()
                origin = inVol.getOrigin(
                    force=True)
                outVol.setSamplingRate(sampling)
                outVol.setOrigin(origin)

                if norVolName.endswith('.mrc'):
                    norVolName = norVolName + ":mrc"
                outFileName = self._getVolumeFileName(norVolName)
                outVol.setFileName(outFileName)
                outputs = {"output3DMap_%04d" % counter: outVol}
                counter += 1
                self._defineOutputs(**outputs)
                self._defineSourceRelation(inVol, outVol)
        else:
            print("skip save normalized vol")

        if os.path.isfile(self._getExtraPath('STOPPROTCOL')):
            self.setStatus(STATUS_FINISHED)
            # NOTE: (ROB) can a dirty way to make an interactive process finish but I do not
            # think there is a clean one
            self._steps[self.step-1].setInteractive(False)

    # --------------------------- INFO functions ---------------------------
    def _validate(self):
        errors = []

        if not validVersion(7, 0.056):
            errors.append("CCP4 version should be at least 7.0.056")

        if self.inputProtocol.get() is not None and \
                self.inputProtocol.get().getClassName().startswith("PhenixProtRunMolprobity"):
                 return errors
        else:
            # # Check that the input volume exist
            # if self.pdbFileToBeRefined.hasValue():
            #     if (not self.pdbFileToBeRefined.get().hasVolume()) \
            #             and self.inputVolumes.isEmpty():
            #         errors.append("Error: You should provide a volume.\n")

            return errors

    @classmethod
    def validateInstallation(cls):

        # Check that the programs exist
        installed, message = Plugin.checkBinaries(cls.COOT)
        if not installed:
            return [message]
        else:
            return []

    def _summary(self):
        #  Think on how to update this summary with created PDB
        summary = []
        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes(EMObject):
                summary.append("*%s:* \n %s " % (key, output.getObjComment()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")

        return methodsMsgs

    def _citations(self):
        return ['Emsley_2004']

    # --------------------------- UTILS functions --------------------------

    def _getVolumesList(self):
        print("_getVolumesList", flush=True)
        # test loop over inputVol
        inVolumes = []
        norVolumesNames = []
        # if self.inputVolumes is None:
        if len(self.inputVolumes) == 0:
            if self.pdbFileToBeRefined.get().getVolume() is not None:
                vol = self.pdbFileToBeRefined.get().getVolume()
                inFileName = vol.getFileName()
                inVolumes.append(vol)
                norVolumesNames.append(self._getVolumeFileName(inFileName))
        else:
            for vol in self.inputVolumes:
                inFileName = vol.get().getFileName()
                inVolumes.append(vol.get())
                norVolumesNames.append(
                    self._getVolumeFileName(inFileName))

        return inVolumes, norVolumesNames

    def _getVolumeFileName(self, inFileName):
        return os.path.join(self._getExtraPath(''),
                            pwutils.replaceBaseExt(inFileName, 'mrc'))

    def replace_at_index(self, tup, ix, val):
        return tup[:ix] + (val,) + tup[ix+1:]

    def getCounter(self):
        template = self._getExtraPath(COOTPDBTEMPLATEFILENAME)
        counter = 1
        while os.path.isfile(template % counter):
            counter += 1
        return counter  # returns next free

cootScriptHeader = '''import ConfigParser
import os
import subprocess
try:
    import coot_python
    has_gui = True
except:
    print("ERROR: coot_python module not found. If you are running in no gui mode this is OK")
    has_gui = False
mydict={{}}
mydict['imol']={imol}
mydict['aa_main_chain']="A"
mydict['aa_auxiliary_chain']="AA"
mydict['aaNumber']=17
mydict['step']=5
mydict['outfile'] = '{templateNameAtomStruct}'
cootPath='{cootFileName}'
editorPath='{editorFileName}'
databasePath='{outpuDataBaseNameWithLabels}'
table_name = '{table_name}'
TYPE_3DMAP = {TYPE_3DMAP}
TYPE_ATOMSTRUCT = {TYPE_ATOMSTRUCT}
protId={protId}
'''

cootScriptBody = '''

def beep(time):
   """I simply do not know how to create a portable beep sound.
      This system call seems to work pretty well if you have sox
      installed"""
   try:
      command = "play --no-show-progress -n synth %f sin 880"%time
      # print command
      os.system(command)
   except:
      pass

def _change_chain_id(signStep):
    """move a few aminoacid between chains"""
    global mydict
    dic = dict(mydict)
    if signStep < 0:
        dic['fromAaNumber'] = mydict['aaNumber'] - dic['step'] +1
        dic['toAaNumber']   = mydict['aaNumber']
        dic['fromAaChain']  = mydict['aa_auxiliary_chain']
        dic['toAaChain']    = mydict['aa_main_chain']
    else:
        dic['fromAaNumber'] = mydict['aaNumber']
        dic['toAaNumber']   = mydict['aaNumber'] + dic['step'] -1
        dic['fromAaChain']  = mydict['aa_main_chain']
        dic['toAaChain']    = mydict['aa_auxiliary_chain']
    mydict['aaNumber'] = mydict['aaNumber'] + (dic['step'] * signStep)
    command = "change_chain_id(%(imol)d, '%(fromAaChain)s', '%(toAaChain)s', 1, %(fromAaNumber)d, %(toAaNumber)d)"%dic

    doIt(command)

def _refine_zone(signStep):
    """Execute the refine command"""
    global  mydict
    dic = dict(mydict)
    if signStep <0:
        dic['fromAaNumber'] = mydict['aaNumber'] - dic['step']
        dic['toAaNumber']   = mydict['aaNumber'] + 2
        mydict['aaNumber']  = mydict['aaNumber'] - dic['step']
    else:
        dic['fromAaNumber'] = mydict['aaNumber'] - 2
        dic['toAaNumber']   = mydict['aaNumber'] + dic['step']
        mydict['aaNumber']  = mydict['aaNumber'] + dic['step']
    command = 'refine_zone(%(imol)s, "%(aa_main_chain)s", %(fromAaNumber)d, %(toAaNumber)d, "")'%dic

    doIt(command)

def _updateMol():
    """update global variable using a file as
    [myvars]
    imol: 0
    aa_main_chain: A
    aa_auxiliary_chain: AA
    aaNumber: 82
    step: 15
    called protocolDirectory/extra/coot.ini"""

    def open_xdg(my_file):
        """open text file with default editor. I tried xdg-open but it
        does not work.
        
        """
        # unset are need otherwise python ccp4 and python scipion are mixed
        os.system("unset PYTHONPATH; unset PYTHONHOME; python %s %s" % (editorPath, my_file))
    open_xdg(os.environ.get('COOT_INI', cootPath))
    global mydict
    config = ConfigParser.ConfigParser()
    config.read(os.environ.get('COOT_INI', cootPath))
    try:
        mydict['imol']               = int(config.get("myvars", "imol"))
        mydict['aa_main_chain']      = config.get("myvars", "aa_main_chain")
        mydict['aa_auxiliary_chain'] = config.get("myvars","aa_auxiliary_chain")
        mydict['aaNumber']           = int(config.get("myvars", "aaNumber"))
        mydict['step']               = int(config.get("myvars", "step"))
        mydict['outfile']            = config.get("myvars", "outfile")
    except ConfigParser.NoOptionError:
        pass
    add_status_bar_text("Global variable updated using coot.ini")
    beep(0.1)


def getOutPutFileName(template, imol):
    """get name based on template that does not exists
    %04d will be incremented untill it does not exists"""
    counter=1
    if "%04d" in template:
        while os.path.isfile(template%(protId, imol, counter)):
             counter += 1

    return template % (protId, imol, counter)

def storeFileNameDataBase(imol, outFileName, outLabel=None, type=TYPE_ATOMSTRUCT):
    import sqlite3
    conn = sqlite3.connect(databasePath)
    c = conn.cursor()
    # create_database if it does not exists
    # in new version database is always created
    # sqlCommand = """create table if not exists %s 
    #                       (id integer primary key AUTOINCREMENT,
    #                        modelId integer,
    #                        fileName text,
    #                        labelName text,
    #                        type int,
    #                        saved integer default 1
    #                        )""" % (DATABASETABLENAME)
    #
    # conn.execute(sqlCommand)
    
    # insert record
    if outLabel is None:
        outLabel = os.path.splitext(os.path.basename(outFileName))[0]

    # sqlCommand = """INSERT INTO %s (modelId, fileName, labelName, type, saved) values
    #                     (%d, '%s', '%s', %d)"""
    saved = 0 # saved = 0 -> This file has not been aaded to scipion
    sql = 'insert into ' + table_name + """ (modelId, fileName, labelName, type, saved) values
                                            (%d, '%s', '%s', %d, %d)""" % (imol, outFileName, outLabel, type, saved)
    c.execute(sql)

    # commit
    conn.commit()
    # close connection
    conn.close()

def _write(imol=-1, outLabel=None):
    """write pdb file, default names
       can be overwritted using coot.ini"""
    global mydict
    aa_imol = imol
    if imol == -1:
        with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
            print("saving last active molecule with #ID", aa_imol)

    dic = dict(mydict)
    dic['imol'] = aa_imol
    outFileName=getOutPutFileName(dic['outfile'], aa_imol)

    if outLabel is None:
        outLabel = os.path.splitext(os.path.basename(outFileName))[0]
    # else:
    #    ext = os.path.splitext(outFileName)[1]
    #    dir = os.path.dirname(outFileName)
    #    basename = os.path.splitext(outLabel)[0]
    #    outFileName = os.path.join(dir, basename + ext)
    
    dic['outfile'] = outFileName
    save_coordinates(dic['imol'], outFileName)
    
    if os.path.isfile(outFileName):
        type = TYPE_ATOMSTRUCT
        storeFileNameDataBase(imol, outFileName, outLabel, type)
        add_status_bar_text("Saved imol: %(imol)s as %(outfile)s" % dic)
    else:
        add_status_bar_text("I do not know how to export a 3D map. File NOT saved.")
        dic['outfile'] = outFileName.replace(".pdb", ".mrc")
        command = "export_map(%(imol)s,'%(outfile)s')" % dic
        # TODO: the file is saved but it is not
        # clear how to handle it
        # No Scipion object will be created
        result = doIt(command)
        type = TYPE_3DMAP
        
    #store information     
    beep(0.1)

def scipion_write(imol=0, outLabel=None):
    """scipion utility for writting files
    args: model number, 0 by default.
    This function is use for python scripts"""
    global mydict
    mydict['imol']=imol
    _write(imol, outLabel)

def doIt(command):
    """launch command"""
    return eval(command)
    #beep(0.1)

def _printEnv():
    for key in os.environ.keys():
       print "%30s %s \\n" % (key,os.environ[key])

def _finishProj():
    global mydict
    filenName = mydict['outfile']%(1,1,1)
    dirPath = os.path.dirname(filenName)
    fileName = os.path.join(dirPath,"STOPPROTCOL")
    open(fileName,"w").close()
    beep(0.1)
    coot_real_exit(0)


# create scipion menu
print("Loading Scipion Coot extensions...")

if has_gui: 
    menubar = coot_python.main_menubar()
    toolbar = coot_python.main_toolbar()
    menu = coot_menubar_menu("Scipion")
    add_simple_coot_menu_menuitem(menu, "write last active model", lambda func: _write(imol = -1))
    add_simple_coot_menu_menuitem(menu, "end protocol", lambda func: _finishProj())
    add_simple_coot_menu_menuitem(menu, "update cootini", lambda func: _updateMol())

#change chain id
add_key_binding("change_chain_id_down","x", lambda: _change_chain_id(-1))
add_key_binding("change_chain_id_down","X", lambda: _change_chain_id(1))

#refine aminoacid segment
add_key_binding("refine zone m","z", lambda: _refine_zone(1))
add_key_binding("refine zone m","Z", lambda: _refine_zone(-1))

#update global variables
add_key_binding("init global variables","U", lambda: _updateMol())

#write file
add_key_binding("write pdb file","w", lambda: _write())

#print environ
add_key_binding("print enviroment","E", lambda: _printEnv())

#finish project
add_key_binding("finish project","e", lambda: _finishProj())

'''
def getModels(outpuDataBaseNameWithLabels, table_name):
    # open database
    conn = sqlite3.connect(outpuDataBaseNameWithLabels)
    if not _checkTableExists(conn, table_name):
        conn.close()
        return

    c = conn.cursor()

    # read filename and label in a loop
    # get id of the last copy per each file
    sqlCommand = """SELECT fileName 
                    FROM %s  NATURAL JOIN lastid
                    WHERE type = %d""" % \
                    (table_name, TYPE_3DMAP)
    c.execute(sqlCommand)

    listOfMaps = []
    for row in c:
        listOfMaps.append(row[0])

    c.execute("""SELECT fileName 
                 FROM %s  NATURAL JOIN lastid
                 WHERE type = %d """ %
              (table_name, TYPE_ATOMSTRUCT))

    listOfAtomStructs = []
    for row in c:
        listOfAtomStructs.append(row[0])
    return listOfMaps, listOfAtomStructs

def createScriptFile(imol,  # problem PDB id
                     scriptFile,  # name of the coot script file
                     templateNameAtomStruct,  # default template name for new files
                     extraCommands='',  # extra commands to add at the
                                       # end of the file
                                       # mainly used for testing
                     cootFileName='/tmp/coot.ini',
                     editorFileName='/tmp/editor.py',
                     outpuDataBaseNameWithLabels='output.db',
                     table_name='pdb',
                     protId=0
                     ):

    listOfMaps, listOfAtomStructs = getModels(outpuDataBaseNameWithLabels,
                                              table_name)

    f = open(scriptFile, "w")
    d = {'imol':imol,
         'templateNameAtomStruct':templateNameAtomStruct,
         'cootFileName':cootFileName,
         'editorFileName':editorFileName,
         'outpuDataBaseNameWithLabels':outpuDataBaseNameWithLabels,
         'table_name':table_name,
         'TYPE_3DMAP':TYPE_3DMAP,
         'TYPE_ATOMSTRUCT':TYPE_ATOMSTRUCT,
         'protId':protId}

    f.write(cootScriptHeader.format(**d))
    f.write(cootScriptBody)

    # load PDB and MAP
    imol_counter = 0
    f.write("\n#load Atomic Structures\n")  # problem atomic structure must be
    for pdb in listOfAtomStructs:
        f.write("read_pdb('%s')\n" % pdb) #
        f.write("set_mol_active(%d, %d)\n" % (imol_counter, imol_counter < 1))
        imol_counter += 1

    f.write("\n#load 3D maps\n")
    f.write("map_colour = (0.0, 0.5, 1.0)\n")
    for vol in listOfMaps:
        f.write("handle_read_ccp4_map('%s', 0)\n" % vol)
        #imol_counter += 1
    # set color of first map
    f.write(f"set_map_colour({imol_counter}, *map_colour)\n")

    f.write("\n#Extra Commands\n")
    f.write(extraCommands)
    f.close()


    # create coot.ini if it does not exist
    if os.path.exists(cootFileName):
        pass
    else:
        f = open(cootFileName,"w")
        f.write("""[myvars]
imol: 0
aa_main_chain: A
aa_auxiliary_chain: AA
aaNumber: 100
step: 10
""")
        f.close()
    # create editor if it does not exist
    if os.path.exists(editorFileName):
        pass
    else:
        f = open(editorFileName,"w")
        f.write("""#
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import sys
    
# save current file
def save_file():
    # try to get current file path
    try:
        # get file path from window title
        path = main_window.title().split('-')[1][1:]
    
    # init path variable to empty string on exception
    except:
        print('Unexpected error:', sys.exc_info()[0])

    # check if file path is available
    if path != '':
        # write file
        with open(path, 'w') as f:
            content = text_area.get('1.0', tk.END)
            f.write(content)
                
    # clear 'is_modified' flag
    text_area.edit_modified(0)


# create main window instance
main_window = tk.Tk()

# configure main window
main_window.title('Notepad')
main_window.geometry('800x600')

# create menu bar instance
menubar = tk.Menu(main_window)

# create 'File' menu items
file_menu = tk.Menu(menubar, tearoff=0)
file_menu.add_command(label="Save", command=save_file)
file_menu.add_command(label="Exit", command=main_window.quit)

# add 'File' menu to the menu bar
menubar.add_cascade(label="File", menu=file_menu)

# create text area to input text
text_area = tk.Text(main_window)
text_area.pack(expand = tk.YES, fill = tk.BOTH, side = tk.LEFT)

# create scrollbar and link it to the text area
scroll_bar = ttk.Scrollbar(main_window, orient=tk.VERTICAL, command=text_area.yview)
scroll_bar.pack(fill=tk.Y, side=tk.RIGHT)
text_area['yscrollcommand'] = scroll_bar.set

# Connect menubar to the window
main_window.config(menu=menubar)

# run main application loop
path = sys.argv[1]
main_window.title('Notepad - ' + path)
with open(path, 'r') as f:
    # clear text area and insert file content
    content = f.read()
    text_area.delete('1.0', tk.END)
    text_area.insert('1.0', content)
    # clear 'is_modified' flag
    text_area.edit_modified(0)

tk.mainloop()
""")
        f.close()
def _checkTableExists(dbcon, tablename):
    dbcur = dbcon.cursor()
    dbcur.execute("""
        SELECT COUNT(*)
        FROM sqlite_master
        WHERE type='table'
         AND name='%s'
        """%tablename)
    if dbcur.fetchone()[0] == 1:
        dbcur.close()
        return True

    dbcur.close()
    return False
