===================
CCP4 scipion plugin
===================

This plugin allows to use CCP4 programs within the Scipion framework. **You need to install the CCP4 suite before installing the plugin**, see section "Binary Files" for details.

CCP4, from Collaborative Computational Project Number 4, is a software suite that allows model building of macromolecule structures obtained by X-ray crystallography, and that has been extended to other techniques like cryo-EM (see `CCP4 home page <http://www.ccp4.ac.uk/>`_ for details).

Programs from CCP4 included in the Scipion framework for model building:

  * coot
  * refmac

===================
Install this plugin
===================

You will need to use `2.0.0 <https://github.com/I2PC/scipion/releases/tag/v2.0>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-em-ccp4
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-ccp4.git

2. Install:

.. code-block::

           scipion installp -p path_to_scipion-em-ccp4 --devel



- **Binary files** 

CCP4 binaries will *NOT* be installed automatically with the plugin. The independent installation of CCP4 software suite by the user is required before running the programs. Default installation path assumed is */usr/local/ccp4-7.0*; this path or any other of your preference has to be set in *CCP4_HOME* in *scipion.conf*. We recommend to install CCP4 version 7.0.056 or higher. (see http://www.ccp4.ac.uk/download/#os=linux)



- **Tests**

Tested with CCP4 versions: 7.0.056 and 7.0.066.

To check the installation, simply run the following Scipion test: 

* scipion test ccp4.tests.test_protocol_coot_refmac



- **Supported versions of CCP4**

7.0.056 or higher.

- **Additional Instruction**

see https://github.com/scipion-em/scipion-em-ccp4/wiki

=========
Protocols
=========

* coot refinement: Molecular interactive graphics application used for flexible fitting, refinement, model completion, and validation of structures of macromolecules regarding electron density maps. See the `details <https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/>`_ of *Coot* utilities. 
* refmac: Automatic refinement program in Fourier space of macromolecule structures regarding electron density maps. See ` <http://www.ccp4.ac.uk/html/refmac5/description.html>`_ of *Refmac* utilities.




========
Examples
========

See `Model Building Tutorial <https://github.com/I2PC/scipion/wiki/tutorials/tutorial_model_building_basic.pdf>`_




===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/ccp4_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/ccp4_prod.svg

