#############################################################################################################
# Global Protein Surface Survey (GPSS) PyMOL Plugin
#
# Author: T.A. Binkowski (abinkowski@anl.gov)
# Date:   July 20, 2006
#
#
# Description: The GPSS PyMOL allows for the visualization of precomputed
#              surface features from protein structures in the Protein
#              Data Bank.  The plugin was developed from a suggestion by
#              Sebestien Moretti.
#
# Acknowledgements: The following publically available sources provide data
#                   for the plugin:
#                    * Computed Atlas of Surface Topography of Proteins (CASTp, http://cast.engr.uic.edu)
#                    * Catalytic Site Atlas (CSA, http://www.ebi.ac.uk/thornton-srv/databases/CSA/)
#                    
#
# Copyright Information:  The Global Protein Surface Survey (GPSS) Pymol
#                         Plugin source code is Copyright (c) 2006 by
#                         T. Andrew Binkowski (abinkowski@anl.gov).
#
# All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
#
#
#
# More Information: Information and documentation can be found
#                   at http://gpss.mcsg.anl.gov.
# 
#
#############################################################################################################

import tkSimpleDialog
import tkMessageBox
import tkFileDialog
from pymol import cmd
import sys, urllib, urllib2, zlib, os

def __init__(self):
    self.menuBar.addcascademenu('Plugin',
                                'GPSS',
                                'GPSS',
                                label = 'Global Protein Surface Survey')
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'GPSS All',
                             label = 'GPSS Sufaces (All)',
                             command = lambda s=self : gpssAllDialog(s))
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'CASTp',
                             label = 'CASTp Surfaces',
                             command = lambda s=self : gpssCastpDialog(s))
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'Ligand',
                             label = 'Ligand Binding Surfaces',
                             command = lambda s=self : gpssLigandDialog(s))
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'Metal',
                             label = 'Metal Binding Surfaces',
                             command = lambda s=self : gpssMetalDialog(s))
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'Peptide',
                             label = 'Peptide Binding Surfaces',
                             command = lambda s=self : gpssPeptideDialog(s))
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'DNA',
                             label = 'DNA Binding Surfaces',
                             command = lambda s=self : gpssDnaDialog(s))
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'Sites',
                             label = 'Functional Sites',
                             command = lambda s=self : gpssSitesDialog(s))
    self.menuBar.addmenuitem('GPSS','separator')
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'File Upload',
                             label = 'File Upload',
                             command = lambda s=self : gpssFileUploadDialog(s))
    self.menuBar.addmenuitem('GPSS','separator')
    self.menuBar.addmenuitem('GPSS',
                             'command',
                             'Help',
                             label = 'Help',
                             command = lambda s=self : gpssHelpDialog(s))

    
    ########################################################################
    # gpssVersionCheck
    #
    ########################################################################
    def _gpssVersionCheck(version):
        versionCheck = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                       'gpssServerVersionCheck.php').read()        
        if versionCheck!=version:
            tkMessageBox.showinfo('GPSS - Help                         ',
                                  'There is a newer version of the GPSS Plugin available at'+
                                  ' http://gpss.mcsg.anl.gov.');
            
    ########################################################################
    # _gpssLoadPdb
    #
    ########################################################################
    def _gpssLoadPdb(pdbCode,reinitialize):

        pdbCode = pdbCode.lower()

        _gpssVersionCheck('0.3');

        if cmd.get_names().count('GPSS-'+pdbCode) < 1:
                
            if reinitialize==1:
                cmd.do('reinitialize')
                
            pdbFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                      'gpssServerPymol.php?pdbId='+
                                      pdbCode + '&mode=pdbFile').read()
            
            
            if pdbFile:
                pdbFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                          'gpssServerPymol.php?pdbId='+
                                          pdbCode + '&mode=pdbFile')
                cmd.read_pdbstr(pdbFile.read(),'GPSS-'+pdbCode)
            
                # Make structure pretty
                cmd.do('set bg_rgb=[1,1,1]')
                cmd.do('hide everything')
                cmd.do('dss')
                cmd.do('show cartoon, GPSS-' + pdbCode)
                cmd.do('color gray')
                return True
            else:
                tkMessageBox.showerror('GPSS                                    ',
                                       pdbCode + ' was not found on the GPSS server.')

            cmd.do('zoom GPSS-'+pdbCode)

        else:

            return True
        
    ########################################################################
    # gpssAll
    #
    ########################################################################
    def gpssAllDialog(app):
        pdbCode = tkSimpleDialog.askstring('GPSS - All                      ',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode:
            gpssAll(pdbCode)

    def gpssAll(pdbCode):
            gpssLigand(pdbCode,1)
            gpssMetal(pdbCode,0)
            gpssPeptide(pdbCode,0)
            gpssDna(pdbCode,0)
            gpssSitesCsa(pdbCode,0)
            gpssCastp(pdbCode,0)

    ########################################################################
    # gpssCastp
    #
    ########################################################################
    def gpssCastpDialog(app):
        pdbCode = tkSimpleDialog.askstring('GPSS - CastP                    ',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode:
            gpssCastp(pdbCode,1)
        
    def gpssCastp(pdbCode,reinitialize):
        if _gpssLoadPdb(pdbCode,reinitialize):

            # Get pocket count
            castpInfo = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                        'gpssServerPymol.php?pdbId='+
                                        pdbCode + '&mode=castpInfo')
            castpLength = int(castpInfo.readlines()[-1][13:16])
            
            # Get pocket atoms
            while castpLength > 0:
                castpSurfaceFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                                   'gpssServerPymol.php?pdbId='+
                                                   pdbCode + '&surfaceId='+str(castpLength)+
                                                   '&mode=castpSurface')
                objectId = 'CASTp-'+str(castpLength)
                cmd.read_pdbstr(castpSurfaceFile.read(),objectId)
                cmd.do('hide everything,' + objectId)
                cmd.do('show surface,' + objectId)
                cmd.do('set transparency=0.2,'+ objectId) 
                castpLength-=1
                
        cmd.do('zoom GPSS-'+pdbCode)

    ########################################################################
    # gpssLigand
    #
    ########################################################################
    def gpssLigandDialog(app):
        pdbCode = tkSimpleDialog.askstring('GPSS - Ligand                   ',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode:
            gpssLigand(pdbCode,0)

    def gpssLigand(pdbCode,reinitialize):
        
        if _gpssLoadPdb(pdbCode,reinitialize):

            # Get ligand surfaces
            ligandInfo = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                          'gpssServerPymol.php?pdbId='+
                                          pdbCode + '&mode=ligandInfo').read()
            if ligandInfo:
                ligandInfoArray = ligandInfo.split(':')

                # Get pocket atoms
                for ligandSurface in ligandInfoArray:
                    ligandSurfaceFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                                        'gpssServerPymol.php?pdbId='+ pdbCode +
                                                        '&surfaceId=' + ligandSurface +
                                                        '&mode=ligandSurface')
                    ligandId = ligandSurface.split('.')
                    objectId = 'Ligand-' + ligandId[1] + '.' + ligandId[2] + '.' + ligandId[3]
                    cmd.read_pdbstr(ligandSurfaceFile.read(),objectId)
                    cmd.do('hide everything,' + objectId)
                    cmd.do('show surface,' + objectId)
                    cmd.do('show sticks,' + objectId + ' and resn ' + ligandId[1])
                    
                cmd.do('zoom GPSS-'+pdbCode)

    ########################################################################
    # gpssMetal
    #
    ########################################################################
    def gpssMetalDialog(app):
        pdbCode = tkSimpleDialog.askstring('GPSS - Metal                   ',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode:
            gpssMetal(pdbCode,0)

    def gpssMetal(pdbCode,reinitialize):
        
        if _gpssLoadPdb(pdbCode,reinitialize):

            # Get ligand surfaces
            metalInfo = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                        'gpssServerPymol.php?pdbId='+
                                        pdbCode + '&mode=metalInfo').read()

            if metalInfo:
                metalInfoArray = metalInfo.split(':')

                # Get pocket atoms
                for metalSurface in metalInfoArray:
                    metalSurfaceFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                                       'gpssServerPymol.php?pdbId='+ pdbCode +
                                                       '&surfaceId=' + metalSurface +
                                                       '&mode=metalSurface')
                    metalId = metalSurface.split('.')
                    objectId = 'Metal-' + metalId[1] + '.' + metalId[2] + '.' + metalId[3]
                    cmd.read_pdbstr(metalSurfaceFile.read(),objectId)
                    cmd.do('hide everything,' + objectId)
                    cmd.do('show surface,' + objectId)
                    cmd.do('show sphere,' + objectId + ' and resn ' + metalId[1])

                    cmd.do('zoom GPSS-'+pdbCode)

    ########################################################################
    # gpssPeptide
    #
    ########################################################################
    def gpssPeptideDialog(app):
        pdbCode = tkSimpleDialog.askstring('GPSS - Peptide                   ',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode:
            gpssPeptide(pdbCode,0)

    def gpssPeptide(pdbCode,reinitialize):
        
        if _gpssLoadPdb(pdbCode,reinitialize):

            # Get ligand surfaces
            peptideInfo = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                          'gpssServerPymol.php?pdbId='+
                                          pdbCode + '&mode=peptideInfo')
            # Get pocket atoms
            for peptideSurface in peptideInfo:
                peptideSurfaceFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                                     'gpssServerPymol.php?pdbId='+ pdbCode +
                                                     '&surfaceId=' + peptideSurface +
                                                     '&mode=peptideSurface')
                peptideId = peptideSurface.split('.')
                objectId = 'Peptide-' + peptideId[1]
                cmd.read_pdbstr(peptideSurfaceFile.read(),objectId)
                cmd.do('hide everything,' + objectId)
                cmd.do('show surface,' + objectId)
                cmd.do('set transparency=0.2,'+ objectId) 
            cmd.do('zoom GPSS-'+pdbCode)

    ########################################################################
    # gpssDna
    #
    ########################################################################
    def gpssDnaDialog(app):
        pdbCode = tkSimpleDialog.askstring('GPSS - DNA                   ',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode:
            gpssDna(pdbCode,0)
            
    def gpssDna(pdbCode,reinitialize):
        
        if _gpssLoadPdb(pdbCode,reinitialize):

            # Get sites
            dna = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                  'gpssServerPymol.php?pdbId='+
                                  pdbCode + '&mode=dnaSurface').read()
            if dna:
                objectId = 'DNA-'+pdbCode;
                cmd.read_pdbstr(dna,objectId)
                cmd.do('hide everything,' + objectId)
                cmd.do('show surface,' + objectId)
                cmd.do('set transparency=0.2,'+ objectId) 
                cmd.do('zoom GPSS-'+pdbCode)

    ########################################################################
    # gpssCsa
    #
    ########################################################################
    def gpssSitesDialog(app):
        pdbCode = tkSimpleDialog.askstring('GPSS - Sites                   ',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode:
            gpssSitesCsa(pdbCode,0)
            
    def gpssSitesCsa(pdbCode,reinitialize):
        
        if _gpssLoadPdb(pdbCode,reinitialize):

            # Get sites
            sites = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                    'gpssServerPymol.php?pdbId='+
                                    pdbCode + '&mode=sitesCsa').read()
            if sites:
                objectId = 'SITES-CSA-'+pdbCode;
                cmd.read_pdbstr(sites,objectId)
                cmd.do('hide everything,' + objectId)
                cmd.do('show spheres,' + objectId)
                cmd.do('color yellow,' + objectId)
                cmd.do('zoom GPSS-'+pdbCode)


    ########################################################################
    # gpssFileUpload
    #
    ########################################################################
    def gpssFileUploadDialog(app):
        file = tkFileDialog.askopenfile(parent=app.root,
                                        mode='rb',
                                        title='GPSS - File Upload')
        if file != None:
            data = file.read()
            file.close()
            email = tkSimpleDialog.askstring('GPSS - File Upload                      ',
                                             'Please enter you email address.\n\nYou will be notified when your structure is available to be viewed.',
                                             parent=app.root)
            
            params = urllib.urlencode({'email': email,'userfile':file.name,'file':data})
            response = urllib.urlopen('http://gpss.mcsg.anl.gov/webservices/gpssServerFileUpload2.php',
                                      params)
            if response:
                tkMessageBox.showinfo('GPSS - File Upload',response.readlines())
            else:
                tkMessageBox.showerror('GPSS - File Upload',"There was a problem connecting to the server, please try again.")

    ########################################################################
    # gpssHelp
    #
    ########################################################################
    def gpssHelpDialog(app):
        tkMessageBox.showinfo('GPSS - Help                         ',
                              'Global Protein Surface Survey\nPymol Plugin v0.1\n\nT. Andrew Binkowski\nArgonne National Laboratory\n\n\nFor more information visit http://gpss.mcsg.anl.gov');


    # Extend pymol functionality
    #cmd.extend('gpssAll', gpssAll)



