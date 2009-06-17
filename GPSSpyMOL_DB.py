#############################################################################################################
# Global Protein Surface Survey (GPSS) PyMOL Plugin
#
# Version: 0.3.1 (14-JAN-2009)
#       David Borhani, 14-JAN-2009
#       Modified gpssCastp extensively to make this function faster, more useful.
#
#
# Version: 0.3 (July 9, 2007)
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
import sys, urllib, urllib2, zlib, os, subprocess, string

def __init__(self):
    
    # David Borhani, 14-JAN-2009
    # Get Xwindow geometry if not on Windows
    if (sys.platform).upper() != 'WIN32':
        window_geom = os.popen('xwininfo -root | grep geometry').read().rstrip().split()[1]
        self.window_width = int(window_geom.split('x')[0])
        self.window_height = int(window_geom.split('x')[1].replace('+','-').split('-')[0])
    
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
        
        GPSSpdbCode = 'GPSS_' + pdbCode

        _gpssVersionCheck('0.3');

        if cmd.get_names().count(GPSSpdbCode) < 1:
                
            if reinitialize==1:
                cmd.do('reinitialize')
                
            pdbFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                      'gpssServerPymol.php?pdbId='+
                                      pdbCode + '&mode=pdbFile').read()
            
            
            if pdbFile:
                pdbFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                          'gpssServerPymol.php?pdbId='+
                                          pdbCode + '&mode=pdbFile')
                cmd.read_pdbstr(pdbFile.read(),GPSSpdbCode)
            
                # Make structure pretty
                # David Borhani, 14-JAN-2009
                # Changed to color by chain, show lines only, keep user's bg color.
                #cmd.do('set bg_rgb=[1,1,1]')
                #cmd.do('hide everything')
                cmd.hide('everything', GPSSpdbCode)
                #cmd.do('dss')
                cmd.dss(GPSSpdbCode)
                #util.cbc(selection='(GPSS-' + pdbCode + ')',first_color=7,quiet=1,legacy=1,_self=cmd)
                #sel_string = '("' + GPSSpdbCode + '")'
                #cmd.do('util.cbc(selection=' + sel_string + ',first_color=9,quiet=1,legacy=1)')
                cmd.do('util.cbc(selection=("' + GPSSpdbCode + '"),first_color=9,quiet=1,legacy=1)')
                #cmd.do('show lines, ' + GPSSpdbCode)
                cmd.show('lines', GPSSpdbCode)
                #cmd.do('show cartoon, ' + GPSSpdbCode)
                #cmd.do('color gray')
                return True
            else:
                tkMessageBox.showerror('GPSS                                    ',
                                       pdbCode + ' was not found on the GPSS server.')

            #cmd.do('zoom GPSS-'+pdbCode)
            cmd.zoom(GPSSpdbCode, animate=0)

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
        pdbCode = tkSimpleDialog.askstring('GPSS - CASTp                    ',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode:
            gpssCastp(pdbCode,1)
        
    def yield_ranges(strlist):
        # David Borhani, 14-JAN-2009
        '''
        Convert a list of (numeric) strings <strlist> to list of 
        (numberic) strings with consecutive ranges for each list element.
        
        E.G.['1924','2136','2138','2139','2140','2146','2148','2151','2152','2153']
        converts to ['1924','2136','2138-2140','2146','2148','2151-2153']
        '''
        from operator import itemgetter
        from itertools import groupby
        
        intlist = []
        strlist.sort()
        for e in strlist:
            intlist.append(int(e))
        
        rangelist = []
        for k, g in groupby(enumerate(intlist), lambda (i,x):i-x):
            temp = map(itemgetter(1), g)
            if len(temp) == 1:
                rangelist.append(str(temp[0]))
            else:
                rangelist.append(str(temp[0]) + '-' + str(temp[-1]))
        
        return rangelist
    
    def gpssCastp(pdbCode,reinitialize):
        # David Borhani, 14-JAN-2009
        # Modified extensively to make this faster, more useful.
        
        GPSSpdbCode = 'GPSS_' + pdbCode
        if _gpssLoadPdb(pdbCode,reinitialize):
            
            # Get pocket count and pocket/cavity summary information
            castpInfo = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                        'gpssServerPymol.php?pdbId='+
                                        pdbCode + '&mode=castpInfo')
            castp_pocinfo = castpInfo.readlines()
            castpLength = int(castp_pocinfo[-1][13:16])
            n_pockets = castpLength
            
            # Write out pocket information file:
            pocinfo = pdbCode + '.pocinfo'
            f_pocinfo = open(pocinfo, 'w')
            for line in castp_pocinfo:
                f_pocinfo.write(line)
            f_pocinfo.close()
            
            # Write out CAVITY ("n_mouth" = 0) information file:
            cavinfo = pdbCode + '.cavinfo'
            f_cavinfo = open(cavinfo, 'w')
            f_cavinfo.write('\t\t\t>>>>>  CAVITIES ONLY  <<<<<\n')
            f_cavinfo.write(castp_pocinfo[0]) # Header line
            for line in castp_pocinfo[1:]:
                if (line[16:20] == '   0'):
                        f_cavinfo.write(line)
            f_cavinfo.close()
            
            # Get ALL pocket atoms, ordered from largest pocket to smallest...
            poc = pdbCode + '.poc'
            f_poc = open(poc, 'w')
            while castpLength > 0:
                castpSurfaceFile = urllib2.urlopen('http://gpss.mcsg.anl.gov/webservices/'+
                                                   'gpssServerPymol.php?pdbId='+
                                                   pdbCode + '&surfaceId='+str(castpLength)+
                                                   '&mode=castpSurface')
                CASTpSelId = 'CASTp_' + pdbCode + '_' + str(castpLength)
                # Get the URL data all as one string (complete w/ \n's)
                castp_surf = castpSurfaceFile.read()
                
                # Write it to the file
                f_poc.write('>>>>> POCKET '+CASTpSelId+'\n')
                f_poc.write(castp_surf)
                
                # Create a PyMOL object from the string...
                # David Borhani, 14-JAN-2009
                # NO, WAY TOO SLOW. Create a SELECTION instead.
                #cmd.read_pdbstr(castp_surf,CASTpSelId)
                #cmd.do('hide everything,' + CASTpSelId)
                
                # List of atom (numbers, as STRING) in this pocket
                poc_atoms = []
                for line in castp_surf.splitlines():
                    atmnum = line[6:11].strip()
                    if(atmnum.isdigit()):
                        if len(poc_atoms) == 0:
                            poc_atoms = [atmnum]
                        else:
                            poc_atoms.append(atmnum)
                
                # Create the PyMOL selection:
                poc_atoms_range = yield_ranges(poc_atoms)
                # Key point from Warren DeLano:
                # Create selection in chunks so that don't exceed ~1024 characters
                chunk = 50
                cmd.select(CASTpSelId, 'none', quiet=0, enable=0)
                while len(poc_atoms_range):
                    id_pat = string.join(poc_atoms_range[0:chunk],"+")
                    cmd.select(CASTpSelId,'?'  +CASTpSelId + ' & ' + GPSSpdbCode + \
                        ' | (id ' + id_pat + ')', quiet=0, enable=0)
                    poc_atoms_range = poc_atoms_range[chunk:]
                
                ### WAY TOO SLOW TO SHOW THE SURFACES OF ALL POCKETS!
                # Also, the surfaces are inaccurate, as only some of the totality of atoms
                # are present, so true cavities may appear to just be invaginations.
                ######cmd.do('show surface,' + CASTpSelId)
                ######cmd.do('set transparency=0.2,'+ CASTpSelId) 
                
                castpLength -= 1
                
            f_poc.close()
            
            # Display the three files, so the user can DETERMINE which
            #  pockets/cavities are of interest.
            file_to_open = [poc, pocinfo, cavinfo]
            if (sys.platform).upper() != 'WIN32':
                # Set up xwindow geometry
                window_offset = 40
                window_x = self.window_width/2
                #window_y = self.window_height
                for (i, win) in enumerate(file_to_open):
                    x = window_x + i*window_offset
                    y = i*window_offset
                    geometry = '-geometry 80x40+%d+%d' % (x,y)
                    xterm_command= 'xterm -j +ls -rightbar ' + geometry + ' -title ' + \
                      win + ' -e vim ' + win
                    task = subprocess.Popen(xterm_command.split())
            else:
                for win in file_to_open:
                    wordpad_command= 'C:\\Program Files\\DeLano Scientific\\PyMOL\\modules\\pmg_tk\\wordpad.exe ' + '.\\' + win
                    task = subprocess.Popen(wordpad_command.split())
                        
        # And shut off all selections
        cmd.deselect()
        cmd.zoom(GPSSpdbCode, animate=0)
        cmd.set('two_sided_lighting', '1', quiet=0)
        print '\nCASTp> Created %d pocket/cavity selections for PDB entry %s\n' \
          % (n_pockets,pdbCode)

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



