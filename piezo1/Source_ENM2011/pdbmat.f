CoM=====================================================================
CoM...Pdbmat: Calculates the Hessian matrix using an ENM.
CoM---------------------------------------------------------------------
CoM
CoM   Goal: To compute the low-frequency (usually collective) normal 
CoM   modes of vibration of the system.
CoM
CoM   Action: Mass-weighted second derivatives energy matrix is computed,
CoM   using Tirion's model, that is, an elastic network model (ENM).
CoM
CoM   In ENMs, pairs of particles (atoms) are linked by springs.
CoM
CoM   In Tirions's model, springs are set between atoms less than 
CoM   CUToff Angstroms away from each others.
CoM
CoM   In HINSEN's version, springs are weighted by exp(-dij/rkh)**2.
CoM   Herein, both kind of typical distances (CUToff & rkh) can be mixed.
CoM
CoM   Next: To obtain the modes, the matrix produced by pdbmat has to be 
CoM   diagonalized. This can be done using programs like DIAGSTD or, 
CoM   for larger systems, DIAGRTB, BLZPACK, etc.
CoM
CoM---------------------------------------------------------------------
CoM   INPUT: 
CoM   ******
CoM   A parameter file named pdbmat.dat   
CoM
CoM   Note that each run of pdbmat produces a pdbmat.dat_run file,
CoM   where parameter values are shown (and shortly commented).
CoM   pdbmat.dat_run can be modified and used as a pdbmat.dat file, 
CoM   for further runs. 
CoM
CoM   So, to begin with, the simplest is to compile and run pdbmat...
CoM   and have a look at the pdbmat.dat_run file so produced.
CoM   If you want to get there the syntax of other available commands,
CoM   just raise the PRINting level.
CoM
CoM  -Among the parameters: 
CoM
CoM  *The name of a file with the coordinates of the system, 
CoM   in FREE or PDB (protein data bank) format.
CoM
CoM   Free format: x, y, z, mass.
CoM   PDB  format: Lignes with ATOM or HETATm keywords are considered.
CoM        Masses can be given in the Bfactors column.
CoM   Note that masses can be read, but they are not required.
CoM   In the later case, they are all set to 1.0 (as often done).
CoM
CoM  *A way to identify pairs of neighbors:
CoM   Either a CUTOFF value (standard Tirion's model)
CoM   Or a file with a list of pairs of neighbors.
CoM   Ligne format of this latter file: atom-number atom-number
CoM
CoM   Alternatively, Hinsen's model can be used.
CoM
CoM   For one-atom-per-residue protein models, typical values are:
CoM   Tirion's distance cutoff = 10-12 Angstroms.
CoM   Hinsen's typical range   = 3 Angstroms.
CoM   For all-atom protein models (like in most PDB cases), smaller
CoM   values have to be used. For Tirion's model, a rule-of-thumb is 
CoM   to take a value as small as possible, but so that the network
CoM   is not split into smaller ones (if it is, you will obtain more
CoM   than six zero-eigenvalues).
CoM
CoM   OUTPUT:
CoM   *******
CoM  -Default output matrix filename:
CoM   Formatted file: pdbmat.sdijf 
CoM   Binary    file: pdbmat.sdijb 
CoM
CoM   This is a matrix in format: i, j, non-zero-i-j-matrix-element
CoM
CoM  -Output coordinate filename (in free format; not mandatory):
CoM   pdbmat.xyzm
CoM
CoM   This is a coordinate file with, for each atom: 
CoM   x, y, z, mass, block-number
CoM   For a pdb file, the block-number is the amino-acid residue number
CoM  (it is of use only in programs like DIAGRTB).
CoM
CoM   More specialized ones:
CoM   ----------------------
CoM  -VMD-command filename: 
CoM   file_produced_by_pdbmat.vmd
CoM
CoM   This is a command file for vizualising the elastic network
CoM   with VMD, the Visual Molecular Dynamics program (checked for v1.8):
CoM   vmd -e file_produced_by_pdbmat.vmd
CoM   But if you just want to vizualise a standard ENM (defined with
CoM   a CUTOff distance), there are more efficient (faster) ways...
CoM
CoM  -Molscript-command-file:
CoM   file_produced_by_pdbmat.in
CoM
CoM   This is a command file for vizualising the elastic network, 
CoM   using Molscript (checked for v2.1), i.e.:
CoM   molscript < file_produced_by_pdbmat.in > file.postscript
CoM   But if you just want to vizualise a standard ENM (defined with a 
CoM   CUTOff distance) with Molscript, there are much faster ways...
CoM
CoM.....................................................................
      program pdbmat
      implicit none
      integer natmax, nresmx, ntopmax, nvoismx

CoM   This is a fortran 77 program (Sorry), so it has predefined:
c   
c     **************
CoM   MEMORY LIMITS:
c     **************

c     NATMAX  :  Maximum number of atoms (particles).
c     NRESMX  :  Maximum number of residues (groups of particles).
c     NVOISMX :  Maximum number of pairs of neighbors.  
c     NTOPMAX :  Maximum number of chemical bonds (of any kind).
   
      parameter( NATMAX=50000 )
      parameter( NRESMX=50000 )

      parameter( NVOISMX=1000000 )
      parameter( NTOPMAX=10*natmax )

CoM   Increase them if needed. To (re)compile pdbmat, type:
CoM   make pdbmat
CoM   or:
CoM   g77 -o pdbmat pdbmat.f
CoM   or use your favorite fortran compiler instead (of g77).
CoM
CoM   To run it in the current directory, type: ./pdbmat
CoM.....................................................................
CoM
CoM   ABOUT Tirion's model (Still in french, sorry; maybe one day...):
CoM   ****************************************************************
CoM   Principe du modele (Tirion, 1996): 
CoM   
CoM   Tous les atomes a moins de "cutoff" les uns des autres 
CoM   sont supposes lies par des ressorts, qui ont tous
CoM   la meme raideur.
CoM   Simplification supplementaire par rapport au modele initial: 
CoM   les atomes sont supposes avoir tous la meme taille
CoM  (le cutoff est le meme pour toutes les paires d'atomes).
CoM   On peut de plus poser qu'ils ont tous la meme masse.
CoM   Sinon, celles-ci sont lues dans la colonne des 
CoM   facteurs B du fichier pdb.
CoM
CoM   Principaux resultats:
CoM
CoM   Les modes de vibration de basse frequence obtenus
CoM   a partir d'un tel modele sont tres voisins de ceux
CoM   obtenus avec un modele beaucoup plus detaille, tels
CoM   ceux utilises lors des etudes de Dynamique Moleculaire.
CoM
CoM   Dans le cas ou le mouvement fonctionnel d'une proteine est un 
CoM   mouvement d'ensemble (collectif), on constate qu'il peut tres 
CoM   souvent etre decrit comme une combinaison lineaire de quelques
CoM   uns de ces modes (de un a trois).
CoM
CoM   Principaux avantages:
CoM
CoM   Pas besoin de prendre en compte tous les atomes.
CoM   Pas besoin de minimisation d'energie prealablement
CoM   au calcul des modes de vibration (E=0 par construction).
CoM
CoM.....................................................................
CoM
CoM   MAIN REFERENCES:
CoM   ****************
CoM   1) M.M. Tirion (1996):
CoM  "Large amplitude elastic motions in proteins from
CoM   a single-parameter, atomic analysis",
CoM   Phys. Rev. letters vol.77(9), p1905-1908.
CoM
CoM   2) K. Hinsen (1998):
CoM  "Analysis of domain motions by approximate normal mode calculations"
CoM   Proteins vol.33, p417-429.
CoM
CoM   3) F. Tama, Y.H. Sanejouand (2001):
CoM  "Conformational change of proteins arising 
CoM   from normal modes calculations"
CoM   Protein Engineering vol.14, p1-6.
CoM
CoM   4) A.R. Atilgan, S.R. Durell, R.L. Jernigan, M.C. Demirel,
CoM   O. Keskin, I. Bahar (2001):
CoM  "Anisotropy of Fluctuation Dynamics of Proteins with an Elastic 
CoM   Network Model"
CoM   Biophys. J. vol.80, p.505-515.
CoM
CoM   5) S. Nicolay, Y.H. Sanejouand (2006):
CoM  "Functional modes of proteins are among the most robust"
CoM   Phys. Rev. letters vol.96, p078104.
CoM
CoM.....................................................................
c     Previous versions:
c     YHS-Nov-1996: Version 1.00 (Toulouse).
c     Version used by the ELNEMO Web site (http://www.elnemo.org):
c     YHS-Feb-2004: Version 3.46 (Lyon).
c     Versions released (http://ecole-modelisation.free.fr/modes.html):
c     YHS-Mar-2001: Version 3.31 (Bordeaux). 
c     YHS-Feb-2004: Version 3.50 (Lyon). 
c     YHS-Feb-2008: Version 3.73 (Lyon). 
 
CoM   Since previously released version (v3.73):
CoM   New features:
CoM  *Anisotropic on-site harmonic potential.
CoM  *Statistics for the number of neighbors.
CoM  *Crude H-bond energy term. 
CoM  *Molscript input is now usable as is.
CoM   Bug fixes:
CoM  *x-y-z input coordinate format should now work.
CoM  *C-character in column 1 (some compilers do mind)
CoM
CoM   In case of problem, feel free to contact: 
CoM   Yves-Henri.Sanejouand@univ-nantes.fr
CoM  (bug reports may help you, but also others)
CoM.....................................................................
ChnG  Test ifort
      integer nmotsmax
      parameter(nmotsmax=100)
      integer fatres(nresmx+1), 
     .        i, iathb(ntopmax), idmax, idres(nresmx), ii, imax, imin, 
     .        ires, iresat(natmax), iseed, ivois(nvoismx),
     .        j, jangle(ntopmax), jat, jathb(ntopmax), jbond(ntopmax), 
     .        jdihe(ntopmax), jj, jsomang(ntopmax), jvois(nvoismx), 
     .        k, kcom, kk, klist, l, ll, lmot, lnom, lnomlst, lnommls, 
     .        lnommtx, lnompdb, lnomvmd, ltypmas,
     .        namax, namin, nangle(ntopmax), nangles, natom, 
     .        nbig, nbmax, nbmin, nbond(ntopmax), nbonds, nc, ndat, 
     .        ndihe(ntopmax), ndihs, nh, nhb, nkij, nl, nmax, nmin, 
     .        nmots, nn, nntr, nnzero, no, noff, nonsit, npair, nres, 
     .        numc(nresmx), numh(nresmx), numn(nresmx), numo(nresmx), 
     .        nunit, nunknown, nvois, nvoisat(natmax), nzero,
     .        prtlev, statvois(natmax), 
     .        uninp, unlst, unmol, unout, unpdb, unrsd, unvmd
      double precision cutbnd, cuthb, cutoff, dcn, dco, ddf, defmas,
     .        der2(3,3*natmax), dist, dist2, dmax, dmin, dmoy, drms, 
     .        elemnt, elmax, fvois(natmax), 
     .        kangle, kbond, kdihe, kfce, kfmax, kfmin, khb, kij, knonb, 
     .        konsite(3*natmax), kvois(natmax), kx, ky, kz,
     .        levelshft, massat(natmax), nmoy, nrms,
     .        random, rave, rbig, rdev, rdmass, rinput, rkh, rmax, rmin, 
     .        rsmall, rx, ry, rz, 
     .        trace, unknown, xat(natmax), yat(natmax), zat(natmax) 
      logical qbinary, qchim, qenm, qerror, qexist, qfread, qhet, 
     .        qinter, qlist, qmasse, qmtx, qonsit, qok, qpdb, 
     .        qvois(natmax)
      character atonam(natmax)*4, cformat*32, codpdb*4, csep*1, 
     .        cstatus*32, 
     .        lign80*80, motinp*80, mots(nmotsmax)*132, 
     .        nomfich*64, nomlst*64, nommls*64, nommtx*64, nomout*66,
     .        nompdb*64, nomvmd*64,
     .        program*8, progrer*11, progrwn*11, 
     .        residus_standards*132, residus_stshort*21, 
     .        resnam(natmax)*4, segid(natmax)*4, ssunam(natmax)*1, 
     .        ssusel*1, typbond*80, typmas*80, typout*80, version*32
      parameter(rbig=1e10,rsmall=1e-10,unknown=9999.d9)
c.......................................................................
      version=' Version 3.89, November 2011.'
c.......................................................................
CoM   Specific details:
CoM   -----------------
      idmax=21
      residus_standards='   ILE PHE TRP LEU CYS VAL MET TYR ALA HIS '//
     .                     'GLY THR SER PRO ARG GLN ASN ASP GLU LYS '
      residus_stshort='IFWLCVMYAHGTSPRQNDEKX'
 
      program=' Pdbmat>'
      progrer='%Pdbmat-Er>'
      progrwn='%Pdbmat-Wn>'

      write(6,'(2A)') program,
     .' Computes the Hessian matrix, using an Elastic Network Model.'
      write(6,'(2A)') program,version

CoM   Fixed parameters for H-bonding:
CoM   Approximate distance between C and N involved in H-bonds: 3 Angs.
      cuthb=4.2
CoM   Approximate distance between bonded C and O atoms: 1.2 Angs.
      dco=1.2
c     Default mass value:
      defmas=1.0

c     ===============================================
c     Ouverture et Lecture du fichier d'instructions:
c     ===============================================
 
c     Default values:
      nompdb='structure.ent'
      cutoff=10.d0
      qenm=.false.
c     Neighbor list:
      nomlst='NONE'
      qlist=.false.
c     Typical distance for Hinsen's weigth:
c    (negative value means 1/rkh=0)
      rkh=-1.d0
      nommls='NONE'
      nomvmd='NONE'
      typmas='CONS'
      ltypmas=4
      qmasse=.false.
      knonb=1.0d0
c     Topological terms:
      typbond='NONE'
c     Distance-cutoff value for defining covalent (chemical) bonds:
      cutbnd=4.0d0
      kangle=0.0d0
      kbond=1000.0d0
      kdihe=0.0d0
      khb=0.0d0
c     Others:
      typout='  FREE'
      qmtx=.false.
      qbinary=.false.
      prtlev=0
      levelshft=rsmall
c     Used only with a levelshift, to add some noise (not useful ?).
      iseed=27041961

      nunit=10
      uninp=nunit
      nunit=nunit+1
      nomfich='pdbmat.dat'
      lnom=10
      cformat="FORMATTED"
      cstatus="old"
      call openam(nomfich,lnom,cformat,cstatus,uninp,.false.,
     .     qinter,qexist)
      if (qinter.or..not.qexist) then 
          write(6,'(/2A/(2A))') progrwn,
     .  ' No pdbmat.dat command file found.',
     .    progrwn,' Defaults assumed for all options. ',
     .    progrwn,' See the pdbmat.dat_run file if you need an example.'
          goto 110
      else
          write(6,'(/2A)') program,
     .  ' Options to be read in pdbmat.dat file.'
      endif
 
 50   continue
      read(uninp,'(A)',end=100) lign80
 
      kcom=index(lign80,'!')
      k=index(lign80,'=') 
      motinp=' '
      if (k.gt.0.and.(kcom.le.0.or.kcom.gt.k)) then 
          motinp=lign80(1:k)
      else
          if (k.le.0.and.kcom.gt.1) then
          write(6,'(/2A/A)') progrwn,
     .  ' No separator (=) in command ligne:',
     .    lign80
          write(6,'(2A)') progrwn,' This ligne is skipped.'
          endif
          goto 50
      endif
      call mintomaj(motinp)

      kcom=index(lign80(k+1:80),'!')
      if (kcom.gt.0) lign80(k+kcom:80)=' '
      klist=index(lign80(k+1:80),'?')

      if (index(motinp,' FILENAME').gt.0.or.
     .    index(motinp,' SCRIPT').gt.0) then 
          if (index(motinp,'MATRI').gt.0) then
              nommtx=lign80(k+1:80)
              qmtx=.true.
          elseif (index(motinp,'LIST ').gt.0.or.
     .            index(motinp,'NEIGH').gt.0) then
              nomlst=lign80(k+1:80)
              qlist=.true.
          elseif (index(motinp,'MOLS').gt.0) then
              nommls=lign80(k+1:80)
          elseif (index(motinp,'VMD').gt.0) then
              nomvmd=lign80(k+1:80)
          else
              nompdb=lign80(k+1:80)
          endif
      else if (index(motinp,' DEFINITION').gt.0) then
          typbond=lign80(k+1:80)
          call mintomaj(typbond)
          call stringcl(typbond,lnom)
          if (typbond(1:3).eq.'ALL') then 
              typbond=' ALL'
          else if (typbond(1:3).eq.'NEI') then 
              typbond=' ALL'
          else if (typbond(1:3).eq.'NON') then 
              typbond='NONE'
          else if (typbond(1:3).eq.'CON') then
              typbond='CONSECUTIF'
          else
              write(6,'(/3A)') progrwn,' Bond definition :',
     .        typbond(1:4)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: NONe, ALL, NEIGhbors, CONsecutive.'
              write(6,'(A)') ' Default assumed.'
              typbond='NONE'
          endif
      else if (index(motinp,'MASS').gt.0) then
          typmas=lign80(k+1:80)
          call mintomaj(typmas)
          call stringcl(typmas,ltypmas)
          if (typmas(1:3).eq.'PDB'.or.typmas(1:3).eq.'COO') then
              qmasse=.true.
              typmas='COOR'
              ltypmas=4
          else if (typmas(1:3).eq.'CON') then
              qmasse=.false.
              typmas='CONS'
              ltypmas=4
          else
              read(typmas(1:ltypmas),*,end=60,err=60) rdmass
              if (rdmass.gt.0) then
                  defmas=rdmass
                  goto 65
              endif
 60           continue
              write(6,'(/3A)') progrwn,' MASS values: ',
     .        typmas(1:3)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: CONstant, COOr, PDB, or a real number.'
              write(6,'(A)') ' Default assumed.'
              qmasse=.false.
              typmas='CONS'
              ltypmas=4
 65           continue
          endif
      else if (index(motinp,'FORMAT').gt.0) then
          typout=lign80(k+1:80)
          call mintomaj(typout)
          call stringcl(typout,lnom)
          if (typout(1:1).eq.'B'.or.typout(1:1).eq.'U') then
              qbinary=.true.
              typout='BINARY'
          else if (typout(1:1).ne.'F') then
              write(6,'(/3A)') progrwn,' Kind of matrix format :',
     .        typout(1:1)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: Free, Binary, Formatted, Unformatted.'
              write(6,'(A)') ' Default assumed.'
              qbinary=.false.
              typout='  FREE'
          else
              qbinary=.false.
              typout='  FREE'
          endif
      else 
          qok=.false.
          read(lign80(k+1:80),*,end=90,err=90) rinput
          if (index(motinp,'SHIFT ').gt.0) then 
               qok=.true.
               levelshft=rinput
          else if (index(motinp,'CUTOF').gt.0.or.
     .             index(motinp,'DISTANCE').gt.0) then
               qok=.true.
               cutoff=rinput
          else if (index(motinp,'HINSEN').gt.0) then
               qok=.true.
               rkh=rinput
          else if (index(motinp,'INTERAC').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   knonb=rinput
               endif
               if (index(motinp,' CUTOF').gt.0.or.
     .             index(motinp,' DIST').gt.0) then
                   qok=.true.
                   cutoff=rinput
               endif 
          else if (index(motinp,'BOND').gt.0.and.
     .        (index(motinp,' FORCE ').gt.0.or.
     .         index(motinp,' CONST').gt.0)) then
               qok=.true.
               kbond=rinput
          else if (index(motinp,' LENGTH').gt.0) then
               qok=.true.
               cutbnd=rinput
          else if (index(motinp,'PRINT').gt.0) then
               qok=.true.
               prtlev=int(rinput)
          else if (index(motinp,'ANGLE').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   kangle=rinput
               endif
          else if (index(motinp,'DIHE').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   kdihe=rinput
               endif
          else if (index(motinp,'HBON').gt.0.or.
     .             index(motinp,'H-BON').gt.0.or.
     .             index(motinp,'HYDRO').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   khb=rinput
               endif
          endif
  90      continue
          if (.not.qok) then
               write(6,'(/2A/A)') progrwn,
     .       ' No known or incomplete set of keywords in ligne:',
     .         motinp
               write(6,'(2A)') progrwn,
     .       ' This command ligne is skipped.'
          endif
      endif
      goto 50
 
 100  continue
      close(uninp)
 110  continue

      call stringcl(nompdb,lnompdb)
      call stringcl(nomlst,lnomlst)
      call stringcl(nommls,lnommls)
      call stringcl(nomvmd,lnomvmd)
      if (nomlst.eq.'none'.or.nomlst.eq.'NONE') qlist=.false.
      if (nommls.eq.'none') nommls='NONE'
      if (nomvmd.eq.'none') nomvmd='NONE'
      
      if (.not.qmtx) then
      if (qbinary) then
        nommtx="matrix.sdijb"
      else
        nommtx="matrix.sdijf"
      endif
      endif
      call stringcl(nommtx,lnommtx)

c     Resume des commandes:
c     ---------------------

      write(6,'(/3A)') program,' Coordinate filename     = ',
     .      nompdb(1:lnompdb)
      if (qlist) then
      write(6,'(3A)') program,' Neighbor-list filename  = ',
     .      nomlst(1:lnomlst)
      else
      write(6,'(/2A,F10.2)') program,
     .        ' Distance cutoff         = ',cutoff
      endif
      if (rkh.gt.0.d0)
     .write(6,'(8X,A,F10.2)') " Hinsen's typical range  = ",rkh
      if (.not.qlist)
     .write(6,'(A,F10.2)') 
     .'         Force constant          = ',knonb
      if (typbond.ne.'NONE') then
      write(6,'(A,6X,A)') 
     .'         Kind of bond definition = ',typbond(1:4)
      write(6,'(A,F10.2)') 
     .'         Maximum bond length     = ',cutbnd,
     .'         Bond force constant     = ',kbond,
     .'         Angle force constant    = ',kangle,
     .'         Dihedral force constant = ',kdihe
      endif
      if (khb.gt.0) 
     .write(6,'(A,F10.2)') 
     .'         Hydrogen force constant = ',khb
      write(6,'(A,6X,A)') 
     .'         Mass values             = ',typmas(1:ltypmas)
      write(6,'(3A)') program,' Matrix filename         = ',
     .      nommtx(1:lnommtx)
      if (prtlev.gt.0) then
      write(6,'(2A,1PG10.1)') program,
     .        ' Levelshift              = ',levelshft
      write(6,'(A,3X,I7)') 
     .'         Printing level          = ',prtlev
      endif
      if (nommls.ne.'NONE')
     .write(6,'(3A)') program,' Molscript filename      = ',
     .      nommls(1:lnommls)
      if (nomvmd.ne.'NONE')
     .write(6,'(3A)') program,' VMD script filename     = ',
     .      nomvmd(1:lnomvmd)

c     Sauvegarde du fichier de commandes complet:
c     -------------------------------------------

      uninp=nunit
      nunit=nunit+1
      nomfich='pdbmat.dat_run'
      lnom=14
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,lnom,cformat,cstatus,uninp,.false.,
     .     qinter,qexist)

c     Pour etre plus clair:
      if (typbond.eq.'NONE') then
          cutbnd=0.d0
          kbond=0.d0
          kangle=0.d0
          kdihe=0.d0
      endif
      if (qlist) cutoff=-1

      write(uninp,'(2A)') 
     .'! This file can be modified and used as a command file',
     .' (named pdbmat.dat) for pdbmat.'
      write(uninp,'(2A)') ' Coordinate FILENAME        = ',
     .      nompdb(1:lnompdb)
      write(uninp,'(2A)') ' MATRIx FILENAME            = ',
     .      nommtx(1:lnommtx)
      write(uninp,'(A,F10.3,A)') ' INTERACtion DISTance CUTOF = ',
     .      cutoff,' ! For defining the list of interacting atoms.'
      write(uninp,'(A,F10.3,A)') ' INTERACtion FORCE CONStant = ',
     .      knonb,' ! For specifying frequency units.'
      if (prtlev.gt.0.or.rkh.gt.0.d0)
     .write(uninp,'(A,F10.3,A)') " HINSEN's typical range     = ",
     .      rkh,' ! Force constant weighting (if negative: none).'
      write(uninp,'(A,6X,2A)') ' MASS values                = ',
     .typmas(1:ltypmas),' ! CONstant, a value, or from COOrdinate file.'
      write(uninp,'(A,8X,I2,A)') ' Output PRINTing level      = ',
     .      prtlev,' ! =1: more detailled. =2: debug level.'
c     Rarely used:
      if (prtlev.gt.0.or.typbond.ne.'NONE') then
      write(uninp,'(A,6X,2A)') ' Bond DEFINITION            = ',
     .      typbond(1:4),
     .  ' ! NONe, ALL, between NEIGhbors or CONsecutive atoms.'
      write(uninp,'(A,F10.3,A)') ' Maximum bond LENGTH        = ',
     .      cutbnd,' ! Atoms closer are assumed to be bonded.'
      write(uninp,'(A,F10.3)') ' BOND FORCE CONStant        = ',kbond
      write(uninp,'(A,F10.3)') ' ANGLE FORCE CONStant       = ',kangle
      write(uninp,'(A,F10.3)') ' DIHEdral FORCE CONStant    = ',kdihe
      endif
      if (prtlev.gt.0.or.khb.gt.0)
     .write(uninp,'(A,F10.3)') ' HYDROgen FORCE CONStant    = ',khb
      if (prtlev.gt.0) then
      write(uninp,'(A,1PG10.1,A)') ' LevelSHIFT                 = ',
     .      levelshft,
     .  ' ! Non-zero value often required (numerical reasons).'
      write(uninp,'(A,4X,2A)') ' Matrix FORMAT              = ',
     .      typout(1:6),' ! Free, or Binary, matrix saved.'
      endif
c     Not often used:
      if (prtlev.gt.0.or.nomlst(1:lnomlst).ne.'NONE')
     .write(uninp,'(3A)') ' NEIGHbor-list FILENAME     = ',
     .      nomlst(1:lnomlst),
     .  ' ! For defining the list of interacting pairs yourself.'
      if (prtlev.gt.0.or.nommls(1:lnommls).ne.'NONE') 
     .write(uninp,'(3A)') ' MOLScript command FILEname = ',
     .      nommls(1:lnommls),
     .  ' ! To draw the network with Molscript.'
      if (prtlev.gt.0.or.nomvmd(1:lnomvmd).ne.'NONE') 
     .write(uninp,'(3A)') ' VMD command FILEname       = ',
     .      nomvmd(1:lnomvmd),
     .  ' ! vmd -e this-file (to visualize the network with VMD).'
      close(uninp)

c     Tests:
      if ((.not.qlist.and.cutoff.lt.0.d0.and.rkh.lt.0.d0).or.
     .    (.not.qlist.and.knonb.le.0.d0).or.
     .   (typbond(1:4).ne.'NONE'.and.(cutbnd.lt.0.d0.or.kbond.lt.0.d0)))
     .    then
          write(6,'(/2A)') progrer,
     .  ' Distances and force constants can not have negative values !' 
          stop '*Commands are not consistent*'
      endif

c     On recherche l'information/sous-unite:
 
      call string_split(nompdb,lnompdb,":",
     .                  mots,nmotsmax,nmots)
      call stringcl(mots(1),lnom)
 
      if (nmots.gt.1) then
          call stringcl(mots(2),lnom)
          ssusel=mots(nmots)
          write(6,'(3A)') program,' Subunit to be selected: ',ssusel
          if (nmots.gt.2) then
              write(6,'(4A)') progrwn,' The end of filename, ',
     .        nompdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel=' '
      endif
      nompdb=mots(1)
      call stringcl(nompdb,lnompdb)
c                                          
c     Lecture du fichier de coordonnees:
c     ==================================

      if (prtlev.gt.0)
     .  write(6,'(/(4A))') program,
     .' Coordinate file ',nompdb(1:lnompdb),' to be opened.'

      unpdb=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam(nompdb,lnompdb,cformat,cstatus,unpdb,.true.,
     .     qinter,qexist)
      if (qinter) stop '*No readable coordinate file found*'

c     Format pdb ?

      nl=0
      qpdb=.false.
 120  continue
      read(unpdb,'(A)',end=130) lign80
      if (lign80(1:5).eq.'ATOM '.or.lign80(1:6).eq.'HETATM ') then
          qpdb=.true.
          goto 130
      else
          nl=nl+1
      endif
      goto 120
 130  continue
      rewind(unpdb)

      do i=1,natmax
         xat(i)=unknown
         yat(i)=unknown
         zat(i)=unknown
         massat(i)=unknown
         iresat(i)=i
      enddo

      if (qpdb) then
          write(6,'(/2A)') program,
     .  ' Coordinate file in PDB format.'

          call rdatompdb(unpdb,ssusel,qhet,xat,yat,zat,massat,
     .         atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .         fatres,nresmx,nres,codpdb,qerror,prtlev)

          if (natom.eq.nres) then
              qenm=.true.
              write(6,'(/2A)') program,' Study of a standard ENM model.'
          else
              write(6,'(/2A)') progrwn,
     .      ' Study of a several-atom-per-residue ENM '//
     .      '(this is not that standard).'
          endif
      else
          if (nl.eq.0) then
              write(6,'(/2A)') progrer,' Empty coordinate file.'
              stop
          endif
          write(6,'(/2A)') program,
     .  ' Coordinate file in Free format.'

          call readxyz(unpdb,xat,yat,zat,massat,iresat,natmax,natom,
     .         ndat,qerror,prtlev)

          if (qmasse.and.ndat.lt.4) then
              write(6,'(/2A)') progrer,
     .      ' Masses were not all found, as expected.'
              qmasse=.false.
          endif

c         Par defaut dans ce cas:
          qenm=.true.
          nres=natom
      endif

c     Tests:
c     ======

      if (qerror) stop '*While reading coordinates**'
      if (natom.le.1) then
          write(6,'(2A)') progrer,
     .  ' Not enough atoms found in file. Nothing done.'
          stop
      endif

c     Identification des residus.
c     Recherche des atomes pouvant etre impliques dans des liaisons H.

      nc=0
      no=0
      nn=0
      nh=0
      do i=1,nres
         numc(i)=-1
         numo(i)=-1
         numn(i)=-1
         numh(i)=-1
      enddo

      if (qpdb) then
      nunknown=0
      do i=1,nres
         ires=fatres(i)
         idres(i)=index(residus_standards,resnam(ires))/4
         if (idres(i).le.0) then
             nunknown=nunknown+1
             if (nunknown.lt.10) then
                 write(6,'(4A)') progrwn," residue:'",
     .           resnam(ires),"' is not a well known amino-acid."
                 idres(i)=idmax
             else if (nunknown.eq.10) then
                 write(6,'(2A)') progrwn,' ........'
                 idres(i)=idmax
             endif
         endif
         if (.not.qenm.and.fatres(i+1)-1.ge.ires) then
         do j=ires,fatres(i+1)-1
            if (atonam(j).eq.'C') then
                nc=nc+1
                numc(i)=j
            elseif (atonam(j).eq.'O') then
                no=no+1
                numo(i)=j
            elseif (atonam(j).eq.'N') then
                nn=nn+1
                numn(i)=j
            elseif (atonam(j).eq.'H') then
                nh=nh+1
                numh(i)=j
            endif
         enddo
         endif
      enddo
      if (nunknown.gt.0) 
     .write(6,'(/A,I6,A)') progrwn,nunknown,
     .    ' unexpected amino-acid residue name(s).'
      endif

      if (.not.qenm.and.(nc.ne.nn.or.prtlev.gt.0)) then
          write(6,'(/A,I6,A)') program,nc,' C-atoms found.'
          write(6,'(A,I6,A)') program,no,' O-atoms found.'
          write(6,'(A,I6,A)') program,nn,' N-atoms found.'
          write(6,'(A,I6,A)') program,nh,' H-atoms found.'
      endif
      if (.not.qenm.and.(nc.eq.0.or.no.eq.0.or.nn.eq.0)) then
          write(6,'(/2A)') progrwn,' No complete peptidic bond found ?'
          qenm=.true.
      endif

CoM   Hydrogen bond detection: through distance-cutoff criteria.
CoM   Requires all-atom models, but not the hydrogens. 

      npair=0
      nhb=0
      if (qenm.and.khb.gt.0.d0) then
          write(6,'(/2A)') progrwn,
     .  ' HYDROgen-bond term not allowed in the case of standard ENMs.' 
          khb=0.d0
      endif
      if (khb.gt.0.d0) then
          do i=1,nres
          if (numc(i).gt.0.and.numo(i).gt.0) then
             do j=1,nres
c            Avec plusieurs chaines, on peut en rater...
             if (abs(j-i).gt.1.and.numn(j).gt.0) then
                ii=numc(i)
                jj=numn(j)
                rx=xat(ii)-xat(jj)
                ry=yat(ii)-yat(jj)
                rz=zat(ii)-zat(jj)
                dist=dsqrt(rx*rx + ry*ry + rz*rz)
                if (dist.gt.cuthb) goto 140
                dcn=dist
                ii=numo(i)
                rx=xat(ii)-xat(jj)
                ry=yat(ii)-yat(jj)
                rz=zat(ii)-zat(jj)
                dist=dsqrt(rx*rx + ry*ry + rz*rz)
                if (dist.lt.cuthb-dco) then
                if (prtlev.gt.0)
     .          write(6,'(2A,I6,A,I6)') program,
     .        ' H-bond between residues ',i,' and ',j
                if (prtlev.gt.1) then
                    write(6,'(2A,F8.2)') program,' d(C-N)= ',dcn
                    write(6,'(2A,F8.2)') program,' d(O-N)= ',dist
                endif
                if (nhb+4.le.ntopmax) then
                    npair=npair+1
                    nhb=nhb+1
                    iathb(nhb)=numc(i)
                    jathb(nhb)=numn(j)
                    nhb=nhb+1
                    iathb(nhb)=numo(i)
                    jathb(nhb)=numn(j)
                    if (numh(j).gt.0) then
                        nhb=nhb+1
                        iathb(nhb)=numo(i)
                        jathb(nhb)=numh(j)
                        nhb=nhb+1
                        iathb(nhb)=numc(i)
                        jathb(nhb)=numh(j)
                    endif
                else
                    write(6,'(/2A,I12)') progrer,
     .            ' Too many H-bonds. Maximum is: ',ntopmax
                    stop '*You have to recompile, Sorry*'
                endif
                endif
 140            continue
             endif
             enddo
          endif
          enddo
          if (nhb.gt.0) then
              write(6,'(/A,I6,A,I6,A)') program,nhb,
     .      ' hydrogen bond terms established between ',
     .        npair,' pairs of residues.'
          else
              write(6,'(/2A)') progrwn,' No hydrogen bond found ?'
              khb=0.d0
          endif
c     Pas un PDB:
      else
          if (ndat.ge.5) then
          nres=1
          fatres(nres)=1
          do i=2,natom
             if (iresat(i).ne.iresat(i-1)) then
                 nres=nres+1
                 if (nres.le.nresmx) fatres(nres)=i
             endif
          enddo
          fatres(nres+1)=natom+1
          if (nres.le.nresmx) then
            write(6,'(A,I8,A)') program,nres,' residues found in file.'
          else
            write(6,'(/A,I8,A)') progrer,nres,' residues found in file.'
            stop '*More than allowed. You have to recompile (sorry).'
          endif
          endif
      endif

      write(6,'(/2A)') program,' Coordinate statistics: '

      call vecstat(xat,natom,rmin,rmax,rave,rdev,.false.)
      write(6,'(4(A,F12.6))') 
     .' <x>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      call vecstat(yat,natom,rmin,rmax,rave,rdev,.false.)
      write(6,'(4(A,F12.6))') 
     .' <y>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      call vecstat(zat,natom,rmin,rmax,rave,rdev,.false.)
      write(6,'(4(A,F12.6))') 
     .' <z>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      if (qmasse) then
          write(6,'(/2A)') program,' Mass statistics: '
          call vecstat(massat,natom,rmin,rmax,rave,rdev,.false.)
          write(6,'(4(A,F12.6))') 
     .  ' <m>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

          if (rmin.le.0.d0) then
              write(6,'(2A)') progrer,
     .      ' Negative or null masses found !'
              qmasse=.false.
          endif
      endif

      if (.not.qmasse) then
          write(6,'(2A,F9.2)') program,
     .  ' Masses are all set to ',defmas
          do i=1,natom
             massat(i)=defmas
          enddo
      endif

c     Constante de force pour un eventuel potentiel on-site:
      nonsit=0
      do i=1,3*natom
         konsite(i)=0.d0
      enddo

CoM   Covalent bond detection: also through a distance criterium.
c     Bonds i-j and j-i are stored, because below the matrix is
c     calculated and saved three lines at a time.
 
      nbonds=0
      nangles=0
      ndihs=0
      if (typbond.ne.' ALL'.and.typbond.ne.'CONSECUTIF') goto 200
 
      if (typbond.eq.' ALL') then
         nbmax=0
         nbmin=999
         imax=-1
         imin=-1
         k=1
         nbond(1)=1
         do i=1,natom
            do j=1,natom
               if (i.ne.j) then
                   rx=xat(i)-xat(j)
                   ry=yat(i)-yat(j)
                   rz=zat(i)-zat(j)
                   dist=dsqrt(rx*rx + ry*ry + rz*rz)
                   if (dist.le.cutbnd) then
                       jbond(k)=j
                       k=k+1
                   endif
               endif
            enddo
            nbond(i+1)=k
            if (nbond(i+1)-nbond(i).gt.nbmax) then
                nbmax=nbond(i+1)-nbond(i)
                imax=i
            endif
            if (nbond(i+1)-nbond(i).lt.nbmin) then
                nbmin=nbond(i+1)-nbond(i)
                imin=i
            endif
            if (k-1.gt.ntopmax) then
                write(6,'(/2A,I12)') progrer,
     .        ' Too many bonds. Maximum is: ',ntopmax
                stop '*You have to recompile, Sorry*'
            endif
         enddo
         nbonds=k-1
      else if (typbond.eq.'CONSECUTIF') then
 
CoM      Even when the CONsecutif keyword is used.
c        On fait attention aux distances...
c        Il peut y avoir plusieurs molecules,
c        plusieurs chaines, dans le systeme.
 
         k=1
         do i=1,natom
            nbond(i)=k
            if (i.gt.1) then
            j=i-1
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist=dsqrt(rx*rx + ry*ry + rz*rz)
            if (dist.le.cutbnd) then
                jbond(k)=j
                k=k+1
            endif
            endif
            if (i.lt.natom) then
            j=i+1
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist=dsqrt(rx*rx + ry*ry + rz*rz)
            if (dist.le.cutbnd) then
                jbond(k)=j
                k=k+1
            endif
            endif
            if (k.gt.ntopmax) then
                write(6,'(/2A,I12)') progrer,
     .        ' Too many bonds. Maximum is: ',ntopmax
                stop '*You have to recompile, Sorry*'
            endif
         enddo
         nbond(natom+1)=k
         imax=2
         imin=1
         nbmin=1
         nbmax=2
         nbonds=k
      endif
 
      if (nbonds.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond found.'
          goto 200
      endif
 
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  nbonds/2,' covalent bonds, i.e.,',
     .  float(nbonds)/float(2*natom),' per atom.'

      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)') 
     .'         Maximum number found =',nbmax,' for atom ',imax,
     .'         Minimum number found =',nbmin,' for atom ',imin
 
CoM   Bond-ANGLes are allowed only if there are BONDs.
 
      if (kangle.le.0.d0) then
          write(6,'(/2A)') program,' BOND but no ANGLe energy term.'
          goto 200
      endif

      namax=0
      namin=9999
      imax=-1
      imin=-1
      ii=1
      nangle(1)=1
      do i=1,natom
         if (nbond(i+1).gt.nbond(i)) then
         do jj=nbond(i),nbond(i+1)-1
            j=jbond(jj)
            if (nbond(j+1).gt.nbond(j)) then
            do kk=nbond(j),nbond(j+1)-1
               k=jbond(kk)
               if (k.ne.i) then
                   jangle(ii)=k
                   jsomang(ii)=j
                   ii=ii+1
               endif
            enddo
            endif
         enddo
         endif
         nangle(i+1)=ii
         if (nangle(i+1)-nangle(i).gt.namax) then
             namax=nangle(i+1)-nangle(i)
             imax=i
         endif
         if (nangle(i+1)-nangle(i).lt.namin) then
             namin=nangle(i+1)-nangle(i)
             imin=i
         endif
         if (ii.gt.ntopmax) then
             write(6,'(/2A,I12)') progrer,
     .     ' Too many angles. Maximum is: ',ntopmax
             stop '*You have to recompile, Sorry*'
         endif
      enddo
      nangles=ii-1
 
      if (nangles.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond-angle found.'
          goto 200
      endif
 
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  nangles/2,' valence angles, i.e.,',
     .  float(nangles)/float(2*natom),' per atom.'

      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)') 
     .'         Maximum number found =',namax,' for atom ',imax,
     .'         Minimum number found =',namin,' for atom ',imin

CoM   Bond-DIHEdrals are allowed only if there are BONDs and ANGLes.
 
      if (kdihe.le.0.d0) then
         write(6,'(/2A)') program,' ANGLes but no DIHEdral energy term.'
         goto 200
      endif

      namax=0
      namin=9999
      imax=-1
      imin=-1
      ii=1
      ndihe(1)=1
c     For each atom:
      do i=1,natom
c        All angles where it is the first atom:
         if (nangle(i+1).gt.nangle(i)) then
         do jj=nangle(i),nangle(i+1)-1
c           For the "top-atom" of this angle:
            j=jsomang(jj)
c           All angles where it is the first atom:
            if (nangle(j+1).gt.nangle(j)) then
            do kk=nangle(j),nangle(j+1)-1
               k=jangle(kk)
               if (k.ne.i.and.j.ne.jsomang(kk).and.i.ne.jsomang(kk))then
                   jdihe(ii)=k
                   ii=ii+1
               endif
            enddo
            endif
         enddo
         endif
         ndihe(i+1)=ii
         if (ndihe(i+1)-ndihe(i).gt.namax) then
             namax=ndihe(i+1)-ndihe(i)
             imax=i
         endif
         if (ndihe(i+1)-ndihe(i).lt.namin) then
             namin=ndihe(i+1)-ndihe(i)
             imin=i
         endif
         if (ii.gt.ntopmax) then
             write(6,'(/2A,I12)') progrer,
     .     ' Too many dihes. Maximum is: ',ntopmax
             stop '*You have to recompile, Sorry*'
         endif
      enddo
      ndihs=ii-1
 
      if (ndihs.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond-dihedral found.'
          goto 200
      endif
 
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  ndihs/2,' dihedral angles, i.e.,',
     .  float(ndihs)/float(2*natom),' per atom.'

      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)') 
     .'         Maximum number found =',namax,' for atom ',imax,
     .'         Minimum number found =',namin,' for atom ',imin

c     ==============================
c     Matrice des derivees secondes:
c     ==============================                                
 200  continue
CoM
CoM   Neighbor list:
CoM   --------------
CoM   Allowed file format:
CoM   atom-number atom-number 
CoM   atom-number atom-number force-constant
CoM   For on-site potentials:
CoM   atom-number x-force-constant y-force-constant z-force-constant

      if (qlist) then
      write(6,'(/2A)') program,' Neighbor list to be read.'
      unlst=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam(nomlst,lnomlst,cformat,cstatus,unlst,.true.,
     .     qinter,qexist)
      if (qinter.or..not.qexist) stop '*Required*'

      nvois=0
      nkij=0
 210  continue
      read(unlst,'(A)',end=230,err=220) lign80

c     Les commentaires ne sont pas pris en compte:
      kcom=index(lign80,'!') 
      if (kcom.le.0) kcom=index(lign80,'#')
      if (kcom.gt.0) then
          if (kcom.eq.1) then
              write(6,'(2A)') lign80(1:50),' ...'
              goto 210
          else
              lign80=lign80(1:kcom-1)
          endif
      endif

c     Plusieurs separateurs sont possibles: , ; ou blanc
      csep=","
      k=index(lign80,csep)
      if (k.le.0) then
          csep=";"
          k=index(lign80,csep)
      endif
      if (k.le.0) csep=" "
      
      qfread=.false.
      qonsit=.false.

      call string_split(lign80,80,csep,
     .     mots,nmotsmax,nmots)

      if (nmots.lt.2) then 
          goto 225 
      else if (nmots.eq.3) then
          qfread=.true.
      else if (nmots.eq.4) then
          qonsit=.true.
      else if (nmots.ne.2) then
          goto 225 
      endif
       
c     Lecture de: i, j (kij, le cas echeant).

      read(mots(1),*,end=225,err=225) ii
      jj=-1
      if (.not.qonsit) 
     .read(mots(2),*,end=225,err=225) jj

      if (qfread) then
          read(mots(3),*,end=215,err=215) kfce
          nkij=nkij+1
          if (nkij.eq.1.or.kfce.lt.kfmin) kfmin=kfce
          if (nkij.eq.1.or.kfce.gt.kfmax) kfmax=kfce
          goto 217
 215      continue
ChnG      A compter:
          qfread=.false.
 217      continue
      endif

c     Valeur par defaut:
      if (.not.qfread) kfce=knonb

      if (ii.gt.natom.or.jj.gt.natom) then
          write(6,'(/2A,I6,A,I6,A)') progrer,
     .  ' Atom number: ',ii,' or ',jj,
     .  ' larger than the number of atoms.'
          goto 225
      endif
      if (ii.le.0.or.(jj.le.0.and..not.qonsit)) then
          write(6,'(/2A)') progrer,
     .  ' Null or negative atom number found.'
          goto 225
      endif
      if (ii.eq.jj) then
          write(6,'(/2A,I6,A)') progrwn,
     .  ' Atom ',ii,' is bonded to itself ? (ignored)'
          goto 210
      endif

      if (qonsit) then
          qok=konsite(3*ii-2).eq.0.d0.and.konsite(3*ii-1).eq.0.d0.and.
     .        konsite(3*ii).eq.0.d0
          if (.not.qok) write(6,'(2A,I6)') 
     .        progrwn,' On-site potential defined twice for atom ',ii

          read(mots(2),*,end=227,err=227) kx
          read(mots(3),*,end=227,err=227) ky
          read(mots(4),*,end=227,err=227) kz

          if (kx.lt.0.d0.or.ky.lt.0.d0.or.kz.lt.0.d0) then
              write(6,'(2A,I6)') progrer,
     .      ' Negative on-site force constant for atom ',ii
              goto 227
          else if (kx.gt.0.d0.or.ky.gt.0.d0.or.kz.gt.0.d0) then
              nonsit=nonsit+1
              konsite(3*ii-2)=2.d0*kx
              konsite(3*ii-1)=2.d0*ky
              konsite(3*ii)=2.d0*kz
          endif
          goto 210
      endif

      nvois=nvois+1
      if (nvois.gt.nvoismx) then
          write(6,'(/2A,I6,A)') progrer,' More than ',nvoismx,
     .  ' pairs of neighbors, the maximum allowed. Sorry.'
          stop '*Recompile with larger array*'
      endif

      ivois(nvois)=ii
      jvois(nvois)=jj
      fvois(nvois)=kfce

c     Ligne suivante:
      goto 210

c     Probleme de lecture:
 220  continue
      write(6,'(2A,I6)') progrer,
     .    ' While reading neigbors pair: ',nvois+1
      stop '*Wrong or corrupted file*'

 225  continue
      write(6,'(2A/A)') progrer,
     .' No (or wrong) pair of atom numbers found in ligne: ',lign80
      stop '*Wrong or corrupted file*'

 227  continue
      write(6,'(2A/A)') progrer,
     .' While reading force constant(s) in ligne: ',lign80
      stop '*Wrong or corrupted file*'

c     Fin du fichier:
 230  continue     
c     En cas de fichier "variable":
      if (nkij.gt.0) qfread=.true.
      if (qfread) then
          write(6,'(2(/A,I6,A),1PG10.4,A,1PG10.4)') 
     .    program,nvois,' pairs of neighbors.',
     .    program,nkij,' force constants were also read,'//
     .  ' ranging between ',kfmin,' and ',kfmax
          if (kfmin.lt.0) stop '*Negative values found (wrong file ?)*'
      else
          write(6,'(/A,I6,A)') program,nvois,' pairs of neighbors.' 
          if (knonb.lt.0) then
          write(6,'(/2A)') progrer,' No value read for FORce constants.'
          stop '*Required somewhere*'
          endif
      endif
      if (nonsit.gt.0) 
     .write(6,'(A,I6,A)') program,nonsit,' on-site terms.'
c     qlist:
      endif
CoM
CoM   Coordinates and masses used can be saved, in FREE format.
CoM   This depends upon the PRINTing level (must be > 0).             

      if (prtlev.gt.0) then
      nomfich="pdbmat.xyzm"
      lnom=11
      if (nomfich(1:lnom).ne.nompdb(1:lnompdb)) then
      unout=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,lnom,cformat,cstatus,unout,.true.,
     .     qinter,qexist)

      write(unout,'(A,4(6X,A,11X),2X,A)') 
     .  '!',' x ',' y ',' z ',' m ',' block '
      do i=1,natom
         write(unout,'(4(1PG20.12),I9)')  
     .   xat(i), yat(i), zat(i), massat(i), iresat(i)
      enddo
      close(unout) 
      write(6,'(2A)') program,
     .    ' Coordinates and masses considered are saved.'
      else
      write(6,'(/3A)') progrwn,
     .    ' Coordinates and masses considered are NOT saved '//
     .    '(Coordinate FILENAME is the default one for output).'
      endif
      endif

c     Fichier de commandes pour VMD:
c     ------------------------------

      unvmd=-1
      if (nomvmd.ne.'NONE') then
      unvmd=nunit
      nunit=nunit+1
      nomfich=nomvmd
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,lnomvmd,cformat,cstatus,unvmd,.true.,
     .     qinter,qexist)

      write(unvmd,'(A)') '#!/usr/local/bin/vmd'
      write(unvmd,'(A)') '# script for VMD (Visual Molecular Dynamics)'
      write(unvmd,'(A)') '# Goal: visualizing the elastic network'
      write(unvmd,'(A)') '# Type: vmd -e this-file'
      write(unvmd,'(A)') 'color Display {Background} white'
      write(unvmd,'(A)') 'mol new'
      write(unvmd,'(A)') 'draw color black'
      endif

c     Fichier de commandes pour Molscript:
c     ------------------------------------

      unmol=-1
      if (nommls.ne.'NONE') then
      unmol=nunit
      nunit=nunit+1
      nomfich=nommls
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,lnommls,cformat,cstatus,unmol,.true.,
     .     qinter,qexist)

      nomout='"'//nompdb(1:lnompdb)//'"'
      write(unmol,'(A)') 
     .   '! Script for Molscript (Kraulis, 1993)',
     .   '! Goal: visualizing the elastic network',
     .   ' plot ',
     .   ' frame off ;'
      write(unmol,'(3A)') 
     .   ' read mol1 ',nomout(1:lnompdb+2),' ;'
      write(unmol,'(A)') 
     .   ' set bonddistance 99.0 ;'
      endif

c     Matrice:
c     --------

      if (qbinary) then
        if (.not.qmtx) nommtx="matrix.sdijb"
        cformat="UNFORMATTED"
      else
        if (.not.qmtx) nommtx="matrix.sdijf"
        cformat="FORMATTED"
      endif
      if (.not.qmtx) lnommtx=12
      unout=nunit
      nunit=nunit+1
      cstatus="ove"
      call openam(nommtx,lnommtx,cformat,cstatus,unout,.true.,
     .     qinter,qexist)

c     ========================================
c     Les atomes sont tous lies deux a deux,
c     par un potentiel "universel" (M.Tirion).
c     ========================================
 
      elmax=0.d0
      trace=0.d0
      dmin=0.d0
      dmax=0.d0
      dmoy=0.d0
      drms=0.d0
      nnzero=0
      nntr=0
      nbig=0
      ll=0

      do i=1,natom
         ii=3*i-2
         nvoisat(i)=0

c        Liste eventuelle des voisins de i:
c        ----------------------------------
         if (qlist) then
             do j=1,natom
                qvois(j)=.false.
                kvois(j)=0.d0
             enddo
             do j=1,nvois
                if (ivois(j).eq.i) then
                    qvois(jvois(j))=.true.
                    if (qfread) kvois(jvois(j))=fvois(j)
                endif
                if (jvois(j).eq.i) then
                    qvois(ivois(j))=.true.
                    if (qfread) kvois(ivois(j))=fvois(j)
                endif
             enddo
         endif
 
c        On calcule trois lignes de la matrice a la fois:
c        ================================================

         do j=1,3*natom
            der2(1,j)=0.d0
            der2(2,j)=0.d0
            der2(3,j)=0.d0
         enddo
 
         do j=1,natom
            if (.not.qlist.or.(qlist.and.qvois(j))) then
            if (i.ne.j) then
            jj=3*j-2
            kij=knonb
            if (qlist.and.qfread) kij=kvois(j) 

            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist2=rx*rx + ry*ry + rz*rz
            dist=dsqrt(dist2)
 
            if (dist.lt.rsmall) then
                write(6,'(/2A,1PG10.4,A/2(I6,2A,I6,A,1X,2A))') 
     .          progrer,' Too small distance = ',dist,
     .        ' between following atoms.',
     .          i,': ',resnam(i),iresat(i),ssunam(i),atonam(i),' and ',
     .          j,': ',resnam(j),iresat(j),ssunam(j),atonam(j)
                stop '*Wrong coordinates*'
            endif
CoM
CoM         Hinsen's version: force constant * exp(-(dist/rkh)**2.d0),
CoM         except for topological bond force constants, e.g. ANGLes.
            if (rkh.gt.0.d0) kij=kij*exp(-(dist/rkh)**2.d0)

c           Constantes de force chimiques.
c           ------------------------------
c           1) i et j peuvent etre en interaction 1-2:
 
            qchim=.false.
            if (nbonds.gt.0) then
            if (nbond(i+1).gt.nbond(i)) then
                do k=nbond(i),nbond(i+1)-1
                   if (jbond(k).eq.j) then
                       kij=kbond
                       qchim=.true.
                       goto 300
                   endif
                enddo
            endif
            else
c           Pas de bond: pas d'angle.
            goto 300
            endif
 
c           2) i et j peuvent etre en interaction 1-3:

            if (nangles.gt.0) then
            if (nangle(i+1).gt.nangle(i)) then
                do k=nangle(i),nangle(i+1)-1
                   if (jangle(k).eq.j) then
                       kij=kangle
                       qchim=.true.
                       goto 300
                   endif
                enddo
            endif
            else
c           Pas d'angle: pas de diedre.
            goto 300
            endif

c           3) i et j peuvent etre en interaction 1-4:

            if (ndihs.gt.0) then
            if (ndihe(i+1).gt.ndihe(i)) then
                do k=ndihe(i),ndihe(i+1)-1
                   if (jdihe(k).eq.j) then
                       kij=kdihe
                       qchim=.true.
                       goto 300
                   endif
                enddo
            endif
            endif
 300        continue

c           4) i et j peuvent faire partie d'une liaison hydrogene:

            if (nhb.gt.0) then
            do k=1,nhb
               if ((iathb(k).eq.i.and.jathb(k).eq.j).or.
     .             (iathb(k).eq.j.and.jathb(k).eq.i)) then
                   kij=khb
                   qchim=.true.
                   goto 310
               endif
            enddo
            endif
 310        continue

c           Calcul des elements: (potentiel harmonique)
c           -------------------------------------------
            if (dist.le.cutoff.or.qlist.or.qchim.or.
     .         (cutoff.le.0.d0.and.rkh.gt.0.d0)) then

                ll=ll+1
                nvoisat(i)=nvoisat(i)+1

                if (j.gt.i) then
                   if (unvmd.gt.0) then
                   write(unvmd,'(A,3F12.4,A,3F12.4,A)') 'draw line {',
     .             xat(i),yat(i),zat(i),'} {',xat(j),yat(j),zat(j),'}'
                   endif
                   if (unmol.gt.0) then
                   write(unmol,'(A,I6,3A)')  
     .           ' bonds require in residue ',iresat(i),', atom ',
     .             atonam(i),' and in molecule mol1 ',
     .           ' require in residue ',iresat(j),', atom ',
     .             atonam(j),' and in molecule mol1 ; '
                   endif
c                  Potentiel harmonique: 1/eval*knonb*(d - rval)**eval
                endif

                if (ll.eq.1.or.dist.lt.dmin) dmin=dist
                if (ll.eq.1.or.dist.gt.dmax) dmax=dist
                dmoy=dmoy+dist
                drms=drms+dist2

c               Elements diagonaux des blocs i et j:
c               -----------------------------------
                ddf=kij/dist2
                elemnt=rx*rx*ddf
                der2(1,ii)=der2(1,ii)+elemnt
                der2(1,jj)=der2(1,jj)-elemnt
                elemnt=ry*ry*ddf
                der2(2,ii+1)=der2(2,ii+1)+elemnt
                der2(2,jj+1)=der2(2,jj+1)-elemnt
                elemnt=rz*rz*ddf
                der2(3,ii+2)=der2(3,ii+2)+elemnt
                der2(3,jj+2)=der2(3,jj+2)-elemnt
 
c               Elements extra-diagonaux des deux blocs:
c               ---------------------------------------
                elemnt=rx*ry*ddf
                der2(1,ii+1)=der2(1,ii+1)+elemnt
                der2(2,ii)=der2(2,ii)+elemnt
                der2(1,jj+1)=der2(1,jj+1)-elemnt
                der2(2,jj)=der2(2,jj)-elemnt
                elemnt=rx*rz*ddf
                der2(1,ii+2)=der2(1,ii+2)+elemnt
                der2(3,ii)=der2(3,ii)+elemnt
                der2(1,jj+2)=der2(1,jj+2)-elemnt
                der2(3,jj)=der2(3,jj)-elemnt
                elemnt=ry*rz*ddf
                der2(2,ii+2)=der2(2,ii+2)+elemnt
                der2(3,ii+1)=der2(3,ii+1)+elemnt
                der2(2,jj+2)=der2(2,jj+2)-elemnt
                der2(3,jj+1)=der2(3,jj+1)-elemnt
            endif
            endif
            endif
         enddo
 
c        Sortie de la matrice-bande calculee:
c        -----------------------------------
c       (Uniquement la demi-matrice superieure)
 
c        Level-shift, pour eviter les zeros numeriques,
c        lors de la diagonalisation a venir 
c       (la minimisation est parfaite, par definition).
c        Le hasard est la pour lever la degenerescence
c        des six valeurs propres nulles, et differencier
c        rotations et translations.
CoM
CoM      Level-shift and on-site potentials are mutually exclusive (why not ?).
CoM      The SHIFt is added as a noise (with a maximum value).
CoM      Purpose: to lift the degeneracy of zero-eigenvalue modes.

         If (nonsit.eq.0) then
            der2(1,ii)  =der2(1,ii)   + levelshft*random(iseed)
            der2(2,ii+1)=der2(2,ii+1) + levelshft*random(iseed)
            der2(3,ii+2)=der2(3,ii+2) + levelshft*random(iseed)
         else
c           Ajout du potentiel on-site:
            der2(1,ii)  =der2(1,ii)   + konsite(ii)
            der2(2,ii+1)=der2(2,ii+1) + konsite(ii+1)
            der2(3,ii+2)=der2(3,ii+2) + konsite(ii+2)
         endif
 
         do j=ii,3*natom
            jat=(j-1)/3+1
            if (der2(1,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii,j,der2(1,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)') 
     .          ii,j,der2(1,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(1,j)).gt.rbig)  nbig=nbig+1
                if (dabs(der2(1,j)).gt.elmax) elmax=dabs(der2(1,j))
            endif
         enddo
         do j=ii+1,3*natom
            jat=(j-1)/3+1
            if (der2(2,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii+1,j,der2(2,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)') 
     .          ii+1,j,der2(2,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(2,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(2,j)).gt.elmax) elmax=dabs(der2(2,j))
            endif
         enddo
         do j=ii+2,3*natom
            jat=(j-1)/3+1
            if (der2(3,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii+2,j,der2(3,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)') 
     .          ii+2,j,der2(3,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(3,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(3,j)).gt.elmax) elmax=dabs(der2(3,j))
            endif
         enddo
         elemnt=(der2(1,ii)+der2(2,ii+1)+der2(3,ii+2))/massat(i)
         if (elemnt.eq.0.d0) then
             write(6,'(2A,I6,A)') progrwn,
     .     ' Atom ',i,' has a null second derivatives...'
         else
             nntr=nntr+1
         endif
         trace=trace+elemnt
      enddo
      close(unout)

      if (unvmd.gt.0) then
          write(unvmd,'(2A)') 'mol load pdb ',nompdb(1:lnompdb)
          write(6,'(/2A)') program,' VMD command file saved.'
          close(unvmd)
      endif
      if (unmol.gt.0) then
          write(unmol,'(A)') ' end_plot'
          write(6,'(/2A)') program,' MOLScript command file saved.'
          close(unmol)
      endif
      if (unrsd.gt.0) close(unrsd)
 
      nmin=natom
      nmoy=0.d0
      nrms=0.d0
      nmax=0
      noff=0
      do i=1,natom
         if (nvoisat(i).gt.nmax) nmax=nvoisat(i)
         if (nvoisat(i).lt.nmin) nmin=nvoisat(i)
         if (nvoisat(i).eq.0) then
         noff=noff+1
         if (noff.eq.1) then 
         write(6,'(/2A,I6,A)') progrwn,' Atom ',i,' has no neighbor.'
         else 
         write(6,'(2A,I6,A)') progrwn,' Atom ',i,' has no neighbor.'
         endif
         endif
         nmoy=nmoy+nvoisat(i)
         nrms=nrms+nvoisat(i)**2.d0
      enddo
      nmoy=nmoy/float(natom)
      nrms=nrms/float(natom)-nmoy**2.d0
      if (nrms.gt.0.d0) nrms=dsqrt(nrms)
 
      if (ll.eq.0) then
          write(6,'(/2A,I12,A)') progrer,
     .  ' No atom-atom interaction found. Too short cutoff ?'
          stop '*Empty matrix*'
      else
          dmoy=dmoy/float(ll)
          drms=drms/float(ll)-dmoy**2.d0
          if (drms.gt.0.d0) drms=dsqrt(drms)
      endif

      if (prtlev.gt.0)
     .write(6,'(/2A)') program,' Matrix statistics:'

      write(6,'(/2A,F8.4,A)') program,' The matrix is ',
     .  100.d0*dfloat(nnzero)/dfloat(3*natom*(3*natom+1)/2),' % Filled.'
      write(6,'(A,I12,A)') program,nnzero,' non-zero elements.'

      if (prtlev.gt.0) then
      write(6,'(A,I12,A)') program,ll/2,' atom-atom interactions.'
      write(6,'(/2A,F9.2,A,F9.2/(A,I6))') program,
     .        ' Number per atom= ',nmoy,' +/- ',nrms,
     .'         Maximum number = ',nmax,
     .'         Minimum number = ',nmin

      if (prtlev.gt.1) then
          do i=1,natmax
             statvois(i)=0
          enddo
          do i=1,natom
ChnG         Pas normal...
             if (nvoisat(i).gt.natmax) nvoisat(i)=natmax
             if (nvoisat(i).gt.0) 
     .       statvois(nvoisat(i))=statvois(nvoisat(i))+1
          enddo
          write(6,'(/2A)') program,' Number of neighbors, count: '
          do i=1,natom
             if (statvois(i).gt.0) write(6,'(10X,2I5)') i, statvois(i)
          enddo
      endif

      write(6,'(/2A,F9.2,A,F9.2/(A,F9.2))') program,
     .        ' Average dist.  = ',dmoy,' +/- ',drms,
     .'         Maximum dist.  = ',dmax,
     .'         Minimum dist.  = ',dmin
      endif
 
      write(6,'(/2A,1PG12.6)') program,' Matrix trace   = ',trace

      if (prtlev.gt.0) then
      write(6,'(2A,1PG12.6)') program,' Larger element = ',elmax
      write(6,'(A,I6,A,1PG8.1)') program,
     .      nbig,' elements larger than +/- ',rbig
      endif
 
      write(6,'(/2A)') program,' Hessian matrix ready.'
      write(6,'(2A)') program,
     .' To diagonalize it (and get the modes),'//
     .' you may use diagstd, diagrtb, blzpack...'

      if (nonsit.eq.0) then
          nzero=6
      else if (nonsit.eq.1) then
          nzero=3
      else if (nonsit.eq.2) then
          nzero=1
      else
          nzero=0
      endif

      if (noff.eq.0) then
          write(6,'(/2A,I1,A)') program,' At least ',
     .    nzero,' zero-eigenvalues expected.'
      else
          write(6,'(/A,I6,A,I6,A)') progrwn,
     .    noff,' atoms without any neighbor: at least ',3*noff+nzero,
     .  ' zero-eigenvalues expected.'
      endif
      if (dabs(levelshft).gt.rsmall) write(6,'(2A)') program,
     .  ' But remember that they were also randomly SHIFted.'

      write(6,'(/2A)') program,' Normal end.'
      stop
      end
c-----------------------------------------------------------------------
      subroutine mintomaj(chaine)
 
c     Les caracteres minuscules sont mis en MAJUSCULES.
c     Les autres ne sont pas touches.
 
c     YHS-Oct-98: Premiere version (Toulouse).
c     YHS-Sep-03: Dernieres modifications (Lyon).
 
      character*(*) chaine
c Local:
      integer icar, ilettre, taille
      character*26  carmaj, carmin
 
      carmin='qwertyuiopasdfghjklzxcvbnm'
      carmaj='QWERTYUIOPASDFGHJKLZXCVBNM'
 
      taille=len(chaine)
      if (taille.le.0) return

      do icar=1,taille
         ilettre=index(carmin,chaine(icar:icar))
         if (ilettre.gt.0) then
             chaine(icar:icar)=carmaj(ilettre:ilettre)
         endif
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine readxyz(uninp,x,y,z,w,ic,nmax,ncoor,ndat,qerror,prtlev)

c     Reads at most NMAX coordinates in free format. 

c     Either:
c     x, y, z
c     or:
c     x, y, z, w
c     or:
c     x, y, z, w, ic

c     If first word in ligne is not a number, the whole ligne is
c     assumed to be a title or a commentary.

c     YHS-Sep-03: First version (Lyon).
c     YHS-Jun-08: v1.1.
cI/O:
      implicit none
      logical qerror
      integer ic(*), ncoor, ndat, nmax, prtlev, uninp
      double precision w(*), x(*), xc, y(*), yc, z(*), zc
cLocal:
      integer nmotsmax
      parameter(nmotsmax=255)

      integer i, k, lmot, nchi, nfound, nlmax, nl, nlu, nmots, 
     .        stats(0:nmotsmax)
      double precision rlu, wc
      character chiffres*15, lignlg*(nmotsmax), 
     .        mots(nmotsmax)*(nmotsmax), program*9, progrer*11, 
     .        progrwn*11

cBegin:
      program=' Readxyz>'
      progrer='%Readxyz-Er>'
      progrwn='%Readxyz-Wn>'

      chiffres='1234567890.eE+-'

      do i=0,nmotsmax
         stats(i)=0
      enddo

c     Lecture ligne a ligne:
c     ----------------------

      if (prtlev.gt.1) write(6,'(/2A)') program,
     .  ' Comments, or lignes with less than three numbers: '

      qerror=.false.
      ncoor=0
      nl=0
 100  continue
      read(uninp,'(A)',end=200) lignlg 
      nl=nl+1
      call string_split(lignlg,nmotsmax," ",mots,nmotsmax,nmots)

      nfound=0
      do i=1,nmots
         call stringcl(mots(i),lmot)
         if (lmot.le.0) goto 150

c        Commentaire ?
         if (mots(i)(1:1).eq.'!') goto 150

c        Chiffre ?

         do k=1,lmot
            if (index(chiffres,mots(i)(k:k)).le.0) goto 110
         enddo

         nfound=nfound+1

         if (nfound.le.4) then
             read(mots(i)(1:lmot),*,err=110) rlu
             if (nfound.eq.1) xc=rlu
             if (nfound.eq.2) yc=rlu
             if (nfound.eq.3) zc=rlu
             if (nfound.eq.4) wc=rlu
         else if (nfound.eq.5) then 
             read(mots(i)(1:lmot),*,err=110) nlu
         endif

c        Mot suivant:
 110     continue

c        Le premier mot n'est pas un chiffre => ligne de commentaires

         if (nfound.eq.0) goto 150
      enddo
 150  continue

c     Stockage des coordonnees:
c     -------------------------

      stats(nfound)=stats(nfound)+1

      if (nfound.ge.3) then
          ncoor=ncoor+1
          if (ncoor.le.nmax) then
              x(ncoor)=xc
              y(ncoor)=yc
              z(ncoor)=zc
              if (nfound.eq.4) then
                  w(ncoor)=wc
              else
                  w(ncoor)=wc
                  ic(ncoor)=nlu 
              endif
          else
              write(6,'(/2A,I9,A)') progrer,' More than ',
     .        nmax,' particles in file.'
              write(6,'(2A)') progrer,
     .      ' Please increase program memory limits (Sorry for that).'
              stop
          endif
      else
          if (prtlev.gt.1) then
              write(6,'(2A)') lignlg(1:72),'...'
          endif
      endif
      
c     Ligne suivante:
      goto 100
 200  continue
      if (prtlev.gt.1.and.nl.eq.ncoor) write(6,'(2A)') program,' None.'

      write(6,'(/2A,I7)') program,
     .' Nb of particles found in file (with x,y,z coordinates): ',ncoor
 
      if (ncoor.eq.0) then
          write(6,'(/2A)') progrer,' No coordinate found in file.'
          qerror=.true.
      endif

      nchi=0
      ndat=0
      nlmax=0
      do i=1,nmotsmax
         if (stats(i).gt.0) nchi=nchi+1
         if (stats(i).gt.nlmax) then
             nlmax=stats(i)
             ndat=i
         endif
      enddo

      do i=0,nmotsmax
         if (stats(i).gt.0.and.(prtlev.gt.1.or.nchi.gt.1)) then
            write(6,'(A,I6,A,I7,A)') program,i,
     .    ' numbers found in ',stats(i),' lignes.'  
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine string_split(chaine,taille,delimiteur,
     .                        souschaine,nbremax,nbre)
c
c     "Chaine" est coupee en "nbre" "souschaine" de part et d'autre du
c     "delimiteur"
c
c     YHS-Sep-1993: Premiere version (Uppsala).
c I/O:
      integer taille, nbremax, nbre
      character*(*) chaine, souschaine(*), delimiteur
c Local:
      integer icar, iprev
c
      nbre=1
      iprev=1
      souschaine(1)=chaine
      do icar=1,taille
         if (chaine(icar:icar).eq.delimiteur) then
            if (icar-1.ge.iprev) then
               souschaine(nbre)=chaine(iprev:icar-1)
               nbre=nbre+1
               if (nbre.le.nbremax) then
                  if (icar+1.le.taille.and.
     .               chaine(icar+1:taille).ne.' ') then
                     souschaine(nbre)=chaine(icar+1:taille) 
                  else
                     nbre=nbre-1
                     return
                  endif
               else
                  write(6,'(A,I6,A/A)') 
     .               ' %String_split-Err: more than ',nbremax,
     .               ' substrings in : ',chaine
                  return
               endif
            endif
            iprev=icar+1
         endif
      enddo
c
      return
      end
C-----------------------------------------------------------------------
      REAL*8 FUNCTION RANDOM(ISEED)
C-----------------------------------------------------------------------
C     RANDOM NUMBER GENERATOR: UNIFORM DISTRIBUTION (0,1)
C     ISEED: SEED FOR GENERATOR. ON THE FIRST CALL THIS HAS TO
C     HAVE A VALUE IN THE EXCLUSIVE RANGE (1, 2147483647)
C     AND WILL BE REPLACED BY A NEW VALUE TO BE USED IN
C     FOLLOWING CALL.
C
C     REF: Lewis, P.A.W., Goodman, A.S. & Miller, J.M. (1969)
C     "Pseudo-random number generator for the System/360", IBM
C     Systems Journal 8, 136.
C
C     This is a "high-quality" machine independent generator.
C     INTEGERS are supposed to be 32 bits or more.
C     The same algorithm is used as the basic IMSL generator.
C
C     Author: Lennart Nilsson
C
      implicit none
      INTEGER ISEED
      REAL*8 DSEED,DIVIS,DENOM,MULTIP
      DATA  DIVIS/2147483647.D0/
      DATA  DENOM /2147483711.D0/
      DATA  MULTIP/16807.D0/
C
      IF(ISEED.LE.1) ISEED=314159
      DSEED=MULTIP*ISEED
      DSEED=MOD(DSEED,DIVIS)
      RANDOM=DSEED/DENOM
      ISEED=DSEED
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine stringcl(chaine,nonblancs)
c
c     Les caracteres "blancs" de la CHAINE sont retires (a gauche et au milieu).
c     L'entier NONBLANCS donne la position du dernier caractere.
c
c     YHS-Jan-1995: Premiere version (Toulouse).
c     YHS-Oct-2000: Derniere modification (Bordeaux).
c I/O:
      integer nonblancs
      character*(*) chaine
c Local:
      integer icar, ncar, taille
c Begin:
      nonblancs=0
      taille=len(chaine)
      if (taille.le.0) return
c
      if (index(chaine(1:taille),' ').le.0) then
          nonblancs=taille
          return
      endif
c
c*****Nettoyage des blancs a gauche.
c     Premier non-blanc:
c
      do icar=1,taille
         if (chaine(icar:icar).ne.' ') goto 150
      enddo
      icar=taille
 150  continue
      chaine=chaine(icar:taille)
c
c*****Nettoyage des blancs au milieu.
c
          icar=1
          ncar=1
 170      continue
          icar=icar+1
          ncar=ncar+1
          if (chaine(icar:icar).eq.' ') then
              chaine=chaine(1:icar-1)//chaine(icar+1:taille) 
              icar=icar-1
          endif
          if (ncar.lt.taille-1) goto 170
c
      nonblancs=index(chaine,' ')-1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine openam(namfil,lnom,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c======================================================================= 

c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT.
 
c======================================================================= 
c     input:
c        namfil: nom du fichier a ouvrir. 
c        lnom: longueur de ce nom.
c       "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.

c     output: 
c        qexist: flag / existence du fichier 
c        qinterr: Pas de nom pour le fichier cherche.
 
c.......................................................................
c     YHS-oct-1993: Premiere version.
c     YHS-jul-2007: v1.5
c.......................................................................
      logical qexist, qinterr, qverbos
      integer lnom, unit
      character cformat*12, cstatus*12, namfil*(*)
c Local
      character ordrunix*132
c Begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
 
      qinterr=.false.
      qexist=.false.
 
      if (namfil.eq.' ') then 
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
 
      if (namfil.eq.'stop'.or.namfil.eq.'end'.or. 
     .    namfil.eq.'fin'.or.namfil.eq.'quit'.or.
     .    namfil.eq.'STOP'.or.namfil.eq.'END'.or.                     
     .    namfil.eq.'FIN'.or.namfil.eq.'QUIT') then 
         write(6,'(A)') 'Openam> Program stopped on user request.'
         stop                                                                   
      endif 
 
c     Checks if filename is consistent with the opening:
 
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          write(6,'(/3A)') '%Openam-Err> File <',
     .         namfil(1:lnom),'> not found.'
          return
      endif
 
      if (qexist.and.cstatus.eq.'NEW') then
         write(6,'(/A)') 
     .      '%Openam-Err> This file exists:',namfil(1:lnom)
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil(1:lnom)
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
                                                                    
      if (qverbos) then
         write(6,'(/A,I2,2A)')
     . ' Openam> File on opening on unit ',unit,': ',namfil(1:lnom)
      endif
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)                
         
      return                                                                       
      end
c-----------------------------------------------------------------------
      subroutine rdatompdb(unpdb,ssusel,qhet,xat,yat,zat,binfo,
     .           atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .           fatres,nresmx,nres,codpdb,qerror,prtlev)
CoS=====================================================================
CoS...Rdatompdb: reads a PDB coordinate file ligne by ligne.
CoS---------------------------------------------------------------------
CoS   Lignes read:
CoS   Only those starting by 'ATOM' (also HETATM if qhet=T).
CoS   Only those from the chosen subunit.
CoS   First MODEL only.
CoS
CoS   Gets PDB code, expected in columns 63:66 from HEADer ligne.
CoS
CoS   A negative sign can be given to a residue identifier.
CoS   This may happen when it can NOT be read as an integer number.
CoS   Purpose: keeping the information that the current residue is
CoS   different from the previously read one.
CoS   Ex: 139A, 139B, 139C become: 139, -139, 139.
CoS---------------------------------------------------------------------
c     fatres(i): numero du premier atome du residu i.
c-----------------------------------------------------------------------
c     YHS-Nov-1996: Premiere version (Toulouse).
c     YHS-Jan-2007: Version 1.32 (Lyon).
c     YHS-Jan-2011: Version 1.42 (Nantes).
 
      implicit none
cI/O:
      integer unpdb, natmax, iresat(*), natom,  lnom,
     .        nresmx, nres, fatres(*), prtlev
      double precision xat(*), yat(*), zat(*), binfo(*)
      logical qerror, qhet
      character atonam(*)*4, codpdb*4, resnam(*)*4, segid(*)*4, 
     .        ssusel*1, ssunam(*)*1
cLocal:
      integer iatom, irs, irsprev, ndiff, nerr, ntit,
     .        i, j, k, ii
      double precision bfact, x, y, z
      character atncur*5, lign80*80, notuse*6, numbers*10,
     .        program*11, progrer*14, progrwn*14, 
     .        ren*4, resid*7, residprev*7, residrd*7, segat*4, ssu*1
cBegin:
      program=' Rdatompdb>'
      if (prtlev.gt.0)
     .write(6,'(/2A)') program,' Reading pdb file.'

      numbers='0123456789'
      progrer='%Rdatompdb-Er>'
      progrwn='%Rdatompdb-Er>'
      qerror=.false.
      nerr=0
 
      codpdb='NONE'
      residprev='X'
      irsprev=-1
      ndiff=0
      ntit=0
      nres=0
      iatom=0
 105  continue   
      read(unpdb,'(A)',end=200,err=110) lign80 
  
      goto 120                                
 110  continue
      nerr=nerr+1                            
 
 120  continue                              
c     Try to catch the PDB code of this structure: 
      if (lign80(1:4).eq.'HEAD') then
c         Usually there:
          codpdb=lign80(63:66)
          call stringcl(codpdb,lnom)
          if (lnom.ne.4) then
c             Often also here:
              codpdb=lign80(73:76)
              call stringcl(codpdb,lnom)
c             Maybe the last word of the ligne ?
              if (lnom.ne.4) then
                  call stringcl(lign80,lnom)
                  codpdb=lign80(lnom-3:lnom)
              endif
          endif
          goto 105
      endif

      if (lign80(1:4).eq.'ATOM'.or.
     .   (qhet.and.lign80(1:6).eq.'HETATM')) then
      segat=' '
c     Residue number (irs) columns are (used to be ?) ill-defined in PDB.
c     Notuse: usually empty columns.
      read(lign80,'(12X,A4,1X,A4,A1,A7,1X,3F8.3,6X,F6.2,A6,A4)',
     .     end=195,err=195) 
     .     atncur, ren, ssu, resid, x, y, z, bfact, notuse, segat 

      if (iatom.lt.natmax) then
          if (ssu.eq.ssusel.or.ssusel.eq.' ') then
          iatom=iatom+1
          xat(iatom)=x
          yat(iatom)=y
          zat(iatom)=z
          binfo(iatom)=bfact
 
          call stringcl(atncur,lnom)
          atonam(iatom)=atncur
          call stringcl(ren,lnom)
          resnam(iatom)=ren
          ssunam(iatom)=ssu
          segid(iatom)=segat
          residrd=resid
          call stringcl(resid,lnom)

          if (resid.eq.residprev) then
              irs=irsprev
          else
              if (index(numbers,resid(1:1)).gt.0.and.
     .            index(numbers,resid(lnom:lnom)).gt.0) then
                  read(resid,*,end=180,err=180) irs
                  goto 185
 180              continue
                  write(6,'(/3A)') progrer,
     .          ' Wrong residue identifier: ',residrd
                  stop '*Can not handle that one (sorry)*'
 185              continue
              else
                  ndiff=ndiff+1
                  irs=-irsprev
                  if (prtlev.gt.0)
     .            write(6,'(4A)') progrwn,' Residue identifier: ',
     .            resid,' is atypical (not an integer number).'
              endif
          endif
          iresat(iatom)=irs
CoS
CoS       Residue numbers found in PDB are not trusted.               
CoS       Differences in residue identifiers or chain name mark a new residue.
CoS       Differences in segment names are not taken into account:
CoS       The information found in last column is not expected to be standard.

          if (iatom.eq.1.or.resid.ne.residprev.or.
     .        ssunam(iatom).ne.ssunam(iatom-1)) then
              nres=nres+1
              if (nres.gt.nresmx) then
                  write(6,'(/2A/A,I6)') progrer,
     .          ' Too many residues in this file.',
     .          ' Maximum allowed is = ',nresmx
                  stop '*Larger arrays required*'
              endif
              residprev=resid
              irsprev=irs
              fatres(nres)=iatom
          else
              if (resnam(iatom).ne.resnam(iatom-1))
     .        write(6,'(4A)') progrer,
     .      ' Several kinds of residues with id: ',resid,ssunam(iatom)  
          endif
          endif
      else
          write(6,'(/2A/A,I6)') progrer,
     .      ' Too many atoms in this file.',
     .      ' Maximum allowed is = ',natmax
          stop '*Larger arrays required*'
      endif
      else if (lign80(1:6).eq.'REMARK'.and.prtlev.gt.0) then
          ntit=ntit+1
          if (ntit.le.10) then
              write(6,'(A)') lign80
          else if (ntit.eq.11) then
              write(6,'(A)') ' .... '
          endif
      else if (lign80(1:6).eq.'ENDMDL') then
          write(6,'(/2A)') progrwn,
     .  ' ENDMDL encountered. Remaining ignored.'
          goto 200
      endif
 
c     2) Ligne suivante du fichier pdb :
 
      goto 105
 
c     3) Fin de la lecture du fichier pdb :

c     Erreur de lecture:
 195  continue
      qerror=.true.
      write(6,'(/2A/A)') progrer,
     .    ' Unable to read coordinate ligne: ',lign80
 
 200  continue 
      if (prtlev.gt.1) then
      write(6,'(2A)') program,' End of file reached.'
      write(6,'(2A,I6)') program,' Number of I/O errors: ',
     .            nerr
      endif
 
      natom=iatom
      fatres(nres+1)=natom+1
      irs=0
      if (natom.gt.0) irs=iresat(natom)
 
      write(6,'(/2A,I6)') program, 
     .           ' Number of residues found = ',nres 
      write(6,'(A,I6)') 
     .'            First residue number     = ',iresat(1),
     .'            Last  residue number     = ',irs,
     .'            Number of atoms found    = ',natom
      if (prtlev.gt.0.and.nres.gt.0)
     .write(6,'(A,F8.1)') 
     .'            Mean number per residue  = ',float(natom)/float(nres)
      if (ndiff.gt.0)
     .write(6,'(/A,I6)') 
     .'            N.of negative residue id.= ',ndiff
 
      if (natom.eq.0) then
          write(6,'(2A)') progrer,' No atom found in file.'
          qerror=.true.
      endif
      if (nres.eq.0) then
          write(6,'(2A)') progrer,' No residue found in file.'
          qerror=.true.
      endif

      close(unpdb)
      return
      end
c-----------------------------------------------------------------------
      subroutine vecstat(vect,nmax,rmin,rmax,rave,rdev,qprint)
CoS---------------------------------------------------------------------
CoS   Subroutine vecstat:
CoS
CoS   Statistics for REAL vector Vect(NMAX):
CoS   minimum, maximum, average (rave) and standard deviation (rdev).
CoS
c     YHS-Sep-03: First version (Lyon).
c     YHS-Jul-06: version 1.02  (Lyon).
c     YHS-Jun-11: version 1.03  (Nantes).

cI/O:
      integer nmax
      logical qprint
      double precision rave, rdev, rmax, rmin, vect(*)
cLocal:
      integer i
      character program*9, progrer*12, progrwn*12

cBegin:
      program=' Vecstat>'
      progrer='%Vecstat-Er>'
      progrwn='%Vecstat-Wn>'

      rave=0.d0
      rdev=0.d0
      rmin=-9999.d0
      rmax=9999.d0

      if (nmax.le.0) then
          write(6,'(2A)') progrer,' Zero-length vector.'
          return
      endif

      do i=1,nmax
         if (vect(i).gt.rmax.or.i.eq.1) rmax=vect(i)
         if (vect(i).lt.rmin.or.i.eq.1) rmin=vect(i)
         rave=rave+vect(i)
         rdev=rdev+vect(i)**2.0
      enddo

      rave=rave/dfloat(nmax) 
      rdev=rdev/dfloat(nmax)-rave*rave
      if (rdev.gt.0.d0) rdev=dsqrt(rdev)

      if (qprint) then
          write(6,'(/A,2(A,1PG14.6))') program,'  Mean = ',
     .    rave,' +/- ',rdev 
          write(6,'(A,2(A,1PG14.6))') program,'  Mini = ',
     .    rmin,' Mx: ',rmax
      endif

      return
      end
