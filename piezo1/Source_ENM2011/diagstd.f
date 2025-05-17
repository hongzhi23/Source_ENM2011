CoM=====================================================================
CoM...Diagstd: Full diagonalization of a matrix, real, symmetrical.
CoM---------------------------------------------------------------------
CoM
CoM   Diagonalization routine: TQLI (EISPACK).
CoM  -simple, public domain, but slow (ALL eigenvectors are computed).
CoM   
CoM---------------------------------------------------------------------
c     Exemple: 
c     Risoul (2008), matrix order=3N=2700: 35 min. 
CoM
CoM   INPUT matrix filename (expected from, e.g., PDBMAT program):
CoM   ************************************************************
CoM   Default (first one found picked):
CoM   CERFACS -Formatted : matrix.sdijf
CoM   CERFBIN -Binary    : matrix.sdijb
CoM   Otherwise, a formatted matrix filename is asked for.
CoM
CoM   Input matrix format: i, j, non-zero-ij-element.
CoM
CoM   It can start with a title, recognized by: 
CoM   !,# in first column, or 'program-name>' as first word.
CoM   ALL the matrix is put into memory (including zeroes).
CoM   As a consequence, only SMALL systems can be handled.
CoM
CoM   OUTPUT:
CoM   *******
CoM   An eigenvector file, in CERFACS format.
CoM   Its filename is obtained by replacing the matrix filename suffix.
CoM.....................................................................
      program diagstd
      implicit none
      integer natmax, ndim
CoM
CoM   This is a fortran 77 program (sorry), so it has predefined:
CoM   ==============
CoM   MEMORY LIMITS:
CoM   ==============
CoM   Modify them if needed, that is, if the program complains or
CoM   if you are studying (too) large systems. 
CoM
c     NATMAX: Maximum number of particles (atoms) allowed.

      parameter( natmax=2000,
     .           ndim=3*natmax )

CoM   Then, to (re)compile this program, type:
CoM   make diagstd
CoM   or:
CoM   g77 -o diagstd diagstd.f
CoM   or use your favorite fortran compiler instead (of g77).
CoM   
CoM   To run it in the current directory, type: ./diagstd
c-----------------------------------------------------------------------
      logical qcrois, qexist, qinfo, qinterr, qok
      integer cmot, evord(ndim), i, ii, j, jj, k, lmot, lnomeig, lnomf, 
     .        natom, nbig, nmots, nmotsmx, nok, nord, nredond, ntit, 
     .        ntrace, nvec, nvecout, nunit, rdunit, 
     .        unmess, unmodes
      double precision amat(ndim,ndim), ev(ndim), evsort(ndim),
     .       matrd, trace, work(ndim)
      parameter(nmotsmx=80)
      character cformat*20, cstatus*20, eige*4, 
     .       lignrd*50, lign80*80, matrice*20, mots(nmotsmx)*80, 
     .       nomfich*128, nommat*128, 
     .       program*9, progrer*12, progrwn*12, version*32
CoM---------------------------------------------------------------------
CoM   In case of problem(s), feel free to tell:
CoM   Yves-Henri.Sanejouand@univ-nantes.fr (bug reports may help others).
CoM---------------------------------------------------------------------
c     YHS-Sep-2002: Premiere version, from Diagijr v1.13 (Bordeaux).
c     YHS-Mar-2008: v1.10 (Lyon).
ChnG  Les subroutines
      version=' Version 1.13, October 2011.'
 
      program=' Diagstd>'
      progrer='%Diagstd-Er:'
      progrwn='%Diagstd-Wn:'
 
c     Sortie standard: 
      unmess=6
      write(unmess,'(2A)') program,version

CoM   By default, eigenvector are ranked by increasing eigenvalues.
c    (LOWE) or by decreasing values (HIGH).
      eige='LOWE'
      nvec=3*natmax

      nunit=10
      rdunit=nunit
      nunit=nunit+1
 
c     Detection de la matrice d'entree:
c     --------------------------------
      cformat='FORMATTED'
      matrice='CERFACS'
      cstatus='OLD'

      nomfich='matrix.sdijf'
      inquire(file=nomfich,exist=qexist)
      if (qexist) then
          nommat=nomfich
          qok=.true.
          nok=1
      else
          qok=.false.
          nok=0
      endif

      nomfich='matrice.sdijf'
      inquire(file=nomfich,exist=qexist)
      if (qexist) then
          nommat=nomfich
          qok=.true.
          nok=nok+1
      endif

      nomfich='pdbmat.sdijf'
      inquire(file=nomfich,exist=qexist)
      if (qexist) then
          nommat=nomfich
          qok=.true.
          nok=nok+1
      endif

      nomfich='matrix.sdijb'
      inquire(file=nomfich,exist=qexist)
      if (qexist) then
          cformat='UNFORMATTED'
          matrice='CERFBIN'
          nommat=nomfich
          qok=.true.
          nok=nok+1
      endif
      
      nomfich='matrice.sdijb'
      inquire(file=nomfich,exist=qexist)
      if (qexist) then
          cformat='UNFORMATTED'
          matrice='CERFBIN'
          nommat=nomfich
          qok=.true.
          nok=nok+1
      endif

      nomfich='pdbmat.sdijb'
      inquire(file=nomfich,exist=qexist)
      if (qexist) then
          cformat='UNFORMATTED'
          matrice='CERFBIN'
          nommat=nomfich
          qok=.true.
          nok=nok+1
      endif

      if (qok.and.nok.gt.1) then
          write(unmess,'(/A,I1,A)') 
     .    progrer,nok,' possible matrix (default) filenames found.'//
     .  ' Please specify which.' 
          qok=.false.  
      endif

      if (qok) then
         call stringcl(nommat,lnomf)
      else
         call getnam("Matrix filename ? (formatted)",nommat,lnomf,qok)
         if (.not.qok) stop '*Required*'
      endif

      call openam(nommat,cformat,cstatus,rdunit,.true.,
     .     qinterr,qexist)
      if (qinterr.or..not.qexist)
     .    stop '*Matrix file could not be opened*'

      call getrep("Eigenvectors rank: lowest eigenvalues first ? (Y/N)",
     .     qinfo,qok)
      if (.not.qinfo) eige='HIGH'

      call getnum("Number of eigenvectors saved ? (<0: all)",nvecout,
     .    -1,ndim,.false.,qok)
      if (.not.qok) nvecout=-1

      call string_split(nommat,lnomf,'.',mots,nmotsmx,nmots)
      if (nmots.eq.1) then
          nomfich=nommat(1:lnomf)//'.eigenfacs'
      else
          call stringcl(mots(1),lmot)
          nomfich=mots(1)(1:lmot)
          cmot=lmot
          if (nmots.gt.2) then
          do i=2,nmots-1
             call stringcl(mots(i),lmot)
             nomfich(cmot+1:cmot+lmot+1)='.'//mots(i)(1:lmot)
             cmot=cmot+lmot+1
          enddo
          endif
          nomfich=nomfich(1:cmot)//'.eigenfacs'
      endif
      lnomeig=lnomf+10
      inquire(file=nomfich,exist=qexist)
      if (qexist) then
      write(unmess,'(4A)') progrer,' File >',nomfich(1:lnomeig),
     .  '< already exists.'
      stop '*Please remove it yourself*'
      endif

      write(unmess,'(3A)') program,
     .    ' Matrix to be read from file: ',nommat(1:lnomeig)
 
c     ============================================
c     Lecture matrice d'entree (CERFACS, CERFBIN):
c     ============================================
 
      write(unmess,'(/2A)') program,
     .    ' Matrix to be read is in CERFACS Format.'

c     1) La matrice a un titre ?
c     --------------------------
      ntit=0
  55  continue
      if (matrice.eq.'CERFACS') then
          read(rdunit,'(A)',end=60,err=900) lign80 
      else
          read(rdunit,end=60,err=900) lign80 
      endif
      if (lign80(1:1).eq.'!'.or.lign80(1:1).eq.'#') then
          ntit=ntit+1
          write(unmess,'(2A)') lign80(1:50),' ...'
          goto 55
      else
          lignrd=lign80
          call string_split(lign80,80,' ',mots,nmotsmx,nmots)
          call stringcl(mots(1),lmot)
          if (mots(1)(lmot:lmot).eq.'>') then
              ntit=ntit+1
              write(unmess,'(2A)') lignrd,' ...'
              goto 55
          endif
      endif

  60  continue
      rewind(rdunit)
      write(unmess,'(2A,I6,A)') program,' It has ',ntit,
     .    ' title ligne(s).'

c     2) Ordre de la matrice, nombre de lignes:
c     -----------------------------------------
      if (ntit.gt.0) then
          do i=1,ntit
             if (matrice.eq.'CERFACS') then
                 read(rdunit,*)
             else
                 read(rdunit)
             endif
          enddo
      endif

      k=0
      nord=0
  90  continue
      if (matrice.eq.'CERFACS') then
      read(rdunit,*,end=100,err=910) i,j
      else
      read(rdunit,end=100,err=910) i,j
      endif
      k=k+1
      if (i.le.0.or.j.le.0) then
          write(unmess,'(/2A,I9,2(A,I6))')
     .    progrer,' in ligne: ',k,' I= ',i,' J= ',j
          stop
      endif
      if (i.gt.nord) nord=i
      if (j.gt.nord) nord=j
      goto 90
 100  continue
 
      write(unmess,'(/2A,I9)')
     .     program,' Matrix dimension  (Nord)  =',nord
      write(unmess,'(2A,I9)')
     .     program,' Number of non-zero elements',k
      if (k.eq.0) stop '*Not enough*'
 
      natom=nord/3
      if (natom.gt.natmax.or.nord.gt.ndim) then
          write(unmess,'(2A)')
     .    progrer,' Matrix can not be read.'
          if (natom.gt.natmax) write(unmess,'(2(A,I9))')
     .   ' Natom= ',natom,' > natmax= ',natmax
          if (nord.gt.ndim) write(unmess,'(2(A,I9))')
     .   ' Nord=  ',nord,' > Ndim=  ',ndim
          stop
      endif
 
c     3) Lecture de la matrice:
c     -------------------------
      rewind(rdunit)
      if (ntit.gt.0) then
          do i=1,ntit
          if (matrice.eq.'CERFACS') then
              read(rdunit,*)
          else
              read(rdunit)
          endif
          enddo
      endif
 
      nredond=0
      ntrace=0
      trace=0.d0
      nbig=0
      do i=1,nord
        do j=1,nord
         amat(i,j)=0.d0
        enddo
      enddo

      do jj=1,k
         if (matrice.eq.'CERFACS') then
         read(rdunit,*,err=95) i,j,matrd
         else
         read(rdunit,err=95) i,j,matrd
         endif
 
         if (amat(i,j).eq.0.d0) then
             amat(i,j)=matrd
             amat(j,i)=matrd
             if (i.eq.j) then 
                trace=trace+matrd
                ntrace=ntrace+1
             endif
             if (matrd.gt.1E+10) then
                 nbig=nbig+1
                 if (nbig.lt.10) then
                     write(unmess,'(2A,2I12,A,G12.3)') 
     .               progrwn,' Element: ',i,j,' = ',matrd
                 else 
                     if (nbig.eq.10) write(unmess,*) '...'
                 endif
             endif
         else
             nredond=nredond+1
         endif
      enddo
      goto 105
  95  continue
      write(unmess,'(/2A,I6)')
     .     progrer,' while reading ligne ',k
      write(unmess,'(2I6,F16.8)') ' i, j, matrd= ',i,j,matrd
      stop '*Wrong matrix format ? (see documentation)*'
 105  continue
 
      write(unmess,'(2A,I9)') program,
     .    ' Nb of elements found twice:',nredond
      if (nredond.gt.0) 
     .write(unmess,'(2A/)') progrwn,' Ok ?'
      write(unmess,'(2A,I9)') program,
     .    ' Nb of elements    > 1E+10 :',nbig
      if (nbig.gt.0) 
     .write(unmess,'(2A/)') progrwn,' Ok ?'
      write(unmess,'(2A,F31.7)') program,
     .    ' Matrix trace:',trace
      if (nord-ntrace.gt.0)
     .write(unmess,'(2A,I7)') progrwn,
     .    ' Nb of zero elements there:',nord-ntrace
      if (ntrace.eq.0) then
      write(unmess,'(2A,I7)') progrwn,
     .    ' No diagonal element found.'
      write(unmess,'(2A)') progrwn,' Ok ?'
      endif
 
c     Diagonalisation:
c     ----------------
      write(unmess,'(/2A)') program,' Diagonalization.'
 
      if (nvec.gt.nord) nvec=nord
      write(unmess,'(A,I6,A)') program,
     .      nvec,' eigenvectors are about to be computed. '
      if (nvecout.lt.0) nvecout=nvec

c     Initialisations:
      do i=1,ndim
         ev(i)=0.d0
      enddo
 
c     Eigenvalues/Matrix Diagonalization
CoM
CoM   TRED2 and TQLI (based on the original EISPACK library) 
CoM   perform a diagonalization of a real symmetric matrix based 
CoM   on the QL algorithm. 

      CALL TRED2(amat,nord,ndim,ev,work)
      CALL TQLI(ev,work,nord,ndim,amat)
 
      trace=0.d0
      do i=1,nvec
         trace=trace+ev(i)
      enddo
      write(unmess,'(/2A,F21.7)') program,
     .     ' Sum of eigenvalues (trace of eigen-matrix) =',trace
 
c     Trier par ordre croissant ou decroissant:
      qcrois=.true.
      if (eige.eq.'HIGH') qcrois=.false.

      call trier(ev,nvec,ndim,evsort,evord,qcrois)

c     Ecriture des modes normaux au format 'CERFACS':
c     -----------------------------------------------
      if (nvecout.gt.0) then
      write(unmess,'(/2A/(5F15.7))') program,
     .    ' Eigenvalues: ',(ev(evord(i)),i=1,nvecout)

ChnG  Guess if it is...
      WRITE(unmess,'(/2A/(5F15.7))') program,
     .    ' Frequencies (cm-1, '//
     .     'if the matrix is a hessien in CHARMM units):',
     .    (sqrt(dabs(ev(evord(i))))*108.591365,i=1,nvecout)
 
      cformat='FORMATTED'
      cstatus='new'
      unmodes=nunit
      nunit=nunit+1
      call openam(nomfich,cformat,cstatus,
     .     unmodes,.true.,
     .     qinterr,qexist)
 
      do j=1,nvecout
         i=evord(j)
         write(unmodes,'(A,I5,7X,A,1PG12.4)') 
     .       ' VECTOR',j,'VALUE',ev(i)
         write(unmodes,'(1X,35(1H-))') 
         write(unmodes,'(3(1PG12.4))') 
     .        (amat(k,i),k=1,nord)
      enddo
      else
      write(unmess,'(/2A/(5F15.7))') program,
     .    ' Eigenvalues: ',(ev(evord(i)),i=1,nvec)
      endif
 
      write(unmess,'(/2A)')
     .      program,' Normal end.'
 
      stop 
 900  continue
      write(unmess,'(/2A)')
     .      progrer,' I/O error while seeking for title lines.'
      write(unmess,'(A,I6,A)')
     .      progrer,ntit,' title line(s) read.'
      stop '*Wrong matrix file*'
 910  continue
      write(unmess,'(/2A)')
     .      progrer,' I/O error while matrix elements.'
      write(unmess,'(A,I6,A)')
     .      progrer,k+1,' matrix element(s) read.'
      stop '*Wrong matrix format ? (a pair of integers was expected)*'
      end
CoS---------------------------------------------------------------------
      subroutine getnam(message,nomlu,lnomlu,qok)
CoS   GETNAM:
CoS
CoS   Input : a MESSAGE (usually: a question).
CoS   Output: a string (NOMLU) and its length (LNOMLU).
CoS   
c     If something has been read, qok=T.
c     YHS-oct-96
      implicit none
cI/O:
      integer lnomlu
      logical qok
      character*(*) message, nomlu
cLocal:
      integer ntry, ntrymx
cBegin:
c     NTRYMX essais en cas de probleme.
      ntrymx=5
 
      qok=.false.
      ntry=0
 
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
 
      write(6,'(A,A)') ' Getnam> ',message
      read(5,'(A)',end=200,err=100) nomlu
 
      call stringcl(nomlu,lnomlu)
      write(6,'(A,A)') ' Getnam> ',nomlu(1:lnomlu)
 
      if (lnomlu.gt.0) qok=.true.
      return
 200  continue
      return
      end 
c-----------------------------------------------------------------------
      SUBROUTINE openam(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT,
c     a priori suite a une interrogation...
c
c     input:
c        namfil: nom du fichier a ouvrir. 
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output: 
c        qexist: flag / existence du fichier 
c        qinterr: Pas de nom pour le fichier cherche.
c
c     YHS-oct-93
c     YHS-jan-95
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*(*) namfil, cformat, cstatus
c Local
      character*132 ordrunix
c begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
c
      qinterr=.false.
      qexist=.false.
c
      if (namfil.eq.' ') then 
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
c
      if (namfil.eq.'stop'.or.namfil.eq.'end'                         
     &    .or.namfil.eq.'fin'.or.namfil.eq.'quit') then 
         write(*,*) 'Openam> Program is stopping on user request.'
         stop                                                                   
      endif 
c
c     Checks if filename is consistent with the opening:
c
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          if (qverbos) write(6,'(A)') '%Openam-Err> File not found.'
          return
      endif
c
      if (qexist.and.cstatus.eq.'NEW') then
         write(*,'(/A/A)') 
     .      '%Openam-Err> This file exists:',namfil
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
c                                                                   
      if (qverbos) then
         write(*,'(/A,I6,A)')
     .           ' Openam> file on opening on unit ',unit,':'
         write(*,*) namfil
      endif
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)                
c        
      return  
      end
c-----------------------------------------------------------------------
      SUBROUTINE TRED2(A,N,NP,D,E)

c     Reduce the matrix to tridiagonal form.

      integer i, j, k, l, n, np
      double precision A(NP,NP), D(NP), E(NP), f, g, h, hh,
     .  scale

      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.
      E(1)=0.
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE TQLI(D,E,N,NP,Z)

c     Finds the eigenvalues and eigenvectors of a tridiagonal matrix:

      integer i, iter, j, k, l, m, n, np
      double precision b, c, D(NP), dd, E(NP), f, g, p, r, s,
     .       Z(NP,NP)
      
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine trier(y,npoint,nmax,ysort,iord,qcrois)
c
c     Tri par ordre croissant (qcrois=T) ou non.
c     Version triviale...
c     YHS-Jun-2002: Premiere version (Bordeaux).
c     YHS-Sep-2002: Derniere version (Bordeaux).

      implicit none
      logical qcrois
      integer i, icur, iord(*), j, nmax, npoint
      double precision y(*), ycur, ysort(*)
      character progrer*10

      progrer='%Trier-Er>'

      if (npoint.gt.nmax) then
          write(6,'(A,I9,A,I9,A)') progrer,npoint,
     .  ' points to be sorted, i.e., more than ',nmax,' Sorry.'
          stop
      endif

      do i=1,npoint
         ysort(i)=y(i)
         iord(i)=i
      enddo

      do i=1,npoint
        do j=1,npoint
          if (qcrois) then
            if (ysort(i).lt.ysort(j)) then
                ycur=ysort(i)
                icur=iord(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
                iord(i)=iord(j)
                iord(j)=icur
            endif
          else
            if (ysort(i).gt.ysort(j)) then
                ycur=ysort(i)
                icur=iord(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
                iord(i)=iord(j)
                iord(j)=icur
            endif
          endif
        enddo
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
      subroutine getrep(message,qinfo,qok)
c
c     qinfo obtenu en reponse au MESSAGE.
c     NTRYMX essais en cas de probleme.
c     YHS-jan-00
c     YHS-oct-00
c
      implicit none
cI/O:
      logical qok, qinfo
      character*(*) message
cLocal:
      integer ntry, ntrymx
      character*1 cread
cBegin:
      ntrymx=2
c
      qinfo=.false.
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) then
          goto 200
          return
      endif
c
      write(6,'(A,A)') ' Getrep> ',message
      read(5,'(A)',end=200,err=100) cread
c
      if (cread.eq.'T'.or.cread.eq.'t'.or.
     .    cread.eq.'Y'.or.cread.eq.'y'.or.
     .    cread.eq.'O'.or.cread.eq.'o') then
          qinfo=.true.
      else
      if (cread.ne.'F'.and.cread.ne.'f'.and.
     .    cread.ne.'N'.and.cread.ne.'n')
     .    write(6,'(3A)') '%Getrep-W> Unexpected answer:',cread,
     . '. Assumed answer is: NO.'
      endif
c
      write(6,*) 'Getrep> ',qinfo
c
      qok=.true.
      return
 200  continue
      write(6,'(A)') '%Getrep-W> No answer.'//
     .    ' Assumed answer is: NO.'
      return
      end 
c-----------------------------------------------------------------------
      subroutine getnum(message,numlu,nummin,nummax,qcheck,qok)
c
c     NUMLU obtenu en reponse au MESSAGE.
c
c     qcheck=.true. =>
c     NUMLU doit etre inferieur a nummax et superieur a nummin.
c
c     NTRYMX essais en cas de probleme.
c     qok=.false. => Probleme a la lecture.
c
c     YHS-oct-96: version 1.0 (Toulouse).
c     YHS-jan-01: version 4.0 (Bordeaux).
ChnG  Attention aux entiers trop grands...
c
      implicit none
cI/O:
      integer numlu, nummax, nummin
      logical qok, qcheck
      character*(*) message
cLocal:
      integer ntry, ntrymx, iread
cBegin:
      ntrymx=5
c
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
c
      write(6,'(A,A)') ' Getnum> ',message
      read(5,*,end=200,err=100) iread
      numlu=iread
c
      write(6,*) 'Getnum> ',numlu
c
      if (qcheck) then
      if (numlu.gt.nummax) then
          write(6,'(A,I6,A)') 
     .  '%Getnum-Err: Number larger than ',nummax,
     . '. This is not allowed, sorry.'
          numlu=nummax
          return
      else if (numlu.lt.nummin) then
          write(6,'(A,I6,A)') 
     .  '%Getnum-Err: Number smaller than ',nummin,
     . '. This is not allowed, sorry.'
          numlu=nummin
          return
      endif
      endif
c
      qok=.true.
      return
 200  continue
      return
      end 
