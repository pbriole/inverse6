
      PROGRAM INV3

C     PROGRAMME D'INVERSION DE DONNEES DE DEFORMATION DU SOL
C     PIERRE BRIOLE
C     30/8/87
C     Sur sun le 11 Mars 1992
C     15/3/14
C     INVERSE6: "The dirty one" (mixing ascending and descending data), 
C     line 697 the number of ascending pickings

C ------------------------------------------------
C     CE PROGRAMME UTILISE LES VARIABLES SUIVANTES
C ------------------------------------------------

C ****************************************
C     FAILLE: PARAMETRES CARACTERISTIQUES:
C ****************************************

C       XO,YO: PROJECTION EN SURFACE DU CENTRE DE L'ARETE SUPERIEURE DE
C              LA FAILLE
C       PHI: AZIMUT DE LA TRACE EN SURFACE DE LA FAILLE COMPTE POSITIF
C            DANS LE SENS HORAIRE A PARTIR DU NORD
C       H: PROFONDEUR DU CENTRE DE L'ARETE SUPERIEURE DE LA FAILLE
C       L,D: DEMI-LONGUEUR DE LA FAILLE, LARGEUR DE LA FAILLE
C       TETA: PENDAGE DE LA FAILLE
C       STR,DIP,OUV: COMPOSANTES DU DEPLACEMENT RELATIF DES COTES DE LA
C                    FAILLE:
C                           STR: COMPTE POSITIF POUR MOUVEMENT SENESTRE
C                           DIP: COMPTE POSITIF POUR FAILLE NORMALE
C                           OUV: OUVERTURE

C ************************************************
C     POINTS: COORDONNEES X,Y DES POINTS DE CALCUL
C ************************************************

C       POINTV(1): DONNEES DEPLACEMENT VERTICAL (1 DONNEES PAR POINT)
C       POINTV(2): DONNEES DEPLACEMENT HORIZONTAL (2 DONNEES PAR POINT)
C       POINTV(3): DONNEES INCLINAISON (2 DONNEES PAR POINT)

C ***********
C     MESURES
C ***********

C     DEPL: DEPLACEMENTS MESURES ET CALCULES, DE LA FORME
C           DZM,DZC
C           DXM,DXC,DYM,DYC
C           DIXM,DIXC,DIYM,DIYC
C           EN FONCTION DE LA NATURE DES DONNEES
C     IRV,IRH: SI DIFFERENTS DE 0, REFERENCE AU PREMIER POINT DE LA
C              LISTE DES DONNEES DE DEPLACEMENT VERTICAL ET
C              HORIZONTAL

C ****************
C     INCERTITUDES
C ****************

C     ERPA: INCERTITUDE SUR LES PARAMETRES DES FAILLES
C           OPTION MODIFICATION DE CES INCERTITUDES A CHAQUE ITERATION
C     ERD: INCERTITUDE SUR LES MESURES

C ********************
C     UNITES EMPLOYEES
C ********************

C     METRE POUR XO,YO,PHI,L,D,X,Y
C     CENTIMETRES POUR STR,DIP,OUV,DXM,DXC,DYM,DYC,DZM,DZC
C     MICRORADIANS POUR DIXM,DIXC,DIYM,DIYC
C     DEGRES POUR PHI,TETA

C ------------------------------------------------------
C     CE PROGRAMME CONTIENT LES SOUS-PROGRAMMES SUIVANTS
C ------------------------------------------------------

C     INI: PROGRAMME D'ENTREE DES DONNEES
C     SORT: SORTIE SUR FICHIER DES RESULTATS
C     DISLO: CALCULE (DEPL) LES DEPLACEMENTS DES POINTS (POINT) AVEC DES
C            PARAMETRES DE FAILLE DONNES (FAILLE)
C     FAULT: SOUS-PROGRAMME APPELE PAR DISLO
C     DISLOK: SOUS-PROGRAMME APPELE PAR FAULT
C     INVERSE: INVERSE LES PARAMETRES DES FAILLES PAR LA METHODE DE
C              TARANTOLA-VALETTE (FAILLE) ET CALCULE LES DEPLACEMENTS
C              CORRESPONDANTS (DEPL)
C     CHOL: SOUS-PROGRAMME D'INVERSION DE MATRICE

C ----------------------
C     REMARQUES DIVERSES
C ----------------------

C      LE NOMBRE MAXIMUM DE POINTS, FAILLES, INVERSIONS EST REGLABLE
C      PAR PARAMETER
C      IL EST REGLE DANS LA VERSION STANDARD A 200 POINTS, 6 FAILLES ET
C      5 INVERSIONS
C      L'INVERSION DES PARAMETRES EST EFFECTUEE UNIQUEMENT SUR CEUX POUR
C      LESQUELS ERPA EST DIFFERENT DE 0, LES AUTRES PARAMETRES SONT
C      CONSIDERES COMME CONSTANTS, LE PROGRAMME NE FONCTIONNE PAS SI
C      TOUS LES ERPA SONT A 0, MEME POUR LE CALCUL DIRECT

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL*8 (M)
      PARAMETER (NPM=2000,NFM=18,NFMX9=180,NINVM=6)
      CHARACTER*12 POINTV(3),FICHV
      COMMON /DIS/ POINT(NPM,2),DEPL(NPM,4),FAILLE(NFM,10),NV(3),IRV,IRH
      COMMON /INV/ ERPA(NFM,10),ERD(NPM,2),MATV(NFMX9,NFMX9),NPAR,IV,IMP
      COMMON /SOR/ DEPLS(NINVM,NPM,2),FAILS(NINVM,NFMX9)
      COMMON /ES/ POINTV,FICHV

C ----------Initialisation
      CALL INI(NPOINTS,NF,NINV)

C ----------Boucle d'inversion (NINV=0 pour le calcul direct)
      DO 1, I=1,NINV+1
      IF(I.EQ.1) THEN
      CALL DISLO(NPOINTS,1,NF)
      ELSE
      WRITE(*,*) 'Inversion numero ',I-1
      CALL INVERSE(NPOINTS,NF)
C ----------Ouvre un fichier pour MATV lorsque IV different de 0
      IF(IV.NE.0.AND.I.EQ.NINV+1) THEN
C ----------Ne garde que le dernier tableau
      READ(16,'(A)')  FICHV
      OPEN(9,STATUS='NEW',FILE=FICHV)
      DO 3 J=1,NPAR
    3 WRITE(9,11) (MATV(J,L),L=1,NPAR)
   11 FORMAT(20F10.5)
      CLOSE(9)
      ENDIF
      ENDIF

C ----------Remplissage du tableau DEPLS
      DO 2 K=1,NPOINTS
      DEPLS(I,K,1)=DEPL(K,2)
    2 DEPLS(I,K,2)=DEPL(K,4)

C ----------Remplissage du tableau FAILS
      NUM=1
      DO 5 K=1,NF
      DO 5 L=1,10
      FAILS(I,NUM)=FAILLE(K,L)
    5 NUM=NUM+1

    1 CONTINUE

      CALL SORT(NINV,NPOINTS,NF)

      END

C -----------------------------------
C     SOUS-PROGRAMME D'INITIALISATION
C -----------------------------------

      SUBROUTINE INI(K,NFT,NINV)

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL*8 (M)
      PARAMETER (NPM=2000,NFM=18,NFMX9=180)
      CHARACTER*12 POINTV(3),FICHV,PARAM
      COMMON /DIS/ POINT(NPM,2),DEPL(NPM,4),FAILLE(NFM,10),NV(3),IRV,IRH
      COMMON /INV/ ERPA(NFM,10),ERD(NPM,2),MATV(NFMX9,NFMX9),NPAR,IV,IMP
      COMMON /ES/ POINTV,FICHV

C ----------Logo
    5 CONTINUE
C ----------Lecture des parametres sur fichier
      WRITE(*,*) 'Nom du fichier de parametres'
      READ(*,'(A)') PARAM
      OPEN(16,STATUS='OLD',FILE=PARAM)

C ----------Lecture du nombre de failles et d'inversions
C           et de IV (pour avoir les variances a-posteriori)
      READ(16,*) NFT,NINV,IV,IMP

C ----------Introduction des 10 parametres de chaque faille dans
C           l'ordre indique dans les formules d'introduction et sous la
C           forme "P EP" (Parametre et ecart type)

      DO 1, NF=1,NFT
      READ(16,*) (FAILLE(NF,K),ERPA(NF,K),K=1,10)
    1 CONTINUE

C ----------Intoduction des fichiers de donnees
C           Les donnees sont contenues dans des fichiers differents
C           selon leur nature (trois natures differentes possibles)
C           Le programme demande successivement le nom des trois
C           fichiers et l'incertitude sur les donnees associee (on
C           suppose que des donnees de meme nature ont la meme
C           incertitude)
C           Si pas de donnees d'un type appeler le fichier NULL.PTS

      K=1

C ----------Lecture des trois fichiers
      DO 2 I=1,3
C ----------Ouverture d'un fichier bidon NULL.PTS
      OPEN(9,FILE='NULL.PTS')
      CLOSE(9)

C ----------Lecture du nom du fichier
      READ(16,'(A)') POINTV(I)

C ----------Lecture des incertitudes sur les donnees
C           Si les incertitudes sur les donnees varient pour un meme
C           type de donnees, il faut modifier le programme et les
C           introduire dans le fichier des points
      IF(I.EQ.1) THEN
      READ(16,*) ERDA,IRV
      ERDB=0
      ELSE IF(I.EQ.2) THEN
      READ(16,*)ERDA,ERDB,IRH
      ELSE
      READ(16,*) ERDA,ERDB
      ENDIF
 
C ----------Lecture des points sur fichier
      OPEN(7,STATUS='OLD',FILE=POINTV(I))
      DO 3, N=0,NPM
      IF(I.NE.1) THEN
      READ(7,*,END=4) X,Y,DAM,DBM
      ELSE
      READ(7,*,END=4) X,Y,DAM
      DBM=0
      ENDIF
      ERD(K,2)=ERDB
      DEPL(K,3)=DBM
      POINT(K,1)=X
      POINT(K,2)=Y
      DEPL(K,1)=DAM
      ERD(K,1)=ERDA
      K=K+1
    3 CONTINUE

    4 CLOSE(7)

      NV(I)=N
    2 CONTINUE
      K=K-1

C ----------Controle des dimensions
      NMES=NV(1)+2*(NV(2)+NV(3))
      IF(NFT.GT.18.OR.NINV.GT.6) THEN
      WRITE(*,*) 'Trop de failles ou d''inversions'
      STOP
      ELSE IF(NMES.GT.2000) THEN
      WRITE(*,*) 'Trop de points'
      STOP
      ENDIF

      WRITE(*,*) K,'    Points lus'
      WRITE(*,*) NV(1),'    Points de nivellement'
      WRITE(*,*) NV(2),'    Points de mesure de distance'
      WRITE(*,*) NV(3),'    Points de mesure d''inclinaison'
      WRITE(*,*) NFT,'    Failles'
      WRITE(*,*) NINV,'    Inversions'

      RETURN
      END

C ------------------------------
C     SOUS-PROGRAMME D'INVERSION
C ------------------------------

      SUBROUTINE INVERSE(NPOINTS,NF)

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL*8 (M)
      INTEGER M
      PARAMETER (NPM=2000,NFM=18,NFMX9=180,NPMX2=4000)
      COMMON /DIS/ POINT(NPM,2),DEPL(NPM,4),FAILLE(NFM,10),NV(3),IRV,IRH
      COMMON /INV/ ERPA(NFM,10),ERD(NPM,2),MATV(NFMX9,NFMX9),NPAR,IV,IMP

C ----------Matrice des derivees partielles, matrice MAT a inverser et
C           vecteurs
      DIMENSION MATRICE(NPMX2,NFMX9),MAT(NFMX9,NFMX9),X(NFMX9)
      NMES=NV(1)+2*(NV(2)+NV(3))

C ----------Vidage des matrices
      DO 11, I=1,NFMX9
      DO 11, J=1,NFMX9
      MAT(I,J)=0
   11 CONTINUE
      DO 12, I=1,NMES
      DO 12, J=1,NFMX9
      MATRICE(I,J)=0
   12 CONTINUE

C ----------Boucle realisant:
C             remplissage de la matrice
C             remplissage de la premiere partie de MAT (termes lies aux
C             incertitudes sur les parametres)
C             calcul du nombre (NPAR) de parametres a inverser
      NPAR=1
      DO 1, I=1,NF
      DO 1, J=1,10

C ----------Ne calcule pas lorsque ERPA=0 (parametre constant)
      IF(ERPA(I,J).NE.0) THEN
      FAILLE(I,J)=FAILLE(I,J)+0.01
      CALL DISLO(NPOINTS,I,I)
      KB=1

      DO 3, K=1,NPOINTS
      MATRICE(KB,NPAR)=DEPL(K,2)
      KB=KB+1
      IF(K.GT.NV(1)) THEN
      MATRICE(KB,NPAR)=DEPL(K,4)
      KB=KB+1
      ENDIF
    3 CONTINUE

      FAILLE(I,J)=FAILLE(I,J)-0.02
      CALL DISLO(NPOINTS,I,I)
      KB=1

      DO 4, K=1,NPOINTS
      MATRICE(KB,NPAR)=(MATRICE(KB,NPAR)-DEPL(K,2))*100
      KB=KB+1
      IF(K.GT.NV(1)) THEN
      MATRICE(KB,NPAR)=(MATRICE(KB,NPAR)-DEPL(K,4))*100
      KB=KB+1
      ENDIF
    4 CONTINUE

      FAILLE(I,J)=FAILLE(I,J)+0.01
      MAT(NPAR,NPAR)=1/ERPA(I,J)**2
      NPAR=NPAR+1
      ENDIF

    1 CONTINUE
      NPAR=NPAR-1

C ----------Remplissage de la seconde partie de MAT (termes lies aux
C           derivees partielles)
      DO 5, I=1,NPAR
      DO 5, L=1,NPAR
      K=1
      DO 5, J=1,NPOINTS
      MAT(I,L)=MAT(I,L)+MATRICE(K,I)*MATRICE(K,L)/(ERD(J,1)**2)
      K=K+1
      IF(J.GT.NV(1)) THEN
      MAT(I,L)=MAT(I,L)+MATRICE(K,I)*MATRICE(K,L)/(ERD(J,2)**2)
      K=K+1
      ENDIF
    5 CONTINUE

C ----------Calcul facultatif des variances a-posteriori
      IF(IV.NE.0) THEN
      DO 9, J=1,NPAR
      DO 10, I=1,NPAR
      X(I)=0.
   10 CONTINUE
      X(J)=1.
      CALL CHOL(MAT,X,NPAR)
      DO 8, L=1,NPAR
      MATV(L,J)=X(L)
    8 CONTINUE
    9 CONTINUE
      ENDIF

C ----------Calcule les valeurs theoriques des deplacements utilises
C           par l'inversion
      CALL DISLO(NPOINTS,1,NF)

C ----------Calcule le vecteur d'entree
      DO 2, I=1,NPAR
      X(I)=0
      K=1
      DO 2, J=1,NPOINTS
      X(I)=X(I)+MATRICE(K,I)*(DEPL(J,1)-DEPL(J,2))/(ERD(J,1)**2)
      K=K+1
      IF(J.GT.NV(1)) THEN
      X(I)=X(I)+MATRICE(K,I)*(DEPL(J,3)-DEPL(J,4))/(ERD(J,2)**2)
      K=K+1
      ENDIF
    2 CONTINUE

C ----------Calcul du vecteur des variations
      CALL CHOL(MAT,X,NPAR)

C ----------Changement des parametres des failles
      NPAR=1
      DO 17, I=1,NF
      DO 17, J=1,10
      IF(ERPA(I,J).NE.0) THEN
      FAILLE(I,J)=FAILLE(I,J)+X(NPAR)
C ----------Changement facultatif des incertitudes sur les parametres
      IF(IV.NE.0.AND.IMP.NE.0) THEN
      ERPA(I,J)=DSQRT(MATV(NPAR,NPAR))
      ENDIF
      NPAR=NPAR+1
      ENDIF
   17 CONTINUE
      NPAR=NPAR-1

C ----------Controle des nouvelles valeurs
      DO 18 I=1,NF

      IF(FAILLE(I,3).GT.180.) THEN
      FAILLE(I,3)=FAILLE(I,3)-360.
      ENDIF
      IF(FAILLE(I,3).LT.-180.) THEN
      FAILLE(I,3)=FAILLE(I,3)+360.
      ENDIF

      DO 20 J=4,6
      IF(FAILLE(I,J).LT.0.) THEN
      WRITE(*,*) 'Faille nø',I
      IF(J.EQ.4) THEN
      WRITE(*,*) 'Probleme: H negatif - Remis a zero arbitrairement'
      ELSE IF(J.EQ.5) THEN
      WRITE(*,*) 'Probleme: L negatif - Remis a zero arbitrairement'
      ELSE
      WRITE(*,*) 'Probleme: D negatif - Remis a zero arbitrairement'
      ENDIF
      FAILLE(I,J)=0.
      ENDIF
   20 CONTINUE

      FI7=FAILLE(I,7)
      IF(FI7.LT.0.OR.FI7.GT.90.) THEN
      FI9=FAILLE(I,9)
      FI3=FAILLE(I,3)
      FAILLE(I,9)=-FI9
      IF(FAILLE(I,3).GT.0.) THEN
      FAILLE(I,3)=FAILLE(I,3)-180.
      ELSE
      FAILLE(I,3)=FAILLE(I,3)+180.
      ENDIF
      FAILLE(I,7)=180.-FI7
      IF(FI7.LT.0.) THEN
      FAILLE(I,7)=-FI7
      FAILLE(I,4)=FAILLE(I,4)-FAILLE(I,6)*DSIN(FI7)
      DLL=FAILLE(I,6)*DCOS(FI7)
      FAILLE(I,1)=FAILLE(I,1)+DLL*DCOS(FI3)
      FAILLE(I,2)=FAILLE(I,2)-DLL*DSIN(FI3)
      IF(FAILLE(I,4).LT.0.) THEN
      FAILLE(I,4)=0.
      ENDIF
      ENDIF
      ENDIF
      IF(FAILLE(I,10).LT.0.) THEN
      WRITE(*,*) 'Faille nø',I
      WRITE(*,*) 'Probleme: OUV negatif - Remis a zero arbitrairement'
      FAILLE(I,10)=0.
      ENDIF
   18 CONTINUE

C ----------Calcul des deplacements pour la nouvelle faille
      CALL DISLO(NPOINTS,1,NF)

      RETURN
      END

C ------------------------------------------------------------------
C     SOUS-PROGRAMME D'INVERSION DE MATRICE (CALCUL DE X AVEX Y=A*X)
C ------------------------------------------------------------------

      SUBROUTINE CHOL(MAT,X,NPAR)

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL*8 (M)
      PARAMETER (NFMX9=180)
      DIMENSION MAT(NFMX9,NFMX9),X(NFMX9),DIAG(NFMX9)

C ----------Inversion de la matrice MAT par la methode de CHOLEWSKI
C           Le vecteur de sortie est place dans X(I)
      DIAG(1)=DSQRT(MAT(1,1))
      DO 6 I=2,NPAR
    6 MAT(1,I)=MAT(I,1)/DIAG(1)
      DO 10 J=2,NPAR
      S=MAT(J,J)
      KMX=J-1
      DO 7 K=1,KMX
    7 S=S-MAT(K,J)**2
      DIAG(J)=DSQRT(S)
      IF(J.EQ.NPAR) GOTO 10
      IMI=J+1
      DO 8 I=IMI,NPAR
      S=MAT(I,J)
      DO 9 K=1,KMX
    9 S=S-MAT(K,I)*MAT(K,J)
      MAT(J,I)=S/DIAG(J)
    8 CONTINUE
   10 CONTINUE
      X(1)=X(1)/DIAG(1)
      DO 14 I=2,NPAR
      S=X(I)
      JMX=I-1
      DO 13 J=1,JMX
   13 S=S-MAT(J,I)*X(J)
   14 X(I)=S/DIAG(I)
      X(NPAR)=X(NPAR)/DIAG(NPAR)
      DO 16 K=2,NPAR
      I=NPAR+1-K
      JMI=I+1
      S=X(I)
      DO 15 J=JMI,NPAR
   15 S=S-MAT(I,J)*X(J)
   16 X(I)=S/DIAG(I)

      RETURN
      END

C -----------------------------------------------------------
      SUBROUTINE FAULT(XX,YY,XP,YP,AL,W,D,AX,AY,AZ,IND)
C -----------------------------------------------------------

C ----------Calcul du deplacement dans le repere geographique
C           J.C.Ruegg-1986 d'apres Okada(1985,BSSA)
C           Modifie le 30/8/87

C           Faille en tension  I= 1
C           Faille en coulissage  I= 2
C           Faille plongeante  I= 3
C           Coordonnes initiales XX,YY dans repere geographique UTM
C           Repere local de faille X1,Y1
C           Repere Okada  X2,Y2
C           Sortie : Deplacements A,B,C surface libre

      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/COMOB/Q,ST,CT,CTT,SF,CF,ELAS,US(3)
      PI= 3.14159265359

C ----------Passage repere geographique au repere local X1,Y1
      X1= (XX-XP)*SF + (YY-YP)*CF
      Y1=-(XX-XP)*CF + (YY-YP)*SF

C ----------Passage au repere Okada
      X2= X1
      Y2= Y1 + D*CTT

C ----------Coins de la faille dans le repere Okada
      QS1= X2+AL
      QS2= X2-AL
      P1= Y2*CT +D*ST
      P2= P1 - W
      Q= Y2*ST - D*CT

C ----------Boucle de calcul aux quatre coins de la faille
      CALL DISLOK (QS1,P1,UX1,UY1,UZ1,IND)
      CALL DISLOK (QS1,P2,UX2,UY2,UZ2,IND)
      CALL DISLOK (QS2,P1,UX3,UY3,UZ3,IND)
      CALL DISLOK (QS2,P2,UX4,UY4,UZ4,IND)

      A2 = UX1-UX2-UX3+UX4
      B2 = UY1-UY2-UY3+UY4
      C2 = UZ1-UZ2-UZ3+UZ4

C ----------Deformation dans le repere geographique
      AX=(A2*SF - B2*CF)/2./PI
      AY=(A2*CF + B2*SF)/2./PI
      AZ=C2/2./PI

      RETURN
      END

C -------------------------
C     SOUS-PROGRAMME DISLOK
C -------------------------

C ----------Calcul des deplacements ux,uy,uz
C           associes a une dislocation
C           FORMULES DE OKADA (1985, BSSA)
C           D'apres J.C. RUEGG, 1986

      SUBROUTINE DISLOK(QSI,ETA,FUX,FUY,FUZ,IND)

      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/COMOB/Q,ST,CT,CTT,SF,CF,ELAS,US(3)

      BY=ETA*CT+Q*ST
      BD=ETA*ST-Q*CT
      R=DSQRT(QSI**2+ETA**2+Q**2)
      XX=DSQRT(QSI**2+Q**2)

      RE=R+ETA
      RQ=R+QSI
      QRE=Q/R/RE
      QRQ=Q/R/RQ
      ATQE=DATAN(QSI*ETA/Q/R)

      AT= ETA*(XX+Q*CT)+XX*(R+XX)*ST
      BT= QSI*(R+XX)*CT

      FUX=0.
      FUY=0.
      FUZ=0.

      AI4= ELAS*(DLOG(R+BD)-ST*DLOG(R+ETA))/CT
      AI5= 2*ELAS*DATAN(AT/BT)/CT

C ----------Choix entre calcul de donnees verticales ou horizontales
      IF(IND.EQ.1) THEN

C ----------Faille en tension
      IF(US(1).NE.0) THEN
      FUZ=(BY*QRQ+CT*QSI*QRE-CT*ATQE-AI5*ST*ST)*US(1)
      ENDIF

C ----------Faille de coulissage
      IF(US(2).NE.0) THEN
      FUZ=FUZ-(BD*QRE+Q*ST/RE+AI4*ST)*US(2)
      ENDIF

C ----------Faille plongeante
      IF(US(3).NE.0) THEN
      FUZ=FUZ+(-BD*QRQ-ST*ATQE+AI5*ST*CT)*US(3)
      ENDIF

      RETURN

      ELSE

C ----------Donnees horizontales
      AI3=ELAS*(BY/CT/(R+BD)-DLOG(R+ETA))+AI4/CTT
      AI1= -ELAS*QSI/CT/(R+BD)-AI5/CTT
      AI2= -ELAS*DLOG(R+ETA)-AI3

C ----------Faille en tension
      IF(US(1).NE.0) THEN
      FUX=(Q*QRE - AI3*ST*ST)*US(1)
      FUY=-(BD*QRQ+ST*QRE*QSI-ST*ATQE+AI1*ST*ST)*US(1)
      ENDIF

C ----------Faille de coulissage pur (strike-slip)
      IF(US(2).NE.0) THEN
      FUX=FUX-(QSI*QRE+ATQE+AI1*ST)*US(2)
      FUY=FUY-(BY*QRE+Q*CT/RE+AI2*ST)*US(2)
      ENDIF

C ----------Faille plongeante pure ( dip-slip)
      IF(US(3).NE.O) THEN
      FUX=FUX+(-Q/R+AI3*ST*CT)*US(3)
      FUY=FUY-(BY*QRQ+CT*ATQE-AI1*ST*CT)*US(3)
      ENDIF

      ENDIF

      RETURN
      END

C -------------------------------------------
C     SOUS-PROGRAMME DE CALCUL DE DISLOCATION
C -------------------------------------------

      SUBROUTINE DISLO(NPOINTS,NFAIL1,NFAIL2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPM=2000,NFM=18)
      COMMON /DIS/ POINT(NPM,2),DEPL(NPM,4),FAILLE(NFM,10),NV(3),IRV,IRH
      COMMON /COMOB/ Q,ST,CT,CTT,SF,CF,ELAS,US(3)
      PI=3.14159265339
C ----------Le coefficient d'elasticit‚ est choisi egal a 0.5
      ELAS=0.5
      NVH=NV(1)+NV(2)

C ----------Vidage de DEPL
      DO 1, I=1,NPOINTS
      DEPL(I,2)=0
      DEPL(I,4)=0
    1 CONTINUE

C ----------Boucle de calcul par faille
      DO 2, I=NFAIL1,NFAIL2

C ----------Transformation des parametres contenus dans FAILLE aux
C           parametres utilises par le sous-programme FAULT
      FA=FAILLE(I,3)*PI/180
      US(1)=FAILLE(I,10)
      US(2)=FAILLE(I,8)
      US(3)=-FAILLE(I,9)
      TA=FAILLE(I,7)*PI/180
      CT=DCOS(TA)
      ST=DSIN(TA)
      CTT=CT/ST
      CF=DCOS(FA)
      SF=DSIN(FA)
      W=FAILLE(I,6)
      RA=FAILLE(I,5)
      H=FAILLE(I,4)
      D=H+W*ST
      XPA=FAILLE(I,1)-H*CF*CTT
      YPA=FAILLE(I,2)+H*SF*CTT

C ----------Boucle de calcul par point
      DO 2, J=1,NPOINTS
      XX=POINT(J,1)
      YY=POINT(J,2)
      IF(J.LE.NVH) THEN

C ----------Calcul pour donnees verticales ou horizontales
      IF(J.LE.NV(1)) THEN
C ----------Donnee verticale
c      CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,1)
c      DEPL(J,2)=DEPL(J,2)+AZ
c --------Donnees radar

C     This is the number of data in the levelling file
      NIA1=28
c     This is the number of data in the ascending file
      NIA2=75
c     This is the number of data in the descending file
      NIA3=94
      NIA4=NIA1+NIA2

      IF(J.LE.NIA1) THEN
      CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,1)
      DEPL(J,2)=DEPL(J,2)+AZ
      ELSE

      IF(J.LE.NIA4) THEN
C       First incidence angle (here is ascending)
        CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,1)
        DEPL(J,2)=DEPL(J,2)+0.81*AZ
        CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,2)
        DEPL(J,2)=DEPL(J,2)-0.57*AX-0.13*AY
      ELSE
C       Second incidence angle (here is descending)
        CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,1)
        DEPL(J,2)=DEPL(J,2)+0.80*AZ
        CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,2)
        DEPL(J,2)=DEPL(J,2)+0.58*AX-0.14*AY
      ENDIF
      ENDIF





      ELSE
C ----------Donnee horizontale
      CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,2)
      DEPL(J,2)=DEPL(J,2)+AX
      DEPL(J,4)=DEPL(J,4)+AY
      ENDIF

      ELSE

C ----------Calcul pour donnees d'inclinaison (derivation numerique)
      XX1=XX+0.001
      XX2=XX-0.001
      CALL FAULT(XX1,YY,XPA,YPA,RA,W,D,AX,AY,AZ1,1)
      CALL FAULT(XX2,YY,XPA,YPA,RA,W,D,AX,AY,AZ2,1)
      DEPL(J,2)=DEPL(J,2)+5000*(AZ1-AZ2)
      YY1=YY+0.001
      YY2=YY-0.001
      CALL FAULT(XX,YY1,XPA,YPA,RA,W,D,AX,AY,AZ1,1)
      CALL FAULT(XX,YY2,XPA,YPA,RA,W,D,AX,AY,AZ2,1)
      DEPL(J,4)=DEPL(J,4)+5000*(AZ1-AZ2)

      ENDIF

    2 CONTINUE

C ----------Reference eventuelle aux premieres donnees des listes de
C           mesures verticales et horizontales
      IF(IRV.NE.0) THEN
      DO 3, I=1,NV(1)
    3 DEPL(I,2)=DEPL(I,2)-DEPL(1,2)
      ENDIF
      IF(IRH.NE.0) THEN
      NV1=NV(1)+1
      DO 4, I=NV1,NVH
      DEPL(I,2)=DEPL(I,2)-DEPL(NV1,2)
    4 DEPL(I,4)=DEPL(I,4)-DEPL(NV1,4)
      ENDIF

      END

C ----------------------------
C     SOUS PROGRAMME DE SORTIE
C ----------------------------

      SUBROUTINE SORT(NINV,NPOINTS,NF)

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPM=2000,NINVM=6,NFM=18,NFMX9=180)
      CHARACTER*12 DEPL1,FAIL,DIFF
      COMMON /DIS/ POINT(NPM,2),DEPL(NPM,4),FAILLE(NFM,10),NV(3),IRV,IRH
      COMMON /SOR/ DEPLS(NINVM,NPM,2),FAILS(NINVM,NFMX9)

C ----------Tableaux
      READ(16,'(A)') DEPL1,FAIL,DIFF
      OPEN(15,FILE=DEPL1)
      OPEN(7,FILE=FAIL)
      OPEN(9,FILE=DIFF)

C ----------Remplissage de FAIL
      NINV=NINV+1
      DO 2, I=1,NF*10
    2  WRITE(7,11) (FAILS(J,I),J=1,NINV)

C ----------Remplissage de DIFF et de DEPL
      DO 1, I=1,NPOINTS
      WRITE(15,11) (DEPLS(K,I,1),K=1,NINV)
      IF(I.GT.NV(1)) THEN
      WRITE(15,11) (DEPLS(K,I,2),K=1,NINV)
      ENDIF
      X=POINT(I,1)
      Y=POINT(I,2)
      IF(I.LE.NV(1)) THEN
      WRITE(9,12) I,X,Y,DEPL(I,1),DEPL(I,2)
      ELSE
      WRITE(9,12) I,X,Y,DEPL(I,1),DEPL(I,2),DEPL(I,3),DEPL(I,4)
      ENDIF
    1 CONTINUE

      CLOSE(15)
      CLOSE(16)
      CLOSE(7)
      CLOSE(9)
   11 FORMAT(10F10.4)
   12 FORMAT(I5,6F10.2)

      END

   
