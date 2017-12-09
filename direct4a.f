      PROGRAM DIRECT4

C     PROGRAMME DE CALCUL DIRECT D'INTERFEROGRAMMES
C     COUPLE AU PROGRAMME INVERSE4
C     PIERRE BRIOLE
C     16/4/99
C     D'apres le programme INVERSE3 du 30/8/87

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

C ********************
C     UNITES EMPLOYEES
C ********************

C     METRE POUR XO,YO,PHI,L,D,X,Y
C     CENTIMETRES POUR STR,DIP,OUV,DXM,DXC,DYM,DYC,DZM,DZC
C     DEGRES POUR PHI,TETA

C **********************
C     REMARQUES DIVERSES
C **********************

C     LE NOMBRE MAXIMUM DE COLONNES EST ICI DE 3000
C     LE NOMBRE MAXIMUM DE FAILLES EST ICI DE 100

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPM=3000,NFM=100)
      CHARACTER*12 PARAM,PARAMRAD
      CHARACTER*40 FICHIERRAW
      CHARACTER*1 TAB(NPM)
      DIMENSION FAILLE(NFM,10),DEPL(NPM)
      COMMON /COMOB/ Q,ST,CT,CTT,SF,CF,ELAS,US(3)

      PI=3.14159265339
C ----------Le coefficient d'elasticite est choisi egal a 0.5
      ELAS=0.5

C ----------Fichiers d'entree
      WRITE(*,*) 'Nom du fichier de parametres'
      READ(*,'(A)') PARAM
      OPEN(16,STATUS='OLD',FILE=PARAM)
      WRITE(*,*) 'Nom du fichier de parametres pour interferogramme'
      READ(*,'(A)') PARAMRAD
      OPEN(17,STATUS='OLD',FILE=PARAMRAD)


C ----------Lecture du nombre de failles et d'inversions
C           et de IV (pour avoir les variances a-posteriori)
C           puis des 10 parametres de chaque faille

      READ(16,*) NFT,NINV,IV,IMP
      DO 11 NF=1,NFT
   11 READ(16,*) (FAILLE(NF,K),DUMMY,K=1,10)

C ----------Lecture des Coordonnees Angle Scene, Pas, Nligne, Ncolonnes
      READ(17,*) XS,YS,PAS,NL,NC
C ----------Lecture des vecteurs d'incidence de la visee radar
      READ(17,*) VIX,VIY,VIZ


C ----------Fichier de sortie
      WRITE(*,*) 'Nom du fichier RAW de sortie'
      READ(*,'(A)') FICHIERRAW
      OPEN(18,FILE=FICHIERRAW,ACCESS='DIRECT',RECL=NC)

C ----------Controle des dimensions

      IF(NC.GT.NPM) THEN
      WRITE(*,*) 'Trop de colonnes'
      STOP
      ENDIF

C ----------Sortie du sous-programme
      WRITE(*,*) NFT,'    Failles'
      WRITE(*,*) NL ,'    Lignes'
      WRITE(*,*) NC ,'    Colonnes'
      CLOSE(16)
      CLOSE(17)

C ----------Boucle de calcul direct par ligne
      DO 1 IL=1,NL

C ----------Vidage de DEPL
      DO 10 ID=1,NC
   10 DEPL(ID)=0.

C ----------Boucle de calcul par faille
      DO 20 I=1,NFT

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
      DO 21 J=1,NC
      XX= XS + PAS * DBLE(J-1)
      YY= YS - PAS * DBLE(IL-1)
      CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,1)
      DEPL(J)=DEPL(J)+VIZ*AZ
      CALL FAULT(XX,YY,XPA,YPA,RA,W,D,AX,AY,AZ,2)
      DEPL(J)=DEPL(J)+VIX*AX+VIY*AY
   21 CONTINUE
   20 CONTINUE


C ----------Remplissage du tableau resultat
      DO 2 K=1,NC

      IFRANGE=(DEPL(K)/28. - DINT(DEPL(K)/28.)) * 256
    2 TAB(K)=CHAR(IFRANGE)

      WRITE(18,REC=IL) (TAB(KK),KK=1,NC)


    1 CONTINUE

        CLOSE(18)
        STOP
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


