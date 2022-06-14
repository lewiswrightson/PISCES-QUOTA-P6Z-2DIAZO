MODULE p6zprod
   !!======================================================================
   !!                         ***  MODULE p6zprod  ***
   !! TOP :  Growth Rate of the four phytoplanktons groups 
   !!        PISCES-QUOTA version of the module with explicit diazotrophy
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p6z_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   p6z_prod_init  :   Initialization of the parameters for growth
   !!   p6z_prod_alloc :   Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim
   USE p6zlim          !  Co-limitations of differents nutrients
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p6z_prod         ! called in p6zbio.F90
   PUBLIC   p6z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p6z_prod_alloc

   !! * Shared module variables
   REAL(wp), PUBLIC ::  pislopen        !: P-I slope of nanophytoplankton
   REAL(wp), PUBLIC ::  pislopep        !: P-I slope of picophytoplankton
   REAL(wp), PUBLIC ::  pisloped        !: P-I slope of diatoms
   REAL(wp), PUBLIC ::  xadap           !: Adaptation factor to low light
   REAL(wp), PUBLIC ::  excretn         !: Excretion ratio of nanophyto
   REAL(wp), PUBLIC ::  excretp         !: Excretion ratio of picophyto
   REAL(wp), PUBLIC ::  excretd         !: Excretion ratio of diatoms
   REAL(wp), PUBLIC ::  bresp           !: Basal respiration rate
   REAL(wp), PUBLIC ::  thetanpm        !: Maximum Chl/N ratio of picophyto
   REAL(wp), PUBLIC ::  thetannm        !: Maximum Chl/N ratio of nanophyto
   REAL(wp), PUBLIC ::  thetandm        !: Maximum Chl/N ratio of diatoms
   REAL(wp), PUBLIC ::  chlcmin         !: Minimum Chl/C ratio of phytoplankton
   REAL(wp), PUBLIC ::  grosip          !: Mean Si/C ratio of diatoms
   REAL(wp), PUBLIC ::  pislopedz       !: P-I slope of diazotrophs
   REAL(wp), PUBLIC ::  excretdz        !: Excretion ratio of diazotrophs
   REAL(wp), PUBLIC ::  thetandzm       !: Maximum Chl/C ratio of diazos
   REAL(wp), PUBLIC ::  Facul_max       !: limit on Diazotrophs Facultative capacity
   REAL(wp), PUBLIC ::  zmxlcoef        !: zmxlcoef

   REAL(wp), PUBLIC ::  pislopecr       !: P-I slope of diazotroph2
   REAL(wp), PUBLIC ::  excretcr        !: Excretion ratio of diazotroph2
   REAL(wp), PUBLIC ::  thetancrm       !: Maximum Chl/C ratio of diazo2
   REAL(wp), PUBLIC ::  photoinhcr      !: Photoinhibition of Croco
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zdaylen ! day length
   
   REAL(wp) :: r1_rday                !: 1 / rday
   REAL(wp) :: texcretn               !: 1 - excretn 
   REAL(wp) :: texcretp               !: 1 - excretp 
   REAL(wp) :: texcretd               !: 1 - excretd        
   REAL(wp) :: texcretdz              !: 1 - excretdz
   REAL(wp) :: totnfix 
   REAL(wp) :: texcretcr              !: 1 - excretcr
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p6zprod.F90 13233 2020-07-02 18:34:16Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p6z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!              Computes also the uptake of nutrients. PISCES-quota
      !!              relies on a full quota formalism
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsilfac, znanotot, zpicotot, zdiattot, zconctemp, zconctemp2
      REAL(wp) ::   zration, zratiop, zratiof, zmax, zsilim, ztn, zadap
      REAL(wp) ::   zpronmax, zpropmax, zprofmax, zrat
      REAL(wp) ::   zlim, zsilfac2, zsiborn, zprod, zprontot, zproptot, zprodtot
      REAL(wp) ::   zprnutmax, zdocprod, zprochln, zprochld, zprochlp
      REAL(wp) ::   zpislopen, zpislopep, zpisloped
      REAL(wp) ::   zrum, zcodel, zargu, zval, zfeup
      REAL(wp) ::   zfact, zrfact2, zmaxsi, zratiosi, zsizetmp, zlimfac
      REAL(wp) ::   zdiaztot, zprodztot, zprochldz, zpislopedz
      REAL(wp) ::   zpronfmax, zpronfmax2, facul
      REAL(wp) ::   zcroctot, zprocrtot, zprochlcr, zpislopecr
      REAL(wp) ::   zpronfmaxcr, facul2
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj    ) :: zmixnano, zmixpico, zmixdiat, zstrn
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpislopeadn, zpislopeadp, zpislopeadd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprnut, zprnutp, zprmaxp, zprmaxn, zprmaxd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprbio, zprpic, zprdia, zysopt
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprchln, zprchlp, zprchld
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprorcan, zprorcap, zprorcad 
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprofed, zprofep, zprofen
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpronewn, zpronewp, zpronewd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zproregn, zproregp, zproregd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpropo4n, zpropo4p, zpropo4d
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprodopn, zprodopp, zprodopd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zrespn, zrespp, zrespd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zcroissn, zcroissp, zcroissd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmxl_fac, zmxl_chl
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpligprod1, zpligprod2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpislopeaddz, zprmaxdz, zprdiaz, zprchldz
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprorcadz, zprofedz, zpronewdz, zproregdz
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpropo4dz, zprodopdz, zrespdz, zcroissdz
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: znfix, tgfuncnfix, zprnfix, zprnutdz, puptkdz
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: qn_diazo, nutlim_dz, Facul_out, znfix_dz, zmxl_fac_tricho
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpislopeadcr, zprmaxcr, zprcroc, zprchlcr
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprorcacr, zprofecr, zpronewcr, zproregcr
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpropo4cr, zprodopcr, zrespcr, zcroisscr
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: tgfuncnfix2, zprnfix2, zprnutcr, puptkcr
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: qn_croco, nutlim_cr, Facul_out2, znfix_cr
      REAL(wp), DIMENSION(jpi,jpj    ) :: zwork
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_prod')

      ! Initialize the local arrays
      zprorcan(:,:,:) = 0._wp ; zprorcap(:,:,:) = 0._wp ; zprorcad(:,:,:) = 0._wp
      zprofed (:,:,:) = 0._wp ; zprofep (:,:,:) = 0._wp ; zprofen (:,:,:) = 0._wp
      zpronewn(:,:,:) = 0._wp ; zpronewp(:,:,:) = 0._wp ; zpronewd(:,:,:) = 0._wp
      zproregn(:,:,:) = 0._wp ; zproregp(:,:,:) = 0._wp ; zproregd(:,:,:) = 0._wp 
      zpropo4n(:,:,:) = 0._wp ; zpropo4p(:,:,:) = 0._wp ; zpropo4d(:,:,:) = 0._wp
      zprdia  (:,:,:) = 0._wp ; zprpic  (:,:,:) = 0._wp ; zprbio  (:,:,:) = 0._wp
      zprodopn(:,:,:) = 0._wp ; zprodopp(:,:,:) = 0._wp ; zprodopd(:,:,:) = 0._wp
      zysopt  (:,:,:) = 0._wp
      zrespn  (:,:,:) = 0._wp ; zrespp  (:,:,:) = 0._wp ; zrespd  (:,:,:) = 0._wp 

      ! diazo
      zprorcadz(:,:,:) = 0._wp ; zprofedz(:,:,:) = 0._wp ; zpronewdz(:,:,:) = 0._wp
      zproregdz(:,:,:) = 0._wp ; zpropo4dz(:,:,:) = 0._wp ; zprdiaz(:,:,:) = 0._wp
      zprodopdz(:,:,:) = 0._wp ; zrespdz (:,:,:) = 0._wp ; znfix(:,:,:) = 0._wp
      zprnfix(:,:,:) = 0._wp ; zwork(:,:) = 0._wp ; Facul_out(:,:,:) = 0._wp
      qn_diazo(:,:,:) =0._wp ; znfix_dz(:,:,:) = 0._wp
      nutlim_dz(:,:,:) = 0._wp

      ! diazo2
      zprorcacr(:,:,:) = 0._wp ; zprofecr(:,:,:) = 0._wp ; zpronewcr(:,:,:) = 0._wp
      zproregcr(:,:,:) = 0._wp ; zpropo4cr(:,:,:) = 0._wp ; zprcroc(:,:,:) = 0._wp
      zprodopcr(:,:,:) = 0._wp ; zrespcr (:,:,:) = 0._wp 
      zprnfix2(:,:,:) = 0._wp ; Facul_out2(:,:,:) = 0._wp
      qn_croco(:,:,:) =0._wp ; znfix_cr(:,:,:) = 0._wp
      nutlim_cr(:,:,:) = 0._wp
      ! Computation of the optimal production rates and nutrient uptake
      ! rates. Based on a Q10 description of the thermal dependency.
      zprnut (:,:,:) = 0.6_wp * (1.0 + zpsino3 * qnnmax ) * r1_rday * tgfunc(:,:,:)
      zprnutp(:,:,:) =  0.6_wp * (1. + zpsino3 * qnpmax ) * r1_rday * tgfunc3(:,:,:)
      zprmaxn(:,:,:) = ( 0.6_wp * (1. + zpsino3 * qnnmax ) ) * r1_rday * tgfunc(:,:,:)
      zprmaxd(:,:,:) = ( 0.6_wp * (1. + zpsino3 * qndmax ) ) * r1_rday * tgfunc(:,:,:)
      zprmaxp(:,:,:) = ( 0.4_wp * (1. + zpsino3 * qnpmax ) ) * r1_rday * tgfunc3(:,:,:)

      ! Computation of diazotroph growth rate
      
      ! Trichodesmium growth based on Jiang et al. (2018)

     ! set diazo growth rate Tricho
      DO jk = 1, jpkm1
       DO jj = 1, jpj
         DO ji = 1, jpi
         zprnfix(ji,jj,jk) = ((-3.99E-04 * (tsn(ji,jj,jk,jp_tem)**3))+(0.02685 * (tsn(ji,jj,jk,jp_tem)**2)) +  &
         &                   (-0.555*tsn(ji,jj,jk,jp_tem)) + 3.633) * r1_rday !Growth of tricho based on Jiang et al 2018
         IF(tsn(ji,jj,jk,jp_tem) .LT. 17.2 .OR. tsn(ji,jj,jk,jp_tem) .GT. 34.9 ) THEN
         zprnfix(ji,jj,jk) = 0.
         ENDIF
         zprmaxdz(ji,jj,jk) = (1. + zpsino3 * qndzmax) * zprnfix(ji,jj,jk)
         zprnutdz(ji,jj,jk) = (1 + zpsino3 * qndzmax) * zprnfix(ji,jj,jk)
         END DO
        END DO
      END DO
     
     ! Crocosphaera growth rate based on Yang et al. (2021)
      DO jk = 1, jpkm1
       DO jj = 1, jpj
         DO ji = 1, jpi
         zprnfix2(ji,jj,jk) = ((-9.097E-05 * (tsn(ji,jj,jk,jp_tem)**3))+(0.00134 * (tsn(ji,jj,jk,jp_tem)**2)) +  &
         &                   (0.1377*tsn(ji,jj,jk,jp_tem)) - 2.561) * r1_rday  !Growth of Croco based on Yang et al 2021 
 
        IF(tsn(ji,jj,jk,jp_tem) .LT. 20 .OR. tsn(ji,jj,jk,jp_tem) .GT. 34.9 ) THEN
         zprnfix2(ji,jj,jk) = 0.
         ENDIF
         zprmaxcr(ji,jj,jk) = (1. + zpsino3 * qncrmax) * zprnfix2(ji,jj,jk)
         zprnutcr(ji,jj,jk) = (1 + zpsino3 * qncrmax) * zprnfix2(ji,jj,jk)
         END DO
        END DO
      END DO


      ! compute the day length depending on latitude and the day
      ! Astronomical parameterization taken from HAMOCC3
      zrum = REAL( nday_year - 80, wp ) / REAL( nyear_len(1), wp )
      zcodel = ASIN(  SIN( zrum * rpi * 2._wp ) * SIN( rad * 23.5_wp )  )

      ! day length in hours
      zstrn(:,:) = 0.
      DO jj = 1, jpj
         DO ji = 1, jpi
            zargu = TAN( zcodel ) * TAN( gphit(ji,jj) * rad )
            zargu = MAX( -1., MIN(  1., zargu ) )
            zstrn(ji,jj) = MAX( 0.0, 24. - 2. * ACOS( zargu ) / rad / 15. )
         END DO
      END DO

      ! Impact of the day duration and light intermittency on phytoplankton growth
      ! Intermittency is supposed to have a similar effect on production as 
      ! day length (Shatwell et al., 2012). The correcting factor is zmxl_fac. 
      ! zmxl_chl is the fractional day length and is used to compute the mean
      ! PAR during daytime. The effect of mixing is computed using the 
      ! absolute light level definition of the euphotic zone
      ! ------------------------------------------------------------------------- 
      DO jk = 1, jpkm1
         DO jj = 1 ,jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  zval = MAX( 1., zstrn(ji,jj) )
                  IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                     zval = zval * MIN(1., heup_01(ji,jj) / ( hmld(ji,jj) + rtrn ))
                  ENDIF
                  zmxl_chl(ji,jj,jk) = zval / 24.
                  zmxl_fac(ji,jj,jk) = 1.0 - exp( -0.26 * zval )
                  zmxl_fac_tricho(ji,jj,jk) = 1.0 - exp( zmxlcoef * zval )
                  !IF ( gdepw_n(ji,jj,jk+1) == hmld(ji,jj) ) THEN
                  !zmxl_fac_tricho(ji,jj,jk) = zmxl_fac(ji,jj,jk)
                  !ENDIF
               ENDIF
            END DO
         END DO
      END DO

      zprbio(:,:,:) = zprmaxn(:,:,:) * zmxl_fac(:,:,:)
      zprdia(:,:,:) = zprmaxd(:,:,:) * zmxl_fac(:,:,:)
      zprpic(:,:,:) = zprmaxp(:,:,:) * zmxl_fac(:,:,:)
      zprdiaz(:,:,:) = zprmaxdz(:,:,:) * zmxl_fac_tricho(:,:,:)
      !zprdiaz(:,:,:) = zprnutdz(:,:,:) * zmxl_fac(:,:,:)
      zprcroc(:,:,:) = zprmaxcr(:,:,:) * zmxl_fac(:,:,:)

      ! Maximum light intensity
      zdaylen(:,:) = MAX(1., zstrn(:,:)) / 24.

      ! Computation of the P-I slope for nanos, picos and diatoms
      ! The formulation proposed by Geider et al. (1997) has been used.
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  ztn         = MAX( 0., tsn(ji,jj,jk,jp_tem) - 15. )
                  zadap       = xadap * ztn / ( 2.+ ztn )
                  ! Nanophytoplankton
                  zpislopeadn(ji,jj,jk) = pislopen * trb(ji,jj,jk,jpnch)    &
                  &                       /( trb(ji,jj,jk,jpphy) * 12. + rtrn)
                  !zpislopeadn(ji,jj,jk) = 1


                  ! Picophytoplankton
                  zpislopeadp(ji,jj,jk) = pislopep * ( 1. + zadap * EXP( -0.25 * epico(ji,jj,jk) ) )   &
                  &                       * trb(ji,jj,jk,jppch) /( trb(ji,jj,jk,jppic) * 12. + rtrn)

                  ! Diatoms
                  zpislopeadd(ji,jj,jk) = pisloped * trb(ji,jj,jk,jpdch)    &
                     &                    /( trb(ji,jj,jk,jpdia) * 12. + rtrn)
                  ! Diazotrophs
                  zpislopeaddz(ji,jj,jk) = pislopedz * trb(ji,jj,jk,jpchd)    &
                  &                       /( trb(ji,jj,jk,jpcdz) * 12. + rtrn)

                  ! Diazotroph2 (croco)
                  zpislopeadcr(ji,jj,jk) = pislopecr * trb(ji,jj,jk,jpchc)    &
                  &                       /( trb(ji,jj,jk,jpccr) * 12. + rtrn) 
                  !
                  zpislopen = zpislopeadn(ji,jj,jk) / ( zprbio(ji,jj,jk) * rday * xlimphy(ji,jj,jk) + rtrn )
                  zpislopep = zpislopeadp(ji,jj,jk) / ( zprpic(ji,jj,jk) * rday * xlimpic(ji,jj,jk) + rtrn )
                  zpisloped = zpislopeadd(ji,jj,jk) / ( zprdia(ji,jj,jk) * rday * xlimdia(ji,jj,jk) + rtrn )
                  zpislopedz = zpislopeaddz(ji,jj,jk) / ( zprdiaz(ji,jj,jk) *rday * xlimdiaz(ji,jj,jk) + rtrn )
                  zpislopecr = zpislopeadcr(ji,jj,jk) / ( zprcroc(ji,jj,jk) *rday * xlimcroc(ji,jj,jk) + rtrn )
                  ! Computation of production function for Carbon
                  ! Actual light levels are used here 
                  !  ---------------------------------------------
                  zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1.- EXP( -zpislopen * enano(ji,jj,jk) / (zmxl_chl(ji,jj,jk)+rtrn) )  )
                  zprpic(ji,jj,jk) = zprpic(ji,jj,jk) * ( 1.- EXP( -zpislopep * epico(ji,jj,jk) / (zmxl_chl(ji,jj,jk)+rtrn))  )
                  zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediat(ji,jj,jk) / (zmxl_chl(ji,jj,jk)+rtrn))  )
                  zprdiaz(ji,jj,jk) = zprdiaz(ji,jj,jk) * ( 1.- EXP( -zpislopedz * ediaz(ji,jj,jk) / (zmxl_chl(ji,jj,jk)+rtrn))  )

                  zprcroc(ji,jj,jk) = zprcroc(ji,jj,jk) * ( 1.- EXP( -zpislopecr * ecroc(ji,jj,jk) / (zmxl_chl(ji,jj,jk)+rtrn))  ) * EXP( -photoinhcr * ecroc(ji,jj,jk) / (zmxl_chl(ji,jj,jk)+rtrn))

                  !  Computation of production function for Chlorophyll
                  !  Mean light level in the mixed layer (when appropriate)
                  !  is used here (acclimation is in general slower than 
                  !  the characteristic time scales of vertical mixing)
                  !  ------------------------------------------------------
                  zpislopen = zpislopen * zmxl_fac(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zpisloped = zpisloped * zmxl_fac(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zpislopep = zpislopep * zmxl_fac(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zpislopedz = zpislopedz * zmxl_fac_tricho(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zpislopecr = zpislopecr * zmxl_fac(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprchln(ji,jj,jk) = zprmaxn(ji,jj,jk) * ( 1.- EXP( -zpislopen * enanom(ji,jj,jk) )  )
                  zprchlp(ji,jj,jk) = zprmaxp(ji,jj,jk) * ( 1.- EXP( -zpislopep * epicom(ji,jj,jk) )  )
                  zprchld(ji,jj,jk) = zprmaxd(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediatm(ji,jj,jk) )  )
                  zprchldz(ji,jj,jk) = zprmaxdz(ji,jj,jk) * ( 1.- EXP( -zpislopedz * ediazm(ji,jj,jk) )  )
                  zprchlcr(ji,jj,jk) = zprmaxcr(ji,jj,jk) * ( 1.- EXP( -zpislopecr * ecrocm(ji,jj,jk) )  ) * EXP(-0.5 * ecrocm(ji,jj,jk))
               ENDIF
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  ! Si/C of diatoms
                  ! ------------------------
                  ! Si/C increases with iron stress and silicate availability (zsilfac)
                  ! Si/C is arbitrariliy increased for very high Si concentrations
                  ! to mimic the very high ratios observed in the Southern Ocean (zsilfac2)
                  ! A parameterization derived from Flynn (2003) is used for the control
                  ! when Si is not limiting which is similar to the parameterisation
                  ! proposed by Gurney and Davidson (1999).
                  ! -----------------------------------------------------------------------
                  zlim  = trb(ji,jj,jk,jpsil) / ( trb(ji,jj,jk,jpsil) + xksi1 )
                  zsilim = MIN( zprdia(ji,jj,jk) / ( zprmaxd(ji,jj,jk) + rtrn ), xlimsi(ji,jj,jk) )
                  zsiborn = trb(ji,jj,jk,jpsil) * trb(ji,jj,jk,jpsil) * trb(ji,jj,jk,jpsil)
                  IF (gphit(ji,jj) < -30 ) THEN
                    zsilfac2 = 1. + 1. * zsiborn / ( zsiborn + xksi2**3 )
                  ELSE
                    zsilfac2 = 1.
                  ENDIF
                  zratiosi = 1.0 - trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn ) / ( zsilfac2 * grosip * 3.0 + rtrn )
                  zratiosi = MAX(0., MIN(1.0, zratiosi) )
                  zmaxsi  = (1.0 + 0.1**4) * zratiosi**4 / ( zratiosi**4 + 0.1**4 )
                  IF ( xlimsi(ji,jj,jk) /= xlimdia(ji,jj,jk) ) THEN
                     zysopt(ji,jj,jk) = zlim * zsilfac2 * grosip * 1.0 * zmaxsi
                  ELSE
                     zysopt(ji,jj,jk) = zlim * zsilfac2 * grosip * 1.0 * zsilim**0.75 * zmaxsi
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO

      !  Sea-ice effect on production
      ! No production is assumed below sea ice
      ! -------------------------------------- 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zprbio(ji,jj,jk)  = zprbio(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprpic(ji,jj,jk)  = zprpic(ji,jj,jk) * ( 1. - fr_i(ji,jj) ) 
               zprdia(ji,jj,jk)  = zprdia(ji,jj,jk) * ( 1. - fr_i(ji,jj) ) 
               zprnut(ji,jj,jk)  = zprnut(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprnutp(ji,jj,jk) = zprnutp(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprdiaz(ji,jj,jk) = zprdiaz(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprnutdz(ji,jj,jk) = zprnutdz(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprcroc(ji,jj,jk) = zprcroc(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprnutcr(ji,jj,jk) = zprnutcr(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
            END DO
         END DO
      END DO

      ! Computation of the various production and uptake terms of nanophytoplankton 
      ! Interactions between N and P are modeled according to the Chain Model 
      ! of Pahlow et al. (2009). Iron uptake is modeled following traditional
      ! Droop kinetics. When the quota is approaching the maximum achievable
      ! quota, uptake is downregulated according to a sigmoidal function 
      ! (power 2), as proposed by Flynn (2003)
      ! ---------------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto.
                  zprorcan(ji,jj,jk) = zprbio(ji,jj,jk)  * xlimphy(ji,jj,jk) * trb(ji,jj,jk,jpphy) * rfact2
                  ! Size computation
                  ! Size is made a function of the limitation of of phytoplankton growth
                  ! Strongly limited cells are supposed to be smaller. sizena is the 
                  ! size at time step t+1 and is thus updated at the end of the 
                  ! current time step
                  ! --------------------------------------------------------------------
                  zlimfac = xlimphys(ji,jj,jk) * zprchln(ji,jj,jk) / ( zprmaxn(ji,jj,jk) + rtrn )
                  zsizetmp = 1.0 + 1.3 * ( xsizern - 1.0 ) * zlimfac**3/(0.3 + zlimfac**3)
                  sizena(ji,jj,jk) = min(xsizern, max( sizena(ji,jj,jk), zsizetmp ) )
                  ! Maximum potential uptake rate
                  zration = trb(ji,jj,jk,jpnph) / ( trb(ji,jj,jk,jpphy) + rtrn )
                  zratiop = trb(ji,jj,jk,jppph) / ( trb(ji,jj,jk,jpphy) + rtrn )
                  zratiof = trb(ji,jj,jk,jpnfe) / ( trb(ji,jj,jk,jpphy) + rtrn )
                  zprnutmax = zprnut(ji,jj,jk) * fvnuptk(ji,jj,jk) / rno3 * trb(ji,jj,jk,jpphy) * rfact2
                  ! Uptake of nitrogen
                  zrat = 1.0 - MIN( 1., zration / (xqnnmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  zpronmax = zprnutmax * zmax * MAX(0., MIN(1., ( zratiop - xqpnmin(ji,jj,jk) )   &
                  &          / ( xqpnmax(ji,jj,jk) - xqpnmin(ji,jj,jk) + rtrn ), xlimnfe(ji,jj,jk) ) )
                  zpronewn(ji,jj,jk) = zpronmax * xnanono3(ji,jj,jk)
                  zproregn(ji,jj,jk) = zpronmax * xnanonh4(ji,jj,jk)
                  ! Uptake of phosphorus and DOP
                  zrat = 1.0 - MIN( 1., zratiop / (xqpnmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  zpropmax = zprnutmax * zmax * xlimnfe(ji,jj,jk) * 16. / 10.
                  zpropo4n(ji,jj,jk) = zpropmax * xnanopo4(ji,jj,jk)
                  zprodopn(ji,jj,jk) = zpropmax * xnanodop(ji,jj,jk)
                  ! Uptake of iron
                  zrat = 1.0 - MIN( 1., zratiof / qfnmax )
                  zmax = MAX(0., MIN(1., zrat**2/ (0.05**2 + zrat**2) ) )
                  zprofmax = zprnutmax * qfnmax * zmax 
                  zprofen(ji,jj,jk) = zprofmax * xnanofer(ji,jj,jk)    &
                  &          * (1. + 0.8 * xnanono3(ji,jj,jk) / ( rtrn  &
                  &          + xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) ) * (1. - xnanofer(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the various production and uptake terms of picophytoplankton 
      ! Interactions between N and P are modeled according to the Chain Model 
      ! of Pahlow et al. (2009). Iron uptake is modeled following traditional
      ! Droop kinetics. When the quota is approaching the maximum achievable
      ! quota, uptake is downregulated according to a sigmoidal function 
      ! (power 2), as proposed by Flynn (2003)
      ! ---------------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for picophyto.
                  zprorcap(ji,jj,jk) = zprpic(ji,jj,jk)  * xlimpic(ji,jj,jk) * trb(ji,jj,jk,jppic) * rfact2
                  ! Size computation
                  ! Size is made a function of the limitation of of phytoplankton growth
                  ! Strongly limited cells are supposed to be smaller. sizepa is
                  ! size at time step t+1 and is thus updated at the end of the 
                  ! current time step
                  ! --------------------------------------------------------------------
                  zlimfac = zprchlp(ji,jj,jk)  * xlimpics(ji,jj,jk) / ( zprmaxp(ji,jj,jk) + rtrn )
                  zsizetmp = 1.0 + 1.3 * ( xsizerp - 1.0 ) * zlimfac**3/(0.3 + zlimfac**3)
                  sizepa(ji,jj,jk) = min(xsizerp, max( sizepa(ji,jj,jk), zsizetmp ) )
                  ! Maximum potential uptake rate of nutrients
                  zration = trb(ji,jj,jk,jpnpi) / ( trb(ji,jj,jk,jppic) + rtrn )
                  zratiop = trb(ji,jj,jk,jpppi) / ( trb(ji,jj,jk,jppic) + rtrn )
                  zratiof = trb(ji,jj,jk,jppfe) / ( trb(ji,jj,jk,jppic) + rtrn )
                  zprnutmax = zprnutp(ji,jj,jk) * fvpuptk(ji,jj,jk) / rno3 * trb(ji,jj,jk,jppic) * rfact2
                  ! Uptake of nitrogen
                  zrat = 1.0 - MIN( 1., zration / (xqnpmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2/ (0.05**2 + zrat**2) ) )
                  zpronmax = zprnutmax * zmax * MAX(0., MIN(1., ( zratiop - xqppmin(ji,jj,jk) )   &
                  &          / ( xqppmax(ji,jj,jk) - xqppmin(ji,jj,jk) + rtrn ), xlimpfe(ji,jj,jk) ) )
                  zpronewp(ji,jj,jk) = zpronmax * xpicono3(ji,jj,jk) 
                  zproregp(ji,jj,jk) = zpronmax * xpiconh4(ji,jj,jk)
                  ! Uptake of phosphorus
                  zrat = 1.0 - MIN( 1., zratiop / (xqppmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  zpropmax = zprnutmax * zmax * xlimpfe(ji,jj,jk) * 16./10.
                  zpropo4p(ji,jj,jk) = zpropmax * xpicopo4(ji,jj,jk)
                  zprodopp(ji,jj,jk) = zpropmax * xpicodop(ji,jj,jk)
                  ! Uptake of iron
                  zrat = 1.0 - MIN( 1., zratiof / qfpmax )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  zprofmax = zprnutmax * qfpmax * zmax
                  zprofep(ji,jj,jk) = zprofmax * xpicofer(ji,jj,jk)  &
                  &          * (1. + 0.8 * xpicono3(ji,jj,jk) / ( rtrn   &
                  &          + xpicono3(ji,jj,jk) + xpiconh4(ji,jj,jk) ) * (1. - xpicofer(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the various production and uptake terms of diatoms
      ! Interactions between N and P are modeled according to the Chain Model 
      ! of Pahlow et al. (2009). Iron uptake is modeled following traditional
      ! Droop kinetics. When the quota is approaching the maximum achievable
      ! quota, uptake is downregulated according to a sigmoidal function 
      ! (power 2), as proposed by Flynn (2003)
      ! ---------------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for diatomees
                  zprorcad(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk) * trb(ji,jj,jk,jpdia) * rfact2
                  ! Size computation
                  ! Size is made a function of the limitation of of phytoplankton growth
                  ! Strongly limited cells are supposed to be smaller. sizeda is
                  ! size at time step t+1 and is thus updated at the end of the 
                  ! current time step. 
                  ! --------------------------------------------------------------------
                  zlimfac = zprchld(ji,jj,jk) * xlimdias(ji,jj,jk) / ( zprmaxd(ji,jj,jk) + rtrn )
                  zsizetmp = 1.0 + 1.3 * ( xsizerd - 1.0 ) * zlimfac**3/(0.3 + zlimfac**3)
                  sizeda(ji,jj,jk) = min(xsizerd, max( sizeda(ji,jj,jk), zsizetmp ) )
                  ! Maximum potential uptake rate of nutrients
                  zration = trb(ji,jj,jk,jpndi) / ( trb(ji,jj,jk,jpdia) + rtrn )
                  zratiop = trb(ji,jj,jk,jppdi) / ( trb(ji,jj,jk,jpdia) + rtrn )
                  zratiof = trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) + rtrn )
                  zprnutmax = zprnut(ji,jj,jk) * fvduptk(ji,jj,jk) / rno3 * trb(ji,jj,jk,jpdia) * rfact2
                  ! Uptake of nitrogen
                  zrat = 1.0 - MIN( 1., zration / (xqndmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  zpronmax = zprnutmax * zmax * MAX(0., MIN(1., ( zratiop - xqpdmin(ji,jj,jk) )   &
                  &          / ( xqpdmax(ji,jj,jk) - xqpdmin(ji,jj,jk) + rtrn ), xlimdfe(ji,jj,jk) ) )
                  zpronewd(ji,jj,jk) = zpronmax * xdiatno3(ji,jj,jk)
                  zproregd(ji,jj,jk) = zpronmax * xdiatnh4(ji,jj,jk)
                  ! Uptake of phosphorus
                  zrat = 1.0 - MIN( 1., zratiop / (xqpdmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2/ (0.05**2 + zrat**2) ) )
                  zpropmax = zprnutmax * zmax * xlimdfe(ji,jj,jk) * 16./10.
                  zpropo4d(ji,jj,jk) = zpropmax * xdiatpo4(ji,jj,jk)
                  zprodopd(ji,jj,jk) = zpropmax * xdiatdop(ji,jj,jk)
                  ! Uptake of iron
                  zrat = 1.0 - MIN( 1., zratiof / qfdmax )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  zprofmax = zprnutmax * qfdmax * zmax
                  zprofed(ji,jj,jk) = zprofmax * xdiatfer(ji,jj,jk)    &
                  &          * (1. + 0.8 * xdiatno3(ji,jj,jk) / ( rtrn   &
                  &          + xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) ) * (1. - xdiatfer(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the various production terms of diazotrophs 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for diazotrophs
                  zprorcadz(ji,jj,jk) = zprdiaz(ji,jj,jk) * xlimdiaz(ji,jj,jk) * trb(ji,jj,jk,jpcdz) * rfact2
                  ! size computation
                  zlimfac = xlimdiazo(ji,jj,jk) * zprchldz(ji,jj,jk) / ( zprnutdz(ji,jj,jk) + rtrn )
                  zsizetmp = 1.0 + 1.3 * ( xsizerdz - 1.0 ) * zlimfac**3/(0.3 + zlimfac**3 + rtrn)
                  sizedza(ji,jj,jk) = min(xsizerdz,max( sizedza(ji,jj,jk), zsizetmp ) )
                  ! nutrient uptake
                  zration = trb(ji,jj,jk,jpndz) / ( trb(ji,jj,jk,jpcdz) + rtrn )
                  zratiop = trb(ji,jj,jk,jppdz) / ( trb(ji,jj,jk,jpcdz) + rtrn )
                  zratiof = trb(ji,jj,jk,jpfed) / ( trb(ji,jj,jk,jpcdz) + rtrn )
                  zprnutmax = zprnutdz(ji,jj,jk) * fvdzuptk(ji,jj,jk) / rno3 * trb(ji,jj,jk,jpcdz) * rfact2     
                  vnfmax_dz(ji,jj,jk) = zprnutmax
                  ! Uptake of nitrogen
                  zrat = 1.0 - MIN( 1., zration / (xqndzmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  qn_diazo(ji,jj,jk) = zmax
                  zpronmax = zprnutmax * zmax * MAX(0., MIN(1., ( zratiop - xqpdzmin(ji,jj,jk) )   &
                  &          / ( xqpdzmax(ji,jj,jk) - xqpdzmin(ji,jj,jk) + rtrn ), xlimdzfe(ji,jj,jk) ) )
                  zpronewdz(ji,jj,jk) = zpronmax * xdiazno3(ji,jj,jk)
                  zproregdz(ji,jj,jk) = zpronmax * xdiaznh4(ji,jj,jk)
                  ! Diazotrophy 
                  !zpronfmax2 = vnfmax_dz(ji,jj,jk) * qn_diazo(ji,jj,jk) * xlimdiaz(ji,jj,jk) * &
                  !&           ( 1. - ( xdiazno3(ji,jj,jk) + xdiaznh4(ji,jj,jk) ))
                  !nfix_iue(ji,jj,jk) = zpronfmax2 / (vnfmax_dz(ji,jj,jk) +rtrn)
                  facul = ( 1. - ( xdiazno3(ji,jj,jk) + xdiaznh4(ji,jj,jk) ))
                  !IF (facul > Facul_max) THEN 
                  !   facul = 1
                  !ENDIF
                  Facul_out(ji,jj,jk) = facul
                  nutlim_dz(ji,jj,jk) = MAX(0., MIN(1., ( zratiop - xqpdzmin(ji,jj,jk) )   &
                  &          / ( xqpdzmax(ji,jj,jk) - xqpdzmin(ji,jj,jk) + rtrn ), xlimdzfe(ji,jj,jk) ) )
                  !Nitrogen Fixation
                  zpronfmax = vnfmax_dz(ji,jj,jk) * qn_diazo(ji,jj,jk) * nutlim_dz(ji,jj,jk) * &
                  &          facul  !( 1. - ( xdiazno3(ji,jj,jk) + xdiaznh4(ji,jj,jk) ) )
                  znfix_dz(ji,jj,jk) = zpronfmax
                  ! Uptake of phosphorus and DOP
                  zrat = 1.0 - MIN( 1., zratiop / (xqpdzmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  zpropmax = zprnutmax * zmax * xlimdzfe(ji,jj,jk) * 16. / 10. + (znfix_dz(ji,jj,jk) / 122)
                  puptkdz(ji,jj,jk) = zpropmax
                  zpropo4dz(ji,jj,jk) = zpropmax * xdiazpo4(ji,jj,jk)
                  zprodopdz(ji,jj,jk) = zpropmax * xdiazdop(ji,jj,jk)
                  ! Uptake of iron
                  zrat = 1.0 - MIN( 1., zratiof / qfdzmax )
                  zmax = MAX(0., MIN(1., zrat**2/ (0.05**2 + zrat**2) ) )
                  zprofmax = zprnutmax * qfdzmax * zmax
                  zprofedz(ji,jj,jk) = zprofmax * xdiazfer(ji,jj,jk)    &
                  &          * (1. + 0.8 * xdiazno3(ji,jj,jk) / ( rtrn   &
                  &          + xdiazno3(ji,jj,jk) + xdiaznh4(ji,jj,jk) ) * (1. - xdiazfer(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the various production terms of diazotrophs 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for diazotrophs
                  zprorcacr(ji,jj,jk) = zprcroc(ji,jj,jk) * xlimcroc(ji,jj,jk) * trb(ji,jj,jk,jpccr) * rfact2
                  ! size computation
                  zlimfac = xlimcroco(ji,jj,jk) * zprchlcr(ji,jj,jk) / ( zprnutcr(ji,jj,jk) + rtrn )
                  zsizetmp = 1.0 + 1.3 * ( xsizercr - 1.0 ) * zlimfac**3/(0.3 + zlimfac**3 + rtrn)
                  sizecra(ji,jj,jk) = min(xsizercr,max( sizecra(ji,jj,jk), zsizetmp ) )
                  ! nutrient uptake
                  zration = trb(ji,jj,jk,jpncr) / ( trb(ji,jj,jk,jpccr) + rtrn )
                  zratiop = trb(ji,jj,jk,jppcr) / ( trb(ji,jj,jk,jpccr) + rtrn )
                  zratiof = trb(ji,jj,jk,jpfec) / ( trb(ji,jj,jk,jpccr) + rtrn )
                  zprnutmax = zprnutcr(ji,jj,jk) * fvcruptk(ji,jj,jk) / rno3 * trb(ji,jj,jk,jpccr) * rfact2
                  vnfmax_cr(ji,jj,jk) = zprnutmax
                  ! Uptake of nitrogen
                  zrat = 1.0 - MIN( 1., zration / (xqncrmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  qn_croco(ji,jj,jk) = zmax
                  zpronmax = zprnutmax * zmax * MAX(0., MIN(1., ( zratiop - xqpcrmin(ji,jj,jk) )   &
                  &          / ( xqpcrmax(ji,jj,jk) - xqpcrmin(ji,jj,jk) + rtrn ), xlimcrfe(ji,jj,jk) ) )
                  zpronewcr(ji,jj,jk) = zpronmax * xcrocno3(ji,jj,jk)
                  zproregcr(ji,jj,jk) = zpronmax * xcrocnh4(ji,jj,jk)
                  ! Diazotrophy 
                  !zpronfmax2 = vnfmax_dz(ji,jj,jk) * qn_diazo(ji,jj,jk) *
                  !xlimdiaz(ji,jj,jk) * &
                  !&           ( 1. - ( xdiazno3(ji,jj,jk) + xdiaznh4(ji,jj,jk)
                  !))
                  !nfix_iue(ji,jj,jk) = zpronfmax2 / (vnfmax_dz(ji,jj,jk) +rtrn)
                  facul2 = ( 1. - ( xcrocno3(ji,jj,jk) + xcrocnh4(ji,jj,jk) ))
                  !IF (facul > Facul_max) THEN 
                  !   facul = 1
                  !ENDIF
                  Facul_out2(ji,jj,jk) = facul2
                  nutlim_cr(ji,jj,jk) = MAX(0., MIN(1., ( zratiop - xqpcrmin(ji,jj,jk) )   &
                  &          / ( xqpcrmax(ji,jj,jk) - xqpcrmin(ji,jj,jk) + rtrn ), xlimcrfe(ji,jj,jk) ) )
                  !Nitrogen Fixation
                  zpronfmaxcr = vnfmax_cr(ji,jj,jk) * qn_croco(ji,jj,jk) * nutlim_cr(ji,jj,jk) * &
                  &          facul2  !( 1. - ( xdiazno3(ji,jj,jk) + xdiaznh4(ji,jj,jk) ) )
                  znfix_cr(ji,jj,jk) = zpronfmaxcr
                  ! Uptake of phosphorus and DOP
                  zrat = 1.0 - MIN( 1., zratiop / (xqpcrmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., zrat**2 / (0.05**2 + zrat**2) ) )
                  zpropmax = zprnutmax * zmax * xlimcrfe(ji,jj,jk) * 16. / 10. + (znfix_cr(ji,jj,jk) / 122)
                  puptkcr(ji,jj,jk) = zpropmax
                  zpropo4cr(ji,jj,jk) = zpropmax * xcrocpo4(ji,jj,jk)
                  zprodopcr(ji,jj,jk) = zpropmax * xcrocdop(ji,jj,jk)
                  ! Uptake of iron
                  zrat = 1.0 - MIN( 1., zratiof / qfcrmax )
                  zmax = MAX(0., MIN(1., zrat**2/ (0.05**2 + zrat**2) ) )
                  zprofmax = zprnutmax * qfcrmax * zmax
                  zprofecr(ji,jj,jk) = zprofmax * xcrocfer(ji,jj,jk)    &
                  &          * (1. + 0.8 * xcrocno3(ji,jj,jk) / ( rtrn   &
                  &          + xcrocno3(ji,jj,jk) + xcrocnh4(ji,jj,jk) ) * (1. - xcrocfer(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO
      
      znfix(:,:,:) = znfix_dz(:,:,:) + znfix_cr(:,:,:)
     
      ! Production of Chlorophyll. The formulation proposed by Geider et al.
      ! is adopted here.
      ! --------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                     !  production terms for nanophyto. ( chlorophyll )
                  !znanotot = enano(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  znanotot = enanom(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  
                  zprod = rday * (zpronewn(ji,jj,jk) + zproregn(ji,jj,jk)) * zprchln(ji,jj,jk) * xlimphy(ji,jj,jk)
                  zprochln = thetannm * zprod / ( zpislopeadn(ji,jj,jk) * znanotot + rtrn )

                  zprochln = MAX(zprochln, chlcmin * 12. * zprorcan (ji,jj,jk) )
                  !  production terms for picophyto. ( chlorophyll )
                  zpicotot = epicom(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod = rday * (zpronewp(ji,jj,jk) + zproregp(ji,jj,jk)) * zprchlp(ji,jj,jk) * xlimpic(ji,jj,jk)
                  zprochlp = thetanpm * zprod / ( zpislopeadp(ji,jj,jk) * zpicotot + rtrn )
                  zprochlp = MAX(zprochlp, chlcmin * 12. * zprorcap(ji,jj,jk) )
                  !  production terms for diatoms ( chlorophyll )
                  zdiattot = ediatm(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod = rday * (zpronewd(ji,jj,jk) + zproregd(ji,jj,jk)) * zprchld(ji,jj,jk) * xlimdia(ji,jj,jk)
                  zprochld = thetandm * zprod / ( zpislopeadd(ji,jj,jk) * zdiattot + rtrn )
                  zprochld = MAX(zprochld, chlcmin * 12. * zprorcad(ji,jj,jk) )
                  !  production terms for diazotrophs ( chlorophyll )
                  zdiaztot = ediazm(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod = rday * (zpronewdz(ji,jj,jk) + zproregdz(ji,jj,jk) + znfix_dz(ji,jj,jk))    &
                  &       * zprchldz(ji,jj,jk) * xlimdiaz(ji,jj,jk)
                  zprochldz = thetandzm * zprod / ( zpislopeaddz(ji,jj,jk) * zdiaztot + rtrn )
                  zprochldz = MAX(zprochldz, chlcmin * 12. * zprorcadz(ji,jj,jk) )
                  !  production terms for diazotroph2 ( chlorophyll )
                  zcroctot = ecrocm(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod = rday * (zpronewcr(ji,jj,jk) + zproregcr(ji,jj,jk) + znfix_cr(ji,jj,jk))    &
                  &       * zprchlcr(ji,jj,jk) * xlimcroc(ji,jj,jk)
                  zprochlcr = thetancrm * zprod / ( zpislopeadcr(ji,jj,jk) * zcroctot + rtrn )
                  zprochlcr = MAX(zprochlcr, chlcmin * 12. * zprorcacr(ji,jj,jk) )      
                  !   Update the arrays TRA which contain the Chla sources and sinks
                  tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln * texcretn
                  !tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + 1 * texcretn
                  tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld * texcretd
                  tra(ji,jj,jk,jppch) = tra(ji,jj,jk,jppch) + zprochlp * texcretp
                  tra(ji,jj,jk,jpchd) = tra(ji,jj,jk,jpchd) + zprochldz * texcretdz
                  tra(ji,jj,jk,jpchc) = tra(ji,jj,jk,jpchc) + zprochlcr * texcretcr
               ENDIF
            END DO
         END DO
      END DO

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              zprontot = zpronewn(ji,jj,jk) + zproregn(ji,jj,jk)
              zproptot = zpronewp(ji,jj,jk) + zproregp(ji,jj,jk)
              zprodtot = zpronewd(ji,jj,jk) + zproregd(ji,jj,jk)
              zprodztot = zpronewdz(ji,jj,jk) + zproregdz(ji,jj,jk) + znfix_dz(ji,jj,jk)
              zprocrtot = zpronewcr(ji,jj,jk) + zproregcr(ji,jj,jk) + znfix_cr(ji,jj,jk)
              zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)  &
              &          + excretp * zprorcap(ji,jj,jk) + excretdz * zprorcadz(ji,jj,jk)  &
              &          + excretcr * zprorcacr(ji,jj,jk)
              tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zpropo4n(ji,jj,jk) - zpropo4d(ji,jj,jk)  &
              &                     - zpropo4p(ji,jj,jk) - zpropo4dz(ji,jj,jk) - zpropo4cr(ji,jj,jk)
              tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronewn(ji,jj,jk) - zpronewd(ji,jj,jk)  &
              &                     - zpronewp(ji,jj,jk) - zpronewdz(ji,jj,jk) - zpronewcr(ji,jj,jk)
              tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproregn(ji,jj,jk) - zproregd(ji,jj,jk)  &
              &                     - zproregp(ji,jj,jk) - zproregdz(ji,jj,jk) - zproregcr(ji,jj,jk)
              tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorcan(ji,jj,jk) * texcretn    &
                 &                  - zpsino3 * zpronewn(ji,jj,jk) - zpsinh4 * zproregn(ji,jj,jk)   &
                 &                  - zrespn(ji,jj,jk) 
              zcroissn(ji,jj,jk) = tra(ji,jj,jk,jpphy) / rfact2/ (trb(ji,jj,jk,jpphy) + rtrn)
              tra(ji,jj,jk,jpnph) = tra(ji,jj,jk,jpnph) + zprontot * texcretn
              tra(ji,jj,jk,jppph) = tra(ji,jj,jk,jppph) + zpropo4n(ji,jj,jk) * texcretn   &
              &                     + zprodopn(ji,jj,jk) * texcretn
              tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk) * texcretn
              tra(ji,jj,jk,jppic) = tra(ji,jj,jk,jppic) + zprorcap(ji,jj,jk) * texcretp     &
                 &                  - zpsino3 * zpronewp(ji,jj,jk) - zpsinh4 * zproregp(ji,jj,jk)   &
                 &                  - zrespp(ji,jj,jk) 
              zcroissp(ji,jj,jk) = tra(ji,jj,jk,jppic) / rfact2/ (trb(ji,jj,jk,jppic) + rtrn)
              tra(ji,jj,jk,jpnpi) = tra(ji,jj,jk,jpnpi) + zproptot * texcretp
              tra(ji,jj,jk,jpppi) = tra(ji,jj,jk,jpppi) + zpropo4p(ji,jj,jk) * texcretp   &
              &                     + zprodopp(ji,jj,jk) * texcretp
              tra(ji,jj,jk,jppfe) = tra(ji,jj,jk,jppfe) + zprofep(ji,jj,jk) * texcretp
              tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk) * texcretd   &
                 &                  - zpsino3 * zpronewd(ji,jj,jk) - zpsinh4 * zproregd(ji,jj,jk)   &
                 &                  - zrespd(ji,jj,jk) 
              zcroissd(ji,jj,jk) = tra(ji,jj,jk,jpdia) / rfact2 / (trb(ji,jj,jk,jpdia) + rtrn)
              tra(ji,jj,jk,jpndi) = tra(ji,jj,jk,jpndi) + zprodtot * texcretd
              tra(ji,jj,jk,jppdi) = tra(ji,jj,jk,jppdi) + zpropo4d(ji,jj,jk) * texcretd   &
              &                     + zprodopd(ji,jj,jk) * texcretd
              tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) * texcretd
              tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zprmaxd(ji,jj,jk) * zysopt(ji,jj,jk) * rfact2 * trb(ji,jj,jk,jpdia)
             ! Diazotrophy
              tra(ji,jj,jk,jpcdz) = tra(ji,jj,jk,jpcdz) + zprorcadz(ji,jj,jk) * texcretdz  &
                 &                  - zpsino3 * zpronewdz(ji,jj,jk) - zpsinh4 * zproregdz(ji,jj,jk)  &
                 &                  - zrespdz(ji,jj,jk) - zpsinfix * znfix_dz(ji,jj,jk)
              zcroissdz(ji,jj,jk) = tra(ji,jj,jk,jpcdz) / rfact2/ (trb(ji,jj,jk,jpcdz) + rtrn)
              tra(ji,jj,jk,jpndz) = tra(ji,jj,jk,jpndz) + zprodztot  * texcretdz
              tra(ji,jj,jk,jppdz) = tra(ji,jj,jk,jppdz) + zpropo4dz(ji,jj,jk) * texcretdz   &
              &                     + zprodopdz(ji,jj,jk) * texcretdz
              tra(ji,jj,jk,jpfed) = tra(ji,jj,jk,jpfed) + zprofedz(ji,jj,jk) * texcretdz

             ! Diazotroph2 (croco)
              tra(ji,jj,jk,jpccr) = tra(ji,jj,jk,jpccr) + zprorcacr(ji,jj,jk) * texcretcr  &
                 &                  - zpsino3 * zpronewcr(ji,jj,jk) - zpsinh4 * zproregcr(ji,jj,jk)  &
                 &                  - zrespcr(ji,jj,jk) - zpsinfix * znfix_cr(ji,jj,jk)
              zcroisscr(ji,jj,jk) = tra(ji,jj,jk,jpccr) / rfact2/ (trb(ji,jj,jk,jpccr) + rtrn)
              tra(ji,jj,jk,jpncr) = tra(ji,jj,jk,jpncr) + zprocrtot  * texcretcr
              tra(ji,jj,jk,jppcr) = tra(ji,jj,jk,jppcr) + zpropo4cr(ji,jj,jk) * texcretcr   &
              &                     + zprodopcr(ji,jj,jk) * texcretcr
              tra(ji,jj,jk,jpfec) = tra(ji,jj,jk,jpfec) + zprofecr(ji,jj,jk) * texcretcr
              !
              
              tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)  &
              &                     + excretp * zprorcap(ji,jj,jk) + excretdz * zprorcadz(ji,jj,jk)  &
              &                     + excretcr * zprorcacr(ji,jj,jk)
              tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + excretd * zprodtot + excretn * zprontot   &
              &                     + excretp * zproptot + excretdz * zprodztot        &
              &                     + excretcr * zprocrtot 
              tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + excretd * zpropo4d(ji,jj,jk) + excretn * zpropo4n(ji,jj,jk)   &
              &    - texcretn * zprodopn(ji,jj,jk) - texcretd * zprodopd(ji,jj,jk) + excretp * zpropo4p(ji,jj,jk)     &
              &    - texcretp * zprodopp(ji,jj,jk) + excretdz * zpropo4dz(ji,jj,jk) - texcretdz * zprodopdz(ji,jj,jk) &
              &    + excretcr * zpropo4cr(ji,jj,jk) - texcretcr * zprodopcr(ji,jj,jk)
              tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2ut * ( zproregn(ji,jj,jk) + zproregd(ji,jj,jk)   &
                 &                + zproregp(ji,jj,jk) + zproregdz(ji,jj,jk) + zproregcr(ji,jj,jk) ) + ( o2ut + o2nit ) &
                 &                * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) + zpronewp(ji,jj,jk) + zpronewdz(ji,jj,jk) &
                 &                + zpronewcr(ji,jj,jk)  )   &
                 &                + ( o2ut + o2nit ) * ( znfix(ji,jj,jk) * 2.0/3.0) + o2nit * ( znfix(ji,jj,jk) / 3.0) &
                 &                - o2ut * ( zrespn(ji,jj,jk) + zrespp(ji,jj,jk) + zrespd(ji,jj,jk) + zrespdz(ji,jj,jk)  &
                 &                + zrespcr(ji,jj,jk)) 
              zfeup = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk) + texcretp * zprofep(ji,jj,jk)   &
                 &                + texcretdz * zprofedz(ji,jj,jk) + texcretcr * zprofecr(ji,jj,jk)
              tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zfeup
              tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - zprmaxd(ji,jj,jk) * zysopt(ji,jj,jk) * rfact2 * trb(ji,jj,jk,jpdia)
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk) - zprorcap(ji,jj,jk) &
              &                     - zprorcadz(ji,jj,jk) - zprorcacr(ji,jj,jk) &
              &                     + zpsino3 * zpronewn(ji,jj,jk) + zpsinh4 * zproregn(ji,jj,jk)   &
              &                     + zpsino3 * zpronewp(ji,jj,jk) + zpsinh4 * zproregp(ji,jj,jk)   &
              &                     + zpsino3 * zpronewd(ji,jj,jk) + zpsinh4 * zproregd(ji,jj,jk)  &
              &                     + zpsino3 * zpronewdz(ji,jj,jk)  + zpsinh4 * zproregdz(ji,jj,jk) &
              &                     + zpsino3 * zpronewcr(ji,jj,jk)  + zpsinh4 * zproregcr(ji,jj,jk) &
              &                     + zpsinfix * znfix(ji,jj,jk)         &
              &                     + zrespn(ji,jj,jk) + zrespd(ji,jj,jk) + zrespp(ji,jj,jk) + zrespdz(ji,jj,jk)  &
              &                     + zrespcr(ji,jj,jk)
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk)  &
              &                     + zpronewp(ji,jj,jk) + zpronewdz(ji,jj,jk) + zpronewcr(ji,jj,jk) )      &
              &                     - rno3 * ( zproregn(ji,jj,jk) + zproregd(ji,jj,jk)     &
              &                     + zproregp(ji,jj,jk) + zproregdz(ji,jj,jk) + zproregcr(ji,jj,jk))
          END DO
        END DO
     END DO
     
     ! Production and uptake of ligands by phytoplankton. This part is activated 
     ! when ln_ligand is set to .true. in the namelist. Ligand uptake is small 
     ! and based on the FeL model by Morel et al. (2008) and on the study of
     ! Shaked and Lis (2012)
     ! -------------------------------------------------------------------------
     IF( ln_ligand ) THEN
         zpligprod1(:,:,:) = 0._wp    ;    zpligprod2(:,:,:) = 0._wp
         DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji =1 ,jpi
                 zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk) + excretp * zprorcap(ji,jj,jk) &
                 &          + excretdz * zprorcadz(ji,jj,jk)
                 zfeup    = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk) + texcretp * zprofep(ji,jj,jk) &
                 &          + texcretdz * zprofedz(ji,jj,jk)
                 tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + zdocprod * ldocp    &
                 &       - zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) + 2.E3 * (1.0 - plig(ji,jj,jk) ) ) * lthet
                 zpligprod1(ji,jj,jk) = zdocprod * ldocp
                 zpligprod2(ji,jj,jk) = zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) &
                 &                      + 2.E3 * (1.0 - plig(ji,jj,jk) ) ) * lthet
              END DO
           END DO
        END DO
     ENDIF

    ! Total primary production per year
    IF( iom_use( "tintpp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
      & tpp = glob_sum( 'p6zprod', ( zprorcan(:,:,:) + zprorcad(:,:,:) + zprorcap(:,:,:) + zprorcadz(:,:,:) &
      &       + zprorcacr(:,:,:)  ) * cvol(:,:,:) )
    IF( iom_use( "tnfix" ) .OR.  ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  ) &
      & totnfix = glob_sum( 'p6zprod', znfix(:,:,:) * cvol(:,:,:) )    

    IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          ALLOCATE( zw2d(jpi,jpj), zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "PPPHYN" ) .OR. iom_use( "PPPHYD" ) .OR. iom_use( "PPPHYP" ) .OR. iom_use( "PPPHYDZ") .OR. iom_use( "PPPHYCR") )  THEN
              zw3d(:,:,:) = zprorcan(:,:,:) * zfact * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "PPPHYN"  , zw3d )
              !
              zw3d(:,:,:) = zprorcap(:,:,:) * zfact * tmask(:,:,:)  ! primary production by picophyto
              CALL iom_put( "PPPHYP"  , zw3d )
              !
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diatoms
              CALL iom_put( "PPPHYD"  , zw3d )
              !
              zw3d(:,:,:) = zprorcadz(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diazotrophs
              CALL iom_put( "PPPHYDZ"  , zw3d )
              !
              zw3d(:,:,:) = zprorcacr(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diazotroph2
              CALL iom_put( "PPPHYCR"  , zw3d )
          ENDIF
          IF( iom_use( "PPNEWN" ) .OR. iom_use( "PPNEWD" ) .OR. iom_use( "PPNEWP" )) THEN ! .OR. iom_use( "PPNEWDZ") )  THEN
              zw3d(:,:,:) = zpronewn(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by nanophyto
              CALL iom_put( "PPNEWN"  , zw3d )
              !
              zw3d(:,:,:) = zpronewp(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by picophyto
              CALL iom_put( "PPNEWP"  , zw3d )
              !
              zw3d(:,:,:) = zpronewd(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diatoms
              CALL iom_put( "PPNEWD"  , zw3d )
              !
              zw3d(:,:,:) = zpronewdz(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diazotrophs
              CALL iom_put( "PPNEWDZ"  , zw3d ) 
              !
              zw3d(:,:,:) = zpronewcr(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diazotroph2
              CALL iom_put( "PPNEWCR"  , zw3d )
          ENDIF
          IF( iom_use( "PBSi" ) )  THEN
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:) * zysopt(:,:,:) ! biogenic silica production
              CALL iom_put( "PBSi"  , zw3d )
          ENDIF
          IF( iom_use( "PFeN" ) .OR. iom_use( "PFeD" ) .OR. iom_use( "PFeP" ) ) THEN !.OR. iom_use( "PFeDZ") )  THEN
              zw3d(:,:,:) = zprofen(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron uptake by nanophyto
              CALL iom_put( "PFeN"  , zw3d )
              !
              zw3d(:,:,:) = zprofep(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron uptake by picophyto
              CALL iom_put( "PFeP"  , zw3d )
              !
              zw3d(:,:,:) = zprofed(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron uptake by  diatoms
              CALL iom_put( "PFeD"  , zw3d )
              !
              zw3d(:,:,:) = zprofedz(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron uptake by  diazotrophs
              CALL iom_put( "PFeDZ"  , zw3d )
              !
              zw3d(:,:,:) = zprofecr(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron uptake by  diazotroph2
              CALL iom_put( "PFeCR"  , zw3d )
          ENDIF
          IF( iom_use( "LPRODP" ) )  THEN
              zw3d(:,:,:) = zpligprod1(:,:,:) * 1e9 * zfact * tmask(:,:,:)
              CALL iom_put( "LPRODP"  , zw3d )  ! Ligand production by phytoplankton
          ENDIF
          IF( iom_use( "LDETP" ) )  THEN
              zw3d(:,:,:) = zpligprod2(:,:,:) * 1e9 * zfact * tmask(:,:,:)
              CALL iom_put( "LDETP"  , zw3d )  ! Uptake of ligands by phytoplankton
          ENDIF
          IF( iom_use( "Mumax" ) )  THEN
              zw3d(:,:,:) = zprmaxn(:,:,:) * tmask(:,:,:)   ! Maximum growth rate
              CALL iom_put( "Mumax"  , zw3d )
          ENDIF
          IF( iom_use( "MuN" ) .OR. iom_use( "MuD" ) .OR. iom_use( "MuP" )) THEN ! .OR. iom_use("MuDz") )  THEN
              zw3d(:,:,:) = zprbio(:,:,:) * xlimphy(:,:,:) * tmask(:,:,:)  ! Realized growth rate for nanophyto
              CALL iom_put( "MuN"  , zw3d )
              !
              zw3d(:,:,:) = zprpic(:,:,:) * xlimpic(:,:,:) * tmask(:,:,:)  ! Realized growth rate for picophyto
              CALL iom_put( "MuP"  , zw3d )
              !
              zw3d(:,:,:) =  zprdia(:,:,:) * xlimdia(:,:,:) * tmask(:,:,:)  ! Realized growth rate for diatoms
              CALL iom_put( "MuD"  , zw3d )
              !
              zw3d(:,:,:) =  zprdiaz(:,:,:) * xlimdiaz(:,:,:) * tmask(:,:,:) !Realized growth rate for diazotrophs
              CALL iom_put( "MuDz"  , zw3d )
              !
              zw3d(:,:,:) =  zprcroc(:,:,:) * xlimcroc(:,:,:) * tmask(:,:,:) !Realized growth rate for diazotroph2
              CALL iom_put( "MuCr"  , zw3d )
          ENDIF
          IF( iom_use( "LNlight" ) .OR. iom_use( "LDlight" ) .OR. iom_use( "LPlight" ) )  THEN
              zw3d(:,:,:) = zprbio (:,:,:) / (zprmaxn(:,:,:) + rtrn) * tmask(:,:,:) ! light limitation term of nanophytoplankton
              CALL iom_put( "LNlight"  , zw3d )
              !
              zw3d(:,:,:) = zprpic (:,:,:) / (zprmaxp(:,:,:) + rtrn) * tmask(:,:,:) ! light limitation term of picophytoplankton
              CALL iom_put( "LPlight"  , zw3d )
              !
              zw3d(:,:,:) =  zprdia (:,:,:) / (zprmaxd(:,:,:) + rtrn) * tmask(:,:,:)  ! light limitation term of diatoms
              CALL iom_put( "LDlight"  , zw3d )

              zw3d(:,:,:) =  zprdiaz (:,:,:) / (zprmaxdz(:,:,:) + rtrn) * tmask(:,:,:)  ! light limitation term of diazo
              CALL iom_put( "LDZlight"  , zw3d )

              zw3d(:,:,:) =  zprcroc (:,:,:) / (zprmaxcr(:,:,:) + rtrn) * tmask(:,:,:)  ! light limitation term of diazo2
              CALL iom_put( "LCRlight"  , zw3d )
          ENDIF
          IF( iom_use( "MunetN" ) .OR. iom_use( "MunetD" ) .OR. iom_use( "MunetP" )) THEN ! .OR. iom_use("MunetDZ") )  THEN
              zw3d(:,:,:) = zcroissn(:,:,:) * tmask(:,:,:) ! ! Realized growth rate for nanophyto
              CALL iom_put( "MunetN"  , zw3d )
              !
              zw3d(:,:,:) = zcroissp(:,:,:) * tmask(:,:,:) ! ! Realized growth rate for picophyto
              CALL iom_put( "MunetP"  , zw3d )
              !
              zw3d(:,:,:) = zcroissd(:,:,:) * tmask(:,:,:) ! ! Realized growth rate for diatoms
              CALL iom_put( "MunetD"  , zw3d )
              !
              zw3d(:,:,:) = zcroissdz(:,:,:) * tmask(:,:,:) ! ! Realized growth rate for diazotroph
              CALL iom_put( "MunetDZ"  , zw3d )
              !
              zw3d(:,:,:) = zcroisscr(:,:,:) * tmask(:,:,:) ! ! Realized growth rate for diazotroph2
              CALL iom_put( "MunetCR"  , zw3d )

          ENDIF
            IF( iom_use("Nfix"   ) ) CALL iom_put( "Nfix", znfix(:,:,:) * rno3 * zfact * tmask(:,:,:) )  ! nitrogen fixation
            IF( iom_use("Nfix_DZ"   ) ) CALL iom_put( "Nfix_DZ", znfix_dz(:,:,:) * rno3 * zfact * tmask(:,:,:) )  ! nitrogen fixation diazo
            IF( iom_use("Nfix_CR"   ) ) CALL iom_put( "Nfix_CR", znfix_cr(:,:,:) * rno3 * zfact * tmask(:,:,:) )  ! nitrogen fixation diazo2

          IF( iom_use( "tintpp" ) )  CALL iom_put( "tintpp" , tpp * zfact )  !  global total integrated primary production molC/s
          !
          IF( iom_use( "tnfix" ) )  CALL iom_put( "tnfix" , totnfix * zfact * rno3 )
         
          IF( iom_use("zprmaxn") ) CALL iom_put( "zprmaxn", zprmaxn(:,:,:) *tmask(:,:,:) )
!          IF( iom_use("zprmaxdz") ) CALL iom_put( "zprmaxdz", zprmaxdz(:,:,:) *tmask(:,:,:) )
          !
            IF( iom_use("INTNFIX") ) THEN   ! nitrogen fixation rate in ocean (vertically integrated )
               zwork(:,:) = 0.
               DO jk = 1, jpkm1
                 zwork(:,:) = zwork(:,:) + znfix(:,:,jk) * rno3 * zfact * e3t_n(:,:,jk) * tmask(:,:,jk)
               ENDDO
               CALL iom_put( "INTNFIX" , zwork )
            ENDIF
            IF( iom_use("INTNFIX_DZ") ) THEN   ! nitrogen fixation rate of diazo in ocean (vertically integrated )
               zwork(:,:) = 0.
               DO jk = 1, jpkm1
                 zwork(:,:) = zwork(:,:) + znfix_dz(:,:,jk) * rno3 * zfact * e3t_n(:,:,jk) * tmask(:,:,jk)
               ENDDO
               CALL iom_put( "INTNFIX_DZ" , zwork )
            ENDIF
            IF( iom_use("INTNFIX_CR") ) THEN   ! nitrogen fixation rate of diazo2 in ocean (vertically integrated )
               zwork(:,:) = 0.
               DO jk = 1, jpkm1
                 zwork(:,:) = zwork(:,:) + znfix_cr(:,:,jk) * rno3 * zfact * e3t_n(:,:,jk) * tmask(:,:,jk)
               ENDDO
               CALL iom_put( "INTNFIX_CR" , zwork )
            ENDIF
          DEALLOCATE( zw2d, zw3d )
       ENDIF
     ENDIF

          IF( iom_use( "PUptkDZ"   ) ) CALL iom_put( "PUptkDZ"  ,puptkdz(:,:,:) * tmask(:,:,:) )  ! P uptake diazotrophs
          IF( iom_use( "PO4UptkDZ"   ) ) CALL iom_put( "PO4UptkDZ" ,zpropo4dz(:,:,:) * tmask(:,:,:) )  ! P uptake diazotrophs
           IF( iom_use( "DOPUptkDZ"   ) ) CALL iom_put( "DOPUptkDZ"  , zprodopdz(:,:,:) * tmask(:,:,:) )  ! DOP uptake diazotrophs
          IF( iom_use("zprnfix") ) CALL iom_put( "zprnfix", zprnfix(:,:,:) *tmask(:,:,:) )
          IF( iom_use("vnfmax_dz") ) CALL iom_put( "vnfmax_dz", vnfmax_dz(:,:,:)* tmask(:,:,:) )
          IF( iom_use("qn_diazo") ) CALL iom_put( "qn_diazo", qn_diazo(:,:,:) *tmask(:,:,:) ) 
          IF( iom_use("nutlim_dz") ) CALL iom_put( "nutlim_dz", nutlim_dz(:,:,:)* tmask(:,:,:) )
          IF( iom_use("cost_nfix") ) CALL iom_put( "cost_nfix",(zpsinfix*znfix(:,:,:)) * tmask(:,:,:) )
          IF( iom_use("cost_nh4") ) CALL iom_put( "cost_nh4",(zpsinh4*zproregn(:,:,:)) * tmask(:,:,:) ) 
          IF( iom_use("zprnutdz") ) CALL iom_put( "zprnutdz",(zprnutdz(:,:,:)) * tmask(:,:,:) )
          IF( iom_use("zprmaxdz") ) CALL iom_put( "zprmaxdz",(zprmaxdz(:,:,:)) * tmask(:,:,:) )
          IF( iom_use("Facul") ) CALL iom_put( "Facul",(Facul_out(:,:,:)) * tmask(:,:,:) )
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_prod')
      !
   END SUBROUTINE p6z_prod


   SUBROUTINE p6z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p6z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the namp6zprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp6zprod
      !!----------------------------------------------------------------------
      INTEGER :: ios    ! Local integer output status for namelist read
      !!
      NAMELIST/namp6zprod/ pislopen, pislopep, pisloped, excretn, excretp, excretd,     &
         &                 thetannm, thetanpm, thetandm, chlcmin, grosip, bresp, xadap,  &
         &                 pislopedz, excretdz, thetandzm, Facul_max, zmxlcoef,                  &
         &                 pislopecr, excretcr, thetancrm, photoinhcr
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist namp6zprod in reference namelist : Pisces phytoplankton production
      READ  ( numnatp_ref, namp6zprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp6zprod in reference namelist' )

      REWIND( numnatp_cfg )              ! Namelist namp6zprod in configuration namelist : Pisces phytoplankton production
      READ  ( numnatp_cfg, namp6zprod, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp6zprod in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, namp6zprod'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    mean Si/C ratio                           grosip       =', grosip
         WRITE(numout,*) '    P-I slope                                 pislopen     =', pislopen
         WRITE(numout,*) '    P-I slope  for diatoms                    pisloped     =', pisloped
         WRITE(numout,*) '    P-I slope  for picophytoplankton          pislopep     =', pislopep
         WRITE(numout,*) '    Acclimation factor to low light           xadap        =', xadap
         WRITE(numout,*) '    excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(numout,*) '    excretion ratio of picophytoplankton      excretp      =', excretp
         WRITE(numout,*) '    excretion ratio of diatoms                excretd      =', excretd
         WRITE(numout,*) '    basal respiration in phytoplankton        bresp        =', bresp
         WRITE(numout,*) '    Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(numout,*) '    Minimum Chl/N in nanophytoplankton        thetannm     =', thetannm
         WRITE(numout,*) '    Minimum Chl/N in picophytoplankton        thetanpm     =', thetanpm
         WRITE(numout,*) '    Minimum Chl/N in diatoms                  thetandm     =', thetandm
         WRITE(numout,*) '    P-I slope diazotrophspislopedz    =', pislopedz
         WRITE(numout,*) '    excretion ratio of diazotrophs            excretdz=', excretdz
         WRITE(numout,*) '    Minimum Chl/N in diazotrophsthetandzm    =', thetandzm
         WRITE(numout,*) '    limit on diazos facultative ability      =', Facul_max
         WRITE(numout,*) '    zmxlcoef                                 =', zmxlcoef
         WRITE(numout,*) '    P-I slope diazotroph2                  pislopecr    =', pislopecr
         WRITE(numout,*) '    excretion ratio of diazotroph3          excretcr=', excretcr
         WRITE(numout,*) '    Minimum Chl/N in diazotroph2           thetancrm    =', thetancrm
         WRITE(numout,*) '    Photoinhibition of Croco               photoinhcr   =', photoinhcr
      ENDIF
      !
      r1_rday   = 1._wp / rday 
      texcretn  = 1._wp - excretn
      texcretp  = 1._wp - excretp
      texcretd  = 1._wp - excretd
      tpp       = 0._wp
      texcretdz = 1._wp - excretdz
      texcretcr = 1._wp - excretcr
      totnfix = 0._wp
      !
   END SUBROUTINE p6z_prod_init


   INTEGER FUNCTION p6z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( zdaylen(jpi,jpj), STAT = p6z_prod_alloc )
      !
      IF( p6z_prod_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p6z_prod_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p6z_prod_alloc
   !!======================================================================
END MODULE p6zprod
