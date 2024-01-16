SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,           &
     DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, &
     PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,  &
     COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,  &
     KSPT, KSTEP, KINC)
  
    
  !ARREGLOS EXTERNOS CODE ASTER
  !----------------------------------------------------------------------------------------------------------------------------------------
  INTEGER                         :: NDI
  INTEGER                         :: NSHR
  INTEGER                         :: NTENS
  INTEGER                         :: NSTATV
  INTEGER                         :: NPROPS
  INTEGER                         :: NOEL
  INTEGER                         :: NPT
  INTEGER                         :: LAYER
  INTEGER                         :: KSPT
  INTEGER                         :: KSTEP
  INTEGER                         :: KINC
  REAL(8)                         :: SSE
  REAL(8)                         :: SPD
  REAL(8)                         :: SCD
  REAL(8)                         :: RPL
  REAL(8)                         :: DTIME
  REAL(8)                         :: TEMP
  REAL(8)                         :: DTEMP
  REAL(8)                         :: PNEWDT
  REAL(8)                         :: DRPLDT
  REAL(8)                         :: CELENT
  REAL(8), DIMENSION(1)           :: PREDEF
  REAL(8), DIMENSION(1)           :: DPRED
  REAL(8), DIMENSION(2)           :: TIME
  REAL(8), DIMENSION(3)           :: COORDS
  REAL(8), DIMENSION(3,3)         :: DROT
  REAL(8), DIMENSION(3,3)         :: DFGRD0
  REAL(8), DIMENSION(3,3)         :: DFGRD1
  REAL(8), DIMENSION(NTENS)       :: STRESS
  REAL(8), DIMENSION(NTENS)       :: DDSDDT
  REAL(8), DIMENSION(NTENS)       :: DRPLDE
  REAL(8), DIMENSION(NTENS)       :: STRAN
  REAL(8), DIMENSION(NTENS)       :: DSTRAN
  REAL(8), DIMENSION(NPROPS)      :: PROPS
  REAL(8), DIMENSION(NSTATV)      :: STATEV
  REAL(8), DIMENSION(NTENS,NTENS) :: DDSDDE
  CHARACTER(LEN = 80)             :: CMNAME    
  

  !ARREGLOS LOCALES
  !----------------------------------------------------------------------------------------------------------------------------------------
  real(8), dimension(6)            :: dstress              !incremento de la tension 
  real(8), dimension(6)            :: dstranE              !incremento de deformacion elastica
  real(8), dimension(6)            :: stranE               !deformacion elastica acumulada
  real(8), dimension(6)            :: stranP               !deformacion plastica acumulada
  real(8), dimension(6)            :: strain               !deformacion total
  real(8), dimension(6)            :: stranCreep = 0.d0    !deformacion plastica debido al creep
  real(8), dimension(6)            :: stranSwell = 0.d0    !deformacion plastica debido al swelling
  real(8), dimension(6)            :: stranPlast = 0.d0    !deformacion plastica propiamente dicha
  real(8), dimension(6)            :: deviatoric           !tensor de tensiones deviatorio
  real(8), dimension(6)            :: n                    !normal a la sup de fluencia vectorizado
  real(8), dimension(6)            :: flow                 !normal a la sup de fluencia vectorizado
  real(8), dimension(6,6)          :: dN                   !derivada de la normal respecto de la tension
  real(8), dimension(6)            :: T2                   !Tensor identidad vectorizado
  real(8), dimension(6,6)          :: T4                   !Tensor identidad de cuarto orden
  real(8), dimension(6)            :: dstranCreep          !incremento de deformacion por irradiacion vector
  real(8), dimension(6)            :: dstranSwell          !incremento de deformacion por sweelling vector
  real(8), dimension(6)            :: dstranP              !incremento de deformacion plastica vector
  real(8),save                     :: dp                   !incremento del multiplicador plastico
  real(8),save                     :: dpi                  !incremento del miltiplicador creep por irradiacion 
  real(8),save                     :: dg                   !incremento del multiplicador sweeling
  real(8),save                     :: deta                 !incremento umbral de irradiacion 
  real(8),save                     :: Ai                   !coeficiente material Ai (eq 4.3-1 de R5.23.03)
  real(8),save                     :: Ag                   !coeficiente material Ag (eq 4.4-3 de R5.23.03)
  real(8),save                     :: Cf                   !funcion zeta f  (eq 4.3-1 de R5.23.03)
  real(8)                          :: R00                  !coeficiente material R (eq 4.4-4 de R5.23.03)
  real(8),save                     :: S                  
  real(8)                          :: sigmap               !pendiente plastica
  real(8),save                     :: eqplas               !incremento plastico
  
  
  
  
  !Parametros para calcular tension de fluencia
  real(8),parameter                :: ZERO = 0.D0
  real(8),parameter                :: ONE = 1.D0
  real(8),parameter                :: TWO = 2.D0
  real(8),parameter                :: THREE = 3.D0
  real(8),parameter                :: SIX = 6.D0
  real(8),parameter                :: CE1 = -8.1220E01
  real(8),parameter                :: CP1 = 7.2000E-05
  real(8),parameter                :: SYSU0 = 2.0000E02
  real(8),parameter                :: SYS1 = 8.0000E02
  real(8),parameter                :: CY0 = 1.5300E-03   
  real(8),parameter                :: CY1 = 1.1440E-06
  real(8),parameter                :: SUTSU0 = 4.5000E02
  real(8),parameter                :: SUTSU1 = 8.1000E02
  real(8),parameter                :: r = 0.01d0
  real(8),parameter                :: CU0 = 1.5300E-03
  real(8),parameter                :: CU1 = 1.1440E-06
  real(8),parameter                :: toler = 1.0D-6
  real(8),parameter                :: E0=2.4868E-01
  real(8),parameter                :: fue0 = 1.2105E0
  real(8),parameter                :: fue1 = -1.2001E0
  real(8),parameter                :: rue = 1.4517E-01
  real(8),parameter                :: gue0 = 1.0000E0
  real(8),parameter                :: gue1 = -9.8750E-01
  real(8),parameter                :: cue0 = 5.5869E-01
  real(8),parameter                :: cue1 = -1.8930E-03
  real(8),parameter                :: cue2 = 5.3656E-06
  real(8),parameter                :: cue3 = -4.5779E-09
  real(8),parameter                :: ftu0 = -9.8861E-01
  real(8),parameter                :: ftu1 = 1.1433E0
  real(8),parameter                :: ftu2 = 1.3991E0
  real(8),parameter                :: RTU1=1.4517E-01
  real(8),parameter                :: RTU2=1.5102E-01
  real(8),parameter                :: GTU0=8.9000E0
  real(8),parameter                :: GTU1=-7.4000E0
  real(8),parameter                :: GTU2=-7.9000E0
  real(8),parameter                :: CTE0=7.1632E-01
  real(8),parameter                :: CTE1=-1.9560E-03
  real(8),parameter                :: CTE2=5.7562E-06
  real(8),parameter                :: CTE3=-7.1266E-09
  real(8),parameter                :: CTE4=3.2172E-012
  
  
  !Logicas
  logical,save                    :: iniciar = .true.            !inicializa variables de estado
  !Transferidos via PROPS
  REAL(8)                         :: E                           !Modulo de elasticidad
  REAL(8)                         :: NU                          !Coeficiente de Poisson
  REAL(8)                         :: IRRA                        !Dosis de irradiacion [dpa]
  REAL                            :: XN                          !Exponente ley de endurecimiento
  
  
  
  
  !Newton-Rapson
  real(8),dimension(10)           :: XresV                       ! Residuo para las variables de estado
  real(8),dimension(6)            :: Xres                        ! Residuo para determinar tension de correccion
  real(8),dimension(10,10)        :: JacobV                      ! Jacobiano para resolver las variables de estado
  real(8),dimension(6,6)          :: Jacob                       ! Jacobiano para determinar la tension de correccion
  real(8),dimension(10)           :: dVs                         ! Incremento de las variables de estado
  real(8),dimension(10)           :: nVs                         ! valor nuevo de las variables de estado
  real(8),dimension(10),save      :: vVs                         ! valor viejo de las variables de estado
  real(8),parameter               :: tolerNR = 1.0d-03           ! Tolerancia admisible para N - R
  real(8),dimension(10)           :: Coef                        ! Coeficiente para NR, busquedad lineal
  real(8),dimension(3,3)          :: X                           ! Auxiliar
  real(8),dimension(6)            :: Xe                          ! Auxiliar residuo elastico
  real(8),dimension(6,6)          :: Xj                          ! Auxiliar jacobiano
  real(8),dimension(6)            :: q                           ! derivada parcial de residuo plastico respecto a incremento elastico
  real(8),dimension(6)            :: reta                        ! derivada parcial de residuo eta respecto a incremento elastico
  real(8),dimension(6)            :: gcreep                      ! derivada parcial de residuo creep respecto a incremento elastico
  integer                         :: iNR                         ! Numero de iteraciones
  integer,parameter               :: iNRmax = 100                ! Numero de iteraciones maximas
  real(8),dimension(1)            :: errNR                       ! Error N-R
  
  
  
  !Dummy
  real(8), dimension(6,6)         :: dum66   = 0.d0
  real(8), dimension(3,3,3,3)     :: dum3333 = 0.d0
  integer,parameter               :: UOC3=2001
  integer,parameter               :: UOC4=2002
  integer                         :: IER
  integer,save                    :: plastize = 0
  real(8),dimension(6,6)          :: dndT
  logical                         :: prim = .true.
  real(8)                         :: nn,inicio_int,fin_int
  real(8)                         :: K
  real(8)                         :: sy
  real(8)                         :: su
  real(8)                         :: epsi
  real(8)                         :: p0
  real(8)                         :: a
  real(8)                         :: spe
  real(8)                         :: n1
  real(8)                         :: pe
  
  
  !INICIALIZACION
  !----------------------------------------------------------------------------------------------------------------------------------------
   DDSDDE = 0.0D0
   
   E = PROPS(1)
   NU = PROPS(2)
   IRRA=PROPS(3)
   Ai = PROPS(4)
   Cf = PROPS(5)
   Cg = PROPS(6)
   R00 = PROPS(7)
   

   EMOD=E+CE1*TEMP
   ENU=NU+CP1*TEMP

   
   EBULK3=EMOD/(ONE-TWO*ENU)
   EG2=EMOD/(ONE+ENU)
   EG=EG2/TWO
   EG3=THREE*EG
   ELAM=(EBULK3-EG2)/THREE
  


! Tensor constitutivo Elstico (3x3) 

  do k1=1 , NDI
   do k2=1 , NDI
     DDSDDE(k2,k1) = ELAM
   end do
     DDSDDE(k1,k1) = EG2+ELAM
  end do
  do k1=NDI+1, NTENS
     DDSDDE(k1,k1) = EG
  end do 
  

!  Evalua una estimacion elastica de la tension    
  dstress = matmul(DDSDDE,dstran)    !Incremento de la tension elastica
  stress = stress + dstress          !Tension de prediccion en tiempo t

  
!!  calcula la tension de fluencia  
  s2 = sysu0 + two*(sys1-sysu0)*r - (sys1 - sysu0)*r**2
  s1 = s2 + (sys1 - s2)*(1 - exp(-irra/3.0000e0))
  sy = s1*exp((-cy0 + cy1*s1)*(temp - 330))


!! calcula maximo valor de tension   
  su2 = SUTSU0 + two*(SUTSU1 - SUTSU0 )*r - (SUTSU1 - SUTSU0) * r**2
  su1 = su2  +(SUTSU1-su2)*(1-exp(-irra/3.0000e0))
  su  = su1*exp((-CU0 + CU1*su1)*(temp-330))

  
!! diferencia entre la tension maxima y tension de fluencia
  dsusy = (su - sy)*exp(0.002/E0)

!! recuperar la tension de fluencia 
  sy = su - dsusy



!!   strain for the ultimate stress
  etaue = fue0 + fue1*(1 - exp(-r/rue))
  xiue= gue0 + gue1*(1-exp(-irra/1.d0))
  stranue1 = cue0 + cue1*temp + cue2*temp**2 + cue3*temp**3         
  stranue = etaue*xiue*stranue1

  
!! ultimate strain   
  etatu = ftu0 + ftu1*(1-exp(-r/rtu1)) + ftu2*exp(-r/rtu2)
  xitu = gtu0 + gtu1*(1-exp(-irra/2.5000e0)) + gtu2*exp(-irra/1.0000e0)
  epsite1 = cte0 + cte1*temp + cte2*temp**2 + cte3*temp**3 + cte4*temp**4
  deltaepsi = etatu*xitu*(epsite1 - stranue1)
  epsi = stranue + deltaepsi

! coeficiente de endurecimiento 
  xn = log(one + epsi)
  
  
! tension equivalente de  von mises
 smises = (stress(1) - stress(2))**2 + (stress(2) - stress(3))**2 +(stress(3) - stress(1))**2
 smises = smises + six * stress(4) ** 2 + six * stress(5) ** 2 + six * stress(6) ** 2
 smises = sqrt(smises/two)
 
  
!Obtenemos los parametros de Plasticidad a partir de ley IRRAD3M
call irrmat(sy,epsi,su,Ai,722.d0,R00,1.d0,1.d0,0.8d0,irra,k,p0,a,spe,pe,n1)

! condicion para determinar si comienza a plastificar
 if(smises > (1.d0+toler)*spe)then
 plastize = 1
 else
 plastize = 0
 end if

! Tensor identidad de segundo orden vectorizado
  T2(1:3) = 1.d0
  T2(4:6) = 0.d0

! Tensor identidad de cuarto orden
  do k1=1,6
  T4(k1,k1) = 1.d0
  end do
  
! Normal a la superficie de fluencia
 shydro = (stress(1) + stress(2) + stress(3)) / 3.D0
 n(1:3) = (stress(1:3)-shydro)/smises
 n(4:6) = stress(4:6) / smises
 flow = n
 n = three/two * n
 

 
 
! tensor dN = dn / dstressC

 do k1 = 1, 6
   do i = 1, 6
      dndT(k1, i) = 0.0
      do j = 1, 6
         if ((k1 .le. 3) .and. (i .le. 3)) then
            if (i == j) then
               if (i == k1) then
                  dndT(k1, i) = dndT(k1, i) + (1.0 - 1.0/3.0) / smises
               end if
            end if
         else
            if (i == j) then
               if (i == k1) then
                  dndT(k1, i) = dndT(k1, i) + 1.0 / smises
               end if
            end if
         end if
      end do
   end do
end do

dN = dndT


! Esquema Newton Raphson
!------------------------------------------------------------------------------------------------------------------------------------
 !! recuperar valores del incremento anterior, sino inicializar 
 vVs(1:6) = statev(1:6)
 vVs(7) = statev(14)
 vVs(8) = statev(15)
 vVs(9) = statev(16)
 vVs(10) = statev(17)

 iNR = 0
 
!- Coeficientes para N-R 
 Coef = 1.d0

 NEWTON: do
 
  iNR = iNR +1
  
   
  
!! inicializar variables de estado
 if (iniciar) then
 iniciar = .false.
 vVs(1) = 0.d0 !incremento de def elastica
 vVs(2) = 0.d0
 vVs(3) = 0.d0
 vVs(4) = 0.d0
 vVs(5) = 0.d0
 vVs(6) = 0.d0
 vVs(7) = 0.d0   !mult plastico
 vVs(8) =  Cf * smises * IRRA  ! incremento eta 
 vVs(9) =  0.d0 ! mult creep
 vVs(10) =  0.d0     ! mult swelling 
 end if


  
! armar el vector Residuo en virtud a tabla 3.1 de integracion implicita de R5.03.23

!-Umbral 
  umbrals = Cf * smises * IRRA
  umbral = 722 ![Mpa*dpa]  
 
  !Cálculo del residuo elástico
IF (plastize == 0 .AND. Ai == 0.0d0 .AND. R00 == 0.0d0) THEN
    Xe = vVs(1:6) - dstran
ELSE IF (plastize == 0 .AND. Ai == 0.0d0 .AND. R00 /= 0.0d0) THEN
    Xe = vVs(1:6) + vVs(10) * T2 - dstran
ELSE IF (plastize == 0 .AND. Ai /= 0.0d0 .AND. R00 == 0.0d0) THEN
    Xe = vVs(1:6) + vVs(9) * n  - dstran
ELSE IF (plastize == 0 .AND. Ai /= 0.0d0 .AND. R00 /= 0.0d0 .AND. umbrals > umbral) THEN
    Xe = vVs(1:6) + vVs(9) * n + vVs(10) * T2 - dstran
ELSE IF (plastize == 0 .AND. Ai /= 0.0d0 .AND. R00 /= 0.0d0 .AND. umbrals < umbral) THEN
    Xe = vVs(1:6) + vVs(10) * T2 - dstran
ELSE IF (plastize == 1 .AND. Ai == 0.0d0 .AND. R00 == 0.0d0) THEN
    Xe = vVs(1:6) + vVs(7) * n - dstran
ELSE IF (plastize == 1 .AND. Ai /= 0.0d0 .AND. R00 == 0.0d0 .AND. umbrals > umbral) THEN
    Xe = vVs(1:6) + vVs(7) * n + vVs(9) * n - dstran
ELSE IF (plastize == 1 .AND. Ai /= 0.0d0 .AND. R00 == 0.0d0 .AND. umbrals < umbral) THEN
    Xe = vVs(1:6) + vVs(7) * n - dstran
ELSE IF (plastize == 1 .AND. Ai /= 0.0d0 .AND. R00 /= 0.0d0 .AND. umbrals > umbral) THEN
    Xe = vVs(1:6) + vVs(7) * n + vVs(9) * n + vVs(10) * T2 - dstran
ELSE IF (plastize == 1 .AND. Ai /= 0.0d0 .AND. R00 /= 0.0d0 .AND. umbrals < umbral) THEN
    Xe = vVs(1:6) + vVs(7) * n + vVs(10) * T2 - dstran
ELSE IF (plastize == 1 .AND. Ai == 0.0d0 .AND. R00 /= 0.0d0) THEN
    Xe = vVs(1:6) + vVs(7) * n + vVs(10) * T2 - dstran
ELSE
    WRITE(*, *) 'Invalid combination of plastize, Ai, and R00'
END IF

  XresV(1:6) = Xe
 
!  residuo plastico
IF (plastize == 1) THEN
    XresV(7) = smises   - k*((vVs(7)+pe+p0)**n1)
ELSE IF (plastize == 0) THEN
    XresV(7) = vVs(7) - 0.0d0
ELSE
    WRITE(*, *) 'Invalid value of plastize'
END IF

  
! residuo eta
IF (umbrals > umbral) THEN
    XresV(8) = vVs(8) - (Cf/2.d0) * smises * IRRA
ELSE IF (umbrals < umbral) THEN
    XresV(8) = vVs(8) - 0.0d0
END IF
  
! residuo creep irradiacion
  if((umbral-umbrals) > 0.d0)then
    H = 0.d0
    else
    H = 1.d0
    end if
  XresV(9) = vVs(9) - H * Ai * smises * IRRA 

  
! residuo swelling
  R0 = R00 * Cg
  Ag = (one/three) * R0 * (one - (exp(one - IRRA)/(1 + exp(one - IRRA)))) * 1e-02
  XresV(10) = vVs(10) - Ag * IRRA 



  
! Matriz Jacobiana de 10x10 ---------------------------------------------------------
  ! derivada parcial de Residuo elastico
IF (plastize .eq. 0 .AND. Ai .eq. 0.0d0) THEN
    Xj = T4
ELSE IF (plastize .eq. 0 .AND. Ai .ne. 0.0d0 .AND. umbrals > umbral) THEN
    Xj = T4 + vVs(9) * matmul(dN, DDSDDE)
ELSE IF (plastize .eq. 0 .AND. Ai .ne. 0.0d0 .AND. umbrals < umbral) THEN
    Xj = T4
ELSE IF (plastize .eq. 1 .AND. Ai .eq. 0.0d0) THEN
    Xj = T4 + vVs(7) * matmul(dN, DDSDDE)
ELSE IF (plastize .eq. 1 .AND. Ai .ne. 0.0d0 .AND. umbrals > umbral) THEN
    Xj = T4 + vVs(7) * matmul(dN, DDSDDE) + vVs(9) * matmul(dN, DDSDDE)
ELSE IF (plastize .eq. 1 .AND. Ai .ne. 0.0d0 .AND. umbrals < umbral) THEN
    Xj = T4 + vVs(7) * matmul(dN, DDSDDE)
END IF

  JacobV(1:6,1:6) = Xj
  JacobV(1:6,7) = n
  JacobV(1:6,8) = 0.d0
  JacobV(1:6,9) = n
  JacobV(1:6,10) = T2
  
  !derivada parcial de Residuo plastico 
  q = matmul(flow,DDSDDE)
IF (plastize .eq. 1) THEN
    JacobV(7, 1:6) = q
    JacobV(7, 7) = - n1*k*((dVs(7)+pe+p0)**(n1-1.d0))
ELSE IF (plastize .eq. 0) THEN
    JacobV(7, 1:6) = 0.0d0
    JacobV(7, 7) = 1.0d0
END IF
  JacobV(7,8:10) = 0.d0
  Sf = k*((dVs(7)+pe+p0)**n1)
   
  !derivada parcial de Residuo eta
IF (umbrals < umbral) THEN
    JacobV(8, 1:6) = 0.0d0
ELSE IF (umbrals > umbral) THEN
    reta = -Cf * IRRA * matmul(n, DDSDDE)
    JacobV(8, 1:6) = reta
END IF
  JacobV(8,7) = 0.d0
  JacobV(8,8) = 1.d0
  JacobV(8,9:10) = 0.d0
  sf = k*((dVs(7)+pe+p0)**n1)
  !derivada parcial de Residuo creep
IF (umbrals < umbral) THEN
    JacobV(9, 1:6) = 0.0d0
ELSE IF (umbrals > umbral) THEN
    gcreep = -Cf * Ai * matmul(n, DDSDDE) * IRRA
    JacobV(9, 1:6) = gcreep
END IF

  JacobV(9,7:8) = 0.d0
  JacobV(9,9) = 1.d0
  JacobV(9,10) = 0.d0
  
  !derivada parcial de Residuo swelling
  JacobV(10,1:9) = 0.d0
  JacobV(10,10) = 1.d0
  
  

  
  !inversa del Jacobiano
  call LU_INVERSE(JacobV,10)
 
  
  !incremento de las variables de estado
  dVs = - Coef * matmul(JacobV,XresV)
  vVs = vVs + dVs


  ! Cálculo del error máximo
    errNR = abs(dVs(1))
    do i = 2, 10
        if (i /= 8) then
            errNR = max(errNR, abs(dVs(i)))
        end if
    end do




  
  if(errNR(1) < tolerNR ) exit NEWTON 

  if(iNR > iNRmax) then

  stop       

  end if
  

  
  enddo NEWTON
  

  
  dstranE = vVs(1:6)
  dp      = vVs(7) ! plastico
  dg      = vVs(10)! swelling
  dpi     = vVs(9) ! creep
  deta    = vVs(8) !umbral irradiacion
    

if(plastize .eq. 1)then  

! actualizo estado tensional con retorno radial
 stress(1:3) = flow(1:3) * Sf + shydro
 stress(4:6) = flow(4:6) * Sf
 
 ! Operador Tangente Consistente
call LU_INVERSE(DDSDDE,6)
if(Ai .ne. 0.d0 .and. umbrals > umbral)then
Jacob = DDSDDE + dVs(7) * dN  + dVs(9) * dN
else
Jacob = DDSDDE + dVs(7) * dN
end if

call LU_INVERSE(Jacob,6)
 DDSDDE = Jacob

end if



 


 
 !almacenar datos en el arreglo STATEV
 statev(1:ntens) = vVs(1:6)
 statev(14) = vVs(7)
 statev(15) = vVs(8)
 statev(16) = vVs(9)
 statev(17) = vVs(10)
 statev(7) = dp     !EPSEQ
 statev(8) = deta   !SEUIL
 statev(9) = dpi    !EPQIRRA
 statev(10) = dg    !GONF
 statev(11) = plastize !INDIPLAS   
 statev(12) = IRRA  !IRRA
 statev(13) = TEMP  !TEMP
 statev(18) = eqplas 
 
           
 
  CONTAINS
 
  INCLUDE 'library7F90.sub'

      SUBROUTINE LU_INVERSE (A,N)

! *** INVERTS A MATRIX USING LU DECOMPOSITION

      INTEGER                        :: N
      REAL(8), DIMENSION(N,N)        :: A,Y
	  INTEGER, DIMENSION(N)          :: INDX      ! MAY CHOKE SOME COMPILERS
!      DIMENSION A(5,5),Y(5,5),INDX(5)
      INTEGER                        :: I,J
	  REAL(8)                        :: AMAX,DUM,D
	  INTEGER                        :: ISINGULAR

!     write(*,*) 'A(i,j) matrix inside lu_inverse'
!     write(*,'(5e12.3)') ((a(i,j),j=1,n),i=1,n)
!     pause

! **************************************************************
! *** BLOCK ADDED 03/DEC/05 TO AVOID NUMERICALLY SINGULAR MATRIX
      AMAX=0.D0
      DO I=1,N
      DO J=1,N
        DUM=ABS(A(I,J))
        IF(DUM .GT. AMAX) AMAX=DUM
      ENDDO
      ENDDO
      DO I=1,N
      DO J=1,N
        A(I,J)=A(I,J)/AMAX      ! normalize the matrix
      ENDDO
      ENDDO
! **************************************************************

      DO I=1,N
        DO J=1,N
          Y(I,J)=0.
        ENDDO
        Y(I,I)=1.
      ENDDO

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)

      IF(ISINGULAR.EQ.1) THEN
        WRITE(*,*) ' *** SINGULAR MATRIX IN LU_INVERSE !!'
!        PAUSE
        STOP
      ENDIF

      DO J=1,N
        CALL LUBKSB(A,N,N,INDX,Y(1,J))
      ENDDO

      DO I=1,N
      DO J=1,N
        A(I,J)=Y(I,J) /AMAX      ! renormalize the inverse
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
   
   




! --------------------------------------------------------------------
! Copyright (C) 1991 - 2017 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine irrmat(r02,eu,rm,ai0,etais,rg0,alpha,phi0,kappa,irrad,k,p0,a,spe,pe,n1)
!
! person_in_charge: jean-luc.flejou at edf.fr
    implicit none
!#include "asterfort/assert.h"
!#include "asterfort/irrnvi.h"
!#include "asterfort/rcvalb.h"
!#include "asterfort/rcvarc.h"
!#include "asterfort/utmess.h"
!include 'utmess.F90'
!include 'rcvalb.F90'

    character(len=8) :: model
    character(len=3) :: matcst
!    character(len=*) :: fami
    integer :: imat, ndt, ndi, nr, nvi, kpg, ksp, iret, itmax
    real(kind=8) :: materd(30, 2), materf(30, 2), rela
!
!     ----------------------------------------------------------------
!     IRRAD3M   : RECUPERATION DU MATERIAU A T(TEMPD) ET T+DT(TEMPF)
!                 NB DE CMP DIRECTES/CISAILLEMENT , NB VAR. INTERNES
!                MATER(*,1) = E , NU , ALPHA
!     ----------------------------------------------------------------
!     IN  FAMI   :  FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!         KPG,KSP:  NUMERO DU (SOUS)POINT DE GAUSS
!         IMAT   :  ADRESSE DU MATERIAU CODE
!         MODEL  :  TYPE DE MODELISATION
!         NMAT   :  DIMENSION  DE MATER
!         ITMAX  :  NOMBRE D ITERATION MAX
!         RELA   :  TOLERANCE RELATIVE DES VALEURS MATERIAUX
!         VIND   :  VARIABLES INTERNES A T
!     OUT MATERD :  COEFFICIENTS MATERIAU A T
!         MATERF :  COEFFICIENTS MATERIAU A T+DT
!                   MATER(*,1) = CARACTERISTIQUES   ELASTIQUES
!                   MATER(*,2) = CARACTERISTIQUES   AUTRE
!         MATCST :  'OUI' SI  MATERIAU A T = MATERIAU A T+DT
!                   'NON' SINON
!         NDT    :  NB TOTAL DE COMPOSANTES TENSEURS
!         NDI    :  NB DE COMPOSANTES DIRECTES  TENSEURS
!         NR     :  NB DE COMPOSANTES SYSTEME NL
!         NVI    :  NB DE VARIABLES INTERNES, DANS LE SYSTEME NL
!     ----------------------------------------------------------------
    integer :: iterat, nbcara
!     NOMBRE DE PARAMETRES DE LA LOI : NBCARA
    parameter   (nbcara = 12 )
    integer :: cerr(nbcara)
    character(len=8) :: nomcir(nbcara)
    real(kind=8) :: mat(nbcara)
    character(len=8) :: nomcel(3)
!
    real(8) :: p0, irrad, irraf, pe, k, a, tempd, tempf
    real(8) :: r02, eu, rm, ai0, etais, rg0, alpha, phi0, kappa, zetag
    real(8) :: zetaf
    real(8) :: n0, n1, f0, f1, fe, pasn, exph, exp0, spe, coeffa
!
    real(kind=8) :: valrm(12)
    integer :: valim(2)
    character(len=10) :: valkm(2)
!
!    data pe    /2.0d-3/
!
    data nomcel /'E       ','NU      ','ALPHA   '/
!
    data nomcir /'R02     ','EPSI_U  ','RM      ','AI0     ',&
     &             'ETAI_S  ','RG0     ','ALPHA   ','PHI0    ',&
     &             'KAPPA   ','ZETA_F  ','ZETA_G  ','TOLER_ET'/
!
!     NOM                         a t-                 a t+ (t-+dt)
!     -------------------------------------------------------------
!     E                           MATERD(1,1)          MATERF(1,1)
!     NU                          MATERD(2,1)          MATERF(2,1)
!     ALPHA                       MATERD(3,1)          MATERF(3,1)
!
!     AI0                         MATERD(4,2)          MATERF(4,2)
!     ETAI_S                      MATERD(5,2)          MATERF(5,2)
!     AG                          MATERD(6,2)          MATERF(6,2)
!     K                           MATERD(7,2)          MATERF(7,2)
!     N                           MATERD(8,2)          MATERF(8,2)
!     P0                          MATERD(9,2)          MATERF(9,2)
!     KAPPA                       MATERD(10,2)         MATERF(10,2)
!     R02                         MATERD(11,2)         MATERF(11,2)
!     ZETAF                       MATERD(12,2)         MATERF(12,2)
!     PENTE EN PE                 MATERD(13,2)         MATERF(13,2)
!     PK                          MATERD(14,2)         MATERF(14,2)
!     PE                          MATERD(15,2)         MATERF(15,2)
!     CONTRAINTE EN PE            MATERD(16,2)         MATERF(16,2)
!     ZETAG                       MATERD(17,2)         MATERF(17,2)
!
!     IRRADIATION                 MATERD(18,2)         MATERF(18,2)
!     AGINT                       MATERD(19,2)         MATERF(19,2)
!
!     TOLER SUR SEUIL             MATERD(20,2)         MATERF(20,2)
!     ERREUR SUR SEUIL            MATERD(21,2)         MATERF(21,2)
!
!     TEMPERATURE                 MATERD(22,2)         MATERF(22,2)
!
!     INCREMENT IRRADIATION       MATERD(23,2)         MATERF(23,2)
!     INCREMENT TEMPERATURE       MATERD(24,2)         MATERF(24,2)
!
!
! -   PROTECTION SUR LA DIMENSION DES TABLEAUX : MATERD MATERF
!    ASSERT(nmat.ge.30)
!
! -   NB DE COMPOSANTES / VARIABLES INTERNES -------------------------
!    call irrnvi(model, ndt, ndi, nr, nvi)
!
! === ================================================
!
!     RECUPERATION MATERIAU A TEMPD ET IRRAD
!
! === ================================================
!     CARACTERISTIQUES ELASTIQUES A TEMP- ET IRRA-
!    call rcvalb(fami, kpg, ksp, '-', imat,&
!                ' ', 'ELAS', 0, ' ', [0.0d0],&
!                3, nomcel, materd(1, 1), cerr, 1)
!
!     TEMPERATURE A T-
!    call rcvarc('F', 'TEMP', '-', fami, kpg,&
!                ksp, tempd, iret)
!     IRRADIATION A T-
!    call rcvarc('F', 'IRRA', '-', fami, kpg,&
!                ksp, irrad, iret)
!     CARACTERISTIQUES MATERIAU A TEMP- ET IRRA-
!    call rcvalb(fami, kpg, ksp, '-', imat,&
!                ' ', 'IRRAD3M', 0, ' ', [0.0d0],&
!                nbcara, nomcir, mat, cerr, 1)
!
    
    mat(1) = r02
    mat(2) = eu 
    mat(3) = rm 
    mat(4) = ai0
    mat(5) = etais
    mat(6) = rg0
    mat(7) = alpha
    mat(8) = phi0 
    mat(9) = kappa
    itmax = 500  
    rela = 0.001
    pe = 0.002
    
    mat(10) =   1!mat(10)
    mat(11) =   1!mat(11)
!     POUR PLUS DE CLARETE, JE RENOMME LES GRANDEURS
    if (cerr(10) .eq. 0) then
        zetaf = mat(10)
    else
        zetaf = 1.0d0
    endif
    if (cerr(11) .eq. 0) then
        zetag = mat(11)
    else
        zetag = 1.0d0
    endif

    !r02 = mat(1)
    !eu = mat(2)
    !rm = mat(3)
    !ai0 = mat(4)
    !etais = mat(5)
    !rg0 = mat(6)
    !alpha = mat(7)
    !phi0 = mat(8)
    !kappa = mat(9)
    !irrad = 3
!
!     CALCUL DE LA PUISSANCE PAR DICHOTOMIE
!       - LA FONCTION EST MONOTONE DECROISSANTE
!       - NORMALISATION PAR R02
    coeffa = rm*exp(eu)/r02
!     F(n) = 1.0 - RM*EXP(EU)*((PE+n-EU)**n)/((n**n)*R02)
!     Finf = Limite F(n)       Fzero = Limite F(n)
!            n->infini                 n->0+
    n0 = eu - pe
    f1 = 1.0d0 - coeffa*exp(-n0)
!     L'équation peut ne pas avoir de solution, pour le vérifier on
!     calcule sa valeur FE à la 1ère borne de recherche +PE/10000.0
!     C'est avec FE que l'on vérifie que l'on à des solutions
    if (n0 .ge. 0.0d0) then
        n1 = n0 + pe/1000.0d0
    else
        n1 = pe/1000.0d0
    endif
    fe = 1.0d0 - coeffa*((n1-n0)**n1)/(n1**n1)
    if ((fe*f1.gt.0.0d0) .or. (n0.eq.0.0d0)) then
!        VALEURS PAR DEFAUT
        n1 = eu
!        VALEUR DE K , N
        if (n1 .gt. 0.0d0) then
            materd(7,2) = rm*exp(eu)/(n1**n1)
            materd(8,2) = n1
        else
            materd(7,2) = rm
            materd(8,2) = 0.0d0
        endif
!        VALEUR DE P0
        materd(9,2) = 0.0d0
!        -----------------
        k = materd(7,2)
        spe = k*(pe**n1)
        a = n1*k*(pe**(n1-1.d0))
    else
        if (n0 .gt. 0.0d0) then
            f0 = 1.0d0
            pasn = n0/10.0d0
            n1 = n0 - (pasn*0.9999d0)
        else
            f0 = 1.0d0 - coeffa
            pasn = pe/10.0d0
            n1 = - (pasn*0.9999d0)
        endif
        iterat = 0
!        WHILE TRUE
10      continue
        n1 = n1 + pasn
        f1 = 1.0d0 - coeffa*((n1-n0)**n1)/(n1**n1)
        if (abs(f1) .le. rela) goto 12
        iterat=iterat+1
        if (iterat .gt. itmax) then
            valkm(1) = 'PREMIERE'
            valim(1) = itmax
            valrm(1) = f0
            valrm(2) = f1
            valrm(3) = n1
            valrm(4) = pasn
            valrm(5) = rm
            valrm(6) = eu
            valrm(7) = r02
            valrm(8) = rela
!              VALEURS INITIALES
            valrm(9) = eu - pe
            valrm(10) = 1.0d0 - coeffa*exp(pe-eu)
            valrm(11) = 1.0d0 - coeffa
            valrm(12) = fe
!            call utmess('F', 'COMPOR1_55', sk=valkm(1), si=valim(1), nr=12,&
!                        valr=valrm)
        endif
        if (f1*f0 .gt. 0.0d0) then
            f0 = f1
        else
            n1 = n1 - pasn
            pasn = pasn * 0.5d0
        endif
        goto 10
12      continue
!        VALEUR DE K
        materd(7,2) = rm*exp(eu)/(n1**n1)
!        VALEUR DE N
        materd(8,2) = n1
!        VALEUR DE P0
        materd(9,2) = n1 - eu
!        ---------------------
        k = materd(7,2)
        p0 = materd(9,2)
        spe = k*((pe+p0)**n1)
        a = n1*k*((pe+p0)**(n1-1.d0))
    endif
    if (a .gt. 0.0d0) then
!        VALEUR DE LA PENTE EN PE
        materd(13,2) = a
!        VALEUR DE PK
        materd(14,2) = pe - (spe - kappa*r02)/a
    else
!        VALEUR DE LA PENTE EN PE
        materd(13,2) = 0.0d0
!        VALEUR DE PK
        materd(14,2) = 0.0d0
    endif
!     VALEUR DE AI0
    materd(4,2) = ai0
!     VALEUR DE ETAI_S
    materd(5,2) = etais
!     VALEUR DE AG
    exph = exp(alpha*(phi0-irrad))
    materd(6,2) = rg0/(1.0d0+exph)/3.0d0
!     VALEUR DE KAPPA
    materd(10,2) = kappa
!     VALEUR DE R02
    materd(11,2) = r02
!     VALEUR DE ZETAF
    materd(12,2) = zetaf
!     VALEUR DE PE
    materd(15,2) = pe
!     VALEUR DE LA CONTRAINTE EN PE
    materd(16,2) = spe
!     VALEUR DE ZETAG
    materd(17,2) = zetag
!     IRRADIATION
    materd(18,2) = irrad
!     VALEUR DE AG DEJA INTEGRE
    if (alpha .gt. 0.0d0) then
        exp0 = exp(alpha*phi0)
        exph = exp(alpha*irrad)
        materd(19,2) = rg0*log((exp0+exph)/(1.0d0+exp0))/(3.0d0*alpha)
    else
        materd(19,2) = 0.0d0
    endif
!     TOLERENCE ET ERREUR SUR LE FRANCHISSEMENT DU SEUIL
    materd(20,2) = mat(12)
    materd(21,2) = 0.0d0
!     TEMPERATURE
    materd(22,2) = tempd
!
! === ================================================
!
!     RECUPERATION MATERIAU A TEMPF ET IRRAF
!
! === ================================================
!     CARACTERISTIQUES ELASTIQUES A TEMP+ ET IRRA+
!    call rcvalb(fami, kpg, ksp, '+', imat,&
!                ' ', 'ELAS', 0, ' ', [0.0d0],&
!                3, nomcel, materf(1, 1), cerr, 1)
!
!     TEMPERATURE A T+
!    call rcvarc('F', 'TEMP', '+', fami, kpg,&
!                ksp, tempf, iret)
!     IRRADIATION A T+
!    call rcvarc('F', 'IRRA', '+', fami, kpg,&
!                ksp, irraf, iret)
!     L'IRRADIATION NE PEUT PAS DECROITRE
    if (irrad .gt. irraf*1.00001d0) then
        valrm(1) = tempd
        valrm(2) = tempf
 !       call utmess('I', 'COMPOR1_57', nr=2, valr=valrm)
        valrm(1) = irrad
        valrm(2) = irraf
 !       call utmess('I', 'COMPOR1_56', nr=2, valr=valrm)
    endif
    if (irrad .gt. irraf) then
        irraf = irrad
    endif
!     CARACTERISTIQUES MATERIAU A TEMP+ ET IRRA+
!    call rcvalb(fami, kpg, ksp, '+', imat,&
!                ' ', 'IRRAD3M', 0, ' ', [0.0d0],&
!                nbcara, nomcir, mat, cerr, 1)
!
!     POUR PLUS DE CLARETE
    if (cerr(10) .eq. 0) then
        zetaf = mat(10)
    else
        zetaf = 1.0d0
    endif
    if (cerr(11) .eq. 0) then
        zetag = mat(11)
    else
        zetag = 1.0d0
    endif
    !r02 = mat(1)
    !eu = mat(2)
    !rm = mat(3)
    !ai0 = mat(4)
    !etais = mat(5)
    !rg0 = mat(6)
    !alpha = mat(7)
    !phi0 = mat(8)
    !kappa = mat(9)
!
!     CALCUL DE LA PUISSANCE PAR DICHOTOMIE
!       - LA FONCTION EST MONOTONE DECROISSANTE
!       - NORMALISATION PAR R02
    coeffa = rm*exp(eu)/r02
!     F(n) = 1.0 - RM*EXP(EU)*((PE+n-EU)**n)/((n**n)*R02)
!     Finf = Limite F(n)       Fzero = Limite F(n)
!            n->infini                 n->0+
    n0 = eu - pe
    f1 = 1.0d0 - coeffa*exp(-n0)
!     L'équation peut ne pas avoir de solution, pour le vérifier on
!     calcule sa valeur FE à la 1ère borne de recherche +PE/1000.0
!     C'est avec FE que l'on vérifie que l'on à des solutions
    if (n0 .ge. 0.0d0) then
        n1 = n0 + pe/1000.0d0
    else
        n1 = pe/1000.0d0
    endif
    fe = 1.0d0 - coeffa*((n1-n0)**n1)/(n1**n1)
    if ((fe*f1.ge.0.0d0) .or. (n0.eq.0.0d0)) then
!        VALEURS PAR DEFAUT
        n1 = eu
!        VALEUR DE K , N
        if (n1 .gt. 0.0d0) then
            materf(7,2) = rm*exp(eu)/(n1**n1)
            materf(8,2) = n1
        else
            materf(7,2) = rm
            materf(8,2) = 0.0d0
        endif
!        VALEUR DE P0
        materf(9,2) = 0.0d0
!        -----------------
        k = materf(7,2)
        spe = k*(pe**n1)
        a = n1*k*(pe**(n1-1.d0))
    else
        if (n0 .gt. 0.0d0) then
            f0 = 1.0d0
            pasn = n0/10.0d0
            n1 = n0 - (pasn*0.9999d0)
        else
            f0 = 1.0d0 - coeffa
            pasn = pe/10.0d0
            n1 = - (pasn*0.9999d0)
        endif
        iterat = 0
!        WHILE TRUE
20      continue
        n1 = n1 + pasn
        f1 = 1.0d0 - coeffa*((n1-n0)**n1)/(n1**n1)
        if (abs(f1) .le. rela) goto 22
        iterat=iterat+1
        if (iterat .gt. itmax) then
            valkm(1) = 'DEUXIEME'
            valim(1) = itmax
            valrm(1) = f0
            valrm(2) = f1
            valrm(3) = n1
            valrm(4) = pasn
            valrm(5) = rm
            valrm(6) = eu
            valrm(7) = r02
            valrm(8) = rela
!              VALEURS INITIALES
            valrm(9) = eu - pe
            valrm(10) = 1.0d0 - coeffa*exp(pe-eu)
            valrm(11) = 1.0d0 - coeffa
            valrm(12) = fe
!            call utmess('F', 'COMPOR1_55', sk=valkm(1), si=valim(1), nr=12,&
!                        valr=valrm)
        endif
        if (f1*f0 .gt. 0.0d0) then
            f0 = f1
        else
            n1 = n1 - pasn
            pasn = pasn * 0.5d0
        endif
        goto 20
22      continue
!        VALEUR DE K
        materf(7,2) = rm*exp(eu)/(n1**n1)
!        VALEUR DE N
        materf(8,2) = n1
!        VALEUR DE P0
        materf(9,2) = n1 - eu
!        ---------------------
        k = materf(7,2)
        p0 = materf(9,2)
        spe = k*((pe+p0)**n1)
        a = n1*k*((pe+p0)**(n1-1.d0))
    endif
    if (a .gt. 0.0d0) then
!        VALEUR DE LA PENTE EN PE
        materf(13,2) = a
!        VALEUR DE PK
        materf(14,2) = pe - (spe - kappa*r02)/a
    else
!        VALEUR DE LA PENTE EN PE
        materf(13,2) = 0.0d0
!        VALEUR DE PK
        materf(14,2) = 0.0d0
    endif
!     VALEUR DE AI0
    materf(4,2) = ai0
!     VALEUR DE ETAI_S
    materf(5,2) = etais
!     VALEUR DE AG
    exph = exp(alpha*(phi0-irraf))
    materf(6,2) = rg0/(1.0d0+exph)/3.0d0
!     VALEUR DE KAPPA
    materf(10,2) = kappa
!     VALEUR DE R02
    materf(11,2) = r02
!     VALEUR DE ZETAF
    materf(12,2) = zetaf
!     VALEUR DE PE
    materf(15,2) = pe
!     VALEUR DE LA CONTRAINTE EN PE
    materf(16,2) = spe
!     VALEUR DE ZETAG
    materf(17,2) = zetag
!     IRRADIATION
    materf(18,2) = irraf
!     VALEUR DE AG DEJA INTEGRE
    if (alpha .gt. 0.0d0) then
        exp0 = exp(alpha*phi0)
        exph = exp(alpha*irraf)
        materf(19,2) = rg0*log((exph+exp0)/(1.0d0+exp0))/(3.0d0*alpha)
    else
        materf(19,2) = 0.0d0
    endif
!     TOLERENCE ET ERREUR SUR LE FRANCHISSEMENT DU SEUIL
    materf(20,2) = mat(12)
    materf(21,2) = 0.0d0
!     TEMPERATURE
    materf(22,2) = tempf
!
!     INCREMENT IRRADIATION
    materd(23,2) = materf(18,2) - materd(18,2)
    materf(23,2) = materd(23,2)
!     INCREMENT TEMPERATURE
    materd(24,2) = materf(22,2) - materd(22,2)
    materf(24,2) = materd(24,2)
!
! -   MATERIAU CONSTANT ?
! -   ON NE PEUT PAS SAVOIR A L AVANCE DONC NON
    matcst = 'NON'
end subroutine
     

      
      END SUBROUTINE UMAT 
