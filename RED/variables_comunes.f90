module variables_comunes
    use Def_prec
    implicit none

    ! Parámetros constantes
    integer (kind=entero), parameter :: npmax = 500   ! número máximo de partículas, por ejemplo
    integer (kind=entero), parameter :: numk  = 5     ! distancia entre particulas
    real    (kind=doblep), parameter :: pi = 3.141592653589d00     ! aproximación de π
    real    (kind=doblep), parameter :: densidad=0.5

    ! Variables de uso compartido
    real (kind=doblep) :: V !Volumen caja V=N/densidad.V=L^3. 
    real (kind=doblep) :: pl, pli, vol, dens, rc, rc2
    real (kind=doblep) :: dt=0.0001d00, dt12, dt2
    integer (kind=entero):: kpaso=1000, ktotal=50000!0

    ! Variables de corrección acumulativa
    real (kind=doblep) :: corr_ener   = 0.d00
    real (kind=doblep) :: corr_sum_rvp = 0.d00
    real (kind=doblep) :: corr_sum_r2vp = 0.d00
end module variables_comunes