module mod_subrutinas
    implicit none
    
contains
    subroutine potencial(rx, ry, rz, ax, ay, az, pl, Ep, rc, W)
        use Def_prec
        ! use variables_comunes
        implicit none

        ! --- Argumentos
        real(kind=doblep), intent(in)            :: rx(:), ry(:), rz(:) !Los dos puntos reservan memoria dinamica del ordenador 
        real(kind=doblep), intent(out)           :: ax(:), ay(:), az(:) !Así la memoria utilizada será la de mis vectores r,v,a definida en el programa principal
        real(kind=doblep), intent(in)            :: pl
        real(kind=doblep), intent(out) :: Ep
        real(kind=doblep), intent(in) :: rc
        real(kind=doblep), intent(out) :: W


        ! --- Locales
        integer(kind=entero) :: i, j, n
        real(kind=doblep)    :: dx, dy, dz, r2, invr2, sr2, sr6, sr12
        real(kind=doblep)    :: v_ij, Ep_local
        real(kind=doblep)    :: rc2
        real(kind=doblep)    :: fmod, fx, fy, fz
        real(kind=doblep)    :: W_local, w_ij

        real(kind=doblep), parameter :: sigma=1.0_doblep, epsilon=1.0_doblep

        ax = 0.0_doblep !!Importante fijar a 0.d00 pq si no en la dinamica la aceleración crece en cada paso hasta infinito y más alla
        ay = 0.0_doblep
        az = 0.0_doblep

        n    = size(rx) !Mide la longitud del array rx, es decir, el numero de particulas
        v_ij = 0.d00
        Ep_local=0.d00
        rc2 = rc*rc       ! precomputo de rc^2 para comparar con r2
        W_local = 0.d00


    
    !$omp parallel do default(shared) &
    !$omp& private(i,j,dx,dy,dz,r2,invr2,sr2,sr6,sr12,v_ij,fmod,fx,fy,fz,w_ij) &
    !$omp& reduction(+:Ep_local, W_local) schedule(static)
    do i = 1, n-1 ! do variable_de_inicio=valor_inicial, valor_final, incremento(opcional)
            do j = i+1, n
            ! Desplazamientos
            dx = rx(i) - rx(j)
            dy = ry(i) - ry(j)
            dz = rz(i) - rz(j)

            ! Condiciones periódicas (imagen mínima)
            dx = dx - pl * anint(dx/pl) !Anint devuelve el entero más proximo
            dy = dy - pl * anint(dy/pl)
            dz = dz - pl * anint(dz/pl)

            r2 = dx*dx + dy*dy + dz*dz
            if (r2 == 0.d00) cycle !Cycle es un elemento de control del bucle do. Si se cumple la condicion, el programa salta al siguiente end do más proximo.  
            if (r2 > rc2) cycle   ! Si estamos fuera del radio de corte del potencial ignoramos esa interaccion. 

            invr2 = 1.0_doblep / r2
            sr2   = (sigma*sigma) * invr2 !sigma=1 lo dejo por si en un futuro me interesa meterlo por alguna razon
            sr6   = sr2*sr2*sr2
            sr12  = sr6*sr6
            

            v_ij = 4.0_doblep * epsilon * (sr12 - sr6)!epsilon vale 1, lo meto por si en un futuro me interesa meterlo
            Ep_local  = Ep_local + v_ij

            ! Calculo fuerzas y aceleraciones
            fmod = 24.0d0*invr2*(2.0d0*sr6*sr6 - sr6)   ! escalar que multiplica r_ij
            fx = fmod*dx ; fy = fmod*dy ; fz = fmod*dz
            w_ij    = dx*fx + dy*fy + dz*fz   ! r_ij · F_ij
            W_local = W_local + w_ij

        !$omp atomic
        ax(i) = ax(i) + fx
        !$omp atomic
        ay(i) = ay(i) + fy
        !$omp atomic
        az(i) = az(i) + fz
        !$omp atomic
        ax(j) = ax(j) - fx
        !$omp atomic
        ay(j) = ay(j) - fy
        !$omp atomic
        az(j) = az(j) - fz
            
            end do
    end do
    !$omp end parallel do

        Ep = Ep_local
        W  = W_local
        !if (present(epot_total)) epot_total = epot_local
        !print *, 'E_pot (LJ, total, PBC) = ', epot_local
    end subroutine


    subroutine ajustar_cinetica(vx,vy,vz,E_total,Ep) !Ajusta las velocidades de mi sistema para tener la energía cinetica que me de la energía total requerida
        use Def_prec
        implicit none

        real(kind=doblep), intent(in)            :: E_total, Ep
        real(kind=doblep), intent(inout) :: vx(:), vy(:), vz(:)
        

        ! Variables locales
        integer(kind=entero) :: i,n
        real(kind=doblep) :: v2, K, alpha, Ec !K es la energía cinetica necesaria para tener mi sistema con E_total
    
        !Calculo la energía cinetica generada aleatoriamente por el programa que asigna velocidades. Aprox N/2
        Ec=0.d00
        n=size(vx)
        do i=1, n 
            v2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
            Ec=Ec+v2
        end do
        Ec=Ec/2

        !Calculo la energía cinetica que deberia tener mi sistema para obtener E_total

        K=E_total-Ep

        !Relacion entre la energía cinetica de mi sistema y la cinetica que deberia tener

        alpha=sqrt(K/Ec)
        
        
        vx=vx*alpha
        vy=vy*alpha
        vz=vz*alpha
        
    end subroutine

    subroutine cinetica(vx,vy,vz,Ec) !Simplemente calcula la energía cinetica
        use Def_prec
        implicit none

        real(kind=doblep), intent(in) :: vx(:), vy(:), vz(:)
        real(kind=doblep), intent(out) :: Ec
        integer(kind=entero) :: i,n
        real(kind=doblep) :: v2

        Ec=0.d00
        n=size(vx)
        do i=1, n 
            v2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
            Ec=Ec+v2
        end do
        Ec=Ec/2

    end subroutine


    subroutine verlet(rx,ry,rz,vx,vy,vz,ax,ay,az,Ep,Ec,W)

    use Def_prec
    use variables_comunes

    implicit none


    !varibles de entrada

    !integer(kind=entero), intent(in)    :: np
    real(kind=doblep), intent(inout)    :: rx(npmax), ry(npmax), rz(npmax)
    real(kind=doblep), intent(inout)    :: vx(npmax), vy(npmax), vz(npmax)
    real(kind=doblep), intent(inout)    :: ax(npmax), ay(npmax), az(npmax)
    !real(kind=doblep), intent(out)      :: 
    real(kind=doblep), intent(inout)    :: Ep, Ec
    real(kind=doblep), intent(out)      :: W

   ! real(kind=doblep) :: dt2, dt12
    dt2  = 0.5d0 * dt * dt
    dt12 = 0.5d0 * dt

    rx= rx + vx*dt + ax*dt2
    ry= ry + vy*dt + ay*dt2
    rz= rz + vz*dt + az*dt2
    vx= vx + ax*dt12
    vy= vy + ay*dt12
    vz= vz + az*dt12


    call potencial(rx,ry,rz,ax,ay,az,pl,Ep,rc,W)


    vx= vx + ax*dt12
    vy= vy + ay*dt12
    vz= vz + az*dt12


    Ec= 0.5d00 * (sum(vx*vx) + sum(vy*vy) + sum(vz*vz))
    ! Instantaneous pressure example (units reduced):
    ! V = pl**3
    ! P = (2.0d0*Ec)/(3.0d0*V) + (1.0d0/(3.0d0*V))*W
    ! If using a cutoff rc with LJ, add long-range correction:
    ! dP_LR = (16.0d0*pi/3.0d0)*( ( (N/V)**2 )*( 2.0d0/(3.0d0*rc**9) - 1.0d0/(rc**3) ) )

    return
end subroutine verlet



end module mod_subrutinas