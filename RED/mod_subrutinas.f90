module mod_subrutinas
    implicit none
    
contains
    subroutine potencial(rx, ry, rz, ax, ay, az, pl, Ep,rc)
        use Def_prec
        ! use variables_comunes
        implicit none

        ! --- Argumentos
        real(kind=doblep), intent(in)            :: rx(:), ry(:), rz(:) !Los dos puntos reservan memoria dinamica del ordenador 
        real(kind=doblep), intent(out)           :: ax(:), ay(:), az(:) !Así la memoria utilizada será la de mis vectores r,v,a definida en el programa principal
        real(kind=doblep), intent(in)            :: pl
        real(kind=doblep), intent(out) :: Ep
        real(kind=doblep), intent(in) :: rc


        ! --- Locales
        integer(kind=entero) :: i, j, n
        real(kind=doblep)    :: dx, dy, dz, r2, invr2, sr2, sr6, sr12
        real(kind=doblep)    :: v_ij, Ep_local
        real(kind=doblep)    :: rc2
        real(kind=doblep)    :: fmod, fx, fy, fz

        real(kind=doblep), parameter :: sigma=1.0_doblep, epsilon=1.0_doblep

        n    = size(rx) !Mide la longitud del array rx, es decir, el numero de particulas
        v_ij = 0.d00
        Ep_local=0.d00
        rc2 = rc*rc       ! precomputo de rc^2 para comparar con r2


        

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

            ax(i) = ax(i) + fx ; ay(i) = ay(i) + fy ; az(i) = az(i) + fz
            ax(j) = ax(j) - fx ; ay(j) = ay(j) - fy ; az(j) = az(j) - fz
            
            end do
        end do

        Ep = Ep_local
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


end module mod_subrutinas