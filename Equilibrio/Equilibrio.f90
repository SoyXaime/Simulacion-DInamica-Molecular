program equilibrio
    use Def_prec
    use variables_comunes
    use mod_subrutinas
    use mod_funciones
    
        implicit none
        
        character(len=25)  :: fname
        character(len=25)  :: gname1
        character(len=25) :: gname2

        real (kind=doblep) :: Ep=0.d00,Ec=0.d00,E_total=0.d00, tiempo=0.d00
        real (kind=doblep) :: xnp, factor
        real (kind=doblep) :: rx(npmax),ry(npmax),rz(npmax),vx(npmax),vy(npmax),vz(npmax),ax(npmax),ay(npmax),az(npmax)
        integer (kind=entero) :: kk, kcuenta, np

        



        write(*,*) 'fichero datos simulacion'
        read(*,9000) fname
9000    format (a25)
        open(10, file=trim(fname),status='old')!Leo el archivo con datos txt
            read(10,*) np,pl,pli,rc,rc2!Leo la primera linea
            read(10,*) vol,dens,ktotal,kpaso,dt!Leo la segunda linea
            read(10,9000) gname1 !Fichero donde estan mis rva
            read(10,9000) gname1 !Vuelvo a leer pq el primero llee la primera linea y la del fichero del profe esta vacia
            read(10,9000) gname2 !Fichero donde grabo 1000 pasos de la simultacion
        close(10)



        open(20, file=trim(gname1), form='unformatted', access='stream', status='old', action='read')
            read(20) rx,ry,rz,vx,vy,vz,ax,ay,az
        close(20)

        !print*, rx

        !Calculo valores necesarios
        xnp=dble(npmax)
        !pli=1.d00/pli
        vol=pl*pl*pl
        !dens=0.5d00
        !rc2=5.d00
        dt12=dt/2.d00
        dt2=dt*dt/2.d00
        factor=pi*xnp*xnp/(vol*rc**3)
        corr_ener=8.d00*factor*(1.d00/(3.d00*rc**6)-1.d00)/3.d00 !correccion de la energia de potencial debido a cortar el radio de accion
        corr_sum_rvp=16.d00*factor*(-2.d00/(3.d00*rc**6)+1.d00)
        corr_sum_r2vp=16.d00*factor*(26.d00/(3.d00*rc**6)-7.d00)

        call potencial(rx,ry,rz,ax,ay,az,pl,Ep,rc)
        call cinetica(vx,vy,vz,Ec)

        Ep=Ep+corr_ener !Comprobación de errores
        print*, 'E=Ep+Ec=',Ep,'+',Ec,'=', Ep+Ec


        !Comienza la dinamica

        !Provisional
        !write(*,*) 'Introduzca el tiempo en el termino la anterior simulacion. Si es la primera simulacion introduzca 0.d00'
        !read(*,*) tiempo

        !Restablezco la energía cinetica para que vuelva a -565.0d00

        if (.false.) then !Esto levanta la energía total de mi sistema hasta -565.d00 si es que se ha caido al iniciar la dinamica
            E_total=-665.d00
            call ajustar_cinetica(vx,vy,vz,E_total,Ep)
            call cinetica(vx,vy,vz,Ec)
        endif

        call potencial(rx,ry,rz,ax,ay,az,pl,Ep,rc)
        call cinetica(vx,vy,vz,Ec)

        Ep=Ep+corr_ener !Comprobación de errores
        print*, 'E=Ep+Ec=',Ep,'+',Ec,'=', Ep+Ec
  
        kcuenta=0
        ktotal=500000
        kpaso=100

        print*, kk,kcuenta,ktotal,kpaso,dt, tiempo !Comprobación de errores
        open(10, file=gname2, position='append', action='write', status='unknown')
        ! Escribe t=0 una vez
        !E_total = Ep + Ec
        !write(10,9001) 0.d0, Ep, Ec, E_total    

        do kk=1,ktotal
            !vamos a ver si grabo o no
            if (mod(kk,kpaso)==0) then !Ej: si kpaso=100 entonces grabo cuando kk={100,200,300...}
                kcuenta=kcuenta+1
                tiempo=dble((kcuenta)*kpaso)*dt
                write(10,9001) tiempo, Ep, Ec, E_total !Escribo en ASCII para luego hacer gráficas con python
            end if

            call verlet(rx,ry,rz,vx,vy,vz,ax,ay,az,Ep,Ec) !Llamamos a verlet para avanzar la simulación un tiempo dt
            Ep=Ep+corr_ener !Verlet llama a potencial pero recuerda que potencial no aplica corr_ener
            E_total=Ec+Ep

            if (mod(kk, ktotal/20) == 0) then !Imprime en pantalla cuando los calculos han avanzado un 5%
                print*, (kk*100)/ktotal, '%'
            end if
        end do
        close(10)

        write(*,*) 'Grabados', kcuenta, 'pasos' !Aviso de del fin de los calculos y avisa del numero de pasos grabados

9001    format (1pe13.6,2x,e13.6,2x,e13.6,2x,e13.6)

        !grabo vectores rva de ultima iteracion
        open(20,file=gname1,form='unformatted',access='stream')
            write(20) rx,ry,rz,vx,vy,vz,ax,ay,az
        close(20)






            !!!Luego la equilibracion son 500mil pasos grabando cada 100

    end program equilibrio !!!!Ver Funcion H en el libro y calcularla y ver que sucede
