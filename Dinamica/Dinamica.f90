program dinamica
    use Def_prec
    use variables_comunes
    use mod_subrutinas
    use mod_funciones
    
        implicit none
        
        character(len=25)  :: fname
        character(len=25)  :: gname1
        character(len=25) :: gname2
        character(len=50) :: codigo
        character(len=100) :: carpeta

        real (kind=doblep) :: Ep=0.d00,Ec=0.d00,E_total=0.d00,tiempo=0.d00
        real (kind=doblep) :: W,T_inst=0.d0
        real (kind=doblep) :: Pdm=0.d0, Pcorr=0.d0, dPLR=0.d0, F2_inst=0.d0
        real (kind=doblep) :: xnp, factor,px,py,pz,p2
        real (kind=doblep) :: rx(npmax),ry(npmax),rz(npmax),vx(npmax),vy(npmax),vz(npmax),ax(npmax),ay(npmax),az(npmax)
        integer (kind=entero) :: kk=1, kcuenta, np

        


        print*, "**************************************************"
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

        print*, "**************************************************"
        !Vamos a crear una carpeta especifica para guardar los datos
        ! Pedir código
        write(*,*) 'Introduce el código de la simulación, algo como: TiempoInicial-TiempoFinal'
        read(*,'(A)') codigo

        ! Construir el nombre de la carpeta
         carpeta = '../resultados/dinamica/simulacion_' // trim(adjustl(codigo))

        ! Crear la carpeta ¡¡¡¡IMPORTANTE!!!!! Compatible solo con sistemas UNIX como Linux/macOS
        call system('mkdir -p ' // trim(carpeta))


            ! Confirmar
        write(*,*) '✅ Carpeta creada:', trim(carpeta)
        print*, "**************************************************"

        open(20, file=trim(gname1), form='unformatted', access='stream', status='old', action='read')
            read(20) rx,ry,rz,vx,vy,vz,ax,ay,az
        close(20)

        !print*, rx

        !Calculo valores necesarios
        xnp=dble(npmax)
        !pli=1.d00/pli
        vol=pl*pl*pl
        dens = xnp/vol
        !dens=0.5d00
        !rc2=5.d00
        dt12=dt/2.d00
        dt2=dt*dt/2.d00
        factor=pi*xnp*xnp/(vol*rc**3)
        corr_ener=8.d00*factor*(1.d00/(3.d00*rc**6)-1.d00)/3.d00 !correccion de la energia de potencial debido a cortar el radio de accion
        corr_sum_rvp=16.d00*factor*(-2.d00/(3.d00*rc**6)+1.d00)
        corr_sum_r2vp=16.d00*factor*(26.d00/(3.d00*rc**6)-7.d00)

        !Calculo momento
        px=sum(vx)
        py=sum(vy)
        pz=sum(vz)
        p2=px*px+py*py+pz*pz
        print*, "**************************************************"
        print*, "DATOS RELEVANTES SOBRE LA SIMULACION AL INICIO"
        print*, "**************************************************"
        print*, 'p2=', p2 !Comprobación de errores

        call potencial(rx,ry,rz,ax,ay,az,pl,Ep,rc,W)
        call cinetica(vx,vy,vz,Ec)

        Ep=Ep+corr_ener !Comprobación de errores
        print*, 'E=Ep+Ec=',Ep,'+',Ec,'=', Ep+Ec
        print*,"**************************************************"
        kcuenta=0
        ktotal=500000
        kpaso=100
        print*, "kk       kcuenta     ktotal       kpaso        dt         tiempo"
        print*, kk,kcuenta,ktotal,kpaso,dt, tiempo !Comprobación de errores
        print*,"**************************************************"

        ! Long-range tail correction for truncated LJ (epsilon=sigma=1):
        dPLR   = (16.d00*pi/3.0d00) * (dens*dens) * ( 2.0d00/(3.0d00*rc**9) - 1.0d00/(rc**3) )
        !Comienza la dinamica

        !Provisional
        !write(*,*) 'Introduzca el tiempo en el termino la anterior simulacion. Si es la primera simulacion introduzca 0.d00'
        !read(*,*) tiempo



        open(10, file=trim(carpeta)//'/'//trim(gname2), position='append', action='write', status='unknown')
        write(10,'(A)') ' t   Ec   Ep   E_tot   W   T   P_corr   F2'
        open(11, file=trim(carpeta)//'/velocidad_x.txt',position='append', action='write', status='unknown')
        open(12, file=trim(carpeta)//'/velocidad_y.txt',position='append', action='write', status='unknown') 
        open(13, file=trim(carpeta)//'/velocidad_z.txt',position='append', action='write', status='unknown')  

        do kk=1,ktotal
            !vamos a ver si grabo o no
            if (mod(kk,kpaso)==0) then !Ej: si kpaso=100 entonces grabo cuando kk={100,200,300...}
                kcuenta=kcuenta+1
                tiempo=dble((kcuenta)*kpaso)*dt
                T_inst = 2.0_doblep*Ec/(3.0d0*xnp)
                ! Virial pressure with your (plus) convention; change sign if your professor’s convention:
                Pdm    = (xnp*T_inst)/vol + W/(3.0d0*vol)
                Pcorr  = Pdm + dPLR
                ! Mean squared force per particle (after forces are updated by verlet):
                F2_inst = ( sum(ax*ax) + sum(ay*ay) + sum(az*az) ) / xnp
                write(10,9001) tiempo, Ec, Ep, E_total, W, T_inst, Pcorr, F2_inst !Escribo en ASCII para luego hacer gráficas con python
                write(11,*) vx
                write(12,*) vy
                write(13,*) vz
            end if

            call verlet(rx,ry,rz,vx,vy,vz,ax,ay,az,Ep,Ec,W) !Llamamos a verlet para avanzar la simulación un tiempo dt
            Ep=Ep+corr_ener !Verlet llama a potencial pero recuerda que potencial no aplica corr_ener
            E_total=Ec+Ep

            if (mod(kk, ktotal/20) == 0) then !Imprime en pantalla cuando los calculos han avanzado un 5%
                print*, (kk*100)/ktotal, '%'
            end if
        end do
        close(10)
        close(11)
        close(12)
        close(13)

        write(*,*) 'Grabados', kcuenta, 'pasos' !Aviso de del fin de los calculos y avisa del numero de pasos grabados

9001    format (1pe13.6,8(2x,e13.6))

        !grabo vectores rva de ultima iteracion
        open(20,file=trim(gname1),form='unformatted',access='stream')
            write(20) rx,ry,rz,vx,vy,vz,ax,ay,az
        close(20)

        ! Copiar el archivo de vectores a la carpeta de la simulación
        call system('cp ' // trim(gname1) // ' ' // trim(carpeta) // '/')
        write(*,*) '📦 Archivo vectores.bin copiado a la carpeta:', trim(carpeta)





            !!!Luego la equilibracion son 500mil pasos grabando cada 100


        print*, "**************************************************"
        print*, "DATOS RELEVANTES SOBRE LA SIMULACION AL FINAL"
        print*, "**************************************************"
        px=sum(vx)
        py=sum(vy)
        pz=sum(vz)
        !Corrijo momento
        p2=px*px+py*py+pz*pz

        print*, 'p2=', p2 !Comprobación de errores

        Ep=Ep+corr_ener !Comprobación de errores
        print*, 'E=Ep+Ec=',Ep,'+',Ec,'=', Ep+Ec
        print*,"**************************************************"
        print*, "kk       kcuenta     ktotal       kpaso        dt         tiempo"
        print*, kk,kcuenta,ktotal,kpaso,dt, tiempo !Comprobación de errores
        print*,"**************************************************"


    end program dinamica !!!!Ver Funcion H en el libro y calcularla y ver que sucede


