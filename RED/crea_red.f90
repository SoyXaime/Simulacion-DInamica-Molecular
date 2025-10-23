program crea_red
    use variables_comunes
    use Def_prec
    use mod_funciones
    use mod_subrutinas

    implicit none

      real (kind=doblep) :: rx(npmax), ry(npmax), rz(npmax)
      real (kind=doblep) :: vx(npmax), vy(npmax), vz(npmax)
      real (kind=doblep) :: px, py, pz, p2
      real (kind=doblep) :: ax(npmax), ay(npmax), az(npmax)
      real (kind=doblep) :: Ep,Ec,E_total
      real (kind=doblep), dimension(3,4) :: base 

      !Variables locales
      real (kind=doblep) :: xnp, pa, pma, factor, E_total_comprobacion
      integer (kind=entero) :: i, ix, iy, iz, n, ib, iflag, unit

    !-----------x-------------x----------------x-------!
      !write(*,*) '¿lado de la caja?'
      !read(*,*) pl
      !write(*,*) 'radio de corte del potencial= (maximo lado/2)'
      !read(*,*) rc

      !write(*,*) 'Energía del sistema'
      !read(*,*) E_total 

      !write(*,*) 'fichero para grabar la simulacion'
      !write(*,*) 'número aleatorio entero positivo para iniciar la simulación'
      !read(*,*) iflag
    !-----------x-------------x----------------x-------!

      pl=10.d00 !Lado de la caja
      rc=5.d00 !Radio de corte del potencial
      E_total=-665.d00 !Energía total de mi sistema
      write(*,*) 'El programa se ha iniciado con los siguientes valores por defecto'
      write(*,*) 'Lado de la caja=', pl
      write(*,*) 'Radio de corte del potencial=', rc
      write(*,*) 'Energía total del sistema=', E_total

      write(*,*) 'Elija un número aleatorio para iniciar la simulación'
      read(*,*) iflag !Semilla para la función de números aleatorios.

      xnp=dble(npmax) !xnp= numero maximo de particulas en doble precision
      pli=1.d00/pl ! Inverso del lado
      vol=pl*pl*pl !volumen de la caja
      dens=xnp/vol !densidad=Numero particulas/Volumen
      rc2=rc*rc !radio de corte del potencial al cuadrado
      pa=pl/dble(numk) !Divido el lado de la caja en k elementos
      pma=pa/2.d00 !la mitad del lado de cada elemento
      factor=pi*xnp*xnp/(vol*rc**3)
      corr_ener=8.d00*factor*(1.d00/(3.d00*rc**6)-1.d00)/3.d00 !correccion de la energia de potencial debido a cortar el radio de accion
      corr_sum_rvp=16.d00*factor*(-2.d00/(3.d00*rc**6)+1.d00)
      corr_sum_r2vp=16.d00*factor*(26.d00/(3.d00*rc**6)-7.d00)

      !print*, 'corr energ=',corr_ener !Comprobación de errores

      !Inicializo variables para evitar errores
      rx = 0.d00; ry = 0.d00; rz = 0.d00
      vx = 0.d00; vy = 0.d00; vz = 0.d00
      ax = 0.d00; ay = 0.d00; az = 0.d00
      Ep=0.d00;   Ec=0.d00



      ! Base FCC (3x4) -> (fila x columna) shape(base)=(3,4) pues es la dimension inicializada al declarar la variable
      base = reshape( [ 0.d00, 0.d00, 0.d00, & 
                        0.5d00, 0.5d00, 0.d00, &
                        0.5d00, 0.d00, 0.5d00, &
                        0.d00, 0.5d00, 0.5d00 ], shape(base) ) !Ojo fortran llena por columnas a diferencia de python que lo hace por filas, (0,0,0) es la primera columna

      ! Coloco N = 4*numk**3 particulas (FCC) recorriendo las celdas
      n = 0
      do iz = 0, numk-1
        do iy = 0, numk-1
          do ix = 0, numk-1
            do ib = 1, 4
              n = n + 1
              rx(n) = (dble(ix) + base(1,ib)) * pa
              ry(n) = (dble(iy) + base(2,ib)) * pa
              rz(n) = (dble(iz) + base(3,ib)) * pa
            end do
          end do
        end do
      end do

      call potencial(rx, ry, rz, ax, ay, az, pl, Ep, rc)! Calculo energía potencial antes de descolocar las particulas
      print*,'La energía potencial de la red fcc con las particulas perfectamente bien colocadas Ep=',Ep+corr_ener

      do i=1, npmax
          rx(i)=rx(i)+(2.d00*mi_random(iflag)-1.d00)*(0.2d00*pma)
          ry(i)=ry(i)+(2.d00*mi_random(iflag)-1.d00)*(0.2d00*pma)
          rz(i)=rz(i)+(2.d00*mi_random(iflag)-1.d00)*(0.2d00*pma)
      end do


      !Comprobación de errores ---------x-------x------    
      !PRINT*, 'Particulas colocadas y descolocadas aleatoriamente en como maximo un 20% del parametro de red'
      !print*, n

      !do i = 1, n
      !   print*, rx(i), ry(i), rz(i)
      !end do
      !Comprobación de errores ---------x-------x------  


      !Asigno velocidades aleatorias
          do i=1,npmax
              vx(i)=2.d00*mi_random(iflag)-1.d00
              vy(i)=2.d00*mi_random(iflag)-1.d00
              vz(i)=2.d00*mi_random(iflag)-1.d00
          end do


      !PRINT*, 'Velocidades aleatorias asignadas entre -1 y 1'
      !print*, 'Corrigiendo posible momento no nulo'

      
      !Calculo momento
      px=sum(vx)
      py=sum(vy)
      pz=sum(vz)


      !Corrijo momento
      px=px/xnp
      py=py/xnp
      pz=pz/xnp

      p2=px*px+py*py+pz*pz
      !print*, 'p2=', p2


      vx=vx-px
      vy=vy-py
      vz=vz-pz


      !Lo vuelvo a calcular y verifico
      px=sum(vx)
      py=sum(vy)
      pz=sum(vz)

      p2=px*px+py*py+pz*pz

      !print*, 'p2=', p2 !Comprobación de errores
  

      !write(*,*) 'Calculando energía potencial...' !Comprobación de errores
      call potencial(rx, ry, rz, ax, ay, az, pl, Ep, rc)
      Ep=Ep+corr_ener !Hay que sumar la correccion por truncamiento. Ya que no estamos teniendo en cuenta el potencial de las particulas r>rc
      !print *, "Energía potencial =", Ep !Comprobación de errores

      !write(*,*) 'Ajustando velocidades para que mi sistema tenga E=', E_total !Comprobación de errores
      call ajustar_cinetica(vx,vy,vz,E_total,Ep)

      PRINT*, 'La energía de mi sistema es E=Ep+Ec= ' !Checkeo que Ep+Ec me da E_total
      call cinetica(vx,vy,vz,Ec) !Esto tiene más sentido hacerlo como funcion en lugar de subrutina
      E_total_comprobacion=Ec+Ep

      print*, Ep, '+', Ec, '=', E_total_comprobacion
      !print *, "Energía cinetica =", Ec



      !Guardado de datos relevantes en fichero ascii

      ! Abrir archivo para escritura (sobrescribe si ya existe)
      unit=10 !Canal de I/O para comunicarse con el archivo
      open(unit=unit, file="../resultados/crea_red/resultados.txt", status="replace", action="write")

      ! Escribir cabecera
      write(10,*) "Resultados de la simulación"
      write(10,*) "Numero de partículas:", n
      write(10,*) 'Volumen de la caja', pl**3
      write(10,*) 'Lado de la caja', pl
      write(10,*) 'Radio de corte del potencial', rc 
      write(10,*) "Energía potencial =", Ep
      write(10,*) "Energía cinética  =", Ec
      write(10,*) "Energía total     =", E_total_comprobacion
      write(10,*) "Momento neto  =", sqrt(p2)

      ! Imprimir arrays, con bucle:
      write(10,*) "Posiciones, velocidades, aceleraciones"
      write(10,'(A)') '     i       rx             ry             rz             vx             vy             vz             ax             ay             az'
      do i = 1, n
        write(10,'(I6,9F15.6)') i, rx(i), ry(i), rz(i), vx(i), vy(i), vz(i), ax(i), ay(i), az(i)   !  posiciones | velocidades | aceleraciones
      end do

      ! Cerrar archivos
      close(10)        
      print*, 'Datos guardados en fichero: resultados.txt'


      open(unit=12, file="../resultados/crea_red/vectores.bin", status="replace", action="write", form="unformatted", access="stream") !form="unformatted" quiere decir que guarde los archivos en binario directamente, formated seria en texto ascii por ejemplo. access="stream" hace que los guarde todos seguidos sin meter nada de info extra
      
      ! Escribir vectores en binario
      write(12) rx, ry, rz, vx, vy, vz, ax, ay, az

      close(12)
      print*, 'Vectores r,v,a guardados en fichero: vectores.bin'


      open(unit=15, file="../resultados/crea_red/datta.txt", status="replace", action="write")
        write(15,*) npmax,pl,pli,rc,rc2
        write(15,* ) vol,dens,ktotal,kpaso,dt
        write(15,*) 
        write(15,*) 'vectores.bin'
        write(15,*) 'gname2.txt'
      close(15)
      stop
end program crea_red
