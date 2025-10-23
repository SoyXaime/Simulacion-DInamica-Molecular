    program propiedades_estaticas

    use Def_prec
    use variables_comunes
    use mod_funciones
    use mod_subrutinas

    implicit none


    !varibles de entrada

    !integer(kind=entero), intent(in)    :: np
    real(kind=doblep), intent(inout)    :: rx(npmax), ry(npmax), rz(npmax)
    real(kind=doblep), intent(inout)    :: vx(npmax), vy(npmax), vz(npmax)
    real(kind=doblep), intent(inout)    :: ax(npmax), ay(npmax), az(npmax)
    !real(kind=doblep), intent(out)      :: 
    real(kind=doblep), intent(inout)    :: Ep, Ec





        P=pre
        P=P_md+P_LR ! Presion calculada por pares + Correccion "Long-Range"

        !P_md -> Bucle do grande
        P_LR= (( -16*pi*dens*dens ) / ( 3*(rc*rc*rc) ))*(1-(2/(3*rc**6)))


        dfiv_med=

        P=xnp*temp*/vol -dfiv_med



    end program propiedades_estaticas
