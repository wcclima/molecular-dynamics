c     Programa que dá a distribuicão de velocidades de um 
c     gás bidimensional

      implicit real*8 (a-h,o-z)

      parameter (Lx = 26) !parâmetro que refere-se ao lado do quadrado
      parameter (Ly = 10)
      parameter (N = 200) !parâmetro que refere-se ao número de partículas
      parameter (dt = 1.d-3) !parâmetro que refere-se ao intervalo de tempo
      parameter (itempo = 40000) !parametro que refere-se ao número de iteracões

      real*8 xn(N), xn1(N), yn(N), yn1(N) !vetores que armazenam as posicoes xy
      real*8 vxn(N), vyn(N),vn(N) !vetores que armazenam as velocidades xy e o módulo
      real*8 Fxn1(N), Fyn1(N) !vetores que armazenam as forcas

******************************************************

      rN = N
      tempo = itempo
      vxq = 0.d0
      vyq = 0.d0

********************************************************

      do i = 1,100
         read(1,*) xn(i),yn(i)
         read(2,*) xn1(i),yn1(i)
         read(3,*) vxn(i),vyn(i)
      end do

      do i = 101,200
         read(7,*) xn(i),yn(i)
         read(8,*) xn1(i),yn1(i)
         read(9,*) vxn(i),vyn(i)
         xn(i) = xn(i) + 13.d0
         xn1(i) = xn1(i) + 13.d0
      end do


*************************************************

c     looping que realiza o movimento das partículas 
      do i = 1,itempo

c        zera-se os vetores para a forca
         do j = 1,N
            
            Fxn1(j) = 0.d0
            Fyn1(j) = 0.d0

         end do

c        looping que calcula a forca de interacão da 
c        partícula j com as partículas k que estejam, 
c        no máximo, a uma distância de 3 da partícula j
         do j = 1,N

c           ***************************
        
            do k = (j + 1),N

c             calculo da menor distância entre as partículas
              dx = xn1(j) - xn1(k)
              dy = yn1(j) - yn1(k)
              sigdx = dx/dabs(dx)
              sigdy = dy/dabs(dy)

              if (dabs(dx).gt.(Lx/2.d0)) then
                 dx = dx - sigdx*Lx
              end if

              if (dabs(dy).gt.(Ly/2.d0)) then
                 dy = dy - sigdy*Ly
              end if

              r  = dsqrt(dx*dx + dy*dy)
              
              if (r.le.3.d0) then
                 f = 24.d0*(2.d0/(r**13) - 1.d0/(r**7)) 
              else
                 f = 0.d0
              end if

c             leva-se em conta, para a economia do programa
c             a terceira lei de Newton
              Fxn1(j) = f*dx/r + Fxn1(j)
              Fxn1(k) = -f*dx/r + Fxn1(k)

              Fyn1(j) = f*dy/r + Fyn1(j)
              Fyn1(k) = -f*dy/r + Fyn1(k)
              
            end do

c           **************************

c           calculo da posicão x da partícula no instante t + 1
            xn2 = 2.d0*xn1(j) - xn(j) + Fxn1(j)*dt*dt

c           calculo da velocidade no instante t + 1
            vxn(j) = (xn2 - xn(j))/(2.d0*dt)

c           atualiza-se as posicões x
            xn(j) = xn1(j)
            xn1(j) = xn2

c           aplica-se as condicões de contorno
            if (xn1(j).gt.Lx) then
               xn1(j) = xn1(j) - Lx
               xn(j) = xn(j) - Lx
            end if

            if (xn1(j).lt.0.d0) then
               xn1(j) = xn1(j) + Lx
               xn(j) = xn(j) + Lx
            end if

c           calcula-se a posicão y da partícula no instante t + 1
            yn2 = 2.d0*yn1(j) - yn(j) + Fyn1(j)*dt*dt

c           calcula-se a velocidade no instante t + 1
            vyn(j) = (yn2 - yn(j))/(2.d0*dt)

c           atualiza-se as posicões y da partícula
            yn(j) = yn1(j)
            yn1(j) = yn2

c           calcula-se o módulo da velocidade
            vn(j) = dsqrt(vxn(j)*vxn(j) + vyn(j)*vyn(j))

c           aplica-se as condicões de contorno
            if (yn1(j).gt.Ly) then
               yn1(j) = yn1(j) - Ly
               yn(j) = yn(j) - Ly
            end if

            if (yn1(j).lt.0.d0) then
               yn1(j) = yn1(j) + Ly
               yn(j) = yn(j) + Ly
            end if

c           calcula-se a velocidade quadrática média 
c           na direcão x e na direcão y
            vxq = vxn(j)*vxn(j)/rN + vxq
            vyq = vyn(j)*vyn(j)/rN + vyq

         end do
      end do


c     Nesta parte do program acumula-se a velocidade das partículas e no final
c     dividi-se a distribuicão pelo número de iteracões. Este procedimento é realizado
c     para sumir com as flutuacões na velocidade.

      vxq = vxq/tempo
      vyq = vyq/tempo

c     calcula-se a energia armazenada nos graus de liberdade
      e = (vxq + vyq)/2.d0
      
c     escreve-se o valor na tela
      write(*,*) e

      stop
      end








