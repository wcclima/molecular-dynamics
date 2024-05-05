c     Programa que simula, via o algoritmo de Verlet, 
c     a dinâmica de um conjunto de N partículas confinadas 
c     em um quadrado de lado L e que interagem com um potencial
c     Liennard-Jonnes. As condicões de contorno são periódicas

      implicit real*8 (a-h,o-z)

      parameter (L = 10) !parâmetro que refere-se ao lado do quadrado
      parameter (N = 65) !parâmetro que refere-se ao número de partículas
      parameter (v = 1.d0) !parâmetro que refere-se ao módulo da velocidade 
c                           inicial das partículas
      parameter (dt = 5.d-3) !parâmetro que refere-se ao intervalo de tempo
      parameter (itempo = 2000) !parametro que refere-se ao número de iteracões

      real*8 x0(N), y0(N), xn(N), xn1(N), yn(N), yn1(N) !vetores que armazenam as posicoes xy
      real*8 vxn(N), vyn(N),vn(N) !vetores que armazenam as velocidades xy e o módulo
      real*8 Fxn1(N), Fyn1(N) !vetores que armazenam as forcas

************************************************
      rN = N

      ifile1 = 1 !arquivo de saída
      ifile2 = 2
      ifile3 = 3

*************************************************

      xmed = 0.d0
      ymed = 0.d0
      x2med = 0.d0
      y2med = 0.d0

      do i = 1,N
         read(ifile1,*) xn(i),yn(i)
         read(ifile2,*) xn1(i),yn1(i)
         x0(i) = xn1(i)
         y0(i) = yn1(i)
         read(ifile3,*) vxn(i),vyn(i)
      end do

      write(7,*) 0.d0,0.d0
      write(8,*) 0.d0,0.d0

c     looping que realiza o movimento das partículas 
      do i = 1,itempo

         ri = i

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

              if (dabs(dx).gt.(L/2.d0)) then
                 dx = dx - sigdx*L
              end if

              if (dabs(dy).gt.(L/2.d0)) then
                 dy = dy - sigdy*L
              end if

c             calcula-se a distancia entre as particulas
              r  = dsqrt(dx*dx + dy*dy)
              
c             calcula-se a forca de interacão
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
            if (xn1(j).gt.L) then
               xn1(j) = xn1(j) - L
               xn(j) = xn(j) - L
            end if

            if (xn1(j).lt.0.d0) then
               xn1(j) = xn1(j) + L
               xn(j) = xn(j) + L
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
            if (yn1(j).gt.L) then
               yn1(j) = yn1(j) - L
               yn(j) = yn(j) - L
            end if

            if (yn1(j).lt.0.d0) then
               yn1(j) = yn1(j) + L
               yn(j) = yn(j) + L
            end if

            deltax = xn1(j) - x0(j)
            deltay = yn1(j) - y0(j)
            sigdx = deltax/dabs(deltax)
            sigdy = deltay/dabs(deltay)
            
            if (dabs(deltax).gt.(L/2.d0)) then
               deltax = deltax - sigdx*L
            end if

            if (dabs(deltay).gt.(L/2.d0)) then
               deltay = deltay - sigdy*L
            end if

            xmed = deltax/rN + xmed
            ymed = deltay/rN + ymed

            x2med = deltax*deltax/rN + x2med
            y2med = deltay*deltay/rN + y2med


         end do

         write(7,*) i*dt,(x2med/ri - xmed*xmed/(ri*ri))
         write(8,*) i*dt,(y2med/ri - ymed*ymed/(ri*ri))


c        *************************************

      end do

      stop
      end








