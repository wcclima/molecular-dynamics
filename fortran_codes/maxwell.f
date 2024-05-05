c     Programa que dá a distribuicão de velocidades de um 
c     gás bidimensional

      implicit real*8 (a-h,o-z)

      parameter (L = 10) !parâmetro que refere-se ao lado do quadrado
      parameter (N = 100) !parâmetro que refere-se ao número de partículas
      parameter (dv = 1.d-2) !parâmetro que refere-se ao intervalo de velocidades
      parameter (irange = 2000) !parâmetro que refere-se ao número de divisões do 
c                              intervalo de velocidades
      parameter (dt = 1.d-3) !parâmetro que refere-se ao intervalo de tempo
      parameter (itempo = 500000) !parametro que refere-se ao número de iteracões
      parameter (brange = 0.d0)

      real*8 xn(N), xn1(N), yn(N), yn1(N) !vetores que armazenam as posicoes xy
      real*8 vxn(N), vyn(N),vn(N) !vetores que armazenam as velocidades xy e o módulo
      real*8 Fxn1(N), Fyn1(N) !vetores que armazenam as forcas
      real*8 cont(irange)

******************************************************

      ifile = 50 !arquivo de saída
      rN = N
      tempo = itempo

      do i = 1,N
         read(7,*) xn(i),yn(i)
         read(8,*) xn1(i),yn1(i)
         read(9,*) vxn(i),vyn(i)
      end do

      do i = 1,irange
         cont(i) = 0.d0
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

              if (dabs(dx).gt.(L/2.d0)) then
                 dx = dx - sigdx*L
              end if

              if (dabs(dy).gt.(L/2.d0)) then
                 dy = dy - sigdy*L
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
         end do

c        inicia-se olhando o intervalo [0,dv]
         a = brange
         b = a + dv
            
         do j = 1,irange
         
            icont = 0
c           conta-se o número de partículas com velocidades 
c           compreendidas no intervalo [a,b]
            do k = 1,N
               if ((vn(k).ge.a).and.(vn(k).le.b)) then
                  icont = icont + 1   
               end if
            end do
c           normaliza-se
            cont(j) = cont(j) + icont/rN
c           passa-se para o próximo intervalo de velocidades
            a = b
            b = b + dv
         
         end do
      end do
      
c     looping que grava a distribuicão no arquivo 
      a = brange
      
      do i = 1,irange
c        calcula-se a velocidade "média" das partículas que 
c        tem velocidade entre v  e v + dv, que é v + dv/2
         vm = a + (i - 1)*dv + dv/2.d0
         write(ifile,*) vm,cont(i)/tempo
      end do

      stop
      end








