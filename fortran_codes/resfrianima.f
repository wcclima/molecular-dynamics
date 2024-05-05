      implicit real*8 (a-h,o-z)

      parameter (L = 40)
      parameter (N = 700)
      parameter (v = 1.d0)
      parameter (dt = 5.d-3)
      parameter (itempo1 = 500)
      parameter (itempo2 = 500)
      parameter (itempo3 = 2000)
      parameter (fator = 5.d-2)

      real*8 xn(N), xn1(N), yn(N), yn1(N)
      real*8 vxn(N), vyn(N),vn(N) 
      real*8 Fxn1(N), Fyn1(N)

5     format(a12,f8.3)
6     format(a18)

*****************************************************

      xmin = 0.d0
      xmax = L
      ymin = 0.d0
      ymax = L
      write(*,2)"@focus off"
2     format(a10)
      write(*,3)"@g0 on"
3     format(a6)
      write(*,4)"@with g0"
4     format(a8)

******************************************************

      rN = N
      pi = dacos(-1.d0)
      dist = L/(dsqrt(rN) + 1.d0)
      linha = int(dsqrt(rN)) + 1
      icoluna = int(dsqrt(rN))

      m = 0

*********************************************************

      do i = 1,linha
         do j = 1,icoluna
            m = m + 1
            xn(m) = j*dist
            yn(m) = i*dist
         end do
      end do

      do i = 1,N
         xn(i) = xn(i) + dist*(rand() - 5.d-1)/2.d0
         yn(i) = yn(i) + dist*(rand() - 5.d-1)/2.d0
      
         if (xn(i).gt.L) then
            xn(i) = xn(i) - L
         end if

	 if (xn(i).lt.0.d0) then
            xn(i) = xn(i) + L
         end if

	 if (yn(i).gt.L) then
            yn(i) = yn(i) - L
         end if

	 if (yn(i).lt.0.d0) then
            yn(i) = yn(i) + L
         end if
      end do

      write(*,5)"@world xmin ",xmin
      write(*,5)"@world xmax ",xmax
      write(*,5)"@world ymin ",ymin
      write(*,5)"@world ymax ",ymax

      do i = 1,N
         write(*,*) xn(i),yn(i)
      end do

      write(*,6)"@S_ line type 0    "
      write(*,6)"@S_ symbol 1       "
      write(*,6)"@S_ symbol color 2 "
      write(*,6)"@S_ symbol size .3 "
      write(*,6)"&                  "
      write(*,6)"@redraw            "
      write(*,6)"@sleep 0.02        "
      write(*,6)"@kill S_           "

********************************************************

      do i = 1,N
         vxn(i) = v*dcos(2.d0*pi*rand())
         vyn(i) = v*dsin(2.d0*pi*rand())
      end do

      do i = 1,N
         xn1(i) = xn(i) + vxn(i)*dt
         yn1(i) = yn(i) + vyn(i)*dt

	 if (xn1(i).gt.L) then
            xn1(i) = xn1(i) - L
            xn(i) = xn(i) - L
         end if

         if (xn1(i).lt.0.d0) then
            xn1(i) = xn1(i) + L
            xn(i) = xn(i) + L
         end if

         if (yn1(i).gt.L) then
            yn1(i) = yn1(i) - L
            yn(i) = yn(i) - L
         end if

         if (yn1(i).lt.0.d0) then
            yn1(i) = yn1(i) + L
            yn(i) = yn(i) + L
         end if
      end do

      write(*,5)"@world xmin ",xmin
      write(*,5)"@world xmax ",xmax
      write(*,5)"@world ymin ",ymin
      write(*,5)"@world ymax ",ymax

      do i = 1,N
         write(*,*) xn1(i),yn1(i)
      end do

      write(*,6)"@S_ line type 0    "
      write(*,6)"@S_ symbol 1       "
      write(*,6)"@S_ symbol color 2 "
      write(*,6)"@S_ symbol size .3 "
      write(*,6)"&                  "
      write(*,6)"@redraw            "
      write(*,6)"@sleep 0.02        "
      write(*,6)"@kill S_           "

**************************************************

      do i = 1,itempo1

         do j = 1,N
            
            Fxn1(j) = 0.d0
            Fyn1(j) = 0.d0

         end do

c     *****************************

         do j = 1,N
        
            do k = (j + 1),N

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
              Fxn1(j) = f*dx/r + Fxn1(j)
              Fxn1(k) = -f*dx/r + Fxn1(k)

              Fyn1(j) = f*dy/r + Fyn1(j)
              Fyn1(k) = -f*dy/r + Fyn1(k)
              
            end do

c     ********************************************

            xn2 = 2.d0*xn1(j) - xn(j) + Fxn1(j)*dt*dt
            
            vxn(j) = (xn2 - xn(j))/(2.d0*dt)

            xn(j) = xn1(j)
            xn1(j) = xn2

            if (xn1(j).gt.L) then
               xn1(j) = xn1(j) - L
               xn(j) = xn(j) - L
            end if

            if (xn1(j).lt.0.d0) then
               xn1(j) = xn1(j) + L
               xn(j) = xn(j) + L
            end if

            yn2 = 2.d0*yn1(j) - yn(j) + Fyn1(j)*dt*dt

            vyn(j) = (yn2 - yn(j))/(2.d0*dt)

            yn(j) = yn1(j)
            yn1(j) = yn2

            vn(j) = dsqrt(vxn(j)*vxn(j) + vyn(j)*vyn(j))

            if (yn1(j).gt.L) then
               yn1(j) = yn1(j) - L
               yn(j) = yn(j) - L
            end if

            if (yn1(j).lt.0.d0) then
               yn1(j) = yn1(j) + L
               yn(j) = yn(j) + L
            end if
         end do

         write(*,5)"@world xmin ",xmin
         write(*,5)"@world xmax ",xmax
         write(*,5)"@world ymin ",ymin
         write(*,5)"@world ymax ",ymax

         do j = 1,N
               write(*,*) xn1(j),yn1(j)
         end do

         write(*,6)"@S_ line type 0    "
         write(*,6)"@S_ symbol 1       "
         write(*,6)"@S_ symbol color 2 "
         write(*,6)"@S_ symbol size .3 "
         write(*,6)"&                  "
         write(*,6)"@redraw            "
         write(*,6)"@sleep 0.02        "
         write(*,6)"@kill S_           "

      end do

******************************************************

      do i = 1,N

         vxn(i) = vxn(i)*fator
         vyn(i) = vyn(i)*fator

         xn(i) = xn1(i) - vxn(i)*dt
         yn(i) = yn1(i) - vyn(i)*dt

      end do

******************************************************

      do i = 1,itempo2

         do j = 1,N
            
            Fxn1(j) = 0.d0
            Fyn1(j) = 0.d0

         end do

         do j = 1,N
        
            do k = (j + 1),N

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
              Fxn1(j) = f*dx/r + Fxn1(j)
              Fxn1(k) = -f*dx/r + Fxn1(k)

              Fyn1(j) = f*dy/r + Fyn1(j)
              Fyn1(k) = -f*dy/r + Fyn1(k)
              
            end do

c     **************************************************

            xn2 = 2.d0*xn1(j) - xn(j) + Fxn1(j)*dt*dt
            
            vxn(j) = (xn2 - xn(j))/(2.d0*dt)

            xn(j) = xn1(j)
            xn1(j) = xn2

            if (xn1(j).gt.L) then
               xn1(j) = xn1(j) - L
               xn(j) = xn(j) - L
            end if

            if (xn1(j).lt.0.d0) then
               xn1(j) = xn1(j) + L
               xn(j) = xn(j) + L
            end if

            yn2 = 2.d0*yn1(j) - yn(j) + Fyn1(j)*dt*dt

            vyn(j) = (yn2 - yn(j))/(2.d0*dt)

            yn(j) = yn1(j)
            yn1(j) = yn2

            vn(j) = dsqrt(vxn(j)*vxn(j) + vyn(j)*vyn(j))

            if (yn1(j).gt.L) then
               yn1(j) = yn1(j) - L
               yn(j) = yn(j) - L
            end if

            if (yn1(j).lt.0.d0) then
               yn1(j) = yn1(j) + L
               yn(j) = yn(j) + L
            end if
         end do

         write(*,5)"@world xmin ",xmin
         write(*,5)"@world xmax ",xmax
         write(*,5)"@world ymin ",ymin
         write(*,5)"@world ymax ",ymax

         do j = 1,N
               write(*,*) xn1(j),yn1(j)
         end do

         write(*,6)"@S_ line type 0    "
         write(*,6)"@S_ symbol 1       "
         write(*,6)"@S_ symbol color 3 "
         write(*,6)"@S_ symbol size .3 "
         write(*,6)"&                  "
         write(*,6)"@redraw            "
         write(*,6)"@sleep 0.02        "
         write(*,6)"@kill S_           "

      end do

***********************************************************

      do i = 1,N

         vxn(i) = vxn(i)*fator
         vyn(i) = vyn(i)*fator

         xn(i) = xn1(i) - vxn(i)*dt
         yn(i) = yn1(i) - vyn(i)*dt

      end do

***********************************************************

      do i = 1,itempo3

         do j = 1,N
            
            Fxn1(j) = 0.d0
            Fyn1(j) = 0.d0

         end do

c     **************************************

         do j = 1,N
        
            do k = (j + 1),N

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
              Fxn1(j) = f*dx/r + Fxn1(j)
              Fxn1(k) = -f*dx/r + Fxn1(k)

              Fyn1(j) = f*dy/r + Fyn1(j)
              Fyn1(k) = -f*dy/r + Fyn1(k)
              
            end do

c     ***************************************

            xn2 = 2.d0*xn1(j) - xn(j) + Fxn1(j)*dt*dt
            
            vxn(j) = (xn2 - xn(j))/(2.d0*dt)

            xn(j) = xn1(j)
            xn1(j) = xn2

            if (xn1(j).gt.L) then
               xn1(j) = xn1(j) - L
               xn(j) = xn(j) - L
            end if

            if (xn1(j).lt.0.d0) then
               xn1(j) = xn1(j) + L
               xn(j) = xn(j) + L
            end if

            yn2 = 2.d0*yn1(j) - yn(j) + Fyn1(j)*dt*dt

            vyn(j) = (yn2 - yn(j))/(2.d0*dt)

            yn(j) = yn1(j)
            yn1(j) = yn2

            vn(j) = dsqrt(vxn(j)*vxn(j) + vyn(j)*vyn(j))

            if (yn1(j).gt.L) then
               yn1(j) = yn1(j) - L
               yn(j) = yn(j) - L
            end if

            if (yn1(j).lt.0.d0) then
               yn1(j) = yn1(j) + L
               yn(j) = yn(j) + L
            end if
         end do

         write(*,5)"@world xmin ",xmin
         write(*,5)"@world xmax ",xmax
         write(*,5)"@world ymin ",ymin
         write(*,5)"@world ymax ",ymax

         do j = 1,N
               write(*,*) xn1(j),yn1(j)
         end do

         write(*,6)"@S_ line type 0    "
         write(*,6)"@S_ symbol 1       "
         write(*,6)"@S_ symbol color 4 "
         write(*,6)"@S_ symbol size .3 "
         write(*,6)"&                  "
         write(*,6)"@redraw            "
         write(*,6)"@sleep 0.02        "
         write(*,6)"@kill S_           "

      end do


      stop
      end








