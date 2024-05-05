      implicit real*8 (a-h,o-z)

      parameter (L = 20)
      parameter (N = 400)
      parameter (v = 1.d-1)
      parameter (dt = 5.d-3)
      parameter (itempo = 20)

      real*8 xn(N), xn1(N), yn(N), yn1(N)
      real*8 vxn(N), vyn(N),vn(N) 
      real*8 Fxn1(N), Fyn1(N)

******************************************************

      rN = N
      pi = dacos(-1.d0)
      dist = L/dsqrt(rN)
      linha = int(dsqrt(rN)) + 1
      icoluna = int(dsqrt(rN))
      ifile = 7

      m = 0

********************************************************

      do i = 1,linha
         do j = 1,icoluna
            m = m + 1
            xn(m) = j*dist
            yn(m) = i*dist
         end do
      end do

      do i = 1,N
      
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

****************************************************

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

      do i = 1,N
         write(49,*) xn(i),yn(i)
      end do

************************************************

      do i = 1,itempo

         c = i/10.d0
         ic = i/10.d0

         c = c - ic

         do j = 1,N
            
            Fxn1(j) = 0.d0
            Fyn1(j) = 0.d0

         end do

         do j = 1,N

c           ***************************
        
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

c           **************************

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

         if (c.eq.0.d0) then

            do j = 1,N
               write(ifile,*) xn1(j),yn1(j)
            end do
            ifile  = ifile + 1

         end if

c        *************************************

      end do

      stop
      end








