      implicit real*8 (a-h,o-z)

      parameter (L = 10)
      parameter (N = 20)
      parameter (v = 1.d0)
      parameter (dt = 2.d-2)
      parameter (itempo = 100000)

      real*8 xn(N), xn1(N), yn(N), yn1(N)
      real*8 vxn(N), vyn(N)

      rN = N
      pi = dacos(-1.d0)
      dist = L/(dsqrt(rN) + 1.d0)
      linha = int(dsqrt(rN)) + 1
      icoluna = int(dsqrt(rN))

      m = 0


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

      do i = 1,N
         write(50,*) xn(i),yn(i)
      end do

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
         write(50,*) xn1(i),yn1(i)
      end do

      do i = 1,itempo

         do j = 1,N

            Fxn1 = 0.d0
            Fyn1 = 0.d0

            do k = 1,N
              if (j.ne.k) then

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

	         r  = dsqrt(dx**2 + dy**2)

                 Fxn1 = 24.d0*(2.d0/(r**14) - 1.d0/(r**8))*dx + Fxn1
                 Fyn1 = 24.d0*(2.d0/(r**14) - 1.d0/(r**8))*dy + Fyn1

              end if
            end do

            xn2 = 2.d0*xn1(j) - xn(j) + Fxn1*(dt**2)
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

            yn2 = 2.d0*yn1(j) - yn(j) + Fyn1*(dt**2)
            yn(j) = yn1(j)
            yn1(j) = yn2

            if (yn1(j).gt.L) then
               yn1(j) = yn1(j) - L
               yn(j) = yn(j) - L
            end if

            if (yn1(j).lt.0.d0) then
               yn1(j) = yn1(j) + L
               yn(j) = yn(j) + L
            end if
         end do

         do j = 1,N
               write(50,*) xn1(j),yn1(j),i
         end do
      end do


      stop
      end








