        implicit real*8 (a-h,o-z)
        xmin = 0
        xmax = 40
        ymin = 0
        ymax = 40
        ifin = 50
        ifiout= 82
        nupart = 200
        nuframe = 60000
        write(ifiout,2)"@focus off"
2       format(a10)
        write(ifiout,3)"@g0 on"
3       format(a6)
        write(ifiout,4)"@with g0"
4       format(a8)
        do nuf =1,nuframe
5         format(a12,f8.3)
          write(ifiout,5)"@world xmin ",xmin
          write(ifiout,5)"@world xmax ",xmax
          write(ifiout,5)"@world ymin ",ymin
          write(ifiout,5)"@world ymax ",ymax
          do i=1,nupart
1           format(f8.6,f10.6)
            read(ifin,*) xi,yi
            write(ifiout,*) xi,yi
          end do
6         format(a18)
          write(ifiout,6)"@S_ line type 0    "
          write(ifiout,6)"@S_ symbol 1       "
          write(ifiout,6)"@S_ symbol color 2 "
          write(ifiout,6)"@S_ symbol size .3 "
          write(ifiout,6)"&                  "
          write(ifiout,6)"@redraw            "
          write(ifiout,6)"@sleep 0.02        "
          write(ifiout,6)"@kill S_           "
        end do

        stop
        end
