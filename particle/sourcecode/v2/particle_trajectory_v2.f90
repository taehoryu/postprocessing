!updated by Sep 15th
      PROGRAM readarray
!      use fundamental_constants_module
      IMPLICIT NONE
      integer, parameter :: dp=kind(1.0d0)
      INTEGER :: row,i,index
      INTEGER :: io,io2,setting,pressure_index,n_lines,read_lines
			INTEGER,ALLOCATABLE,DIMENSION(:)::particle_index, process_number,particle_index_00,n_lines_P
      INTEGER,ALLOCATABLE,DIMENSION(:)::index_P
      !real(dp), dimension(:), allocatable :: xi,t,y,yy,x,xx
      real(DP),ALLOCATABLE,DIMENSION(:)::range_x,range_y,range_vx,range_vy,&
																				 range_rho,range_T
      real(DP),ALLOCATABLE,DIMENSION(:,:)::range_x_P,range_y_P,range_vx_P,range_vy_P
      real(DP),ALLOCATABLE,DIMENSION(:,:)::range_rho_P,range_T_P,range_P_P,average_particle_array
      real(DP),ALLOCATABLE,DIMENSION(:)::average_x_P,average_y_P,average_vx_P,average_vy_P
      real(DP),ALLOCATABLE,DIMENSION(:)::average_rho_P,average_T_P,average_P_P



      real(DP) :: x,y,time,vx,vy,density,temperature,average_x,average_y,average_vx,average_vy
			real(DP) :: time_00,x_00,y_00,vx_00,vy_00
			real(DP) :: time_record,maximum_time,time_difference,initial_t,pressure
			real(DP) :: average_rho,average_T,old_time_particle,maximum_particle_number,initial_value
			CHARACTER(LEN=256)::F1
			CHARACTER(LEN=11)::filename
			CHARACTER(LEN=1)::filenumber,charI3
			CHARACTER(LEN=10)::filename2
      CHARACTER(LEN=18)::filename3
			CHARACTER(LEN=2)::filenumber2
      CHARACTER(LEN=9),ALLOCATABLE,DIMENSION(:) ::file_path
			CHARACTER(LEN=4)::charI,charI2
			LOGICAl,ALLOCATABLE,DIMENSION(:) :: file_exists
			logical ::other_files,file_exists_parsed


      TYPE mytype
        INTEGER   :: ints
        REAL(DP)     :: doubles
      ENDTYPE mytype

      TYPE(mytype),DIMENSION(:,:),ALLOCATABLE :: particle_array
      call system("mkdir P001mbar")
      call system("mkdir P010mbar")
      call system("mkdir P100mbar")
      call system("mkdir P0001bar")
      call system("mkdir P0010bar")

      ALLOCATE(particle_array(1000,20))

      F1="(2I15,9f18.10)"
      !print *, 'enter the number of row'
      !read *, i
      !print *, 'enter the number of every extracted data'
      !read *, k
      ! kth data will be extracted
      ! the number of row
      allocate(particle_index(2),range_x(2),range_y(2),range_vx(2),range_vy(2),&
			particle_index_00(2),range_rho(2),range_T(2),process_number(20),file_exists(20))
      allocate(n_lines_P(20),index_P(20),file_path(5))!particle_array(100,6),
      allocate(range_x_P(2,20),range_y_P(2,20),range_vx_P(2,20),range_vy_P(2,20))
      allocate(average_x_P(20),average_y_P(20),average_vx_P(20),average_vy_P(20))
      allocate(average_rho_P(20),average_T_P(20),average_P_P(20))
      allocate(range_rho_P(2,20),range_T_P(2,20),range_P_P(2,20),average_particle_array(7,20))
      file_path(1)='P001mbar/'
      file_path(2)='P010mbar/'
      file_path(3)='P100mbar/'
      file_path(4)='P0001bar/'
      file_path(5)='P0010bar/'

      !F1="(ES24.15E2,A5,ES24.15E2,f15.12,A5,4ES24.15E2)"
      !targetgalaxy=0
			old_time_particle=0.0d0
      index=0
			average_x=0.0d0
			average_y=0.0d0
			average_vx=0.0d0
			average_vy=0.0d0
			average_rho=0.0d0
			average_T=0.0d0
			maximum_particle_number=0
			maximum_time=0.0d0
			file_exists=.false.
			other_files=.true.
			do i=1,2
        if(i==1)then
          initial_value=0.0d0
        else
          initial_value=1D60
        end if
        range_x(i)=initial_value
        range_y(i)=initial_value
        range_vx(i)=initial_value
        range_vy(i)=initial_value
        range_rho(i)=initial_value
        range_T(i)=initial_value
			enddo
			filename="Timestamp_0"
			filename2="Timestamp_"
			time_record=0.0d0
			INQUIRE(file="Timestamp_parsed.txt",exist=file_exists_parsed)
			if (file_exists_parsed) then
        print *, "file [Timestamp_parsed.txt] exists"
        open(1, FILE="Timestamp_parsed.txt")
        do
        READ(1,*,IOSTAT=io) particle_index_00(1),x_00,y_00,time_00
          if(time_record<time_00 .and. time_record/=0.0d0) exit
          time_record=time_00
        !  if(all(particle_index_00(1)/=particle_array))then
        !    particle_array(int(maximum_particle_number+1))=particle_index_00(1)
        !    maximum_particle_number=maximum_particle_number+1
        !  end if
        end do
			close(1)
			print *, 'maximum particle number',int(maximum_particle_number)
			go to 120
			end if
			open(UNIT=92, FILE="Timestamp_parsed.txt", ACTION="write", STATUS="replace")
			maximum_particle_number=0
			
			time_difference=500.0d0
			time_record=-time_difference

			io=0
			open(1, FILE="Timestamp_00")
				DO
				READ(1,*,IOSTAT=io) particle_index_00(1),particle_index_00(2),x_00,y_00,time_00,vx_00,vy_00,density,temperature
				maximum_particle_number=max(dble(particle_index_00(1)+particle_index_00(2)),dble(maximum_particle_number))

          if(time_00>time_record+time_difference )then
            write(92,"(2I3,50ES30.5E4)")particle_index_00(1),particle_index_00(2),x_00,y_00,time_00,vx_00,vy_00,density,temperature
            print *, "t_from 00 file (first)",time_00
            time_record=time_00
            other_files=.false.

          elseif(time_00==time_record)then
            write(92,"(2I3,50ES30.5E4)")particle_index_00(1),particle_index_00(2),x_00,y_00,time_00,vx_00,vy_00,density,temperature
										if(other_files .eqv. .false.)then
					
            do i=1,20
              if(i<10)then
                write(charI,"(i0)")i
                filenumber=trim(charI)
                INQUIRE(file=filename//filenumber,exist=file_exists(i))
              else
                write(charI2,"(i0)")i
                filenumber2=trim(charI2)
                INQUIRE(file=filename2//filenumber2,exist=file_exists(i))
              endif
              if (file_exists(i))then
                if(i<10)then
                  print  *, 'file_exists and open',i, filename//filenumber
                  open (i+10, FILE=filename//filenumber)
                else
                   print  *, 'file_exists and open',i, filename//filenumber2
                  open (i+10, FILE=filename2//filenumber2)
               
                end if
                  do
                    READ(i+10,*,IOSTAT=io2) particle_index(1),particle_index(2),x,y,time,vx,vy,density,temperature
						
                    if(time==time_00)then
                      write(92,"(2I3,50ES30.5E4)")particle_index(1),particle_index(2),x,y,time,vx,vy,density,temperature
                      if(i<10)then
                      print *, 'record from ', filename//filenumber, time
                      else
                      print *, 'record from ', filename2//filenumber2, time
                      end if
                    end if
								
                    if(time>time_00 .or. io2/=0) exit
                  end do
                  close(i+10)
                
              end if

            end do
            other_files=.true.

					end if
				
				
				
				endif
					
				if (io<0) exit
      end do
			
      close(92)
      close(1)

120			old_time_particle=0.0d0
			
			time=0.0d0
      initial_t=0.0d0
      n_lines=0
      n_lines_P=0
			io=0
      particle_array(:,1:3)%ints=0
      particle_array(:,4:6)%doubles=0.0d0
      average_particle_array=0.0d0
      open (1, FILE="Timestamp_parsed.txt")
     ! id=2083
      DO while (io>=0)
      	READ(1,"(2I3,50ES30.5E4)",IOSTAT=io) particle_index(1),particle_index(2),x,y,time,vx,vy,density,temperature

!        n_lines=n_lines+1
        if(initial_t==0.0d0)then
          initial_t=time
        end if
!write(*,"(2I3,50ES12.5E2)"),particle_index(1),particle_index(2),x,y,time,vx,vy,density,temperature
        if(time==initial_t)then
          pressure=density*temperature* (1.3806488D-16)/2.34d0*6.02214129D23
          n_lines=n_lines+1
          if(1D-3<pressure/1D6 .and. pressure/1D6 <=1D-2)then
            particle_array(n_lines,3)%ints=1
          else if(1D-2<pressure/1D6 .and. pressure/1D6<=1D-1)then
            particle_array(n_lines,3)%ints=2
          else if(1D-1<pressure/1D6 .and. pressure/1D6<=1.0d0)then
            particle_array(n_lines,3)%ints=3
          else if(1.0d0<pressure/1D6 .and. pressure/1D6<=1D1)then
            particle_array(n_lines,3)%ints=4
          else if(1D1<pressure/1D6 .and. pressure/1D6<=1D2)then
            particle_array(n_lines,3)%ints=5
          end if

            n_lines_P(particle_array(n_lines,3)%ints)=n_lines_P(particle_array(n_lines,3)%ints)+1
            particle_array(n_lines,1)%ints=particle_index(1)
            particle_array(n_lines,2)%ints=particle_index(2)
          !  particle_array(n_lines,3)%ints=particle_array(n_lines,3)
            particle_array(n_lines,4)%doubles=density
            particle_array(n_lines,5)%doubles=temperature
            particle_array(n_lines,6)%doubles=pressure
            particle_array(n_lines,7)%doubles=x
            particle_array(n_lines,8)%doubles=y
            particle_array(n_lines,9)%doubles=vx
            particle_array(n_lines,10)%doubles=vy
            average_particle_array(1,particle_array(n_lines,3)%ints)=average_particle_array(1,particle_array(n_lines,3)%ints)&
                                                                    +density
            average_particle_array(2,particle_array(n_lines,3)%ints)=average_particle_array(2,particle_array(n_lines,3)%ints)&
                                                                    +temperature
            average_particle_array(3,particle_array(n_lines,3)%ints)=average_particle_array(3,particle_array(n_lines,3)%ints)&
                                                                    +pressure
            average_particle_array(4,particle_array(n_lines,3)%ints)=average_particle_array(4,particle_array(n_lines,3)%ints)+x
            average_particle_array(5,particle_array(n_lines,3)%ints)=average_particle_array(5,particle_array(n_lines,3)%ints)+y
            average_particle_array(6,particle_array(n_lines,3)%ints)=average_particle_array(6,particle_array(n_lines,3)%ints)+vx
            average_particle_array(7,particle_array(n_lines,3)%ints)=average_particle_array(7,particle_array(n_lines,3)%ints)+vy


        else
        exit
        end if



      end do
      do i=1,5
        average_particle_array(:,i)=average_particle_array(:,i)/n_lines_P(i)
      end do

print *, "Total number of particle"
print *, "P~1mbar", n_lines_P(1)
print *, "P~10mbar", n_lines_P(2)
print *, "P~100mbar", n_lines_P(3)
print *, "P~1bar", n_lines_P(4)
print *, "P~10bar", n_lines_P(5)
print *, "Initial average values"

write(*,"(A25,5ES15.5E2)") "P[bar]~",(10.0**(i-4),i=1,5)
print *, "density [g/cm^3] ", (average_particle_array(1,i),i=1,5)
print *, "temperature  [K] ", (average_particle_array(2,i),i=1,5)
print *, "pressure   [bar] ", (average_particle_array(3,i)/1D6,i=1,5)
print *, "x           [km] ", (average_particle_array(4,i)/1D5,i=1,5)
print *, "y           [km] ", (average_particle_array(5,i)/1D5,i=1,5)
print *, "vx        [km/s] ", (average_particle_array(6,i)/1D5,i=1,5)
print *, "vy        [km/s] ", (average_particle_array(7,i)/1D5,i=1,5)
      close(1)

write(*,*) "Data parsed! Continue?"
read(*,*) i

      range_x_P=0.0d0
      range_y_P=0.0d0
      range_vx_P=0.0d0
      range_vy_P=0.0d0
      index_P=0

      read_lines=0
      DO i=1,n_lines
        open (1, FILE="Timestamp_parsed.txt")
        write(charI,"(i0)") i
        write(charI2,"(i0)") particle_array(i,2)%ints
        write(charI3,"(i0)") particle_array(i,3)%ints
        filename3='trajectory_'//trim(charI)//'.txt'

        do while (io>=0)
          READ(1,"(2I3,50ES30.5E4)",IOSTAT=io) particle_index(1),particle_index(2),x,y,time,vx,vy,density,temperature
          if(particle_array(i,1)%ints==particle_index(1).and. particle_array(i,2)%ints==particle_index(2))then
            print *,i,time,particle_array(i,1)%ints,particle_array(i,2)%ints,particle_array(i,3)%ints

            open(UNIT=10, FILE=file_path(particle_array(i,3)%ints)//filename3, ACTION="write",STATUS="unknown",access='append')
            write(10,"(50ES30.5E4)") time,x,y,vx,vy

          end if

        end do
        close(10)
        close(1)
        io=0
      end do



	126 write(*,"(A5)") "info"
		  write(*,"(A30,I5)") "maximum_particle_number",int(maximum_particle_number-1)
			write(*,"(A30,f10.3)") "maximum time",maximum_time


556      PRINT *, 'done'
      END PROGRAM readarray

