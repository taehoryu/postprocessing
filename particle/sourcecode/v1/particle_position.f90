!1000rho
!updated by Aug 14th
      PROGRAM readarray
!      use fundamental_constants_module
      IMPLICIT NONE
      integer, parameter :: dp=kind(1.0d0)
      INTEGER :: row,i,index
      INTEGER :: io,nlines,setting
			INTEGER,ALLOCATABLE,DIMENSION(:):: particle_array
      !real(dp), dimension(:,:), allocatable ::x
      !real(dp), dimension(:), allocatable :: xi,t,y,yy,x,xx
      real(DP),ALLOCATABLE,DIMENSION(:)::particle_index,range_x,range_y,range_vx,range_vy,&
																				 range_rho,range_T
      real(DP) :: x,y,time,vx,vy,density,temperature,average_x,average_y,average_vx,average_vy
			real(DP) :: time_record,maximum_time
			real(DP) :: average_rho,average_T,old_time_particle,maximum_particle_number,initial_value
			CHARACTER(LEN=256)::F1
			CHARACTER(LEN=11)::filename
			CHARACTER(LEN=1)::filenumber
			CHARACTER(LEN=4)::charI
			LOGICAl :: file_exists
      F1="(2I15,9f18.10)"
      !print *, 'enter the number of row'
      !read *, i
      !print *, 'enter the number of every extracted data'
      !read *, k
      ! kth data will be extracted
      ! the number of row
      allocate(particle_index(2),range_x(2),range_y(2),range_vx(2),range_vy(2),&
																				 range_rho(2),range_T(2))
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
			maximum_particle_number=20
			maximum_time=0.0d0
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
			
			INQUIRE(file="Timestamp_average",exist=file_exists)
			if (file_exists) then
			print *, "file [Timestamp_average] exists"
			go to 120
			end if
			open(UNIT=13, FILE="Timestamp_average", ACTION="write", STATUS="replace")
			open(UNIT=12, FILE="Timestamp_parsed", ACTION="write", STATUS="replace")
			
			maximum_particle_number=0
			
			do i=1, 10
				write(charI,"(i0)") i-1
				filenumber=trim(charI)
				INQUIRE(file=filename//filenumber,exist=file_exists)
				if (file_exists)then
					print  *, 'file_exists', filename//filenumber
					open (1, FILE=filename//filenumber)

			io=0
		
				DO
					READ(1,*,IOSTAT=io) particle_index(1),particle_index(2),x,y,time,vx,vy!,density,temperature
					write(12,"(50f20.5)")particle_index(1)+particle_index(2),x,y,time,vx,vy!,density,temperature
					maximum_particle_number=max(maximum_particle_number,particle_index(1)+particle_index(2))
					maximum_time=max(maximum_time,time)
					if (io<0) exit
				end do

				end if
			end do
			close(12)
			close(1)
			
120  		allocate(particle_array(int(maximum_particle_number)))	
			particle_array=0							
			old_time_particle=0.0d0
121			io=0
			nlines=0
			setting=0
				open (1, FILE="Timestamp_parsed")	
				DO 
      	READ(1,"(50f20.5)",IOSTAT=io) particle_index(1),x,y,time,vx,vy!,density,temperature
				maximum_time=max(maximum_time,time)
				nlines=nlines+1
					if(time>old_time_particle .and. setting==0)then
					old_time_particle=time
					time_record=time
					range_x(1)=x
					range_x(2)=x
					range_y(1)=y
					range_y(2)=y
					range_vx(1)=vx
					range_vx(2)=vx
					range_vy(1)=vy
					range_vy(2)=vy
					average_x=0.0d0
					average_y=0.0d0
					average_vx=0.0d0
					average_vy=0.0d0
					
					
					setting=1
					end if

					if(time==old_time_particle .and. setting==1 .and. all(particle_array/=particle_index(1)))then
					particle_array(index+1)=int(particle_index(1))
!					print*, 'here',time,particle_index(1)
					  
						range_x(1)=max(range_x(1),x)
						range_x(2)=min(range_x(2),x)
						range_y(1)=max(range_y(1),y)
						range_y(2)=min(range_y(2),y)
						range_vx(1)=max(range_vx(1),vx)
						range_vx(2)=min(range_vx(2),vx)
						range_vy(1)=max(range_vy(1),vy)
						range_vy(2)=min(range_vy(2),vy)
					!	range_rho(1)=max(range_rho(1),density)
				!		range_rho(2)=min(range_rho(2),density)
			!			range_T(1)=max(range_T(1),temperature)
		!				range_T(2)=min(range_T(2),temperature)

						average_x=(average_x*dble(index)+x)/dble(index+1)
						average_y=(average_y*dble(index)+y)/dble(index+1)
						average_vx=(average_vx*dble(index)+vx)/dble(index+1)
						average_vy=(average_vy*dble(index)+vy)/dble(index+1)
		!				average_rho=(average_rho*dble(index)+density)/dble(index+1)
	!					average_T=(average_T*dble(index)+temperature)/dble(index+1)
	
					index=index+1	
									
					end if
					if(io<0) then
					write(13,"(I5,50ES20.7E3)")int(maximum_particle_number),time_record,average_x/1D5,&
					range_x/1D5,average_y/1D5,range_y/1D5,average_vx/1D5,range_vx/1D5,&
					average_vy/1D5,range_vy/1D5!,average_rho,range_rho,average_T,range_T
					print *, 'record',old_time_particle
					setting=0
					index=0
					particle_array=0
					close(1)
					if(time_record==maximum_time) goto 126
					go to 121
					end if
		end do
		close(13)
					
					
					


	126 write(*,"(A5)") "info"
		  write(*,"(A30,I5)") "maximum_particle_number",int(maximum_particle_number)
			write(*,"(A30,f10.3)") "maximum time",maximum_time


556      PRINT *, 'done'
      END PROGRAM readarray

