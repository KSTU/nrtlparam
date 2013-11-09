module GlobalVar
	integer(4) SubNumber
	character(20),allocatable:: SubName(:)
	integer(4) MaxPoint
	real(8),allocatable:: Temp(:,:)
	real(8),allocatable:: Frac(:,:)
	real(8),allocatable:: LnGam(:,:)
	integer(4),allocatable:: PointNumber(:)
	integer(4),allocatable:: TempCount(:)
	real(8),allocatable:: TempAll(:,:)
end module

program NRTLparam
	use GlobalVar
	implicit none
	integer(4) i,j,k
	integer(4) TempInt
	integer(4) FileEnd
	integer(4) check
	
	print *, 'NRTL'
	MaxPoint=20000
	open(11,file='main.in')	!read input file
		read(11,*) SubNumber
		allocate(SubName(SubNumber+1))
		do i=1,SubNumber
			read(11,*) SubName(i)
		enddo
	close(11)
	print*,'Numbers of substances:',SubNumber
	do i=1,SubNumber
		print*, '    ', trim(SubName(i))
	enddo
	allocate(Temp(SubNumber+1,MaxPoint))	!allocating arrays
	allocate(Frac(SubNumber+1,MaxPoint))
	allocate(LnGam(SubNumber+1,MaxPoint))
	allocate(TempCount(SubNumber+1))
	allocate(TempAll(SubNumber+1,MaxPoint))
	allocate(PointNumber(SubNumber+1))
	
	!read experimental points
	do i=1,SubNumber
		TempInt=0
		
		open(12,file=SubName(i))
			TempInt=TempInt+1
			read(12,*,IOSTAT=FileEnd) Temp(i,TempInt),Frac(i,TempInt),LnGam(i,TempInt)
			do while (.not. IS_IOSTAT_END(FileEnd))
				TempInt=TempInt+1
				read(12,*,IOSTAT=FileEnd) Temp(i,TempInt),Frac(i,TempInt),LnGam(i,TempInt)
			enddo
		close(12)
		PointNumber(i)=TempInt
		print *, 'SubName ', SubName, ' readed ', PointNumber(i), ' points'
	enddo
	!detect Temerature counts
	do i=1,SubNumber
		TempCount(i)=1
		TempAll(i,1)=Temp(i,1)
		do j=2,PointNumber(i)
			TempInt=TempCount(i)
			check=1
			do k=1,TempInt
				if((Temp(i,j)-TempAll(i,k))<0.01) then
					check=0
				endif
				
			
			enddo
			if(check==1) then
				TempCount(i)=TempCount(i)+1
				TempAll(i,TempCount(i))=Temp(i,j)
			endif
		enddo
		print*,'Subroutine ', SubName(i),' have ', TempCount(i), ' Different teperatures'
		do j=1,TempCount(i)
			print*,'      ', i, '  ', TempAll(i,j)
		enddo
	enddo
	
	
	
	
	
	
end program
