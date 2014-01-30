module GlobalVar
	real(8),allocatable:: DbPar(:,:)
	character(20),allocatable:: SubName(:)
	integer(4),allocatable:: PoiNum(:)
	real(8),allocatable:: Temp(:,:)
	real(8),allocatable:: Frac(:,:)
	real(8),allocatable:: LnG(:,:)
	real(8),allocatable:: Jac(:,:)
	real(8),allocatable:: TJac(:,:)
	real(8),allocatable:: MatA(:,:)
	real(8),allocatable:: MatG(:,:)
	real(8),allocatable:: MatR(:,:)
	real(8),allocatable:: MatE(:,:)
	real(8),allocatable:: MatO(:,:)
	real(8),allocatable:: MatDBet(:,:)
	real(8),allocatable:: MatDP(:,:)
	real(8),allocatable:: MatDPN(:,:)
	
	real(8),allocatable:: DevAbsGB(:,:)
	real(8),allocatable:: DevAbsGA(:,:)
	real(8),allocatable:: DevRelGB(:,:)
	real(8),allocatable:: DevRelGA(:,:)
	real(8),allocatable:: LnGB(:,:)
	real(8),allocatable:: LnGA(:,:)
	
	integer(4) MetK
	real(8) MetV
	real(8) MetT
	real(8) MetE1,MetE2
	real(8) MetMu
	integer(4) error
	real(8) Lout
	real(8),allocatable:: SumDAA(:)
	real(8),allocatable:: SumDAB(:)
	real(8) DAA,DAB
	real(8) MyP1,MyP2
	real(8),allocatable:: LTemp(:,:)
	
	real(8),allocatable:: OptMemPar(:,:)
	real(8) SumDL
	real(8) CurDL
	real(8),allocatable:: OptPoi(:,:)
	real(8) MuMin
	integer(4) StepNum
	
	integer(4) AppType
end module

program NRTLparam
	use GlobalVar
	implicit none
	integer(4) TempInt
	integer(4) ArrSub
	integer(4) ArrPoi
	integer(4) ArrPar1,ArrPar2
	integer(4) ArrMP
	integer(4) i,j,k
	integer(4) FileEnd
	integer(4) UseParA,UseParT,SumPoi,ArrPar
	real(8) MetQ
	integer(4) CheckOut
	integer(4) Init
	integer(4) DebEn
	character(120) TempStr
	integer :: seed
	
	seed=31345
	call srand(seed)
	ArrSub=2
	ArrPar1=4
	ArrPar2=7
	allocate(SubName(ArrSub))
	allocate(PoiNum(ArrSub))
	allocate(SumDAA(ArrSub))
	allocate(SumDAB(ArrSub))
	allocate(LTemp(1,1))
	
	open(11,file='main.in')	!read input
		read(11,'(a)') TempStr
		read(11,'(a)') SubName(1)
		read(11,'(a)') TempStr
		read(11,'(a)') SubName(2)
		read(11,'(a)') TempStr
		read(11,*) AppType
		read(11,'(a)') TempStr
		read(11,*) MuMin
		read(11,'(a)') TempStr
		read(11,*) StepNum
		read(11,'(a)') TempStr
		read(11,*) UseParA
		if (UseParA>2) then
			UseParA=2
		endif
		read(11,'(a)') TempStr
		read(11,*) UseParT
		if (UseParT>4) then
			UseParT=4
		endif
	close(11)
	
	do i=1,2	!Check For 
		TempInt=0
		if (trim(SubName(i))/='null') then
			open(12,file=SubName(i))
				TempInt=TempInt+1
				read(12,*,IOSTAT=FileEnd) TempInt !Temp(i,TempInt),Frac(i,TempInt),LnG(i,TempInt)
				do while (.not. IS_IOSTAT_END(FileEnd))
					TempInt=TempInt+1
					read(12,*,IOSTAT=FileEnd) TempInt !Temp(i,TempInt),Frac(i,TempInt),LnG(i,TempInt)
				enddo
			close(12)
			PoiNum(i)=TempInt-1
			!print *, 'SubName ', SubName(i), ' readed ', PoiNum(i), ' points'
		else
			PoiNum(i)=0
		endif
	enddo
	print *, 'prepearing arrays'
	ArrPoi=PoiNum(1)+PoiNum(2)
	ArrPar=UseParA+2*UseParT
	
	allocate(DbPar(ArrPar1,ArrPar2))	!allocating parameters
	allocate(MatDP(ArrPar,1))
	allocate(MatDPN(ArrPar,1))
	allocate(OptMemPar(ArrPar,1))
	
	allocate(Temp(ArrSub,ArrPoi))
	allocate(Frac(ArrSub,ArrPoi))
	allocate(LnG(ArrSub,ArrPoi))
	
	allocate(Jac(ArrPoi,ArrPar))
	allocate(TJac(ArrPar,ArrPoi))
	allocate(MatA(ArrPar,ArrPar))
	allocate(MatR(ArrPoi,1))
	allocate(MatG(ArrPar,1))
	allocate(MAtE(ArrPar,ArrPar))
	allocate(MatDBet(ArrPar,1))
	allocate(MatO(ArrPar,ArrPar))
	
	allocate(DevAbsGB(ArrSub,ArrPoi))
	allocate(DevAbsGA(ArrSub,ArrPoi))
	allocate(DevRelGB(ArrSub,ArrPoi))
	allocate(DevRelGA(ArrSub,ArrPoi))
	allocate(LnGA(ArrSub,ArrPoi))
	allocate(LnGB(ArrSub,ArrPoi))
	
	do i=1,2	!Reading Data
		open(12,file=SubName(i))
		if (trim(SubName(i))/='null') then
			TempInt=1
			read(12,*,IOSTAT=FileEnd) Temp(i,TempInt),Frac(i,TempInt),LnG(i,TempInt)
			LnG(i,TempInt)=log(LnG(i,TempInt))
			do while (.not. IS_IOSTAT_END(FileEnd))
				TempInt=TempInt+1
				read(12,*,IOSTAT=FileEnd) Temp(i,TempInt),Frac(i,TempInt),LnG(i,TempInt)
				LnG(i,TempInt)=log(LnG(i,TempInt))
			enddo
			close(12)
			print *, 'SubName ', SubName(i), ' readed ', PoiNum(i), ' points'
		endif
	enddo
	call ParInit(MatDP,ArrPar,UseParA,UseParT)
	print *, 'Initial parameters readed'
	CheckOut=0
	Init=1
	DebEn=0
	MetK=1
	MetV=2.0
	MetE1=0.005
	MetE2=0.005
	MetT=0.01
	MetMu=0.1
	MyP1=1.3
	MyP2=1.2
	SumDL=1000000.0
	do while (CheckOut==0)
		!call GetYac(DbPar,ArrPar1,ArrPar2,Jac,ArrPoi,ArrPar,UseParA,UseParT)
		!subroutine GetYacNRTL(ParKol,PoiKol,SKol,Param,MatYac,Tem,Frac,PSKol,AType,TType)
		if (AppType==1) then
			call GetYacNRTL(ArrPar,ArrPoi,ArrSub,MatDP,Jac,Temp,frac,PoiNum,UseParA,UseParT)
		else if (AppType==2) then
			call GetYacPower(ArrPar,ArrPoi,ArrSub,MatDP,Jac,Temp,frac,PoiNum,UseParA,UseParT)
		endif
		TJac=transpose(Jac)
		MatA=matmul(TJac,Jac)
		if (AppType==1) then
			call GaNRTL(ArrPoi,ArrPar,ArrSub,MatDP,PoiNum,LnGB,LnG,DevAbsGB,DevRelGB,Temp,&
		&frac,UseParA,UseParT)
		else if (AppType==2) then
			call GaPower(ArrPoi,ArrPar,ArrSub,MatDP,PoiNum,LnGB,LnG,DevAbsGB,DevRelGB,Temp,&
		&frac,UseParA,UseParT)
		endif
		TempInt=1
		do i=1,ArrSub
			do j=1,PoiNum(i)
				MatR(TempInt,1)=DevAbsGB(i,j)
				TempInt=TempInt+1
			enddo
		enddo
		MatG=matmul(TJac,MatR)
		if (DebEn==1) then
			call debug(Jac,ArrPoi,ArrPar,'jac.dbg',len('jac.dbg'))
			call debug(TJac,ArrPar,ArrPoi,'tjac.dbg',len('tjac.dbg'))
			call debug(MatR,ArrPoi,1,'matr.dbg',len('matr.dbg')) 
			call debug(MatA,ArrPar,ArrPar,'mata.dbg',len('mata.dbg'))
			call debug(MatG,ArrPar,1,'matg.dbg',len('matg.dbg'))
		endif
		if (Init==1) then
			init=0
			MetMu=MatA(1,1)
			do i=2,ArrPar
				if (MatA(i,i)>MetMu) then
					MetMu=MatA(i,i)
				endif
			enddo
			MetMu=MetMu*MetT
		endif
		call GetMatE(MatE,ArrPar,MetMu)
		do i=1,ArrPar
			do j=1,ArrPar
				MatO(i,j)=MatA(i,j)
				if (i==j) then
					MatO(i,j)=MatO(i,j)+MetMu*MatE(i,j)	!MatA(i,j) 
				endif
			enddo
		enddo
		call GetDBet(MatG,MatO,MatDBet,ArrPar)
		if (DebEn==1) then
			call debug(MatO,ArrPar,ArrPar,'mato.dbg',len('mato.dbg'))
		endif
		MatDPN=MatDP+MatDBet
		if (DebEn==1) then
			call debug(MatDBet,ArrPar,1,'db.dbg',len('db.dbg'))
			call debug(MatDP,ArrPar,1,'dp.dbg',len('dp.dbg'))
			call debug(MatDPN,ArrPar,1,'dpn.dbg',len('dpn.dbg'))
		endif
		if (AppType==1) then
			call GaNRTL(ArrPoi,ArrPar,ArrSub,MatDPN,PoiNum,LnGA,LnG,DevAbsGA,DevRelGA,Temp,&
			&frac,UseParA,UseParT)
		else if (AppType==2) then
			call GaPower(ArrPoi,ArrPar,ArrSub,MatDPN,PoiNum,LnGA,LnG,DevAbsGA,DevRelGA,Temp,&
			&frac,UseParA,UseParT)
		endif
		!call GetL(ArrPar,MatDBet,MatG,MetMu,LOut)
		MetQ=0.0
		DAA=0.0
		DAB=0.0
		do i=1,ArrSub
			SumDAA(i)=0.0
			SumDAB(i)=0.0
			do j=1,PoiNum(i)
				SumDAB(i)=SumDAB(i)+DevRelGB(i,j)*DevRelGB(i,j)
				SumDAA(i)=SumDAA(i)+DevRelGA(i,j)*DevRelGA(i,j)
				MetQ=MetQ+(DevAbsGB(i,j)*DevAbsGB(i,j)-DevAbsGA(i,j)*DevAbsGA(i,j))
			enddo
			DAA=DAA+SumDAA(i)
			DAB=DAB+SumDAB(i)
		enddo
		DAA=sqrt(DAA)/ArrPoi
		DAB=sqrt(DAB)/ArrPoi
		
		if (mod(MetK,100)==0) then
			write(6,'(a,i5,a,f15.10,a,f10.4)') 'step out', MetK, ' μ ', MetMu, ' minimized function ' ,DAB
!			write(6,'(a,e10.4,a,e10.4,a,e10.4)') 'first ',abs(maxval(MatDBet)), ' > ',MetE2*(maxval(DbPar)+MetE2)&
!			& ,' second ', maxval(MatG)
!			write(6,'(a,e10.4,a,e10.4)') 'DAB  ' ,DAB, '  - DAA  ', DAA
			!pause
		endif
!!!!!-------------------------------------------
		if (DAA*0.99>DAB) then
			MetMu=MetMu*MyP1
		else
			if ((DAB-DAA)/DAB>0.05) then
				MetMu=MetMu*MyP1
			endif
			if ((DAB-DAA)/DAB<0.0001) then
				MetMu=MetMu/MyP2
			endif
		endif
		if (MetK==500) then
			MyP2=MyP1*1.1
		endif
		if (mod(MetK,10000)==0) then
			MetE1=MetE1*1.1
			print *, '!!!!! Errors Up'
			print *, ' New E1 = ', MetE1
		endif
		if (MetMu<MuMin) then
			MetMu=MuMin
		endif
		if (MetMu>1.1) then
			MetMu=1.1
		endif
		MatDP=MatDPN
		MetK=MetK+1
!!!!-----------------------------------------------
!		if (DAA<0.05) then
!			MetMu=MetMu*2.5
!		endif
!-art
!		if ((abs(maxval(MatG))<MetE1).or.(abs(maxval(MatDBet))<MetE2*(abs(maxval(MatDP))+MetE2))) then
!			CheckOut=1
!		else
!			CheckOut=0
!		endif
		if (DAA>MetE1) then
			CheckOut=0
		else
			CheckOut=1
		endif
		if (MetK>StepNum) then
			CheckOut=1
		endif
!		LTemp=0.5*matmul(transpose(MatDBet),MetMu*MatDBet-MatG)
!		Lout=LTemp(1,1)
!		MetQ=MetQ/Lout
!		if (MetQ<0.0) then
!			MetMu=MetMu*MetV
!			MetV=MetV*1.02
!		else
!			MetMu=MetMu/MetV
!			MetV=1.02
!			MatDP=MatDPN
!			MetK=MetK+1
!		endif
!-art
		CurDL=0.0
		do i=1,2
			do j=1,PoiNum(i)
				CurDL=CurDL+abs((exp(LnG(i,j))-exp(LnGA(i,j)))/exp(LnG(i,j)))
			enddo
		enddo
		if (CurDL<SumDL) then
			SumDL=CurDL
			OptMemPar=MatDP
		endif
		!pause
	enddo
	if (AppType==1) then
		call GaNRTL(ArrPoi,ArrPar,ArrSub,OptMemPar,PoiNum,LnGB,LnG,DevAbsGB,DevRelGB,Temp,&
		&frac,UseParA,UseParT)
	else if (AppType==2) then
		call GaPower(ArrPoi,ArrPar,ArrSub,OptMemPar,PoiNum,LnGB,LnG,DevAbsGB,DevRelGB,Temp,&
		&frac,UseParA,UseParT)
	endif
	call PointOut(ArrPar,ArrPoi,ArrSub,LnG,LnGB,frac,Temp,DevAbsGB,DevRelGB,1,SubName(1),PoiNum&
	&,UseParA,UseParT,OptMemPar)
	call PointOut(ArrPar,ArrPoi,ArrSub,LnG,LnGB,frac,Temp,DevAbsGB,DevRelGB,2,SubName(2),PoiNum&
	&,UseParA,UseParT,OptMemPar)
	print *, 'Sucsesfully done'
	
	print *, 'program end'
	
end program

subroutine TauFunc(TauOut,TempIn,TypeIn,Param,n,m,nst)
	implicit none
	real(8) TauOut
	real(8) TempIn
	integer(4) TypeIn
	real(8) Param(n,m)
	integer(4) n,m,nst
	
	if (TypeIn==1) then
		TauOut=Param(nst,1)
	endif
	if (TypeIn==2) then
		TauOut=Param(nst,1)+Param(nst,2)/TempIn
	endif
	if (TypeIn==3) then
		TauOut=Param(nst,1)+Param(nst,2)/TempIn+Param(nst,3)/TempIn/TempIn
		!print *,Param(nst,1),Param(nst,2),Param(nst,3),TauOut,TempIn
		!pause
	endif
	if (TypeIn==4) then
		TauOut=Param(nst,1)+Param(nst,2)/TempIn+Param(nst,3)/TempIn/TempIn+Param(nst,4)*log(TempIn)
	endif
!	if (TypeIn==5) then
!		TauOut=Param(nst,1)+Param(nst,2)/TempIn+Param(nst,3)/TempIn/TempIn+Param(nst,4)&
!		&*log(TempIn)+Param(nst,5)*TempIn**Param(nst,6)
!	endif
!	if (TypeIn==6) then
!		TauOut=Param(nst,1)+Param(nst,2)/TempIn+Param(nst,3)/TempIn/TempIn+Param(nst,4)&
!		&*log(TempIn)+Param(nst,5)*TempIn**Param(nst,6)
!	endif
end subroutine

subroutine DTauFunc(DTauOut,TempIn,TypeIn,Param,n,m,nst)
	implicit none
	real(8) DTauOut
	real(8) TempIn
	integer(4) TypeIn
	real(8) Param(n,m)
	integer(4) n,m,nst
	
	if (TypeIn==1) then
		DTauOut=1.0
	endif
	if (TypeIn==2) then
		DTauOut=1.0/TempIn
	endif
	if (TypeIn==3) then
		DTauOut=1.0/TempIn/TempIn
	endif
	if (TypeIn==4) then
		DTauOut=log(TempIn)
	endif
!	if (TypeIn==5) then
!		DTauOut=TempIn**Param(nst,6)
!	endif
!	if (TypeIn==6) then
!		DTauOut=Param(nst,5)*TempIn**Param(nst,6)*log(TempIn)
!	endif
end subroutine

subroutine AlphFunc(AlphOut,TempIn,TypeIn,Param,n,m,nst)
	implicit none
	real(8) AlphOut
	real(8) TempIn
	integer(4) TypeIn
	real(8) Param(n,m)
	integer(4) n,m,nst
	
	if (TypeIn==1) then
		AlphOut=Param(nst,1)
	endif
	if (TypeIn==2) then
		AlphOut=Param(nst,1)+Param(nst,2)*TempIn
	endif
	
end subroutine

subroutine DAlphFunc(DAlphOut,TempIn,TypeIn,Param,n,m,nst)
	implicit none
	real(8) DAlphOut
	real(8) TempIn
	integer(4) TypeIn
	real(8) Param(n,m)
	integer(4) n,m,nst
	
	if (TypeIn==1) then
		DAlphOut=1.0
	endif
	if (TypeIn==2) then
		DAlphOut=TempIn
	endif
	
end subroutine

subroutine GiFunc(Gi,a,t)
	implicit none
	real(8) Gi,a,t
	
	Gi=exp(-a*t)
	if (Gi<0.0000001) then
		Gi=0.0000001
	endif
end subroutine

subroutine DGiFunc(DGi,a,t,d)
	implicit none
	real(8) DGi,a,t,d
	DGi=-d*exp(-a*t)
	
end subroutine

subroutine GetYacNRTL(ParKol,PoiKol,SKol,Param,MatYac,Tem,Frac,PSKol,AType,TType)
	implicit none
	integer(4) ParKol,PoiKol,SKol
	real(8) Param(ParKol,1)
	real(8) MatYac(PoiKol,ParKol)
	real(8) Tem(SKol,PoiKol)
	real(8) Frac(SKol,PoiKol)
	integer(4) PSKol(SKol)
	integer(4) AType,TType
	
	integer(4) i,j,k
	real(8) t12,t21,a
	real(8) g12,g21
	real(8) dadp,dtdp
	real(8) dgda12,dgda21,dgdt12,dgdt21
	real(8) x2,x1
	integer(4) TempInt,TempInt1,TempInt2
	!integer(4) fArrPar
	!real(8),allocatable:: 
	!fArrPar=nint((SKol+(Skol-1)*(SKol-1)-1)/2)+SKol*Skol-Skol
	real(8) fPar(3,10)
	integer(4) n,m
	
	n=3
	m=10
	TempInt=1
	TempInt1=1
	do i=1,nint((SKol+(SKol-1)*(SKol-1)-1)/2.0)	!Changing to internal parameters
		do j=1,AType
			fPar(TempInt1,j)=Param(TempInt,1)
			TempInt=TempInt+1
		enddo
		TempInt1=TempInt1+1
	enddo
	do i=1,SKol*SKol-SKol
		do j=1,TType
			fPar(TempInt1,j)=Param(TempInt,1)
			TempInt=TempInt+1
		enddo
		TempInt1=TempInt1+1
	enddo
	!debug
!	print *, nint((SKol+(SKol-1)*(SKol-1)-1)/2.0),SKol*SKol-SKol
!	do i=1,ParKol
!		write(6,'(f10.6)') Param(i,1)
!	enddo
!	write(6,'(a)') 'next'
!	do i=1,3
!		do j=1,10
!			write(6,'(f10.6,a,$)') fPar(i,j),' '
!		enddo
!		write(6,'(a)') ' '
!	enddo
!	pause
	TempInt1=1
	do i=1,Skol
		do j=1,PSKol(i)
		call TauFunc(t12,Tem(i,j),TType,fPar,n,m,2)	!second string t12
		call TauFunc(t21,Tem(i,j),TType,fPar,n,m,3)	!
		call AlphFunc(a,Tem(i,j),AType,fPar,n,m,1)	
		call GiFunc(g12,a,t12)
		call GiFunc(g21,a,t21)
		call DGiFunc(dgdt12,a,t12,a)
		call DGiFunc(dgdt21,a,t21,a)
		call DGiFunc(dgda12,a,t12,t12)
		call DGiFunc(dgda21,a,t21,t21)
		
		if (i==1) then
			x2=(1.0-frac(i,j))
			x1=frac(i,j)
			TempInt2=1
!			print *,Tem(i,j), t12,t21,a,g12,g21,x1,x2
!			pause
			do k=1,AType	!a diagonal
				call DAlphFunc(dadp,Tem(i,j),k,fPar,n,m,1)
!				print *, ' a ', dadp
				MatYac(TempInt1,TempInt2)=x2*x2*(2.0*t21*g21*dgda21/(g21*x2+x1)/(g21*x2+x1)-&
				&2.0*t21*g21*g21*dgda21*x2/(g21*x2+x1)/(g21*x2+x1)/(g21*x2+x1)+&
				&t12*dgda12/(x2+g12*x1)/(x2+g12*x1)-&
				&2.0*t12*g12*dgda12*x1/(x2+g12*x1)/(x2+g12*x1)/(x2+g12*x1))*dadp
!				print *, MatYac(TempInt1,TempInt2),g12,2.0*t12*g12*dgda12*x1/(x2+g12*x1)/(x2+g12*x1)/(x2+g12*x1)
				TempInt2=TempInt2+1
			enddo
			do k=1,TType	!t12
				call DTauFunc(dtdp,Tem(i,j),k,fPar,n,m,2)
			!	print *, ' t12 ', dtdp
				MatYac(TempInt1,TempInt2)=x2*x2*(t12*dgdt12/(x2+g12*x1)/(x2+g12*x1)+&
				&g12/(x2+g12*x1)/(x2+g12*x1)-&
				&2.0*t12*g12*dgdt12*x1/(x2+g12*x1)/(x2+g12*x1)/(x2+g12*x1))*dtdp
				TempInt2=TempInt2+1
			enddo
			do k=1,TType	!t21
				call DTauFunc(dtdp,Tem(i,j),k,fPar,n,m,3)
			!	print *, 't21 ', dtdp
				MatYac(TempInt1,TempInt2)=x2*x2*(2.0*t21*g21*dgdt21/(g21*x2+x1)/(g21*x2+x1)+&
				&g21*g21/(g21*x2+x1)/(g21*x2+x1)-&
				&2.0*t21*g21*g21*dgdt21*x2/(g21*x2+x1)/(g21*x2+x1)/(g21*x2+x1))*dtdp
				TempInt2=TempInt2+1
			enddo
		endif
		if (i==2) then
			x2=(frac(i,j))
			x1=(1.0-frac(i,j))
			TempInt2=1
			do k=1,AType	!a diagonal
				call DAlphFunc(dadp,Tem(i,j),k,fPar,n,m,1)
			!	print *, dadp, 'a',Tem(i,j),k,n,m
				MatYac(TempInt1,TempInt2)=x1*x1*(t21*dgda21/(g21*x2+x1)/(g21*x2+x1)-&
				&2.0*t21*g21*dgda21*x2/(g21*x2+x1)/(g21*x2+x1)/(g21*x2+x1)+&
				&2.0*t12*g12*dgda12/(x2+g12*x1)/(x2+g12*x1)-&
				&2.0*t12*g12*g12*dgda12*x1/(x2+g12*x1)/(x2+g12*x1)/(x2+g12*x1))*dadp
				TempInt2=TempInt2+1
			enddo
			do k=1,TType	!t12
				call DTauFunc(dtdp,Tem(i,j),k,fPar,n,m,2)
			!	print *, dtdp, 't12',Tem(i,j),k
				MatYac(TempInt1,TempInt2)=x1*x1*(2.0*t12*g12*dgdt12/(x2+g12*x1)/(x2+g12*x1)+&
				&g12*g12/(x2+g12*x1)/(x2+g12*x1)-&
				&2.0*t12*g12*g12*dgdt12*x1/(x2+g12*x1)/(x2+g12*x1)/(x2+g12*x1))*dtdp
				TempInt2=TempInt2+1
			enddo
			do k=1,TType	!t21
				call DTauFunc(dtdp,Tem(i,j),k,fPar,n,m,3)
			!	print *, dtdp, 't21'
				MatYac(TempInt1,TempInt2)=x1*x1*(t21*dgdt21/(g21*x2+x1)/(g21*x2+x1)+&
				&g21/(g21*x2+x1)/(g21*x2+x1)-&
				&2.0*t21*g21*dgdt21*x2/(g21*x2+x1)/(g21*x2+x1)/(g21*x2+x1))*dtdp
				TempInt2=TempInt2+1
			enddo
		endif
		TempInt1=TempInt1+1
		enddo
	enddo
end subroutine

subroutine ParInit(Param,n,a,t)
	implicit none
	integer(4) n,a,t
	real(8) Param(n,1)
	
	integer(4) i,TempInt
	do i=1,n
		Param(i,1)=0.0
	enddo
	TempInt=1
	if (a==1) then
		Param(TempInt,1)=0.5
		TempInt=TempInt+1
	endif
	if (a==2) then
		Param(TempInt,1)=0.5+0.3*rand()
		TempInt=TempInt+1
		Param(TempInt,1)=0.0001*rand()
		TempInt=TempInt+1
	endif
	do i=1,t
		Param(TempInt,1)=15.0*float(i)   !+15.0*rand()
		if (i==4) then
			Param(TempInt,1)=0.03
		endif
		TempInt=TempInt+1
	enddo
	do i=1,t
		Param(TempInt,1)=15.0*float(i)   !+15.0*rand()
		if (i==4) then
			Param(TempInt,1)=0.03
		endif
		TempInt=TempInt+1
	enddo
	
end subroutine

subroutine MatTrans(MatIn,MatOut,n,m)
	implicit none
	integer(4) n,m
	real(8) MatIn(n,m)
	real(8) MatOut(m,n)
	integer(4) i,j
	
	MatOut=0.0
	do i=1,n
		do j=1,m
			MatOut(j,i)=MatIn(i,j)
		enddo
	enddo
	
end subroutine

subroutine MatPro(MatIn1,n1,m1,MatIn2,m2,n2,MatOut,nn,rr,mm)
	implicit none
	integer(4) n1,m1,n2,m2,nn1,nn2
	real(8) MatIn1(n1,m1)
	real(8) MatIn2(m2,n2)
	real(8) MatOut(n1,n2)
	integer(4) nn,mm,rr
	integer(4) i,j,k
	real(8) Msum
	
	!print *, n1,n2,nn,mm,rr
	do i=1,n1
		do j=1,n2
			!print *, i,j,MatOut(i,j)
			MatOut(i,j)=0.0
		enddo
	enddo
	
	do i=1,nn
		do j=1,mm
			Msum=0.0
			do k=1,rr
				MatOut(i,j)=Msum+MatIn1(i,k)*MatIn2(k,j)
			enddo
		enddo
	enddo
	
end subroutine

subroutine MuMax(MatIn,n,m,MuOut)
	implicit none
	integer(4) n,m
	real(8) MatIn(n,n)
	real(8) MuOut
	real(8) Mu
	integer(4) i
	
	Mu=MatIn(1,1)
	do i=1,m
		if(MatIn(i,i)>Mu) then
			Mu=MatIn(i,i)
		endif
	enddo
	MuOut=Mu
end subroutine

subroutine  GetMatE(MatOut,m,mu)
	implicit none
	integer(4) m
	real(8) mu
	real(8) MatOut(m,m)
	
	integer(4) i,j
	
	do i=1,m
		do j=1,m
			MatOut(i,j)=0.0
		enddo
		MatOut(i,i)=mu
	enddo
	
end subroutine

subroutine GetInMat(pw,n,nm,j1)
integer(4)  nm,n,j1,m,j,i,k,L
real(8) pw(nm,nm)
real(8) qw(100,100),jh(100)
real(8)  P,Q
	j1=1
	m=0
	do 100 j=1,n
 100    jh(j)=1
	do 110 i=1,n
	do 110 j=1,n
	qw(i,j)=0.
	if (i.eq.j) qw(i,j)=1.
 110    continue
 120    p=0.
	do 130 i=1,n
	do 130 j=1,n
	if((jh(i).lt.0).or.(jh(j).lt.0)) goto 130
	if(abs(pw(i,j)).le.p) goto 130
	p=abs(pw(i,j))
	k=i
	l=j
 130    continue
	if (m.eq.0) q=p
	if ((p/q).ge.1.e-20) goto 140
	j1=-1
	return
 140    m=m+1
	do 150 j=1,n
	if (j.ne.l) pw(k,j)=pw(k,j)/pw(k,l)
 150    qw(k,j)=qw(k,j)/pw(k,l)
	do 160 i=1,n
	do 160 j=1,n
	if(i.eq.k) goto 160
	if(j.ne.l) pw(i,j)=pw(i,j)-pw(i,l)*pw(k,j)
	qw(i,j)=qw(i,j)-pw(i,l)*qw(k,j)
 160    continue
	if(k.eq.l) goto 190
	do 170 j=1,n
	p=pw(l,j)
	pw(l,j)=pw(k,j)
 170    pw(k,j)=p
	do 180 j=1,n
	p=qw(l,j)
	qw(l,j)=qw(k,j)
 180    qw(k,j)=p
 190    jh(l)=-1
	if (m.lt.n) goto 120
	do 200 i=1,n
	do 200 j=1,n
 200    pw(i,j)=qw(i,j)
	return
	end
	
subroutine GetDBet(MG,MO,MDB,n)
	implicit none
	integer(4) n
	real(8) MG(n,1)
	real(8) MTemp(n,1)
	real(8) MO(n,n)
	real(8) MDB(n,1)
	integer(4) i,j,error
	do i=1,n
		do j=1,1
			MTemp(i,j)=MG(i,j)
		enddo
	enddo
!	call debug(MO,n,n,'MatO.dbg',len('MatO.dbg'))
	call GetInMat(MO,n,n,error)
!	call debug(MO,n,n,'MatO-1.dbg',len('MatO-1.dbg'))
	MDB=matmul(MO,MTemp)
!	call debug(MDB,n,1,'DP.dbg',len('DP.dbg'))
end subroutine

subroutine Par2Mat(DP,n1,m1,MDP,n2,PA,PT)
	implicit none
	integer(4) n1,m1,n2,PA,PT
	real(8) DP(n1,m1)
	real(8) MDP(n2,1)
	integer(4) i,Tempint
	
	Tempint=1
	do i=1,PA
		MDP(TempInt,1)=DP(1,i)
		Tempint=TempInt+1
	enddo
	do i=1,PT
		MDP(TempInt,1)=DP(2,i)
		Tempint=Tempint+1
	enddo
	do i=1,PT
		MDP(TempInt,1)=DP(3,i)
		TempInt=TempInt+1
	enddo
	
end subroutine

subroutine Mat2Par(DP,n1,m1,MDP,n2,PA,PT)
	implicit none
	integer(4) n1,m1,n2,PA,PT
	real(8) DP(n1,m1)
	real(8) MDP(n2,1)
	integer(4) i,Tempint
	
	Tempint=1
	do i=1,PA
		DP(1,i)=MDP(TempInt,1)
		TempInt=TempInt+1
	enddo
	do i=1,PT
		DP(2,i)=MDP(TempInt,1)
		TempInt=TempInt+1
	enddo
	do i=1,PT
		DP(3,i)=MDP(TempInt,1)
		TempInt=TempInt+1
	enddo
	
end subroutine

subroutine debug(arr,n,m,str,k)
	implicit none
	integer(4) n,m,k
	character(k) str
	real(8) arr(n,m)
	
	integer(4) i,j
	
	open(32,file=trim(str))
		do i=1,n
			do j=1,m
				write(32,'(e16.6,$)') arr(i,j)
			enddo
			write(32,'(a)') 
		enddo
	close(32)
	
end subroutine

subroutine GaNRTL(fAPoi,fAPar,fASub,Param,NPoi,gnrtl,gex,DevAbs,DevRel,Tem,frac,AType,TType)
	implicit none
	integer(4) fAPoi,fAPar,fASub
	real(8) Param(fAPar,1)
	integer(4) NPoi(fASub)
	real(8) gnrtl(fASub,fAPoi)
	real(8) gex(fASub,fAPoi)
	real(8) Tem(fASub,fAPoi)
	real(8) frac(fASub,fAPoi)
	real(8) DevAbs(fASub,fAPoi)
	real(8) DevRel(fASub,fAPoi)
	integer(4) AType,TType
	
	integer(4) i,j,k
	real(8) t12,t21,a
	real(8) g12,g21
	real(8) fPar(10,10)
	integer(4) n,m
	integer(4) TempInt
	real(8) x1,x2
	n=10
	m=10
	do i=1,n
		do j=1,10
			fPar(i,j)=0.0
		enddo
	enddo
	TempInt=1
	do i=1,AType
		fPar(1,i)=Param(TempInt,1)
		TempInt=TempInt+1
	enddo
	do i=1,TType
		fPar(2,i)=Param(TempInt,1)
		TempInt=TempInt+1
	enddo
	do i=1,TType
		fPar(3,i)=Param(TempInt,1)
		TempInt=TempInt+1
	enddo
	!debug
!	do i=1,3
!		do j=1,5
!			write(6,'(f10.6,a,$)') fPar(i,j), ' '
!		enddo
!		write(6, '(a)') ' '
!	enddo
	
	do i=1,fASub
		do j=1,NPoi(i)
			call TauFunc(t12,Tem(i,j),TType,fPar,n,m,2)	!second string t12
			call TauFunc(t21,Tem(i,j),TType,fPar,n,m,3)	!
			call AlphFunc(a,Tem(i,j),AType,fPar,n,m,1)	
			call GiFunc(g12,a,t12)
			call GiFunc(g21,a,t21)
			!print *,Tem(i,j),' parameters ','t12', t12,'t21',t21,'a',a,'g12',g12,'g21',g21
			
			if (i==1) then
				x2=1.0-frac(i,j)
				x1=frac(i,j)
				gnrtl(i,j)=x2*x2*(t21*(g21/(x1+x2*g21))*(g21/(x1+x2*g21))+&
				&t12*g12/(x2+x1*g12)/(x2+x1*g12))
			else if (i==2) then
				x2=frac(i,j)
				x1=1.0-frac(i,j)
				gnrtl(i,j)=x1*x1*(t12*(g12/(x2+x1*g12))*(g12/(x2+x1*g12))+&
				&t21*g21/(x1+x2*g21)/(x1+x2*g21))
			endif
			DevAbs(i,j)=gex(i,j)-gnrtl(i,j)
			if (abs(gex(i,j))>0.001) then
				DevRel(i,j)=DevAbs(i,j)/gex(i,j)
			else
				DevRel(i,j)=0.0
			endif
			!print *,'temp',Tem(i,j),'x1',x1,'x2',x2,'gex', gex(i,j),'nrtl',gnrtl(i,j)
		enddo
	enddo
	!pause
end subroutine

subroutine GetL(fAPar,DB,G,Mu,L)
	implicit none
	integer(4) fAPar
	real(8) DB(fAPar,1)
	real(8) G(fAPar,1)
	real(8) Mu
	real(8) L
	
	real(8) DBT(1,fAPar)
	real(8) Lmat(1,1)
	integer(4) i
	real(8) sumL
	
	DBT=transpose(DB)
	sumL=0.0
!	do i=1,fAPar
!		sumL=sumL+0.5*DBT(1,i)*(Mu*DB(i,1)-G(i,1))
!	enddo
!	L=sumL
	Lmat= matmul(0.5*DBT,Mu*DB-G)
	L=Lmat(1,1)
end subroutine

subroutine PointOut(fAPar,fAPoi,fASub,Gex,Gn,frac,Tem,DAb,DRe,nm,FileOut,fASubPoi,AType,TType,Par)
	implicit none
	integer(4) fAPar,fAPoi,fASub
	real(8) Gex(fASub,fAPoi)
	real(8) frac(fASub,fAPoi)
	real(8) Tem(fASub,fAPoi)
	real(8) Gn(fASub,fAPoi)
	real(8) DAb(fASub,fAPoi)
	real(8) DRe(fASub,fAPoi)
	integer(4) nm
	character(20) FileOut
	character(20) FileOut2
	integer(4) i,j,k
	integer(4) fASubPoi(fASub)
	integer(4) AType,TType
	real(8) Par(fAPar,1)
	
	real(8) AvAbsLnG,MaxAbsLnG
	real(8) AvRelLnG,MaxRelLnG
	real(8) AvAbsG,MaxAbsG
	real(8) AvRelG,MaxRelG
	
	integer,dimension(8) :: values
	integer(4) TempInt
	
	write(FileOut2,'(a,a)') trim(FileOut),'.out'
	
	
	
	open(22,file=trim(adjustl(FileOut2)))
		write(22,'(a)') ' Frac() Temp ln(γe) ln(γNRTL)   lnγe-ln(γNRTL)   abs(ln(γe)-ln(γNRTL)) &
		& (ln(γex)-ln(γNRTL))/ln(γex)  abs((ln(γe)-ln(γNRTL))/ln(γe))  γe   γNRTl   γex-γNRTL  &
		& abs(γe-γNRTL) (γe-γNRTL)/γe  abs((γe-γNRTL)/γe)'
		do i=1,fASubPoi(nm)
			write(22,'(13(f16.10,a1),e16.10,$ )') frac(nm,i), ' ',Tem(nm,i),' ',Gex(nm,i),&
			&' ',Gn(nm,i),' ',DAb(nm,i),' ',abs(DAb(nm,i)),' ',DRe(nm,i),' ',&
			&abs(DRe(nm,i)),' ', exp(Gex(nm,i)),' ', exp(Gn(nm,i)),' ', exp(Gex(nm,i))- exp(Gn(nm,i)),' ',&
			& abs(exp(Gex(nm,i))- exp(Gn(nm,i))),' ',(exp(Gex(nm,i))- exp(Gn(nm,i)))/exp(Gex(nm,i)),' ',&
			& abs((exp(Gex(nm,i))- exp(Gn(nm,i)))/exp(Gex(nm,i)))
			write(22,'(a)') ' '
			
		enddo
	close(22)
	MaxAbsG=-10000.0
	MaxRelG=-10000.0
	MaxAbsLnG=-10000.0
	MaxRelLnG=-10000.0
	AvAbsG=0.0
	AvRelG=0.0
	AvAbsLnG=0.0
	AvRelLnG=0.0
	do i=1,2
		do j=1,fASubPoi(i)
			AvAbsG=AvAbsG+abs(exp(Gex(i,j))- exp(Gn(i,j)))
			AvRelG=AvRelG+abs((exp(Gex(i,j))- exp(Gn(i,j)))/exp(Gex(i,j)))
			if (abs(exp(Gex(i,j))- exp(Gn(i,j)))>MaxAbsG) then
				MaxAbsG=abs(exp(Gex(i,j))- exp(Gn(i,j)))
			endif
			if (abs((exp(Gex(i,j))- exp(Gn(i,j)))/exp(Gex(i,j)))>MaxRelG) then
				MaxRelG=abs((exp(Gex(i,j))- exp(Gn(i,j)))/exp(Gex(i,j)))
			endif
			AvAbsLnG=AvAbsLnG+abs(Gex(i,j)-Gn(i,j))
			if(abs(Gex(i,j)-Gn(i,j))>MaxAbsLnG) then
				MaxAbsLnG=abs(Gex(i,j)-Gn(i,j))
			endif
			if (abs(Gex(i,j))>0.000001) then
				AvRelLnG=AvAbsLnG+abs((Gex(i,j)-Gn(i,j))/Gex(i,j))
				if(abs((Gex(i,j)-Gn(i,j))/Gex(i,j))>MaxRelLnG) then
					MaxRelLnG=abs((Gex(i,j)-Gn(i,j))/Gex(i,j))
				endif
			endif
		enddo
	enddo
	AvAbsG=AvAbsG/(fASubPoi(1)+fASubPoi(2))
	AvAbsLnG=AvAbsLnG/(fASubPoi(1)+fASubPoi(2))
	AvRelG=AvRelG/(fASubPoi(1)+fASubPoi(2))
	AvRelLnG=AvRelLnG/(fASubPoi(1)+fASubPoi(2))
	open(23,file='results.out')
		call date_and_time(VALUES=values)
		write(23,'((i4,x,i2,x,i2,a,i2,a,i2,a,i2))') values(1),values(2),values(3) , ' - ', values(5),&
		':',values(6),':',values(7)
		write(23,'(a)') ' '
		write(23,'(a)') ' α approximation : '
		if (AType==1) then
			write(23,'(a)') 'α=A'
		endif
		if (AType==2) then
			write(23,'(a)') 'α=A+B*T'
		endif
		write(23,'(a)') ' Parameters:'
		if (AType==1) then
			write(23,'(a,f20.10)') 'A = ', Par(1,1)
			TempInt=1
		endif
		if (AType==2) then
			write(23,'(a,f20.10)') 'A = ', Par(1,1)
			write(23,'(a,f20.10)') 'B = ', Par(2,1)
			TempInt=2
		endif
		write(23,'(a)') ' τ approximation : '
		if (TType==1) then
			write(23,'(a)') 'τ(ij)=A(ij)'
		endif
		if (TType==2) then
			write(23,'(a)') 'τ(ij)=A(ij)+B(ij)/T'
		endif
		if (TType==3) then
			write(23,'(a)') 'τ(ij)=A(ij)+B(ij)/T+C(ij)/T^2'
		endif
		if (TType==4) then
			write(23,'(a)') 'τ(ij)=A(ij)+B(ij)/T+C(ij)/T^2+D(ij)*ln(T)'
		endif
		write(23,'(a)') ' Parameters:'
		do i=1,2
			do j=1,2
				if (i/=j) then
					write(23,'(a,i1,i1)') ' τ',i,j
					do k=1,TType
						TempInt=TempInt+1
						if (k==1) then
							write(23,'(a,i1,i1,a,f20.10)') 'A(', i,j,') = ', Par(TempInt,1)
						endif
						if (k==2) then
							write(23,'(a,i1,i1,a,f20.10)') 'B(', i,j,') = ', Par(TempInt,1)
						endif
						if (k==3) then
							write(23,'(a,i1,i1,a,f20.10)') 'C(', i,j,') = ', Par(TempInt,1)
						endif
						if (k==4) then
							write(23,'(a,i1,i1,a,f20.10)') 'D(', i,j,') = ', Par(TempInt,1)
						endif
					enddo
				endif
			enddo
		enddo
		write(23,'(a)') '                                     ln(γ)                   γ'
		write(23,'(a,2f20.10)') 'average absolute error',AvAbsLnG,AvAbsLnG
		write(23,'(a,2f20.10)') 'maximum absolute error',MaxAbsLnG,MaxAbsLnG
		write(23,'(a,2f20.10)') 'average relative error',AvRelLnG,AvRelG
		write(23,'(a,2f20.10)') 'maximum relative error',MaxRelLnG,MaxRelG
	close(23)
end subroutine

subroutine GetYacPower(ParKol,PoiKol,SKol,Param,MatYac,Tem,Frac,PSKol,AType,TType)
	implicit none
	integer(4) ParKol,PoiKol,SKol
	real(8) Param(ParKol,1)
	real(8) MatYac(PoiKol,ParKol)
	real(8) Tem(SKol,PoiKol)
	real(8) Frac(SKol,PoiKol)
	integer(4) PSKol(SKol)
	integer(4) AType,TType
	
	integer(4) i,j,k
	integer(4) TempInt
	
	TempInt=1
	do i=1,SKol
		do j=1,PSKol(i)
			do k=1,ParKol
				MatYac(TempInt,k)=frac(i,j)**(k-1)
			enddo
			TempInt=TempInt+1
		enddo
	enddo
	
end subroutine

subroutine GaPower(fAPoi,fAPar,fASub,Param,NPoi,gnrtl,gex,DevAbs,DevRel,Tem,frac,AType,TType)
	implicit none
	integer(4) fAPoi,fAPar,fASub
	real(8) Param(fAPar,1)
	integer(4) NPoi(fASub)
	real(8) gnrtl(fASub,fAPoi)
	real(8) gex(fASub,fAPoi)
	real(8) Tem(fASub,fAPoi)
	real(8) frac(fASub,fAPoi)
	real(8) DevAbs(fASub,fAPoi)
	real(8) DevRel(fASub,fAPoi)
	integer(4) AType,TType
	
	integer(4) i,j,k
	integer(4) TempInt
	real(8) TempReal
	
	do i=1,fASub
		do j=1,NPoi(i)
			TempReal=0.0
			do k=1,fAPar
				TempReal=TempReal+Param(k,1)*frac(i,j)**(k-1)
			enddo
			gnrtl(i,j)=TempReal
			DevAbs(i,j)=gex(i,j)-gnrtl(i,j)
			if(abs(gex(i,j))>0.0001) then
				DevRel(i,j)=DevAbs(i,j)/gex(i,j)
			else
				DevRel(i,j)=0.0
			endif
		enddo
	enddo
end subroutine

subroutine OptPar()
	implicit none
	
end subroutine
