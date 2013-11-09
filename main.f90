!дома ночь 7-8 12
module generator    !модуль генератора
 integer(4) n1, n2, n3 !случаные числа
 integer(4) m1,a1,b1,a2,b2,m3,m2 !параметры генератора
 integer(4) max1, max2,max3 !максимумы генераторов
 real(8) randmass(500) !массив случайных чисел
 real(8) outrand
contains
subroutine randomn
!use generator
    implicit none

integer(4) i
real(8) sseed1,sseed2,sseed3
integer(4) seed(20)
character(20) a,b,c
real(8) n11,n21,n31
call date_and_time(a,b,c,seed)
sseed1=0.0
sseed2=0.0
sseed3=0.0
!print *,seed
do i=1,20,1 !здесь шарной бред
 sseed1=sseed1+seed(i)
 if (mod(i,2)==0) then
  sseed2=sseed2+seed(i)
 else
  sseed2=sseed2/2.0
 endif
 if (mod(i,3)==0) then
  sseed3=sseed3+seed(i)
 else
  sseed3=sseed3-seed(i)/2.2
 endif
enddo

!pause
!проверка  тригонометрических функций
n11=abs(cos(sseed1/1000))
n21=abs(sin(sseed2/1000))
n31=abs(cos(sseed3/1000))
print *,'- Random seeds -', n11,n21,n31
!call random_seed(int(sseed))
!call random_number(n11)
!call random_number(n21)
!call random_number(n31)
m1=2147483563
a1=40014
b1=12345
m2=2147483399
a2=40692
b2=54321
max3=2147483647
max1=m1-1
max2=m2-1
n1=ceiling((n11*max1)-1)
n2=ceiling((n21*max2)-1)
n3=ceiling((n31*max3)-1)
!   open(15,file='rand.txt')
 do i=1,500
 call randstart()
 randmass(i)=outrand
 enddo
print *,'--- Generator Start ---'
end subroutine
!*******************************************************
!Рандомная функция начало
subroutine  randstart()
!use generator
    implicit none

n1=abs(mod(a1*n1+b1,m1))
n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
if (float(n3)/float(max3)<0.5) then
 outrand=float(n1)/float(max1)
else
 outrand=float(n2)/float(max2)
endif
end subroutine
!*******************************************************
!Рандомная функция
function getrand()
!use generator
    implicit none
real(8) getrand
integer(4) repick
!n1=abs(mod(a1*n1+b1,m1))
!n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
repick=ceiling((float(n3)/float(max3))*500)
getrand=randmass(repick)
call randstart()
randmass(repick)=outrand
end function

end module


program dcvchem
	use generator
	implicit none
	real(8) BoxL1,BoxL2
	real(8) BoxG1,BoxG2
	real(8) BoxH
	character(8) TempString,TempString1
	integer(4) TotalAtom
	integer(4) i,j,k,m
	integer(4) k2,i2,j2
	integer(4) ii,jj,kk
	integer(4) SubNum
	character(5),allocatable:: SubName(:)
	integer(4),allocatable:: NLiq(:)
	integer(4),allocatable:: SubAtomNum(:)
	character(5),allocatable:: SubAtomName(:,:,:)
	real(8),allocatable:: SubAtomX(:,:,:),SubAtomY(:,:,:),SubAtomZ(:,:,:)
	real(8),allocatable:: SubAtomVX(:,:,:),SubAtomVY(:,:,:),SubAtomVZ(:,:,:)
	integer(4) MemAtom
	real(8),allocatable:: MemX(:),MemY(:),MemZ(:)
	character(5), allocatable:: MemName(:)
	integer(4) TempInt,TempInt1,TempInt2
	real(8) TempR1,TempR2,TempR3
	real(8) Sigma
	real(8) MemDelta
	integer(4) Mem1HW
	integer(4) Mem1Len
	real(8) RoL
	integer(4) TotalMol
	
	real(8),allocatable:: MolX(:,:),MolY(:,:),MolZ(:,:)
	integer(4),allocatable:: InOutV1(:,:),InOutV2(:,:)
	integer(4) MemType
	character(20),allocatable:: SubFile(:),SubFileINP(:)
	real(8),allocatable:: SubSigma(:,:)
	real(8),allocatable:: SubEps(:,:)
	real(8),allocatable:: SubMass(:,:)
	real(8),allocatable:: SubCha(:,:)
	
	real(8) SumX,SumY,SumZ,SumMass
	real(8),allocatable:: Vol1X(:,:,:),Vol1Y(:,:,:),Vol1Z(:,:,:)
	real(8),allocatable:: Vol2X(:,:,:),Vol2Y(:,:,:),Vol2Z(:,:,:)
	real(8),allocatable:: Vol3X(:,:,:),Vol3Y(:,:,:),Vol3Z(:,:,:)
	real(8),allocatable:: Vol1VX(:,:,:),Vol1VY(:,:,:),Vol1VZ(:,:,:)
	real(8),allocatable:: Vol2VX(:,:,:),Vol2VY(:,:,:),Vol2VZ(:,:,:)
	real(8),allocatable:: Vol3VX(:,:,:),Vol3VY(:,:,:),Vol3VZ(:,:,:)
	integer(4),allocatable:: SumInV1(:),SumInV2(:),SumOut(:)
	real(8),allocatable:: CenAtomX(:,:),CenAtomY(:,:),CenAtomZ(:,:)
	real(8),allocatable:: CheckX(:),CheckY(:),CheckZ(:)
	real(8) RandX,RandY,RandZ
	
	integer(4) NStep, step
	real(8) Vol1Ak,Vol2Ak
	integer(4) ch, CheckType
	real(8) DeltaEn,drx,dry,drz,rx,ry,rz,r,LJ
	real(8),allocatable:: MixSig(:,:,:,:),MixEps(:,:,:,:),MixCha(:,:,:,:)
	real(8) BoxVol,Temp,Prob
	real(8),allocatable:: AkL(:),AkG(:)
	real(8),allocatable:: Vol1MX(:,:),Vol1MY(:,:),Vol1MZ(:,:)
	real(8),allocatable:: Vol2MX(:,:),Vol2MY(:,:),Vol2MZ(:,:)
	real(8),allocatable:: Vol3MX(:,:),Vol3MY(:,:),Vol3MZ(:,:)
	
	integer(4) StepType,CheckMol
	real(8) rmx,rmy,rmz
	real(8) TotBoxLen
	integer(4) Ntot
	real(8),allocatable:: frac(:)
	character(20) time1,time2,time3
	integer(4) seed(20)
	integer(4) rseed
	integer(4),allocatable:: NMolLiq(:)
	real(8),allocatable:: ToRoLIq(:)
	character(100) LongTempString
	integer(4) TopInt
	
	real(8),allocatable:: Vol1MolX(:,:),Vol1MolY(:,:),Vol1MolZ(:,:)
	real(8),allocatable:: Vol2MolX(:,:),Vol2MolY(:,:),Vol2MolZ(:,:)
	real(8) Mem1B,Mem1E,Mem2B,Mem2E
	real(8) CurTime
	integer(4) First
	real(8) DTime
	
	integer(4) Nslice
	real(8),allocatable:: MolNumSl(:,:)
	integer(4) hist
	
	real(8),allocatable:: InMem1(:),MemIn1(:),MemOut1(:)
	real(8),allocatable:: InMem2(:),MemIn2(:),MemOut2(:)
	real(8),allocatable:: Vol1MolNum(:),Vol2MolNum(:)
	real(8),allocatable:: NearInVol(:),NearOutVol(:)
	
	real(8),allocatable:: ExpSumSl(:,:)
	
	integer(4) SampleNum
	integer(4),allocatable:: InitA(:)
	real(8),allocatable:: InitX(:,:),InitY(:,:),InitZ(:,:)
	real(8),allocatable:: InitS(:,:),InitE(:,:),InitM(:,:),InitQ(:,:)
	real(8),allocatable:: Sigma2(:,:,:,:),Epsi2(:,:,:,:),Qi2(:,:,:,:)
	
	real(8) dr, dxr,dyr,dzr
	real(8) dar,daxr, dayr,dazr,rz2,rz6,sumpot
	integer(4) try,xob,yob,zob
	
	real(8) MemSigma,MemEpsi
	
	call randomn
	
	open(7,file='test.temp')
		read(7,'(f20.10)') BoxL1
		read(7,'(f20.10)') BoxL2
		read(7,'(f20.10)') BoxG1
		read(7,'(f20.10)') BoxG2
		read(7,'(f20.10)') BoxH
		read(7,'(i6)') SubNum
		read(7,'(i6)') TotalMol
		allocate(SubName(SubNum+1))
		allocate(NLiq(SubNum+1))
		allocate(SubAtomNum(SubNum+1))
		allocate(SubFile(SubNum+1))
		allocate(SubFileINP(SubNum+1))
		allocate(frac(SubNum+1))
		allocate(AkG(SubNum+1))
		allocate(AkL(SubNum+1))
		allocate(ToRoLiq(SubNum+1))
		do i=1,SubNum
			read(7,'(a5,2i6)') SubName(i),NLiq(i),SubAtomNum(i)
		enddo
		read(7,'(i6)') MemAtom
		read(7,'(f20.10)') TotBoxLen
		read(7,'(f20.10)') Mem1B
		read(7,'(f20.10)') Mem1E
		read(7,'(f20.10)') Mem2B
		read(7,'(f20.10)') Mem2E
		read(7,'(f20.10)') CurTime
		print *,CurTime
		read(7,'(i6)') First
		allocate(SubAtomName(SubNum+1,TotalMol,30))
		allocate(SubAtomX(SubNum+1,TotalMol,30))
		allocate(SubAtomY(SubNum+1,TotalMol,30))
		allocate(SubAtomZ(SubNum+1,TotalMol,30))
		allocate(SubAtomVX(SubNum+1,TotalMol,30))
		allocate(SubAtomVY(SubNum+1,TotalMol,30))
		allocate(SubAtomVZ(SubNum+1,TotalMol,30))
		allocate(MemName(MemAtom))
		allocate(MemX(MemAtom))
		allocate(MemY(MemAtom))
		allocate(MemZ(MemAtom))
	close(7)
	print *, 'MAIN IN FILE ReAD DONE'
	
	open(9,file='main.in')
		read(9,*) TempString
		read(9,*) SubNum
		do i=1,SubNum
			read(9,*) SubFile(i),frac(i), AkL(i),AkG(i),ToRoLIq(i)
			write(SubFileINP(i),'(2a)') trim(adjustl(SubFile(i))),'.inp'
		enddo
		read(9,*) TempString
		read(9,*) MemType
		if (MemType==1) then
		read(9,*) TempString
			read(9,*) Sigma
			read(9,*) TempString
			read(9,*) MemDelta
			read(9,*) TempString
			read(9,*) Mem1HW
			read(9,*) TempString
			read(9,*) Mem1Len
			read(9,*) TempString
			read(9,*) RoL
		endif
		read(9,*) TempString
		read(9,*) NStep
		read(9,*) TempString
		read(9,*) Temp
		read(9,*) TempString
		read(9,*) TopInt
		read(9,*) TempString
		read(9,*) DTime
	close(9)
	print *, 'MEMIN File read DONE'
		
	allocate(SubSigma(SubNum+1,30))
	allocate(SubEps(SubNum+1,30))
	allocate(SubMass(SubNum+1,30))
	allocate(SubCha(SubNum+1,30))
	
	do i=1,SubNum
		print *,SubFile(i),SubAtomNum(i)
		open(41,file=SubFileINP(i))
		do j=1,SubAtomNum(i)
			read(41,*) SubSigma(i,j),SubEps(i,j),SubMass(i,j),SubCha(i,j)
			print *, i,j,SubSigma(i,j),SubEps(i,j),SubMass(i,j),SubCha(i,j)
		enddo
		close(41)
	enddo
	open(8,file='test2.gro')
		read(8,'(a)') TempString
		read(8,'(i6)') TotalAtom
		do i=1,SubNum
			do j=1,NLiq(i)
				do k=1,SubAtomNum(i)
					read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, TempString,SubAtomName(i,j,k),&
					&TempInt1,SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
					&SubAtomVX(i,j,k),SubAtomVY(i,j,k),SubAtomVZ(i,j,k)
!					print *,i,j,k, TempInt, TempString,SubAtomName(i,j,k),&
!					&TempInt1,SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
!					&SubAtomVX(i,j,k),SubAtomVY(i,j,k),SubAtomVZ(i,j,k)
				enddo
			enddo
		enddo
		do i=1,MemAtom
			read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt,TempString,MemName(i),&
			&TempInt1,MemX(i),MemY(i),MemZ(i),&
			&TempR1,TempR2,TempR3
		enddo
	close(8)
	print *, ' GRO FILE READ DONE '
	
	!замена координат
	allocate(Vol1MolX(SubNum+1,TotalAtom))
	allocate(Vol1MolY(SubNum+1,TotalAtom))
	allocate(Vol1MolZ(SubNum+1,TotalAtom))
	allocate(Vol2MolX(SubNum+1,TotalAtom))
	allocate(Vol2MolY(SubNum+1,TotalAtom))
	allocate(Vol2MolZ(SubNum+1,TotalAtom))
	allocate(MolX(SubNum+1,TotalAtom))
	allocate(MolY(SubNum+1,TotalAtom))
	allocate(MolZ(SubNUm+1,TotalAtom))
	
	do i=1,SubNum
		!print *, NLiq(i)
		do j=1,NLiq(i)
			!print *, i
			SumX=0.0
			SumY=0.0
			SumZ=0.0
			SumMass=0.0
			!print *, SubAtomNum(i)
			do k=1,SubAtomNum(i)
				SumX=SumX+SubAtomX(i,j,k)*SubMass(i,k)
				!print *, SubAtomX(i,j,k)*SubMass(i,k)
				SumY=SumY+SubAtomY(i,j,k)*SubMass(i,k)
				!print *, SubAtomY(i,j,k)*SubMass(i,k)
				SumZ=SumZ+SubAtomZ(i,j,k)*SubMass(i,k)
				!print *, SubAtomZ(i,j,k)*SubMass(i,k)
				SumMass=SumMass+SubMass(i,k)
			enddo
			MolX(i,j)=SumX/SumMass
			MolY(i,j)=SumY/SumMass
			MolZ(i,j)=SumZ/SumMass
			!print *,' test 22'
		enddo
	enddo
	
	do i=1,SubNum
		do j=1,NLiq(i)
			Vol2MolX(i,j)=MolX(i,j)
			Vol2MolY(i,j)=MolY(i,j)
			Vol2MolZ(i,j)=MolZ(i,j)
		enddo
	enddo
	
		print *, ' test 1'
	
	Nslice=ceiling(TotBoxLen/(MemDelta*Sigma/2.0))
	open(56,file='nslice.temp',access='append')
		write(56,'(i10)') Nslice
	close(56)
	allocate(MolNumSl(SubNum+1,Nslice+1))
	allocate(Vol1MolNum(SubNum+1))
	allocate(Vol2MolNum(SubNum+1))
	allocate(MemIn1(SubNum+1))
	allocate(MemOut1(SubNum+1))
	allocate(InMem1(SubNum+1))
	allocate(MemIn2(SubNum+1))
	allocate(MemOut2(SubNum+1))
	allocate(InMem2(SubNum+1))
	
	allocate(ExpSumSl(SubNum+1,Nslice+1))
	
	print *, 'test 2'
	if (First==1) then
		!First=0
		do i=1,SubNum
			do j=1,Nslice
				ExpSumSl(i,j)=0.0
			enddo
		enddo
		SampleNum=0
	else
		
		open(29,file='chem.temp')
		
		do i=1,SubNum
			do j=1,Nslice
				read(29,'(f20.5)') ExpSumSl(i,j)
			enddo
		enddo
		read(29,*) SampleNum
		close(29)
	endif
!			open(28,file='Mu.out')
!			write(28,'(a,$)')  'r,(nm)   '
!			do i=1,SubNum
!				write(28,'(a,$)')  'chem   '
!			enddo
!			write(28,'(a)')  '   '
!			close(28)
	SampleNum=SampleNum+1
	
	call date_and_time(time1,time2,time3,seed)
	rseed=ceiling(abs(cos(float(seed(8))/1000.0))*20000.0)
	print *,' Random integer ', rseed
	call srand(rseed)
	
	allocate(InitA(SubNum+1))
	allocate(InitX(SubNum+1,30))
	allocate(InitY(SubNum+1,30))
	allocate(InitZ(SubNum+1,30))
	allocate(InitS(SubNum+1,30))
	allocate(InitE(SubNum+1,30))
	allocate(InitM(SubNum+1,30))
	allocate(InitQ(SubNum+1,30))
	
	do i=1,SubNum
		write(SubFile(i),'(2a)') trim(adjustl(SubFile(i))),'.gro'
		open(41,file=SubFile(i))
			read(41,*) TempString
			read(41,*) InitA(i)
			do j=1,InitA(i)
				read(41,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, TempString, TempString,&
					&TempInt1,InitX(i,j),InitY(i,j),InitZ(i,j),&
					&TempR1,TempR2,TempR3
			enddo
		close(41)
		!write(SubFileINP(i),'(2a)') trim(adjustl(SubFile(i))),'.gro'
		print *, SubFileINP(i)
		open(42,file=SubFileINP(i))
			do j=1,InitA(i)
				read(42,*) InitS(i,j),InitE(i,j),InitM(i,j),InitQ(i,j)
			enddo
		close(42)
	enddo
	print *,'test 3'
	allocate(Sigma2(30,30,30,30))
	allocate(Epsi2(30,30,30,30))
	allocate(Qi2(30,30,30,30))
	do i=1,SubNum
		do j=1,SubNum
			do k=1,InitA(i)
				do m=1,InitA(j)
					Sigma2(i,j,k,m)=(InitS(i,k)+InitS(j,m))/2.0
					Epsi2(i,j,k,m)=sqrt(InitE(i,k)*InitE(j,m))
					Qi2(i,j,k,m)=InitQ(i,k)*InitQ(j,m)
				enddo
			enddo
		enddo
	enddo
	
	print *, 'test 4'
		
	do i=1,SubNum
		Sumx=0.0
		SumY=0.0
		SumZ=0.0
		SumMass=0.0
		do j=1,InitA(i)
			SumX=SumX+InitX(i,j)*InitM(i,j)
			SumY=SumY+InitY(i,j)*InitM(i,j)
			SumZ=SumZ+InitZ(i,j)*InitM(i,j)
			SumMass=SumMass+InitM(i,j)
		enddo
		InitX(i,j)=SumX/SumMass
		InitY(i,j)=SumY/SumMass
		InitZ(i,j)=SumZ/SumMass
	enddo
	print *, 'Initial start'
	print *,' Temperature', Temp
	
	do i=1,SubNum
		do j=1,Nslice
			print *, i,j,' of  ', Nslice
			do try=1,30000
				RandX=getrand()*BoxH
				RandY=getrand()*BoxH
				RandZ=getrand()*TotBoxLen/float(Nslice)+TotBoxLen/float(Nslice)*float(j)
				!добавить вращение
				sumpot=0.0
				do i2=1,SubNum
					do j2=1,NLiq(i2)
						dxr=abs(RandX-MolX(i2,j2))
						dyr=abs(RandY-MolY(i2,j2))
						dzr=abs(RandZ-MolZ(i2,j2))
						if (dxr>0.5*BoxH) then
							dxr=BoxH-dxr
							xob=1
						else
							xob=0
						endif
						if (dyr>0.5*BoxH) then
							dyr=BoxH-dyr
							yob=1
						else
							yob=0
						endif
						if (dzr>0.5*TotBoxLen) then
							dzr=TotBoxLen-dzr
							zob=1
						else
							zob=0
						endif
						dr=sqrt(dxr*dxr+dyr*dyr+dzr*dzr)
						if (dr<BoxH/2.0) then
							do ii=1,InitA(i)
								do jj=1,InitA(i2)
									if (xob==0) then
										daxr=abs(RandX+InitX(i,ii)-SubAtomX(i2,j2,jj))
									else
										daxr=BoxH-abs(RandX+InitX(i,ii)-SubAtomX(i2,j2,jj))
									endif
									!print *, xob,daxr,dxr
									if (yob==0) then
										dayr=abs(RandY+InitY(i,ii)-SubAtomY(i2,j2,jj))
									else
										dayr=BoxH-abs(RandY+InitY(i,ii)-SubAtomY(i2,j2,jj))
									endif
									!print *, yob,dayr,dyr
									if (zob==0) then
										dazr=abs(RandZ+InitZ(i,ii)-SubAtomZ(i2,j2,jj))
									else
										dazr=TotBoxLen-abs(RandZ+InitZ(i,ii)-SubAtomZ(i2,j2,jj))
									endif
									!print *, zob,dazr,dzr
									dar=sqrt(daxr*daxr+dayr*dayr+dazr*dazr)
									rz=Sigma2(i,i2,ii,jj)/dar
									rz2=rz*rz
									rz6=rz2*rz2*rz2
									sumpot=sumpot+4.0*(rz6*rz6-rz6)*Epsi2(i,i2,ii,jj)  !+Qi2(i,i2,ii,jj)/dar
									!print *,dr, dar, sumpot, 4.0*(rz6*rz6-rz6)*Epsi2(i,i2,ii,jj)
									!print *, Sigma2(i,i2,ii,jj), Epsi2(i,i2,ii,jj)
									!pause
								enddo
							enddo
						endif
					enddo
				enddo
				do i2=1,MemAtom
					dxr=abs(RandX-MemX(i2))
					dyr=abs(RandY-MemY(i2))
					dzr=abs(RandZ-MemZ(i2))
					if (dxr>0.5*BoxH) then
						dxr=BoxH-dxr
						xob=1
					else
						xob=0
					endif
					if (dyr>0.5*BoxH) then
						dyr=BoxH-dyr
						yob=1
					else
						yob=0
					endif
					if (dzr>0.5*TotBoxLen) then
						dzr=TotBoxLen-dzr
						zob=1
					else
						zob=0
					endif
					dr=sqrt(dxr*dxr+dyr*dyr+dzr*dzr)
					if (dr<BoxH/2.0) then
						do ii=1,InitA(i)
							if (xob==0) then
								daxr=abs(RandX+InitX(i,ii)-MemX(i2))
							else
								daxr=BoxH-abs(RandX+InitX(i,ii)-MemX(i2))
							endif
							!print *, xob,daxr,dxr
							if (yob==0) then
								dayr=abs(RandY+InitY(i,ii)-MemY(i2))
							else
								dayr=BoxH-abs(RandY+InitY(i,ii)-MemY(i2))
							endif
							!print *, yob,dayr,dyr
							if (zob==0) then
								dazr=abs(RandZ+InitZ(i,ii)-MemZ(i2))
							else
								dazr=TotBoxLen-abs(RandZ+InitZ(i,ii)-MemZ(i2))
							endif
							!print *, zob,dazr,dzr
							dar=sqrt(daxr*daxr+dayr*dayr+dazr*dazr)
							if (SubNum==1) then
								MemSigma=0.3
								MemEpsi=1.0
							endif
							if (SubNum==2) then
								MemSigma=0.3
								MemEpsi=2.0
							endif
							rz=MemSigma/dar
							rz2=rz*rz
							rz6=rz2*rz2*rz2
							sumpot=sumpot+4.0*(rz6*rz6-rz6)*MemEpsi  !+Qi2(i,i2,ii,jj)/dar
							!print *,dr, dar, sumpot, 4.0*(rz6*rz6-rz6)*Epsi2(i,i2,ii,jj)
							!pause
						enddo
					endif
				enddo
				!if (sumpot<0) then
				!	print *, sumpot,exp(-sumpot/Temp),temp,ExpSumSl(i,j)
					ExpSumSl(i,j)=ExpSumSl(i,j)+exp(-sumpot/Temp)
				!	pause
				!endif
			enddo
		enddo
	enddo
	
	open(29,file='chem.temp')
		do i=1,SubNum
			do j=1,Nslice
				write(29,'(f40.20)') ExpSumSl(i,j)
			enddo
		enddo
		write(29,*) SampleNum
	close(29)
	
	open(28,file='Mu.out')
			write(28,'(a,$)')  'r,(nm)   '
			do i=1,SubNum
				write(28,'(a)')  'chem   '
			enddo
		close(28)
	
	open(18,file="Mu.out",access="append")
		do i=1,Nslice
			write(18,'(f20.10,$)') TotBoxLen*float(i)/float(Nslice)
			do j=1,SubNum
				write(18,'(f20.10,a,$)') -log(ExpSumSl(j,i)/float(SampleNum)/30000.0),' '
			enddo
			write(18,'(a)') ' '
		enddo
		
	close(18)
	
end program





