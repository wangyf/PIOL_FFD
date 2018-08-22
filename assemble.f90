module m_assemble
IMPLICIT  NONE
contains

! generate timeseries from P, pP, sP
subroutine assembly
use m_globals
use m_collective
use m_fieldio

integer :: rootcor(3),i
real,allocatable,dimension(:,:) :: f
character(256) :: tmpstr

rootcor = 0

! appendix about it range
write(tmpstr,'(A,I7.7,A,I7.7)') '-',its,'-',ite
!output summed moment rate tensor
if (master) call mstrwrite2d(gsummrij,trim(oufile)//'/summr'//trim(tmpstr))

! output components of phases or together
write(str,'(I1)') comp
! P wave
if (nt0 > 1 .and. t0s > 0.) then
	allocate(f(nt0,nrec)) !P wave
	if(comp == 0) then
		f = sqrt(sum(P * P, 3))
	elseif(comp<4) then
		f = P(:,:,comp)
	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	call rreduce2(recP,f, 'allsum', rootcor)
	deallocate(f)
	if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave'//trim(str)//trim(tmpstr))
end if

!S wave
if (nt3 > 1 .and. t3s > 0.) then
	allocate(f(nt3,nrec)) 
	
	! S wave
	if(comp == 0) then
		f = sqrt(sum(S * S, 3))
	elseif(comp<4) then
		f = S(:,:,comp)
	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	call rreduce2(recS,f, 'allsum', rootcor)
	if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave'//trim(str)//trim(tmpstr))
	
	! SV wave
	if(comp == 0) then
		f = sqrt(sum(SV * SV, 3))
	elseif(comp<4) then
		f = SV(:,:,comp)
	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	call rreduce2(recSV,f, 'allsum', rootcor)
	if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave'//trim(str)//trim(tmpstr))
	
	! SH wave
	if(comp == 0) then
		f = sqrt(sum(SH * SH, 3))
	elseif(comp<4) then
		f = SH(:,:,comp)
	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	call rreduce2(recSH,f, 'allsum', rootcor)
	if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave'//trim(str)//trim(tmpstr))
	
	deallocate(f)
end if


! pP wave
if (nt1 > 1 .and. t1s > 0.) then
	allocate(f(nt1,nrec)) !P wave
	if(comp == 0) then
		f = sqrt(sum(PP * PP, 3))
	elseif(comp<4) then
		f = PP(:,:,comp)
	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	call rreduce2(recPP,f, 'allsum', rootcor)
	deallocate(f)
	if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave'//trim(str)//trim(tmpstr))
end if

! sP wave
if (nt2 > 1 .and. t2s > 0.) then
	allocate(f(nt2,nrec)) !P wave
	if(comp == 0) then
		f = sqrt(sum(SP * SP, 3))
	elseif(comp<4) then
		f = SP(:,:,comp)
	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	call rreduce2(recSP,f, 'allsum', rootcor)
	deallocate(f)
	if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave'//trim(str)//trim(tmpstr))
end if

end subroutine

end module