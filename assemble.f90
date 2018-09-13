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
!write(str,'(I1)') comp
! P wave
if (nt0 > 1 .and. t0s > 0.) then
	allocate(f(nt0,nrec)) !P wave
	if(comp == 0) then
		f = sqrt(sum(P * P, 3))
		call rreduce2(recP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave0'//trim(tmpstr))
	elseif(comp==1) then

		! 1st component
		f = P(:,:,1)
		call rreduce2(recP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave1'//trim(tmpstr))

		! 2nd component
		f = P(:,:,2)
		call rreduce2(recP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave2'//trim(tmpstr))

		! 3rd component
		f = P(:,:,3)
		call rreduce2(recP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave3'//trim(tmpstr))
	else
		write(0,*) 'comp should be 0 or 1';stop
	end if 
	deallocate(f)
end if

!S wave
if (nt3 > 1 .and. t3s > 0.) then
	allocate(f(nt3,nrec)) 
	
	! S wave
	if(comp == 0) then
		f = sqrt(sum(S * S, 3))
		call rreduce2(recS,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave0'//trim(tmpstr))
	elseif(comp==1) then

		! 1st component
		f = S(:,:,1)
		call rreduce2(recS,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave1'//trim(tmpstr))

		! 2nd component
		f = S(:,:,2)
		call rreduce2(recS,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave2'//trim(tmpstr))

		! 3rd component
		f = S(:,:,3)
		call rreduce2(recS,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave3'//trim(tmpstr))
	else
		write(0,*) 'comp should be 0 or 1';stop
	end if 

	! SV wave
	if(comp == 0) then
		f = sqrt(sum(SV * SV, 3))
		call rreduce2(recSV,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave0'//trim(tmpstr))
	elseif(comp==1) then

		! 1st component
		f = SV(:,:,1)
		call rreduce2(recSV,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave1'//trim(tmpstr))

		! 2nd component
		f = SV(:,:,2)
		call rreduce2(recSV,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave2'//trim(tmpstr))

		! 3rd component
		f = SV(:,:,3)
		call rreduce2(recSV,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave3'//trim(tmpstr))
	else
		write(0,*) 'comp should be 0 or 1';stop
	end if 

	! SH wave
	if(comp == 0) then
		f = sqrt(sum(SH * SH, 3))
		call rreduce2(recSH,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave0'//trim(tmpstr))
	elseif(comp==1) then

		! 1st component
		f = SH(:,:,1)
		call rreduce2(recSH,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave1'//trim(tmpstr))

		! 2nd component
		f = SH(:,:,2)
		call rreduce2(recSH,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave2'//trim(tmpstr))

		! 3rd component
		f = SH(:,:,3)
		call rreduce2(recSH,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave3'//trim(tmpstr))
	else
		write(0,*) 'comp should be 0 or 1';stop
	end if 	
	
	deallocate(f)
end if

!pP wave
if (nt1 > 1 .and. t1s > 0.) then
	allocate(f(nt1,nrec)) !P wave
	if(comp == 0) then
		f = sqrt(sum(PP * PP, 3))
		call rreduce2(recPP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave0'//trim(tmpstr))
	elseif(comp==1) then

		! 1st component
		f = PP(:,:,1)
		call rreduce2(recPP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave1'//trim(tmpstr))

		! 2nd component
		f = PP(:,:,2)
		call rreduce2(recPP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave2'//trim(tmpstr))

		! 3rd component
		f = PP(:,:,3)
		call rreduce2(recPP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave3'//trim(tmpstr))
	else
		write(0,*) 'comp should be 0 or 1';stop
	end if 
	deallocate(f)
end if

!sP wave
if (nt2 > 1 .and. t2s > 0.) then
	allocate(f(nt2,nrec)) !P wave
	if(comp == 0) then
		f = sqrt(sum(SP * SP, 3))
		call rreduce2(recSP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave0'//trim(tmpstr))
	elseif(comp==1) then

		! 1st component
		f = SP(:,:,1)
		call rreduce2(recSP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave1'//trim(tmpstr))

		! 2nd component
		f = SP(:,:,2)
		call rreduce2(recSP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave2'//trim(tmpstr))

		! 3rd component
		f = SP(:,:,3)
		call rreduce2(recSP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave3'//trim(tmpstr))
	else
		write(0,*) 'comp should be 0 or 1';stop
	end if 
	deallocate(f)
end if

end subroutine

end module