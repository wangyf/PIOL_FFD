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
! P wave
if (nt0 > 1 .and. t0s > 0.) then
	allocate(f(nt0,nrec)) !P wave
	if(comp == 0) then
		write(str,'(I1)') comp
		f = sqrt(sum(P * P, 3))
		call rreduce2(recP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave'//trim(str)//trim(tmpstr))
	elseif(comp==1) then
		write(str,'(I1)') 1
		f = P(:,:,1)
		call rreduce2(recP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave'//trim(str)//trim(tmpstr))
		
		write(str,'(I1)') 2
		f = P(:,:,2)
		call rreduce2(recP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave'//trim(str)//trim(tmpstr))
		
		write(str,'(I1)') 3
		f = P(:,:,3)
		call rreduce2(recP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recP,trim(oufile)//'/Pwave'//trim(str)//trim(tmpstr))

	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	deallocate(f)
end if

!S wave
if (nt3 > 1 .and. t3s > 0.) then
	allocate(f(nt3,nrec))  !S wave
	if(comp == 0) then
		write(str,'(I1)') comp
		f = sqrt(sum(S * S, 3))
		call rreduce2(recS,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave'//trim(str)//trim(tmpstr))
		
		f = sqrt(sum(SV * SV, 3))
		call rreduce2(recSV,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave'//trim(str)//trim(tmpstr))
		
		f = sqrt(sum(SH * SH, 3))
		call rreduce2(recSH,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave'//trim(str)//trim(tmpstr))
		
	elseif(comp==1) then
		write(str,'(I1)') 1
		f = S(:,:,1)
		call rreduce2(recS,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave'//trim(str)//trim(tmpstr))
		
		f = SV(:,:,1)
		call rreduce2(recSV,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave'//trim(str)//trim(tmpstr))
		
		f = SH(:,:,1)
		call rreduce2(recSH,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave'//trim(str)//trim(tmpstr))
		
		write(str,'(I1)') 2
		f = S(:,:,2)
		call rreduce2(recS,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave'//trim(str)//trim(tmpstr))
		
		f = SV(:,:,2)
		call rreduce2(recSV,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave'//trim(str)//trim(tmpstr))
		
		f = SH(:,:,2)
		call rreduce2(recSH,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave'//trim(str)//trim(tmpstr))
		
		write(str,'(I1)') 3
		f = S(:,:,3)
		call rreduce2(recS,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recS,trim(oufile)//'/Swave'//trim(str)//trim(tmpstr))
		
		f = SV(:,:,3)
		call rreduce2(recSV,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSV,trim(oufile)//'/SVwave'//trim(str)//trim(tmpstr))
		
		f = SH(:,:,3)
		call rreduce2(recSH,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSH,trim(oufile)//'/SHwave'//trim(str)//trim(tmpstr))

	else
		write(0,*) 'comp should be 0 ~ 1';stop
	end if 
	deallocate(f)
end if

! pP wave
if (nt1 > 1 .and. t1s > 0.) then
	allocate(f(nt1,nrec)) !P wave
	if(comp == 0) then
		write(str,'(I1)') comp
		f = sqrt(sum(PP * PP, 3))
		call rreduce2(recPP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave'//trim(str)//trim(tmpstr))
	elseif(comp==1) then
		write(str,'(I1)') 1
		f = PP(:,:,1)
		call rreduce2(recPP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave'//trim(str)//trim(tmpstr))
		
		write(str,'(I1)') 2
		f = PP(:,:,2)
		call rreduce2(recPP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave'//trim(str)//trim(tmpstr))
		
		write(str,'(I1)') 3
		f = PP(:,:,3)
		call rreduce2(recPP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recPP,trim(oufile)//'/pPwave'//trim(str)//trim(tmpstr))

	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	deallocate(f)
end if

! sP wave
if (nt2 > 1 .and. t2s > 0.) then
	allocate(f(nt2,nrec)) !P wave
	if(comp == 0) then
		write(str,'(I1)') comp
		f = sqrt(sum(SP * SP, 3))
		call rreduce2(recSP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave'//trim(str)//trim(tmpstr))
	elseif(comp==1) then
		write(str,'(I1)') 1
		f = SP(:,:,1)
		call rreduce2(recSP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave'//trim(str)//trim(tmpstr))
		
		write(str,'(I1)') 2
		f = SP(:,:,2)
		call rreduce2(recSP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave'//trim(str)//trim(tmpstr))
		
		write(str,'(I1)') 3
		f = SP(:,:,3)
		call rreduce2(recSP,f, 'allsum', rootcor)
		if (master) call mstrwrite2d(recSP,trim(oufile)//'/sPwave'//trim(str)//trim(tmpstr))

	else
		write(0,*) 'comp should be 0 ~ 3';stop
	end if 
	deallocate(f)
end if

end subroutine

end module