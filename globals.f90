! Global variables
module m_globals
implicit none

! Input parameters, see default-prm.py for documentation
integer, dimension(3) :: np3, nn
real,dimension(3) :: dx
integer :: nt,debug, itio,itstats,mpin,mpout
real :: dt,tm0,tm

! Error indicator
integer,parameter :: ierror = -9999
real,parameter :: rerror = -9999.9999
integer :: &
    it,             & ! current time step
    its,			& ! start it index
    ite,			& ! end it index
    ditse,			& ! interval it index
    ifn,            & ! fault normal component=abs(faultnormal)
    ip,             & ! process rank
    ipid,           & ! processor ID
    np0               ! number of processes available
integer, dimension(3) :: &
    nl3,            & ! number of mesh nodes per process
    nm,             & ! size of local 3D arrays
    nhalo,          & ! number of ghost nodes
    ip3,            & ! 3D process rank
    ip3root,        & ! 3D root process rank
    ip2root,        & ! 2D root process rank
    i1core, i2core, & ! core region
    nnoff             ! offset between local and global indices
logical :: &
    sync,           & ! synchronize processes
    verb,           & ! print messages
    bigendian,		& ! indicate big-endian
    master            ! master process flag
character(256) :: &
    str,		    & ! string for storing file names
    infile,			& ! input data folder
    oufile,         & ! output data folder
    lattice           ! method for station on focal sphere (Fibonacci or lat-lon)
real, allocatable, target, dimension(:,:,:) :: &
    s1, s2		    ! temporary storage  
real, allocatable, target, dimension(:,:,:,:) :: &
    w1, w2

! user customized
!**********************************************
integer :: & !receivers
	nt0, &	 ! P wave time nodes
	nt1, &   ! pP wave time nodes
	nt2, &   ! sP wave time nodes
	nt3, &   ! S wave time nodes
	nrec,&   ! number of total receivers 
	ntheta,& ! number of receivers along theta
	nphi,&   ! number of receivers along phi
    nfibonacci, & !number of receivers for Fibonacci spiral grids
	comp,&   ! output component (0 is amplitude, 1 is x, 2 is y, 3 is z)
	endian,& ! endianness (1 is big-endian input; 0 is little-endian input)
	usrorig,&!# 0-> origin is at origin in SORD coordinate 
			 !#1-> origin on free surface but horizontal coordinate is adjusted to center of fault
             !# 2-> origin is average to fault
	nsourceoffset !how many identical source are read in only changing locations
	
real ::  & !receivers
	dphi,dtheta,ophi,otheta,t0s,t0e,t1s,t1e,t2s, &
        t2e,t3s,t3e,radius,trvmin,trvmax
real :: &
    rho,           & !density kg/m^3
    vp,            & !P wave velocity m/s
    vs               !S wave velocity m/s
real :: &
    x0,     & !user define origin (average of fault)
    y0,     &
    z0

real,allocatable,dimension(:) :: &
    sourceoffset    !this vector will move source when computing waveforms but the receivers
    				! are kept the original based on x0,y0,z0. New source center is
    				! (x0,y0,z0)+sourceoffset
real,allocatable,dimension(:,:) :: xyzoffset	!xoffset yoffset zoffset timeoffset			
!real :: eachsourceoffset(3)    	
    				
real,allocatable,dimension(:,:,:) :: xx,yy,zz
real,allocatable,dimension(:) :: rx,ry,rz
real,allocatable,dimension(:,:,:,:) :: mrij
real,allocatable,dimension(:,:) :: lsummrij,gsummrij
real,allocatable,dimension(:,:,:) :: P,PP,SP,SV,SH,S !3 components
real,allocatable,dimension(:,:) :: recP,recPP,recSP,recS, recSH, recSV

real,allocatable,dimension(:,:,:,:,:) :: &
    tp,           & !P wave travelling time on each subfault
    tpp,          & !pP wave travelling time
    tsp,          & !sP wave travelling time
    ts,		& !S wave travelling time on each subfault
    ap,           & !P wave amplitude
    app,          & !pP wave amplitude
    asp,          & !sP wave amplitude
    as,		& !S wave amplitude
    rayp,         & !ray parameter for P wave
    raypp,        & !ray parameter for pP wave
    raysp,		& !ray parameter for sP wave
    rays		  !ray parameter for S wave
!**********************************************

end module
