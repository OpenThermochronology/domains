! program domains v 1.30 (based on 1.21)

! 9 April 2025: v. 1.21  Patched odd bug found by Brenhin Keller --
!   Some runs, especially those where modeling pushed into skunky data past Kspar breakdown,
!   yielded domain.in files having incorrect Do/a2 values for each domain. This also manifested
!   as the component domain trends on the Arrhenius plot being way off, shoved to the right.
!   Generally the domains values were 1 to 3 units too large (too diffusive). Oddly, this only
!   impacted occasional model runs (but not only those with obviously iffy modeling choices).
! 	I have been unable to pinpoint what data or parameters are causing this.
!
!   I noticed that the calculated Tc values remain correct given the sample's reported Do/ao2 value 
!   and reported sizes. What seemed to be failing were Lovera's ordj values calculated using the ord,
!   xlog, and rpmax variables. Based on the correct Tc values, I used just the value of ord
!   and the size values in the odd members of array a1() to get the needed info. Modest testing
!   seems to show that all isnow  ok using this approach. Each of the changes are documented with
!   this version's 1.30 number (see lines ~650 to 680 of this version). 

!   The original code is too complex to follow easily but I suspect that some inversion issues
! 	are impacting what becomes the rpmax variable. This appears to be part of the original
!   code, not related to changes related to the Domains wrapper I added. CAVEAT EMPTOR!

! April 13, 2025. I decided to simplify code operation by hard-coding a few options, and add this
! 	to the above bug fix to create 1.30. The manual and an alert in this code will note that
!	the min and max number of domains and regression-weighting options can be changed in code
!	with a re-compile.
!   
!   I could have put everything into a text input file but that might be slower and less flexible
! 	for routine operation compared to a few command-line inputs. This splits the difference.
!
!	I have a left in a check that compares the max domain size as previously calculated with
!	the value obtained from the correct a1() array. I also added warnings to the console and report
!   outputs to flag cases where the mismatch has happened.


! ***********************************************************************************************
! PURPOSE: Determine kinetic and domain-distribution parameters for use in MDD diffusion modeling
!          Versions 1.11 and higher handle spherical diffusion geometry

! USAGE: ./domains [SAMPLENAME] [SIZEFILENAME] [PLOTFLAG]
!     [SAMPLENAME] is a string of 10 or fewer characters (extra will be ignored)
!       [PLOTFLAG] is an integer: 1 results in use of gmt plotting commands for display output
!
! AUTHOR: Core routines and algorithm by Oscar Lovera at UCLA (his original autoarr program)
!
!         These revisions, alterations, updates, and clean-ups by:
!
!                  Peter Zeitler
!		Department of Earth and Environmental Sciences
!              1 West Packer Avenue
!               Lehigh University
!              Bethlehem, PA 18015
!                +1-610-758-3671
!             peter.zeitler@lehigh.edu
!
! Send any complaints or bug reports to Zeitler, DO NOT blame or hassle Oscar, who is innocent.
!
! ***** COMPILATION *****
!
! Compile domains like this (you MUST do it this way, using the 'fno-automatic' flag):
!
! gfortran domains-130.f90 -o domains13 -fno-automatic -O2  -fallow-argument-mismatch -w; rm -f *.mod

! (substitute your local source file name for 'domains-130.f90', and your preferred executable name for 'domains13').
! You MUST use the -fno-automatic compiler flag (you did read this, right?). Using O2 optimization then seems to work reliably.
! 
! Depending on your gfortran vintage and how it's configured, domains may compile smoothly
! or the compiler may throw a series of warnings. If you choose the flag -fallow-argument-mismatch
! these should not appear.
!
! If you place the location of the binary in your PATH, you should be able to call the code directly from any directory.
! 
! For Mac OS X, a good source for gcc and gfortran installer packages remains (see top of page):
! 
! http://hpc.sourceforge.net/
! 
! ***** INPUT FILE *****
! 
! Remember that your input file must be in UNIX file format. Data format:
! 
! #-of-steps (note: one less than total since we can't deal with f=1.000)
! temp(degC)   f   time(minutes)
! etc.
! 
! *********************************************************************************************************************
! Changes:
! 27 July 2012: v. 0.985 Mostly just cleaned up introductory comments esp. discussion of weighting. Removed the logD weighting
! option because the simple regression by its nature will weight the early (larger) ABS(logD) values slightly more, which is what we'd prefer.
!
! 2 August 2012: v. 986 Fixed a stupidity: was writing out logRRo values as spectra but kfcorr needs this
! as a simple list. Doh. Because it might be useful to have these values as plottable spectra, changed name
! and added a new output file for the kfcorr program.
!
! 20-21 December 2012: v. 990 Added code to change working directory to the one in which program is launched. This
! allows exectuable to go anywhere in user's path (type echo $PATH to see path; see a UNIX guide for how to make changes),
! and removes need for multiple copies.
! Also, you can now specify a size file having any name as long as it is in .size format
!
! 9 January 2013. Version 1.00. Things seem to be working, so removed the 'beta' status (everything is in beta anyway, eh?).
! Added Tc10 calculation for each domain, and reported this into to output plot and also the report text file. Moved plots up
! on page to more comfortably fit on laptop screen. Changed Ro to ro in labels.
!
! 6 August 2013. Version 1.10b. Started attempt to provide option for spherical geometry. It's a matter of figuring out how Oscar
! used various coefficients. 
!
! 13 August, 2013. Version 1.10b. Domains appears to be working for spherical diffusion geometry, based on quick testing 
! using one sample. Output for slab agrees with previous version and looks good in Size Extractor. Output for sphere looks
! good compared to size extractor. Not yet checked: uncertainties (should be ok).
!
! Also in version 1.10: added comment in summary output file and on ps plot, giving the diffusion geometry used. And fixed glitch
! where domains larger in relative size than 99.9 would not print.

! Version 1.10gmt5 - changed system calls to work with new syntax for gmt5 (basically, put "gmt " in front of all commands
! Version 1.11
! This version actually finished 12/25/14: working with gmt5. Moved temperature labels slightly. Code is now throwing an
! exception (that seems harmless (??): "Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL")
! Need to test if this is from domains code (probably) or gmt (less likely). Due to compiler version??
!
! Version 1.12 was a local offshoot of 1.11 intended for use of analyzing 4He CRH release data.

! Version 1.13 (September 2023) - Just changes to plotting. Used psconvert to create plotfile as pdf.
! Minor changes to wording a few input prompts.

! Version 1.20 (September 2023) - Improved performance for gmt5/6, including plot in pdf form. Clarified wording for some inputs.
! Add option to force activation energy through a selected point.

! Version 1.30 (April 2025) - See above. Bug fix, small simplification to user interface.

! 
! ***********************************************************************

! $$$ *********** Main program *********** 

! Someday someone should explictly declare each variable and then axe the following line (without breaking things!)

module global

integer geometry, option

end module global

! GENERAL WARNINGS for those new to dabbling in fortran. !. Beware the implicit typing of variable
! by their first letter! And also, beware - there is a limit to line lengths for code. For me, this
! is 132 characters beyond that you need to append an & and continue on the next line!!!!!

program domains

	use global
	implicit double precision (a-h,o-z)
	
	parameter(nca=50,ns=200)
	
	logical :: dir_e
	
	dimension a1(nca),a2(nca),auxal(nca,nca),auxa2(nca)
	dimension telab(ns),tilab(ns),auxa1(nca),dyda(nca)
	dimension zi(0:ns),wf(0:ns),xlogd(0:ns),xlogr(0:ns)
	dimension covar(nca,nca),alpha(nca,nca),lista(nca),da(nca)
	dimension f39(0:ns),f(0:ns),wksp(nca),iwksp(nca)
	
	dimension tcten(15)

	external funcs
	
	character  tab1*9, yes*1, sizefile*100, filename*100, eselect*1, systemstring*255, stringfrag*255, gmtfrag*10
	character  directory*50, dirpath*255, morestring*255	
	character samplename*100, command*255
	
    character(len=255) :: path

	double precision imp,maxtemp,siga,sigb,siglogd,sige
	integer gmt, nisample, tempflag, nbreakdown,steplo,stephi,inputs
	integer auxna,guesscount,wopt
	
	tab1=char(09)
	gmt = 0   ! initialize this with a value; real value should supercede this from command line
	geometry = 1  ! 0 = slab 1 = sphere -- give this a value; real value should supercede this from terminal input
	
! get samplename and value of gmt option from commandline using gfortran intrinsic functions

	inputs = command_argument_count()

	if (inputs.lt.3) then
		print *, ' '
		print *, ' domains version 1.11:'
		print *, ' USAGE: domains [samplename] [sizefile_name] [plotflag 1=yes]'
		print *, ' '
		STOP
	end if
	
	if (inputs.gt.3) then
		print *, ' WARNING: extra command-line inputs. Program will continue but be careful.'		
	end if	

	call get_command_argument(1, command)
	read(command,'(A100)') samplename
	
	call get_command_argument(2, command)
	read(command,'(A100)') sizefile
	
	call get_command_argument(3, command)
	read(command,'(I1)') gmt
	
! make sure we are operating in user's working directory

	call getcwd(path)
	call chdir(TRIM(path))

! open input file
	filename=TRIM(sizefile)
	write(*,*) sizefile
	open(unit=10,file=filename,status='old')
	
! initialize some things	
	pi = 3.141592654
	R = 1.98588E-3
	ee = dexp(1.d00) ! I have no immediate idea why Oscar did this!!
	f(0) = 0.
	wf(0) = 0.
	zi(0) = 0.
	b = 8.
	imp = 2.
!	acut = .60
	dchmin = 0.01
	ncons = 0
	ndom = 10   ! program default for maximum number of domains RECOMPILE TO CHANGE
	mdom = 3    ! program default for minimum number of domains RECOMPILE TO CHANGE
!	ndom = 15   ! program seems to work if this is as high as 20; but runs slower
!	mdom = 3
	wopt = 0  ! program default for weighting of observed data (no weighting, errors good) 
				!  RECOMPILE TO CHANGE
	iquality = 0 ! new to version 1.30  - flag that model ran ok - 1 means max domain issue, 2 means 30+ fitting attempts

! greet the user and demonstrate program's start	
	write(*,*) ' '	
	write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'	    
	write(*,*) '                   PROGRAM domains 1.30'
	write(*,*) ' '
	print *,'   Kinetic and domain parameters from Ar-Ar stpheating data'
	print *,' '
	print *,'            Original code autoarr by Oscar Lovera'
	print *,'                 modified by Peter Zeitler'
	print *,' '    
	print *,'                Last updated April 2025'
	print *,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	print *,' '

	write (*,'("First, some necessary questions: ")')
	write (*,*) ' '

	geometry = -1
	do while ((geometry.lt.0).or.(geometry.gt.1))
		write (*,'("Choose diffusion geometry (0 = slab, 1 = sphere): ")',advance='no')
		read *,geometry
	end do

	write (*,'("")')

	maxtemp = -100.
	do while ((maxtemp.gt.2000).or.(maxtemp.lt.800.))
		write (*,'("Enter upper temperature cut-off for modeling (˚C): ")',advance='no')
		read *,maxtemp
	end do
	maxtemp = maxtemp + 273.15
	
	write (*,'("")')

	write (*,'(100a)',advance='no') 'Enter maximum value for log(r/ro) (default=free (0)): '
	read *,auxro
	if (auxro.gt.0) then  ! need to retain and set this legacy flag for use deep in the bowels...
		ncons = 1
	end if

	print *,' '
	print *,'Data from input file - please check:'
	print *,'----------------------------------------'
	print *,'  Step   Temp(C)   Loss     Time(min)'
	print *,'----------------------------------------'	

! get the sample data from file stepheat.size

	read(10,'(I5)') ni
	nimax = ni

! We want to put results in a directory, which means we have to make it, 
! but Fortran can't tell if a directory exists. Because this code compiles under gfortran
! and is intended for UNIX systems, we can use a (non-portable) trick: we look for the '.' file
! that will be located in any directory that exists.

! Note that I decided that since individual result files would get overwritten anyway
! if we worked in a flat-file mode, this code wipes any pre-existing same directory it finds.
! So the users will need to make copies of any previous data if they want to save it. But this will
! be much easier since they just need to change the name of the results directory.

	filename = "./results-"//TRIM(samplename)//"/."

	inquire( file = filename, exist=dir_e )
	if ( dir_e ) then
!		write (*,*) 'directory already exists - wiping it'
		stringfrag = 'rm -R results-'//TRIM(samplename)
		call system(stringfrag)
	end if
	
	stringfrag = "mkdir results-"//TRIM(samplename)
	call system(stringfrag)
	dirpath = "./results-"//TRIM(samplename)//"/"
	
	filename = TRIM(dirpath)//'domains-'//TRIM(samplename)//".in"	
	open(unit=14,file=filename,status='unknown')

	filename = TRIM(dirpath)//'report-'//TRIM(samplename)//".out"	
	open(unit=28,file=filename,status='unknown')	

	tempflag = 0
	do nt = 1,ni
		read(10,*) telab(nt),f39(nt),tilab(nt)
		write(*,'(2X,I3,4X,f7.2,3X,f7.5,3x,f6.2)'),nt,telab(nt),f39(nt),tilab(nt)
		telab(nt) = telab(nt) + 273.15  ! convert temperature to Kelvins
		tilab(nt) = 60.*tilab(nt)         ! convert time to seconds
		if(ni.eq.nimax.and.telab(nt).gt.maxtemp) then   ! was 1423.15 in autoarr - let user determine this
			nimax = nt-1  ! records which step is last one at or below maxtemp 
		end if

		if (tempflag.eq.0.and.telab(nt).gt.maxtemp) then
			tempflag = 1
			nbreakdown = nt - 1
		end if
	end do
	
	close(10)
	
	print *,'----------------------------------------'	
	print *,' '

! First, give user some options about how we get Ea and the Ro reference line
	write (*,'(100a)',advance='no') 'To select E & Ro reference, enter "y". Or enter "n" to proceed: '
	read *,eselect
	print *, ' '
	
! Ok, input done. Now we try start calculating.

! Calculate E,Do, observed Arrhenius and logRRo values.

! Use diff() to get the observed kinetic parameters and calculated some needed information
	call diff(ord,E,f39,telab,tilab,xlogr,xlogd,wf,ni,xro,eselect,steplo,stephi,samplename,dirpath,wopt,siglogd,sige)

! based on first call to diff and user preferences, modify xro to user preference if not flagged as free (0)
	if(auxro.ne.0)xro = auxro

	write (*,'(1X,A4,f7.3,A3,F6.3,A20,f6.3,A4,F5.3)'),'E = ',e,' ±  ',sige,'    Log10(D/Ro^2) = ',ord,' ±  ',siglogd
	print *,' '
	print *,'All the parameters are now set, so Oscar says relax...'	
	print *,' '
		
! Readjust (lower) the number of steps for use in the inversion routines	
! for those steps that are likely to have been recorded after sample breakdown
	nisample = ni
 	ni = nimax

! Get a table of progressive cumulative Dt(i)
	call zita(nisample,zi,e,ord,tilab,telab)
	
	ckchisq = 1.0e30
	ckchin = 1.0e30
	mct = 0
	iseed = xro
	guesscount = 0

! Guess a domain distribution
70	call guess(ndom,a1,a2,xro,iseed)
	
	guesscount = guesscount + 1
	
	if (mct.gt.30) then 
		print *, 'Warning: More than 30 fitting attempts. Adjusting parameters and retrying...'
		iquality = 2  ! flag that model struggled a bit
		if (ncicle.gt.0) ncicle = 4
		amax = 0.
		mct = 0
		chisq = 1.0e30
		goto 54
	endif
	
	nc = 0
	na = 2.*ndom -1
	  
	do j = 1,na
	   lista(j) = j
	end do

	mfit = na
	
	if (ncons.eq.1) mfit = na-1

	alamda = -1.
	kc = 0.
	ch = -1.
	alam = 0.001

! call LM fitting routine
26	call mrqmin(zi,f39,wf,ni,a2,na,lista,mfit,covar,alpha,nca,chisq,funcs,alamda,amax)

	do 52 j = 1,na,2
		if (a2(j+1).lt.-14) amax = -1.
	do 52 k = 1,na,2
		if (j.eq.k) goto 52
		if(a2(j).eq.a2(k)) amax = -1
52	continue

	if(amax.eq.-1.) then
		mct = mct + 1
	  	goto 70
	endif
	
	if (alam.gt.alamda) then
		nc = 0
	else
		nc = nc + 1
		if (nc.le.50) goto 38
		mct = mct + 1
		goto 70
	endif

	chisqn = chisq

	if (chisq.gt.1.) chisqn = 1.
	dchisq = abs((chisq - ch)/chisqn)
	kc = kc + 1
	
	if(dchisq.ge.dchmin.and.kc.le.100.or.kc.lt.5) goto 38
	
!84	write(28, '("    # domains = ",I2,"   Isteps = ",I2,"   nc = ",I2,"   chiq = ",G13.4)'),ndom,kc,nc,chisq
84	continue

	goto 54

72	alamda = 0.
	ndom = (na+1)/2

! Final wrapup call of LM routine (see Press et al.)
	call mrqmin(zi,f39,wf,ni,a2,na,lista,mfit,covar,alpha,nca,chisq,funcs,alamda,amax)
	              
	if (amax.eq.-1) stop 'stop 1: Consult your vodoo (says Oscar)'
	do nt = 1,nisample
		call funcs(zi(nt),a2,y,dyda,na,a1,amax)
		f(nt) = y
		if(amax.eq.-1)stop 'stop 3: Consult your vodoo (says Oscar)'
	end do

	call sort3(2*ndom,a1,a2,da,wksp,iwksp)

	rpmax = a1(na)
	xlog = ord - 2.* dlog10(rpmax)

! Calculate 10C/my closure temperature for each domains	

! loop through domains
	do j = 1,na+1,2   ! borrow data-writing code from above
		tcten(j) = 473.
		
		
! need if statement to toggle TC coefficient		
		
		
		if (geometry.eq.0) then		
			arg0 = (8.7 * 10**ord/a1(j)/a1(j)) * 1e6 * 365.25  *24. * 3600. / (e*1000./1.98588 * 10.)
		else
			arg0 = (55. * 10**ord/a1(j)/a1(j)) * 1e6 * 365.25  *24. * 3600. / (e*1000./1.98588 * 10.)
		end if
		
! do a wildly safe 20 iterations of Tc eqn			
		do jj = 1,20
			arg = arg0*tcten(j)*tcten(j)
			tcten(j) = e*1000/1.98588/log(arg)
		end do
!			write (*,*) j, tcten(j) - 273.15
	end do

! Write a report about the run
	write(14,'(I5)') ndom
! Report regressions results to report-XXX.out
	write(28,'("****** SAMPLE: ",A10," ******")') TRIM(samplename)
	if (geometry.eq.0) then
		write (28,*) '    Infinite-slab diffusion geometry'
	else
		write (28,*) '    Spherical diffusion geometry'	
	end if
	if (auxro.eq.0) then
		write (28,*) '    Maximum value for log(r/ro) free'	
	else
		write (28,'("     Maximum value for log(r/ro) fixed: ",F4.2)') auxro
	end if
	if (iquality.eq.1) then
		write (28,*) '    Mismatch in rpmax and max domain size - results should still be ok'	
	end if
	if (iquality.eq.2) then
		write (28,*) '    More than 30 fitting attempts - results should still be ok'	
	end if
	
	write (28,*) ' '
	write (28,*) '*** Arrhenius regression results: ***'
	if (option.eq.3) then
		write(28,'("    Using step: ",I2," (E forced through point)")') steplo
	else
		write(28,'("    Using steps ",I2," to ",I2,":")') steplo,stephi	
	end if

	if (wopt.eq.0) then
		if (option.eq.3) then
			write(28,'("    E = ",F5.2," forced - error could not be calculated")') e	
		else
			write(28,'("    E = ",F5.2," ± ",F5.2," kcal/mol")') e,sige		
		end if
	else
		write(28,'("    E = ",F5.2," error could not be calculated (weightings)")') e
	end if

	if (wopt.eq.0) then
		if (option.eq.3) then
			write(28,'("    Log10D/Ro^2 ",F7.3)') ord	
		else
			write(28,'("    Log10D/Ro^2 = ",F7.3," ± ",F5.3)') ord,siglogd		
		end if
	else
		write(28,'("    Log10D/Ro^2 ",F7.3," error could not be calculated (weightings)")') ord
	end if		

	select case (wopt)
	case (0)
		write(28,'("    Unweighted regression")')

	case (1)
		write(28,'("    Regression weighted by 1/(Fn - Fn-1).")')
	case (2)
		write(28,'("    Regression weighted by 1/logD.")')
	end select

	write(28,*)' '		
	write(28,*)'******** Domain Parameters for Arvert ********'
	write(28,*)'------------------------------------------------------------------------   ------------'
	write(28,'(1X,"Domain",4x,"  E  ",4x,"Log10(Do/a2)",4X,"Vol. Frac.",4x,"Size(Ro)",5x,"Size(Rmax)",3x,"Tc(10C/m.y.)")')
	write(28,*)'------------------------------------------------------------------------   ------------'
	sumc = 0.

! Do most of the plotting here (with a bit of reporting-to-file embedded)
	if (gmt.eq.1) then
		gmtfrag = "gmt "
! set some gmt parameters
		call system("gmt gmtset FONT_LABEL 14p") 
		call system("gmt gmtset FONT_ANNOT_PRIMARY 12p")
		call system("gmt gmtset MAP_GRID_PEN_PRIMARY 0.20,100/100/100,2_1:0p")
		call system("gmt gmtset PS_CHAR_ENCODING ISOLatin1+")
		call system("gmt gmtset MAP_TICK_LENGTH_PRIMARY 0.1i")
		call system("gmt gmtset MAP_TICK_PEN_PRIMARY 0.75p")
		call system("gmt gmtset MAP_FRAME_PEN 1p")
		call system("gmt gmtset MAP_ORIGIN_Y 0i")
		
! first arrhenius plot	
!   plot arrhenius axes
		stringfrag = "gmt psbasemap -Xc -Yr7.8i -R5.0/15.0/-25.0/-5.0 -Bf1a2g1:'10 000/K':/a5.0f1.0g5.0:'Ln(D/r@+2@+)':WSe "
		systemstring = TRIM(stringfrag)//" -JX3.5i/2.75i -K -P > "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))
	
!   draw top border for plot	   
		stringfrag = "gmt psbasemap -R5.0/15.0/-25.0/-5.0 -Bf0::/::n -JX3.5i/2.75i -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	   
		   
! series of plotting calls to add temperature ticks to top border
		open(unit=99,file='line.scr',status='UNKNOWN')
		write (99,*)'13.83 -5'
		write (99,*)'13.83 -5.5'
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -W0.5p,black -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))		
		   
		open(unit=99,file='line.scr',status='UNKNOWN')
		write (99,*)'11.45 -5'
		write (99,*)'11.45 -5.5'
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -W0.5p,black -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))
		
		open(unit=99,file='line.scr',status='UNKNOWN')
		write (99,*)'9.32 -5'
		write (99,*)'9.32 -5.5'
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -W0.5p,black -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))
		
		open(unit=99,file='line.scr',status='UNKNOWN')
		write (99,*)'7.45 -5'
		write (99,*)'7.45 -5.5'
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -W0.5p,black -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))
		
! annotate the temperature ticks	
		stringfrag = "echo 14.06 -4.75 450 \\260C | gmt pstext -D0/0.075 -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -N -F+f9,0+jCB  "
		systemstring = TRIM(stringfrag)//" -Wwhite  -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))		
	
		stringfrag = "echo 11.68 -4.75 600 \\260C | gmt pstext -D0/0.075 -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -N -F+f9,0+jCB  "
		systemstring = TRIM(stringfrag)//" -Wwhite  -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	
	
		stringfrag = "echo 9.58 -4.75 800 \\260C | gmt pstext -D0/0.075 -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -N -F+f9,0+jCB  "
		systemstring = TRIM(stringfrag)//" -Wwhite  -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	
	
		stringfrag = "echo 7.69 -4.75 1050 \\260C | gmt pstext -D0/0.075 -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -N -F+f9,0+jCB  "
		systemstring = TRIM(stringfrag)//" -Wwhite  -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	
	end if

! the following loop does some output work and we also use it for some plotting
! NOTE: array a1() holds both volume and size values for domains (odd indices = size, even indices = volume fraction)
	
	write(*,*)' '
	write(*,*)'******** CHECK rpmax = largest? ********'
	write(*,'(15X,F6.2,2X,F6.2)')rpmax,a1(2*ndom-1)
	if (abs(rpmax - a1(2*ndom-1)).ge.0.01) then
		write(*,*)'     >>> WARNING: large-domain mismatch (reported results should be ok)'
		iquality = 1  ! flag that non-fatal glitch has occurred
	end if

	do j = 1,na+1,2
		write(14,'(F8.4)') e
		ordj = xlog - 2.*dlog10(a1(j)/rpmax)
! version 1.30 workaround aimed at domains.in bug and also plotting of single-domain lines
!		write(14,'(F8.4)') ordj
		write(14,'(F8.4)') ord - dlog10(a1(j)*a1(j))
		write(14,'(F8.6)') a1(j+1)

		write(28,'(1X,I4,5X,F7.3,6X,F7.4,6X,f8.4,4x,f9.4,5X,f9.4,8X,f5.1)')(j+1)/2,E,ord - dlog10(a1(j)*a1(j)), & 
		a1(j+1),a1(j),a1(j)/a1(2*ndom-1),tcten(j)-273.15		

		sumc = sumc + a1(j+1)
		
! superimpose individual domain trends on arrhenius plot
		if (gmt.eq.1) then
			xpoint = 14.5
! version 1.30 workaround aimed at domains.in bug and also plotting of single-domain lines
!			ypoint = 2.302585*(ordj) - E/10000./R*xpoint			
			ypoint = 2.302585*(ord - dlog10(a1(j)*a1(j))) - E/10000./R*xpoint
			if (ypoint.le.-25.0) then
				ypoint = -24.9
!				xpoint = (2.302585*ordj - ypoint)/(E/10000./R)
! version 1.30 workaround aimed at domains.in bug and also plotting of single-domain lines
!				xpoint = (2.302585*ordj - ypoint)/(E/10000./R)	
				xpoint = (2.302585*(ord - dlog10(a1(j)*a1(j))) - ypoint)/(E/10000./R)				
			end if
			if (ypoint.ge.-5.0) then
				ypoint = -5.1
! version 1.30 workaround aimed at domains.in bug and also plotting of single-domain lines
!				xpoint = (2.302585*(ordj) - ypoint)/(E/10000./R)
				xpoint = (2.302585*(ord - dlog10(a1(j)*a1(j))) - ypoint)/(E/10000./R)
			end if			
			open(unit=99,file='line.scr',status='UNKNOWN')
			write (99,*)xpoint,ypoint
			xpoint = 6.0	
! version 1.30 workaround aimed at domains.in bug and also plotting of single-domain lines		
!			ypoint = 2.302585*ordj - E/10000./R*xpoint
			ypoint = 2.302585*(ord - dlog10(a1(j)*a1(j))) - E/10000./R*xpoint	
			if (ypoint.le.-25.0) then
				ypoint = -24.9			
! version 1.30 workaround aimed at domains.in bug and also plotting of single-domain lines					
!				xpoint = (2.302585*ordj - ypoint)/(E/10000./R)
				xpoint = (2.302585*(ord - dlog10(a1(j)*a1(j))) - ypoint)/(E/10000./R)				
			end if
			if (ypoint.ge.-5.0) then
				ypoint = -5.1
! version 1.30 workaround aimed at domains.in bug and also plotting of single-domain lines	
!				xpoint = (2.302585*ordj - ypoint)/(E/10000./R)
				xpoint = (2.302585*(ord - dlog10(a1(j)*a1(j))) - ypoint)/(E/10000./R)
			end if			
			write (99,*)xpoint,ypoint
			close(99)
			
			! plot each domain's arrhenius line
			stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -W1.5p,gray90 -O -K -P  >> "
			systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
			call system(TRIM(systemstring))			
		end if		
	end do	

! now plot observed data on top of base plot (we're keeping stacking of elements in mind)
	if (gmt.eq.1) then
	! first overlay the observed arrhenius data
		filename = TRIM(dirpath)//"arr-observed-"//TRIM(samplename)//".dat"	
		stringfrag = "gmt psxy "//TRIM(filename)//" -A -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -Sc0.1i -W0.4p,gray83 "
		systemstring = TRIM(stringfrag)//"  -Gblue -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	   

	! Finally, overlay the Ro reference line; remember we are plotting in ln() space, not log10
		xpoint1 = 14.5
		xpoint2 = 7.5
		ypoint1 = 2.302585*ord - E/10000./R*xpoint1
		ypoint2 = 2.302585*ord - E/10000./R*xpoint2	
		open(unit=99,file='line.scr',status='UNKNOWN')
		write (99,*)xpoint1,ypoint1
		write (99,*)xpoint2,ypoint2
		close(99)
		
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -W2p,gray30 -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))				

! next, the RRo plot
	! plot frame, right-side logRRo axis
		stringfrag = "gmt psbasemap -Xr0i -Yr-3.550i -R0/1.0/-1.0/3.0 -Bf0.05a0.2g0.1:'Fractional @+39@+Ar Loss':/f0.5a0.5:"
		morestring = TRIM(stringfrag)//"'log10(r/r@-o@-)':ES -JX3.5i/2.75i -O -K -P  >> "//TRIM(dirpath)//"kinetics-"
		systemstring = TRIM(morestring)//TRIM(samplename)//".ps"		
		call system(TRIM(systemstring))			

! Plot left hand size axis (log scaling)
		stringfrag = "gmt psbasemap -R0/1.0/0.1/1000. -B:/:/f2a1g2:'Size Relative to r@-o@-':W "
		systemstring = TRIM(stringfrag)//"  -JX3.5i/2.75il -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	

! Plot top of frame
		stringfrag = "gmt psbasemap -R0/1.0/0.1/1000. -Bf0::/::n  -JX3.5i/2.75il -O -K -P  >> "//TRIM(dirpath)//"kinetics-"
		systemstring = TRIM(stringfrag)//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	

! Plot (as underlay) depiction of domain distribution as "size spectrum"
! Plot and label high-temperature region where data weren't modeled
		open(unit=99,file='line.scr',status='UNKNOWN')
		write (99,*)f39(nbreakdown),2.15
		write (99,*)1.0,2.15
		close(99)
	! thick line to show region
		stringfrag = "gmt psxy line.scr -A -Yr0i -JX3.5i/2.75i -R0/1.0/-01.0/3.0 -W10p,230/215/215 -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))			
	! label for region	
		write(morestring,'(I4)') INT(maxtemp - 273.15)
		stringfrag = "echo 0.808 2.55 not modeled | gmt pstext -JX3.5i/2.75i -R0/1.0/-1.0/3.0 -N -F+f7,0+jLB "
		systemstring = TRIM(stringfrag)//" -W,255/255/255  -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))			
		stringfrag = "echo 0.808 2.40 above "//trim(morestring)//" \\260C | gmt pstext -JX3.5i/2.75i -R0/1.0/-1.0/3.0 -N -F+f7,0+jLB "
		systemstring = TRIM(stringfrag)//" -W,255/255/255  -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))
		
! Report the fit (chi-sqared) value, for what it's worth
		write(morestring,'(ES10.2)') ckchisq
		stringfrag = "echo 0.74 -0.56 fit = "//trim(morestring)//" | gmt pstext -JX3.5i/2.75i -R0/1.0/-1.0/3.0 -N -F+f7,0+jLB "
		systemstring = TRIM(stringfrag)//" -W,255/255/255 -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))		
		
! Plot (as underlay) depiction of domain distribution as "size spectrum"
		open(unit=99,file='line.scr',status='UNKNOWN')
		totalloss = 0.0
		
		do j = 1,na+1,2
			write (99,*) totalloss,dlog10(a1(j))   
			totalloss = totalloss + a1(j+1)
			write (99,*) totalloss,dlog10(a1(j))
		end do
		close(99)
		
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R0/1.0/-1.0/3.0 -W5p,200/200/200 -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))			
		
! Overlay the observed logRro data first

		filename = TRIM(dirpath)//"logr-observed-"//TRIM(samplename)//".dat"
		stringfrag = "gmt psxy "//filename//" -A -Yr0i -JX3.5i/2.75i -R0/1.0/-1.0/3.0 -W2.5p,blue "
		systemstring = TRIM(stringfrag)//" -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	

! Finally, add to the plot some text information about the sample and the run
! The text file needs to use gmt's arcane paragraph syntax for pstext with the -m option
! We double space to force breaks and white space. We need to use octal code /040 to write out code for spaces
! because pstex in paragraph mode ignores extra whitespace.
! We also have to go through some shenanigans to get chosen regression hi and lo steps saved, and passed back to main program

		slop = e*dlog10(ee)/10000./r

! Write the sample and run info to the source text file
		open(unit=99,file='text.scr',status='UNKNOWN')
		write(99,'("> 15.8 -5.0 12p 2i l")')
		
		write(99,'("SAMPLE: ",A10)') TRIM(samplename)
		write(99,'("")')
!		write(99,'("")')
		
		if (option.eq.3) then
			write(99,'("Using step ",I2,":")') steplo
		else
			write(99,'("Using steps ",I2," to ",I2,":")') steplo,stephi		
		end if
		write(99,'("")')
!		write(99,'("")')
		
		if (wopt.eq.0) then	
			if (option.eq.3) then
				write(99,'("E = ",F5.2," \261 no error - forced")') e
			else
				write(99,'("E = ",F5.2," \261 ",F5.2," kcal/mol")') e,sige		
			end if			
		else
			write(99,'("E = ",F5.2," \261 error N/A")') e
		end if
!		write(99,'("")')
		write(99,'("")')
		
		if (wopt.eq.0) then
			if (option.eq.3) then
				write(99,'("Log@-10@-(D/r@-o@-@+2@+) = ",F7.3)') ord
			else
				write(99,'("Log@-10@-(D/r@-o@-@+2@+) = ",F7.3," \261 ",F5.3)') ord,siglogd	
			end if	
		else
			write(99,'("Log@-10@-(D/r@-o@-@+2@+) = ",F7.3," \261 error N/A")') ord
		end if		
!		write(99,'("")')
		write(99,'("")')

		select case (wopt)
		case (0)
			write(99,'("Unweighted regression")')
!			write(99,'("")')
			write(99,'("")')
		case (1)
			write(99,'("Weighted by step size")')
!			write(99,'("")')
			write(99,'("")')
		end select

		select case (geometry)
		case (0)
			write(99,'("Infinite-slab geometry")')
!			write(99,'("")')
			write(99,'("")')
		case (1)
			write(99,'("Spherical geometry")')
			write(99,'("")')
!			write(99,'("")')		
		end select

		if (auxro.eq.0) then
				write(99,'("Maximum log(r/ro) free")')
				write(99,'("")')
		else
				write(99,'("Maximum log(r/ro) fixed: ",F4.2)') xro  !was auxro, shouldn't matter
				write(99,'("")')
		end if

		select case (iquality)
		case (0)
			write(99,'("No model-run issues")')
!			write(99,'("")')
			write(99,'("")')
		case (1)
			write(99,'("Mismatch in max domain size")')
			write(99,'("")')
!			write(99,'("")')
		case (2)
			write(99,'("More than 30 fitting attempts")')
			write(99,'("")')
!			write(99,'("")')	
		end select		
	
		write(99,'("\040")')
		write(99,'("")')
!		write(99,'("")')
		
		write(99,'(I2,2x"domains:")') ndom
		write(99,'("")')
!		write(99,'("")')
		
		write(99,'("N\040\040\040\040\040\040Vol.\040\040\040\040\040\040Vol@-T@-\040\040\040\040Rel. Size\040\040\040\040\040Tc10")')
		write(99,'("")')
!		write(99,'("")')
		
		totalloss = 0.0
		do j = 1,na+1,2   ! borrow data-writing code from above
		
			totalloss = totalloss + a1(j+1)
			if (a1(j)<10.0) then
				write (99,'(I2," \040\040\040 ",F5.3," \040\040\040 ",F4.2," \040\040\040\040\040 ",F5.3,&
					" \040\040\040\040\040\040\040 ",F4.0)')&
			        (j+1)/2,a1(j+1),totalloss,a1(j),tcten(j)-273.15
			else if (a1(j)<100.0) then
					write (99,'(I2," \040\040\040 ",F5.3," \040\040\040 ",F4.2," \040\040\040\040\040 ",F6.2,&
						" \040\040\040\040\040\040\040 ",F4.0)')&
						(j+1)/2,a1(j+1),totalloss,a1(j),tcten(j)-273.15			
				else
					write (99,'(I2," \040\040\040 ",F5.3," \040\040\040 ",F4.2," \040\040\040\040\040 ",F7.1,&
						" \040\040\040\040\040\040\040 ",F4.0)')&
						(j+1)/2,a1(j+1),totalloss,a1(j),tcten(j)-273.15			
			end if

		write(99,'("")')
!		write(99,'("")')
		end do
		close(99)
		
	! Now draw the sample info (have to shift back in Y since we want text next to Arrhenius plot)
		stringfrag = "gmt pstext text.scr -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -Yr+3.67i -M -F+f8+a0+jLT -N -O -K -P  >> "
		systemstring = TRIM(stringfrag)//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))				
	end if 


! get the model kinetic and RRo data
	call arr(f,tilab,telab,ni,e,ord,nisample,samplename,dirpath)
	      
	if (gmt.eq.1) then	    
		
	! Finally, overlay the chosen model's arrhenius data. Had to re-nudge this with -Y since data table was moved from 3.55 to 3.67
	! in previous plotting call.
		filename = TRIM(dirpath)//"arr-model-"//TRIM(samplename)//".dat"	
		stringfrag = "gmt psxy "//TRIM(filename)//" -A -JX3.5i/2.75i -R5.0/15.0/-25.0/-5.0 -Yr-0.12i -Sa0.15i -W0.4p,black "
		systemstring = TRIM(stringfrag)//" -Gorangered1 -O -K -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	
		
	! And as a last thing, overlay the chosen model's logRro data onto the RRo plot

		filename = TRIM(dirpath)//"logr-model-"//TRIM(samplename)//".dat"
		stringfrag = "gmt psxy "//filename//" -A -JX3.5i/2.75i -Yr-3.55i -R0/1.0/-1.0/3.0 -W2p,orangered1 "
		systemstring = TRIM(stringfrag)//" -O -P  >> "//TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		call system(TRIM(systemstring))	
		
	! clean up scratch files
		call system('rm text.scr')
		call system('rm line.scr')
	end if	
		
! try opening kinetics plot
	if (gmt.eq.1) then
		filename = TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".ps"
		systemstring = "gmt psconvert "//filename//" -Tf -Z "	
		call system(TRIM(systemstring))	
		filename = TRIM(dirpath)//"kinetics-"//TRIM(samplename)//".pdf"
		systemstring = 'open '//filename
		call system(TRIM(systemstring))	
	end if
	
	write(*,'("")')
	
	if (iquality.eq.0) then
		write(*,*) 'Model finished with no worries!'
	else
		write(*,*) 'Model finished but with some struggles.'  ! lots of attempts, or large-domain mismatch
	end if
	write(*,'("")')
	stop

54	continue

	if (ckchisq.gt.chisq) then
		do 64 j = 1,na+1
			auxa2(j) = a2(j)
		do 64 k = 1,na+1
			auxal(j,k) = alpha(j,k)
64		continue
		auxna = na
		ckchisq =  chisq
	endif

	if (ckchin.gt.chisq) then
		call funcs(zi(1),a2,y,dyda,na,auxa1,amax)
		if(amax.eq.-1)stop 'stop 2: Consult your vodoo'
		ckchin = chisq
	endif

	if (ncicle.lt.4)then
	    ncicle = ncicle + 1
	else
		call sort3(2*ndom,auxa1,a2,da,wksp,iwksp)
		ndom = ndom -1
		mct = 0
		ncicle = 0
		ckchin = 1.0e30
		sumc = 0.
		do j = 1,na,2
!			write(32,100)sumc,log(auxa1(j))
			sumc = sumc + auxa1(j+1)
!			write(32,*)sumc,log(auxa1(j))
		end do
!		  write(32,150)
		if(ndom.eq.mdom-1)then
			do 66 j = 1,auxna+1
				a2(j) = auxa2(j)
		    do 66 k = 1,auxna+1
		        alpha(j,k) = auxal(j,k)
66		    continue
			na = auxna
			mfit = na
			if(ncons.eq.1)mfit = na-1
		        
			write(*,'(" Number of domains = ",I2," with ",I4," calls to guess()")') (na +1)/2,guesscount
			amax = 0. 
		    goto 72
		endif
	endif

	goto 70
			
38	ch = chisq
	alam = alamda
	goto 26

end program domains ! main program

! $$$ *********** SUBROUTINE ARR

!	"PARAMETERS"

!	R= GAS CONSTANT KCAL/(MOL-K)
!	D0= FREQUENCY FACTOR (1/SEC) 

!	"INPUT"

!	NUSA=# OF SAMPLES 
!	NSAMP=# OF DIFFERENT DIFF. DOMAINs
!       E= ACTIVATION ENERGY (KCAL/MOL)
!	ORD = LOG (Doi/Ro**2)
!       C(J)= VOL. FRAC. OF Jth DOMAIN
!	RP(J)= PLATEAU SIZE (LOG(R/Ro) PLOT)
!	NI=# OF HEATING STEPS
!	TELAB = TEMP. STEP (K)
!	TILAB= STEP HEATING DURATION (MIN)

!	"OUTPUT"

!	TINV = INVERSE LAB. TEMP. (10000/K)
!	DZX = - LOG(D/R**2) (1/SEC)
!	F(J)*100 = CUMULATIVE % 39-AR RELEASED 
!	AVT = (F(J)+F(J-1))*50
!	XLOGR= LOG(R/Ro)
!	RAD(J) = SIZE OF THE Jth DIFF. DOMAIN
!       C(J) = VOL. FRAC. OF Jth DOMAIN


! ***********************************************************************
subroutine arr(f,tilab,telab,ni,e,ord,nisample,samplename,dirpath)
	use global
	implicit double precision (a-h,o-z)
	parameter(ns=200)
	dimension avt(ns)
	dimension telab(ni),f(0:ni),zx(ns),tilab(ni)
	character  tab1*9, mo*3, filename*70, samplename*10, dirpath*50
	
	tab1 = char(09)
	
	filename = TRIM(dirpath)//"arr-model-"//TRIM(samplename)//".dat"

	open(unit=16,file=filename,status='unknown')
	
	filename = TRIM(dirpath)//"logr-model-"//TRIM(samplename)//".dat"	
	open(unit=18,file=filename,status='unknown')	

	pi = 3.14159265
	ee = dexp(1.d00)
	R = 1.98588E-3

!	INVERSION OF 39-F
!    the following would nice to clean up but it works and I don't
!    want to spend the time to modernize it and then rigorously test it!

!    Also, to facilitate writing out full model, introduced a variable nisample
!    and passed it to several routines. This didn't seem to damage the inversion
!    which is still based on data from below 1150 C. 

	do k = 1,nisample
		select case (geometry)
		case (0)  ! slab
			if (f(k).gt.0.469) then
				zx(k+1) = -dlog(pi**2/8.*(1.-f(k)))*4./pi**2
			else
				zx(k+1) = pi/4.*(f(k))**2
			end if
		
		case (1)  ! sphere
			if (f(k).gt.0.8712) then
				zx(k+1) =  -dlog(pi**2/6.*(1.-f(k)))/pi**2
			else
				zx(k+1) = (2.-pi/3.*f(k)-2.*dsqrt(1.-pi/3.*f(k)))/pi
			end if
		end select
	end do

	zx(1) = 0.
	slop = E*dlog10(ee)/10000./R
!	write(*,*) 'model intercept slope',ord,slop
	do k = 1,nisample
		avt(k) = (f(k)+f(k-1))/2.*100.
		dzx = dlog10((zx(k+1)-zx(k))/tilab(k))  ! removed times 4.... why times four?? need to track this done if it solves the problem
		dzxplot = dlog((zx(k+1)-zx(k))/tilab(k))
		tinv = 1./telab(k)*10000.
		xlogr = (ord-slop*tinv-dzx)/2.  ! this hsould stay because we need to turn r2 to r
		write(18,'(1X,5(F12.8,A1))') f(k-1),tab1,xlogr
		write(18,'(1X,5(F12.8,A1))') f(k),tab1,xlogr
		write(16,'(1X,5(F12.8,A1))') tinv,tab1,dzxplot
	end do

2	continue

	return
end

! $$$ *********** SUBROUTINE DIFF
!   args in are telab, tilab, ni, and yes

subroutine diff(ord,E,f,telab,tilab,xlogr,xlogd,wt,ni,xro,eselect,steplo,stephi,samplename,dirpath,wopt,siglogd,sige)
	use global
	implicit double precision (a-h,o-z)
	parameter(ns=200)
	dimension f(0:ni),telab(ni),tilab(ni),xlogr(0:ni)
	dimension xlogd(0:ni),tinv(ns),wt(ni)
	
	double precision imp,siglogd,sige
	character tab1*9, eselect*1,samplename*10, filename*70, dirpath*50
	integer steplo,stephi,wopt
	
	tab1 = char(09)

	filename = TRIM(dirpath)//"arr-observed-"//TRIM(samplename)//".dat"
	open(unit=20,file=filename,status='unknown')
	filename = TRIM(dirpath)//"logr-observed-"//TRIM(samplename)//".dat"	
	open(unit=22,file=filename,status='unknown')
	filename = TRIM(dirpath)//"logr-"//TRIM(samplename)//".dat"	
	open(unit=24,file=filename,status='unknown')
	
	imp = 2.  !implicit declaration made this real*4 in original code
	b = 8.
	xlogr(0) = 0.
	pi = 3.141592654
	ee = dlog10(dexp(1.d00))  !ln() to log10() factor
	r = 1.98588E-3
	
!	CALCULATION OF LOG(D/R^2)   NOTE: could add better code here rather than approximations with stitch!
	
	do k = 1,ni       ! here xlogr() is cum Dt/a2 NOT LOG(R/Ro)!!!!
		select case (geometry)
		case (0)  ! slab
			if (f(k).gt.0.469) then
				xlogr(k) = -dlog(pi**2/8.*(1.-f(k)))*4./pi**2
			else
				xlogr(k) = pi/4.*(f(k))**2
			end if
		
		case (1)  ! sphere
			if (f(k).gt.0.8712) then
				xlogr(k) =  -dlog(pi**2/6.*(1.-f(k)))/pi**2
			else if (f(k).gt.1) then
					xlogr(k) = (2.-pi/3.*f(k)-2.*dsqrt(1.-pi/3.*f(k)))/pi			
				else
					xlogr(k) = (2.-pi/3.*f(k)-2.*dsqrt(1.-pi/3.*f(k)))/pi			
			end if
		end select
	end do
	if (geometry.eq.0) then
		write(20,*)'GEOMETRY: slab'
	else
		write(20,*)'GEOMETRY: sphere'
	end if
	
	sumwt = 0.
	nix = ni
	do k = 1,ni  !calculate observed arrhenius data and also prepare some weights for later fitting routines
		if(nix.eq.ni.and.telab(k).gt.1423) nix = k-1
		
		xlogd(k) = dlog((xlogr(k)-xlogr(k-1))/tilab(k))  ! report values in ln()
		tinv(k) = 1./telab(k)*10000.
		write(20,'(F7.4,A1,F8.4)')tinv(k),tab1,xlogd(k)
		xlogd(k) = dlog10((xlogr(k)-xlogr(k-1))/tilab(k))
		
		select case (wopt)
			case (0)	! no weighting but following autoarr we calculate them here
				wt(k) = 1./dsqrt(f(k)-f(k-1))			
			case (1)  ! legacy autoarr
	   			wt(k) = 1./dsqrt(f(k)-f(k-1))
	   		case (2)  ! NEW: weight by log D values to avoid regression bias to large values
	   			wt(k) = dsqrt(xlogd(k))
			case default  ! default to no weighting but as in case(0) we calculate legacy weights
				wt(k) = 1./dsqrt(f(k)-f(k-1))
		end select		
		sumwt = sumwt + wt(k)
		
	end do
	close(20)
	
	do k = 1,ni
		wt(k) = wt(k)/sumwt
	end do
	
!	CALCULATION OF E AND Do/Ro^2

	if(eselect.eq.'y') then
		write (*,*) 'Select an option for activation energy determination:'
		write (*,*) '   0 -- use the values chosen by this program'
		write (*,*) '   1 -- specify a range of heating steps to regress'	
		write (*,*) '   2 -- select from a table of values'
		write (*,*) '   3 -- force a value of E through a point'
		write (*,*) ' '		
		write (*,'(20a)',advance='no') ' Enter option: '
		read *,option
		write (*,*) ' '	
	else
		option = 0
	end if

	call param(ni,tinv,xlogd,wt,e,ord,option,steplo,stephi,samplename,dirpath,wopt,siglogd,sige)
	
	slop = e*ee/(r*10000)
	xro = (ord-slop*tinv(nix)-xlogd(nix))/2.*(1.+ (1.-f(nix))/2.)


! calculate and write out observed logR/Ro values (log 10)
! note how xlogr() is being reused here!

	xlogr(1) = (ord-slop*tinv(1)-xlogd(1))/2.
	write(22,'(F7.4,A1,F8.4)')f(0),tab1,xlogr(1)
	write(22,'(F7.4,A1,F8.4)')f(1),tab1,xlogr(1)
	
	do k = 2,ni  !calculate and write out observed logR/Ro values (log 10)
	   xlogr(k) = (ord-slop*tinv(k)-xlogd(k))/2.
	   write(24,'(F7.4,A1,F8.4)')f(k-1),tab1,xlogr(k)
	   write(22,'(F7.4,A1,F8.4)')f(k-1),tab1,xlogr(k)
	   write(22,'(F7.4,A1,F8.4)')f(k),tab1,xlogr(k)
	end do
	close(22)
	close(24)
	return
end

! $$$ *********** SUBROUTINE PARAM
!   this routine sniffs out the best E by regressing progressively more points
!   starting with the first three, and going up to nstop (set as 30 as parameter)
!   The best value is derived by detecting rollover and then backing off a bit,
!   looking at the series of fits.

!	This procedure has its benefits but makes it hard to graft in uncertainties in the recession
!   since this code uses the average of two E's and intercepts.

!   One could substitute, as an option, regression through specified steps
!   Or one could offer up a table of regressed intervals and let the user select
!   their preferred value of E.

subroutine param(ni,tinv,xlogd,wt,e,ord,option,steplo,stephi,samplename,dirpath,wopt,siglogd,sige)
	
	implicit double precision (a-h,o-z)
	
	parameter (ns=500, nstop = 30)
	dimension tinv(ni),xlogd(0:ni),wt(ni),y(ns),alog(ns),sigy(ns),siglogval(ns)
	character tab1*9, filename*70, dirpath*50,samplename*10
	integer option, kstart, variation(ns),varstart(ns),varend(ns),count

	double precision linefit(ns),eact(ns),ordo(ns),esig(ns),logdsig(ns),siglogd,sige,siga,sigb
	double precision temp,sigetemp,siglogdtemp
	integer itemp, model,i,ii,best,minsteps,steplo,stephi,mwt,wopt
	
	! mwt = 0, no weighting  = 1 weighting
	
	select case (wopt)
		case (0)   ! no weighting for linefit
			mwt = 0	
		case (1)
			mwt = 1  ! weighting, nature of weighting determined by wopt option specified by user		
		case default
			mwt = 0	  ! default to no weighting
	end select
	
	tab1 = char(09)
	filename = TRIM(dirpath)//"ener-"//TRIM(samplename)//".out"
!	open(unit=30,file=filename,status='unknown')
	
	ee = dlog10(dexp(1.d00))  !ln() to log10() factor
	r = 1.98588E-3

	select case (option)
		case (3)   ! user forces E through a point
			write (*,*) "You've chosen to force E through a point."
			write (*,'("    Enter E (kcal): ")',advance='no')
			read *,e
			write (*,'("    Enter point: ")',advance='no')
			read *,k
			write (*,*) ' '
			
			sige = 0.000
			slop = e / 10000. / r / 2.302585

			ord = slop * tinv(k) + xlogd(k)
			steplo = k
			stephi = k
			siglogd = 0.0000

!			write(30,'(F7.3)') e
!			write(30,'("User forced E through step: ",I3)') k
			return	
	
		case (1)   ! user enters range to regress
			write (*,*) "You've chosen to regress a single set of contiguous points."
			write (*,*) ' '
			write (*,'("    Enter step number, start of regression interval: ")',advance='no')
			read *,kstart
			write (*,'("                         end of regression interval: ")',advance='no')
			read *,k			
			write (*,*) ' '
			call fit(tinv,xlogd,kstart,k,wt,mwt,a,b,siga,sigb,chi2,q)
			e = -r*b*10000./ee
			sige = -sigb/b*e
			ord = a
			siglogd = siga
!			write(30,'(F7.3,A1,F7.3)') e,tab1,ord
!			write(30,'("User regressed steps: ",I3," to ",I3)') kstart,k
			steplo = kstart
			stephi = k
			return
			
		case (2)   ! user inspects table of choices and then chooses
			write (*,'("You''ve elected to choose E and the Ro reference from a table of models:")') 

			minsteps = -1
			do while ((minsteps.gt.ni).or.(minsteps.lt.3))
				write (*,'("    Enter minimum steps for regression: ")',advance='no')
				read *,minsteps
			end do			
			
			write (*,*) ' '
			
			y(2) = 0.
			
			nst = nstop
			if(ni.lt.nstop) then
				nst = ni
			endif
! following error should never happen 	

			if (nst.le.2) then
				write (*,*) '*** SERIOUS ERROR: less than three steps!!'
				write (*,*) '*** Program halted so you can figure out how you have failed again.'				
				STOP
			end if
			
			count = 0
			do kstart = 1,nst-(minsteps - 1)
				do k = kstart+(minsteps-1),nst
					call fit(tinv,xlogd,kstart,k,wt,mwt,a,b,siga,sigb,chi2,q)
					count = count + 1
!					y(count) = -r*b*10000./ee
!					alog(count) = a
					varstart(count) = kstart
					varend(count) = k		
					variation(count) = count
					linefit(count) = chi2
					eact(count) = -r*b*10000./ee
					ordo(count) = a
					esig(count) = sigb/b*eact(count)
					logdsig(count) = siga
!					a = a/ee
!					write(30,'(I2,A1,I2,A1,F7.3,A1,F6.3,A1,F7.3,A1,F6.3)')kstart,tab1,k,tab1,eact(count),tab1,-esig(count),tab1,&
!					a,tab1,logdsig(count)
				end do
			end do
!			close(30)

! use bubble sort to sort data according to maximum E
			do ii = 1,count-3
				do i = 1,count-1-ii
					if(eact(i+1).ge.eact(i)) then
						temp = eact(i)
						eact(i) = eact(i+1)
						eact(i+1) = temp
						temp = ordo(i)
						ordo(i) = ordo(i+1)
						ordo(i+1) = temp
						temp = esig(i)
						esig(i) = esig(i+1)
						esig(i+1) = temp
						temp = logdsig(i)
						logdsig(i) = logdsig(i+1)
						logdsig(i+1) = temp						
						temp = linefit(i)
						linefit(i) = linefit(i+1)
						linefit(i+1) = temp		
						itemp = variation(i)
						variation(i) = variation(i+1)
						variation(i+1) = itemp						
						itemp = varstart(i)
						varstart(i) = varstart(i+1)
						varstart(i+1) = itemp						
						itemp = varend(i)
						varend(i) = varend(i+1)
						varend(i+1) = itemp						
					end if
				end do
			end do

			best = 20
			write (*,*) 'Model  Start   End     Ea     Ea-err   Log10D/Ro^2  LineFit'
			write (*,*) '------------------------------------------------------------'
			write (*,*) '>>>>> BY ACTIVATION ENERGY:'
			do i = 1,best
				write(*,'(1X,I3,5X,I2,5X,I2,4X,F6.3,3X,F6.3,4X,F7.3,5X,F6.3)') variation(i),varstart(i),varend(i),&
				eact(i),-esig(i),ordo(i),linefit(i)
			end do
			
! use bubble sort to sort data according to minimum LineFit
			do ii = 1,count-3
				do i = 1,count-1-ii
					if(linefit(i+1).le.linefit(i)) then
						temp = eact(i)
						eact(i) = eact(i+1)
						eact(i+1) = temp
						temp = ordo(i)
						ordo(i) = ordo(i+1)
						ordo(i+1) = temp
						temp = esig(i)
						esig(i) = esig(i+1)
						esig(i+1) = temp
						temp = logdsig(i)
						logdsig(i) = logdsig(i+1)
						logdsig(i+1) = temp	
						temp = linefit(i)
						linefit(i) = linefit(i+1)
						linefit(i+1) = temp		
						itemp = variation(i)
						variation(i) = variation(i+1)
						variation(i+1) = itemp						
						itemp = varstart(i)
						varstart(i) = varstart(i+1)
						varstart(i+1) = itemp						
						itemp = varend(i)
						varend(i) = varend(i+1)
						varend(i+1) = itemp						
					end if
				end do
			end do	
			write (*,*) ''
			write (*,*) '>>>>> BY FIT (Ea > 30):'	
			
			i = 0
			j = 0
			do while (j.le.best)  ! we report 20 best fitting lines.. with eact over 30 
				i = i + 1
				if (eact(i).gt.30.) then
					write(*,'(1X,I3,5X,I2,5X,I2,4X,F6.3,3X,F6.3,4X,F7.3,5X,F6.3)') variation(i),varstart(i),varend(i),&
					eact(i),-esig(i),ordo(i),linefit(i)
					j = j + 1
				end if
			end do			
			write (*,*) '---------------------------------------','-------------------------------------------'

! also calculate auto selection, for reference or bailout:

			nst = nstop
			kmax = 2
			dymin = 100.
			y(2) = 0.
			
			if(ni.lt.nstop) then
				nst = ni
			endif
			
			kstart = 1
			do k = 3,nst
				call fit(tinv,xlogd,kstart,k,wt,mwt,a,b,siga,sigb,chi2,q)
				y(k)  =  -r*b*10000./ee
				sigy(k) = sigb/b*y(k)
				alog(k) = a
				siglogval(k) = siga
			end do
			
			do k = 3,nst
				dy = y(k+1)-y(k)
				ddy = y(k-1)-2.*y(k)+ y(k+1)
				if(y(k).gt.y(kmax)) kmax = k
				
				if(abs(dy).le.dymin.and.ddy.le.0.) then  ! look at change in E to detect rollover
					dymin = abs(dy)
					kmin = k
				endif
				
				if(dy.lt.0.) then   !count E's that drop... reset if E goes back up
					ndec = ndec + 1
				else
					ndec = 0.
				endif
				
				if(ndec.gt.4) then  ! if four E's dropped successively, declare rollover and be done
						! This averaging of values looks suspect (to PZ). Certainly the averaging of errors (by PZ) is a fugly
					etemp = (y(kmin) + y(kmin+1))/2.
					ordtemp = (alog(kmin) + alog(kmin+1))/2
					sigetemp = (sigy(kmin) + sigy(kmin+1))/2
					siglogdtemp = (siglogval(kmin) + siglogval(kmin+1))/2
					kmax = 0
				endif
			end do			
			
			write(*,'(3X,A1,5X,I2,5X,I2,4X,F6.3,3X,F6.3,4X,F7.3,5X,F6.3," < determined by program")') '0',kstart,kmin,etemp,&
			-sigetemp,ordtemp,chi2
			write (*,*) '---------------------------------------','-------------------------------------------'
			write (*,*) ' '			
			
			model = -1
			do while ((model.gt.count).or.(model.lt.0))
				write (*,'("Enter preferred model: ")',advance='no')
				read *,model
			end do
			write (*,*) ' '	
			if (model.gt.0) then
				! first have to get index of model in most recently sorted list
				do i = 1,count
					if(model.eq.variation(i)) then
						e = eact(i)
						ord = ordo(i)
						sige = -esig(i)
						siglogd = logdsig(i)
						steplo = varstart(i)
						stephi = varend(i)						
						exit
					end if
				end do
			else
				e = etemp
				ord = ordtemp
				sige = -sigetemp
				siglogd = siglogdtemp
				steplo = kstart
				stephi = kmin
			end if
			
		case default  ! usual autoarr values are used
			if (option.ne.0) then
				write (*,*) "You didn't choose a valid option. Defaulting to this program's  E and Ro reference."
			else
				write (*,*) "You've chosen to use program-determined values for E and the Ro reference."
			end if
			write (*,*) ' '
			nst = nstop
			kmax = 2
			dymin = 100.
			y(2) = 0.
			
			if(ni.lt.nstop) then
				nst = ni
			endif
			
			kstart = 1
			do k=3,nst
				call fit(tinv,xlogd,kstart,k,wt,mwt,a,b,siga,sigb,chi2,q)
				y(k)  =  -r*b*10000./ee
				sigy(k) = sigb/b*y(k)
				alog(k) = a
				siglogval(k) = siga
				write(30,'(I2,A1,F7.3,A1,F6.3,A1,F7.3,A1,F6.3)')k,tab1,y(k),tab1,sigy(k),tab1,a,siglogval(k)
			end do
			
			do k = 3,nst
				dy = y(k+1)-y(k)
				ddy = y(k-1)-2.*y(k)+ y(k+1)
				if(y(k).gt.y(kmax)) kmax=k
				
				if(abs(dy).le.dymin.and.ddy.le.0.) then  ! look at change in E to detect rollover
					dymin = abs(dy)
					kmin = k
				endif
				
				if(dy.lt.0.) then   !count E's that drop... reset if E goes back up
					ndec = ndec + 1
				else
					ndec = 0.
				endif
				
				if(ndec.gt.4) then  ! if four E's dropped successively, declare rollover and be done
					e = (y(kmin) + y(kmin+1))/2.
					ord = (alog(kmin) + alog(kmin+1))/2
					sige = -(sigy(kmin) + sigy(kmin+1))/2
					siglogd = (siglogval(kmin) + siglogval(kmin+1))/2
					kmax = 0
					close(30)
					steplo = kstart
					stephi = kmin
					return
				endif
			end do
			
			if(kmax.gt.0)then
				e = y(kmax)
				ord = alog(kmax)
				sige = -sigy(kmax)
				siglogd = siglogval(kmax)
				print *, ' '
				print *, 'Warning: did not find a real maximum for E.'
				print *, 'You should check the ener.out output file and calculate E manually if necessary.'
				print *, ' '
			endif
			close(30)
			
			steplo = kstart
			stephi = kmin
			
			return
	end select	
end


! $$$ *********** SUBROUTINES FIT, GSER, and GCF  - FUNCTIONS GAMMQ and GAMMLN

subroutine fit(x,y,kstart,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
	implicit double precision (a-h,o-z)
	dimension x(ndata),y(0:ndata),sig(ndata)
	double precision gamser,siga,sigb
	integer kstart
	
	sx = 0.
	sy = 0.
	st2 = 0.
	b = 0.
	
	if(mwt.ne.0) then
		ss = 0.
		do i = kstart,ndata
	  		wt = 1./(sig(i)**2)
	  		ss = ss+wt
	  		sx = sx+x(i)*wt
	  		sy = sy+y(i)*wt
		end do
	else
		do i = kstart,ndata
			sx = sx+x(i)
			sy = sy+y(i)
		end do
		
        ss = float(ndata-kstart+1)
    endif
    
	sxoss = sx/ss

	if (mwt.ne.0) then
		do i = kstart,ndata
			t = (x(i)-sxoss)/sig(i)
			st2 = st2+t*t
			b = b+t*y(i)/sig(i)
		end do
	else
		do i = kstart,ndata
			t = x(i)-sxoss
			st2 = st2+t*t
			b = b+t*y(i)
		end do
	endif
	
	b = b/st2
	a = (sy-sx*b)/ss
	siga = dsqrt((1.+sx*sx/(ss*st2))/ss)
	sigb = dsqrt(1./st2)
	chi2 = 0.
	
	if (mwt.eq.0) then
		do i = kstart,ndata
			chi2 = chi2+(y(i)-a-b*x(i))**2
		end do
		
		q = 1.
		sigdat = dsqrt(chi2/((ndata-kstart+1)-2))
		siga = siga*sigdat
		sigb = sigb*sigdat
	else
		do i = kstart,ndata
			chi2 = chi2+((y(i)-a-b*x(i))/sig(i))**2
		end do
		
		chi2 = chi2/(ndata-kstart - 1)  ! calculate a reduced chi square (PZ)
!		q = gammq(0.5*(ndata-2),0.5*chi2)
	endif
	
    return
end

! $$$ *********** GAMMQ
double precision function gammq(a,x)

	double precision a, x, gamser, gln, gammcf
	if(x.lt.0..or.a.le.0.) stop 'ERROR(GAMMQ): (x.lt.0..or.a.le.0.)'
	if(x.lt.a+1.)then
		call gser(gamser,a,x,gln)
		gammq = 1.-gamser
	else
		call gcf(gammcf,a,x,gln)
		gammq = gammcf
	endif
	return
end

! $$$ *********** GAMSER
subroutine gser(gamser,a,x,gln)

	implicit double precision (a-h,o-z)
	double precision gamser
	parameter (itmax = 100,eps=3.e-7)

	gln = gammln(a)
	
	if(x.le.0.)then
		if(x.lt.0.)stop 'ERROR(GSER): (x.lt.0.)'
		gamser = 0.
		return
	endif
	
	ap = a
	sum = 1./a
	del = sum
	
	do n = 1,itmax
		ap = ap+1.
		del = del*x/ap
		sum = sum+del
		if (abs(del).lt.abs(sum)*eps) goto 1
	end do
	
	stop 'ERROR(GSER): a too large, itmax too small'
	
1	gamser = sum*dexp(-x+a*dlog(x)-gln)
	return
end

! $$$ *********** GAMMCF
subroutine gcf(gammcf,a,x,gln)
         
	implicit double precision (a-h,o-z)
	double precision gammacf,a,x,gln
	parameter (itmax=100,eps=3.e-7)
	
	gln = gammln(a)
	gold = 0.
	a0 = 1.
	a1 = x
	b0 = 0.
	b1 = 1.
	fac = 1.
	
	do n = 1,itmax
		an = float(n)
		ana = an-a
		a0 = (a1+a0*ana)*fac
		b0 = (b1+b0*ana)*fac
		anf = an*fac
		a1 = x*a0+anf*a1
		b1 = x*b0+anf*b1
		
		if (a1.ne.0.)then
			fac = 1./a1
			g = b1*fac
			if(abs((g-gold)/g).lt.eps) goto 1
			gold = g
		endif
	end do
	
	stop 'ERROR(GCF): a too large, itmax too small'
	
1	gammcf = dexp(-x+a*dlog(x)-gln)*g
      return
end

! $$$ *********** GAMMLN
double precision function gammln(xx)

	double precision cof(6),stp,half,one,fpf,x,tmp,ser,xx
	data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,-1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      
	data half,one,fpf/0.5d0,1.0d0,5.5d0/

	x = xx-one
	tmp = x+fpf
	tmp = (x+half)*dlog(tmp)-tmp
	ser = one

	do j = 1,6
		x = x+one
		ser = ser+cof(j)/x
	end do
	
	gammln = tmp+dlog(stp*ser)
	return
end

! $$$ *********** FUNCS

Subroutine funcs(x,b,y,dyda,na,a,amax)
	use global
	implicit double precision (a-h,o-z)
         
    double precision x,b,y,dyda,a,amax
	parameter (nmax = 21)
	dimension a(na+1),dyda(na+1),b(na+1),csh(nmax)

!    the multiplication by 4 stand for the divition of Do ???? related to odd omission of 4x in arr-calcs??

	if(na.eq.0) return

	pi = 3.141592654
	y = 0. 
	as = 1.
	do j = 1,na,2
		if(b(j).lt.-50.)then
			a(j) = 0
		else
			a(j) = dexp(b(j))
		endif
		
		if(b(j+1).lt.-20.)then
			a(j+1) = 0.
			csh(j+1) = 0.
		else
			a(j+1) = (1. + dtanh(b(j+1)))/2.
			csh(j+1) = 0.5/dcosh(b(j+1))**2
		endif
		
		as = as - a(j+1)
	end do
	
	a(na+1)  = as + a(na+1)
	
	if(a(na+1).eq.0.)then
		amax = -1.
	  	return
	endif
	
	b(na+1) = dlog(a(na+1))

	do i = 1,na,2
		   arg = x/a(i)*4.  !CHANGED TEST --This fixes slab option (but why??)
		   ! see what happens with sphere  (still working, or need to make *4 only for slab?)
		   ! works with sphere too. So..... we havea working piece of code
		   ! but we don't know why it works correctly. Not the best of situations....
		   !
		   ! Lame idea at moment: 4.0 arrives from something in this routine, NOT due to
		   ! diffusion geometry (obviously, since this works for slabs and spheres).
		   !
		   ! Maybe the thing to do is ask Oscar!

		select case (geometry)
		case (0)   ! slab
			if (arg.le.0.173) then  
				gf = 2. * dsqrt(arg/pi)
			else
				if ((pi/2)**2*arg.gt.80) then
					gf = 1.
				else
					gf = 1. - 8./pi**2 * dexp(-(pi/2)**2*arg)
				endif
			endif
		
			dgf = 0 
			do j = 1,50000,2
				arg1 = (j*pi/2.)**2*arg
				if (arg1.gt.25.) goto 21
				dgf = dgf + 2. *  dexp(-arg1)
			end do
		case (1) ! sphere
			if (arg.le.0.1575) then  
				gf = 6. * dsqrt(arg/pi) - 3.* arg
			else
				if ((pi)**2*arg.gt.80) then
					gf = 1.
				else
					gf = 1. - 6./pi**2 * dexp(-(pi)**2*arg)
				endif
			endif
		
			dgf = 0 
			do j = 1,50000,1
				arg1 = (j*pi)**2*arg
				if (arg1.gt.25.) goto 21
				dgf = dgf + 6. *  dexp(-arg1)
			end do		
		end select
   		
21		y = y + a(i+1) * gf 
		dyda(i+1)= gf  * csh(i+1)
		dyda(i) = -a(i+1) * dgf * arg
		a(i) = dsqrt(a(i))
	end do
	return
end

! $$$ *********** GUESS 

subroutine guess(ndom,a1,a2,xro,iseed)

	implicit double precision (a-h,o-z)
	dimension a1(2*ndom),a2(2*ndom)

	na = 2*ndom
	sum = 0.
	do j = 1,na,2
		a1(j+1) =  1. + 10. * ran(iseed)
		sum = sum + a1(j+1)
	end do
	
	do j = 1,na,2
		a1(j+1) =  a1(j+1)/sum
	end do

	sum = 0 
	do j = 1,na-2,2
		sum =  1 + 10. * ran(iseed) + sum
		a1(na-2-j) = sum
	end do
	
	do j = 1,na-2,2
		a1(j) = a1(j)/sum * xro
	end do
	
	sum = a1(na)+a1(na-2)
	do j = 3,na-3,2
		a1(j) = a1(j) - dlog10(sum)
		sum = sum + a1(na-(j+1))
	end do

!	'ro' IS THE INVERSE OF Ro (1/Ro)
	
	ro = 10.**(a1(1))
	do j = 3,na-3,2
		a2(j) = a1(1)-a1(na-j)
	end do
	
	do j = 3,na-3,2
		a1(j) = a2(j)
	end do

	a1(na-1) = dlog10(a1(na))
	do j = 1,na-3,2
		a1(j) = a1(j+1)/(10.**a1(j)-10.**a1(j+2))
	end do
	
	a1(na-1) = 1.

	nloop = 0
	
29	continue
	ncont = 0
	nloop = nloop + 1
	
	do j = 1,na-3,2
		rom = 0.
		if (a1(j).gt.1.) stop 'a1(j) > 1.'
		if (a1(j+2).lt.a1(j)) then
			ncont = 0
			do k = j,j+2,2
				rom = rom + a1(k+1)/a1(k)
			end do
			
			a1(j+2) = a1(j)
			a1(j) = a1(j+1)/(rom - a1(j+3)/a1(j+2))
		else
			ncont = ncont + 1
		endif
	end do
	if (nloop.gt.30) stop 'nloop greater than 30 on guess' 
	if (ncont.lt.(ndom-1)) goto 29

	sumro = 0.
	do j = 1,na-1,2
	     sumro = sumro + a1(j+1)/a1(j)
	end do

!	CALCULATION OF A2

	do j = 1,na-1,2
	   a2(j) = 2.*dlog(a1(j)*ro)
	   z = 2. * a1(j+1) - 1.
	   a2(j+1) = 0.5 * dlog((z+1)/abs(z-1)) 
	end do

	return
end

! $$$ *********** RAN

Double precision function ran(iseed)

	parameter(ia=7141,ic=54773,im=259200)
	iseed = mod(iseed*ia+ic,im)
	ran = float(iseed)/float(im)
	return
end

! $$$ *********** INDEXX

subroutine indexx(n,arrin,indx)

	implicit double precision (a-h,o-z)
	dimension arrin(n),indx(n)
	
	do j = 1,n,2
        indx(j) = j
	end do
	
	l = n/4*2+ 1
	ir = n-1
	
10	continue

	if (l.gt.1) then
		l = l-2
		indxt = indx(l)
		q = arrin(indxt)
	else
		indxt = indx(ir)
		q = arrin(indxt)
		indx(ir) = indx(1)
		ir = ir-2
		
		if (ir.eq.1) then
			indx(1) = indxt
			return
		endif
	endif
	i = l
	j = l+l+1
	
20	if (j.le.ir) then
		if(j.lt.ir)then
			if(arrin(indx(j)).lt.arrin(indx(j+2)))j = j+2
		endif

		if (q.lt.arrin(indx(j))) then
			indx(i) = indx(j)
			i = j
			j = j+j+1
		else
			j = ir+2
		endif
		
        goto 20
	endif
	
	indx(i) = indxt
	goto 10
	
end

! $$$ *********** SUBROUTINES MRQMIN, MRQCOF, GAUSSJ, COVSRT
! $$$ *********** MRQMIN

subroutine mrqmin(x,y,sig,ndata,a,ma,lista,mfit,covar,alpha,nca,chisq,funcs,alamda,amax)
	use global
	implicit double precision (a-h,o-z)
!	parameter (mmax = 20)
parameter (mmax = 50)
	dimension x(0:ndata),y(0:ndata),sig(ndata),a(ma+1),lista(ma),covar(nca,nca),alpha(nca,nca),atry(mmax),beta(mmax),da(mmax)
	
	if (alamda.lt.0.) then
		amax = 0.
		do j = 1,ma,2
			if(a(j)*3..gt.amax)amax = a(j)*3.
		end do
	
		kk = mfit+1
		do j = 1,ma
	
			ihit = 0
			do k = 1,mfit
				if (lista(k).eq.j) ihit = ihit+1
			end do
		
			if (ihit.eq.0) then
				lista(kk) = j
				kk = kk+1
			else
				if (ihit.gt.1) then
					stop 'ERROR(MRQMIN): improper permutation in lista'
				end if
			endif
		end do
	
		if (kk.ne.(ma+1)) stop 'ERROR(MRQMIN): improper perm. in lista'
	
		alamda = 0.001
	
		call mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,nca,chisq,funcs,amax)
        
		if (amax.eq.-1) return
        ochisq = chisq
		do j = 1,ma
			atry(j) = a(j)
		end do
	endif
	
	do j = 1,mfit
		do k = 1,mfit
			covar(j,k) = alpha(j,k)
		end do
		
        covar(j,j) = alpha(j,j)*(1.+alamda)
        da(j) = beta(j)
	end do

	call gaussj(covar,mfit,nca,da,1,1,amax)
	
	if (amax.eq.-1.) return

	if (alamda.eq.0.) then
		call covsrt(covar,nca,ma,lista,mfit)
		return
	endif
	
21	sum = 0.

	do j = 1,mfit
		atry(lista(j)) = a(lista(j))+da(j)
	end do
	
	if (ma.ne.mfit) then
		do k = 1,mfit-1,2
			if(atry(k).ge.atry(ma).or.dabs(atry(k)).gt.amax) then
				da(k) = da(k)/2.
				goto 21
			endif
		end do
	else
		do j = 1,mfit,2
			if(dabs(atry(j)).gt.amax)then
				da(j) = da(j)/2.
				goto 21
			endif
		end do
	endif
	 
	do k = 1,mfit-1,2
		sum = sum + (1. + dtanh(atry(k+1)))/2.
	end do

	if (sum.ge.1.)then
		do k = 1,mfit,2
			da(k+1) = da(k+1)/2.
		end do
		
		goto 21
		
	endif
	
	call mrqcof(x,y,sig,ndata,atry,ma,lista,mfit,covar,da,nca,chisq,funcs,amax)
      
	if(amax.eq.-1) return

	if (chisq.lt.ochisq) then
		alamda = 0.1*alamda
		ochisq = chisq
		
		do j = 1,mfit
			do k = 1,mfit
				alpha(j,k) = covar(j,k)
			end do
			
			beta(j) = da(j)
			a(lista(j)) = atry(lista(j))
		end do
	else
		alamda = 10.*alamda
		chisq = ochisq
	endif
  
	return
      
end

! $$$ *********** MRQCOF
subroutine mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,nalp,chisq,funcs,amax)
	use global    
	implicit double precision (a-h,o-z)
!	parameter (mmax=20)
	parameter (mmax=50)	
	dimension x(0:ndata),y(0:ndata),sig(ndata),alpha(nalp,nalp),beta(ma),dyda(mmax),lista(mfit),a(ma+1),a1(mmax)
      
	do j = 1,mfit
		do k = 1,j
			alpha(j,k) = 0.
		end do
		
        beta(j) = 0.
	end do
	
      chisq = 0.
	do i = 1,ndata
	
		call funcs(x(i),a,ymod,dyda,ma,a1,amax)
		
		if (amax.eq.-1) then
			chisq = 1000000.
			return
		endif
		
			sig2i = 1./(sig(i)*sig(i))
			dy = y(i)-ymod
			
		do j = 1,mfit
			wt = dyda(lista(j))*sig2i
			
			do k = 1,j
				alpha(j,k) = alpha(j,k)+wt*dyda(lista(k))
			end do
			
		  	beta(j) = beta(j)+dy*wt
		end do
		chisq = chisq+dy*dy*sig2i
	end do
	
	do j = 2,mfit
		do k = 1,j-1
          alpha(k,j) = alpha(j,k)
		end do
	end do

	return

	end

! $$$ *********** GAUSSJ
subroutine gaussj(a,n,np,b,m,mp,amax)
	
	implicit double precision (a-h,o-z)
	parameter (nmax=50)
	dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
	do j = 1,n
		ipiv(j) = 0
	end do

	do i = 1,n
		big = 0.
        do j = 1,n
			if (ipiv(j).ne.1) then
				do k = 1,n
					if (ipiv(k).eq.0) then
						if (abs(a(j,k)).ge.big) then
							big = abs(a(j,k))
							irow = j
							icol = k
						endif
					else
						if (ipiv(k).gt.1) then
							amax  =  -1.
							return
						endif
					end if
				end do
			endif
		end do
		
		ipiv(icol) = ipiv(icol)+1
        
		if (irow.ne.icol) then	
			do l = 1,n
				dum = a(irow,l)
				a(irow,l) = a(icol,l)
				a(icol,l) = dum
			end do
			
			do l = 1,m
				dum = b(irow,l)
				b(irow,l) = b(icol,l)
				b(icol,l) = dum
			end do
        endif
        
		indxr(i) = irow
		indxc(i) = icol
        
		if (a(icol,icol).eq.0.) then
			amax = -1.
			return
		endif

		pivinv = 1./a(icol,icol)
		a(icol,icol) = 1.
        
		do l = 1,n
			a(icol,l) = a(icol,l)*pivinv
		end do

        do l = 1,m
          b(icol,l) = b(icol,l)*pivinv
		end do
		
		do ll = 1,n
			if (ll.ne.icol) then
				dum = a(ll,icol)
				a(ll,icol) = 0.
				do l = 1,n
					a(ll,l) = a(ll,l)-a(icol,l)*dum
				end do
			
				do l = 1,m
					b(ll,l) = b(ll,l)-b(icol,l)*dum
				end do
			endif
		end do

	end do

	do l = n,1,-1
	
		if (indxr(l).ne.indxc(l)) then
			do k = 1,n
				dum = a(k,indxr(l))
				a(k,indxr(l)) = a(k,indxc(l))
				a(k,indxc(l)) = dum
			end do
		endif
	end do
	
      return
      
end

! $$$ *********** COVSRT

subroutine covsrt(covar,ncvm,ma,lista,mfit)

	implicit double precision (a-h,o-z)
	dimension covar(ncvm,ncvm),lista(mfit)
      
	do j = 1,ma-1
		do i = j+1,ma
			covar(i,j) = 0.
		end do
	end do
	
	do i = 1,mfit-1
		do j = i+1,mfit
			if(lista(j).gt.lista(i)) then
				covar(lista(j),lista(i)) = covar(i,j)
			else
				covar(lista(i),lista(j)) = covar(i,j)
			endif
		end do
	end do
	
	swap = covar(1,1)
	
	do j = 1,ma
		covar(1,j) = covar(j,j)
		covar(j,j) = 0.
	end do

	covar(lista(1),lista(1)) = swap
	
	do j = 2,mfit
		covar(lista(j),lista(j)) = covar(1,j)
	end do
	
	do j = 2,ma
		do i = 1,j-1
			covar(i,j) = covar(j,i)
		end do
	end do

	return

end

! $$$ *********** SORT3

subroutine sort3(n,ra,rb,rc,wksp,iwksp)

	implicit double precision (a-h,o-z)
	dimension ra(n),rb(n),rc(n),wksp(n),iwksp(n)
	call indexx(n,ra,iwksp)

	do j = 1,n
		wksp(j) = ra(j)
	end do

	do j = 1,n,2
		ra(j) = wksp(iwksp(j))
		ra(j+1) = wksp(iwksp(j)+1)
	end do
	
	do j = 1,n
		wksp(j) = rb(j)
	end do
	
	do j = 1,n,2
		rb(j) = wksp(iwksp(j))
		rb(j+1) = wksp(iwksp(j)+1)
	end do

	do j = 1,n
		wksp(j) = rc(j)
	end do

	do j = 1,n,2
		rc(j) = wksp(iwksp(j))
		rc(j+1) = wksp(iwksp(j)+1)
	end do
	
	return
	
end

! $$$ *********** ZITA

subroutine zita(ni,zi,e,ord,tilab,telab)

	implicit double precision (a-h,o-z)
	dimension zi(0:ni)
	dimension telab(ni),tilab(ni)

	pi = 3.141592654d00
	R = 1.98588d-3
	zi(0) = 0.d00
	d0 = 10.d00**ord/4.d00

	do nt = 1,ni
		zi(nt) = d0*tilab(nt)*dexp(-E/R/telab(nt))+zi(nt-1)
	end do
	
	return

end
