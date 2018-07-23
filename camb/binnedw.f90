module binnedw
use precision
use constants
use ModelParams

      implicit none
      logical, parameter :: debugging = .false.
      real(dl), dimension(:),allocatable :: binned_a
      real(dl), dimension(:),allocatable :: binned_z, binned_w,binned_red, rhodeint, redint !output arrays of GP reconstruction
      real(dl), dimension(:),allocatable :: b1, c1, d1                   !coefficients for interpolation
      real(dl), dimension(:),allocatable :: b2, c2, d2                   !coefficients for interpolation
      real(dl)    :: multitheta !double theta function for binning
      integer  :: theta_bin=1, smooth_bin=2, GP=3
   
      !initializing global ODE solver parameters from CAMB
      real(dl), parameter :: initial_z = 0._dl
      real(dl), parameter :: initial_a = 1._dl
      real(dl) :: final_z
      real(dl) :: final_a
      integer, parameter  :: nsteps = 100
      real(dl) :: red_interval

      contains

      subroutine get_wofz(CP, a, wde)
      Type(CAMBparams) CP
      real(dl), intent(in)  :: a
      real(dl)              :: z
      real(dl), intent(out) :: wde
      integer               :: i,j,k

      real(dl), parameter   :: eps=1.e-12 !avoids 1/0

      if (a.gt.0._dl) then
         z = -1+1._dl/a
      else
         z = -1+1._dl/(a+eps)
      end if

      if (CP%mode.eq.theta_bin) then        
         if (z.ge.CP%zb(CP%nb)) then
            wde = CP%wb(CP%nb)
         else
            wde = CP%wb(1)
            do i=1,CP%nb-1
               multitheta = (sign(1d0,(z-CP%zb(i)))+1)/2 - (sign(1d0,(z-CP%zb(i+1)))+1)/2
               wde = wde + (CP%wb(i+1)-CP%wb(1))*multitheta
            end do 
         
         end if
         

      else if (CP%mode.eq.smooth_bin) then
         if (z.ge.CP%zb(CP%nb)) then
            wde = CP%wb(CP%nb)
         else
            wde = CP%wb(1)
            do i=1,CP%nb-1
               if (i.eq.1) then
                  wde = wde + (CP%wb(i+1)-CP%wb(i))/2 * (1+tanh( CP%s*(z-CP%zb(i))/((CP%zb(i))/2)  ) )
               else
                  wde = wde + (CP%wb(i+1)-CP%wb(i))/2 * (1+tanh( CP%s*(z-CP%zb(i))/((CP%zb(i)-CP%zb(i-1))/2)  ) )
               end if
            end do

         end if
      else if (CP%mode.eq.GP) then
         if ((z.ge.binned_z(1)).and.(z.le.binned_z(nsteps))) then
            wde = ispline(z, binned_z, binned_w, b1, c1, d1, nsteps)
         else
            wde = binned_w(nsteps)
         end if
      end if

      end subroutine get_wofz

      subroutine get_integral_rhode(CP)
      !This computes only the time dependent part
      !The initial value of rhode is added later
      Type(CAMBparams) CP
      real(dl)              :: wde, rhode0, integral, wplus, wminus
      integer               :: i,j,k


!Binning for integration------------------
      final_z=CP%zb(CP%nb)

      if (debugging) write(*,*) final_z

      do i=1,nsteps
         redint(i)=(i-1)*final_z/(nsteps-1)
      end do
      red_interval=final_z/(nsteps-1)     

      do i=1,nsteps-1
         binned_red(i)=(i+1-1)*final_z/(nsteps-1)
      end do

!-----------------------------------------

      do j=1, nsteps-1   
         call get_wofz(CP, 1._dl/(1+redint(j)), wminus)
         call get_wofz(CP, 1._dl/(1+redint(j+1)), wplus)
         integral = integral + 0.5*((1+wplus)/(1+redint(j+1))+(1+wminus)/(1+redint(j)))*(final_z/(nsteps-1))
         rhodeint(j) = exp(3._dl*integral)
      end do

      call newspline(binned_red,rhodeint, b2, c2, d2, nsteps-1)

      end subroutine get_integral_rhode

      subroutine get_rhode(a,rhode)
      !This gets the whole rhode, but it's used only for testing purposes
      real(dl), intent(in)  :: a
      real(dl), intent(out) :: rhode
      real(dl)              :: z,lastw, rhode0
      real(dl), parameter   :: eps=1.e-12 !avoids 1/0


      if (a.gt.0._dl) then
         z = -1+1._dl/a
      else
         z = -1+1._dl/(a+eps)
      end if
      if (z.le.redint(nsteps)) then
         rhode = ispline(z, binned_red, rhodeint, b2, c2, d2, nsteps-1)
      else
         call get_wofz(CP, 1/(1+redint(nsteps)), lastw)
         rhode = ((1+z)/(1+redint(nsteps)))**(3*(1+lastw))*ispline(redint(nsteps), binned_red, rhodeint, b2, c2, d2, nsteps-1)
      end if

      end subroutine get_rhode

      subroutine calc_w_de(CP)

      Type(CAMBparams) CP

      !Interface with GP python script
      character(LEN= 1000)                :: abins
      character(LEN= 1000)                :: wbins
      character(LEN= 1000)                :: steps_de
      character(LEN= 1000)                :: a_ini
      character(LEN= 1000)                :: wn
      character(LEN= 1000)                :: a_end
      character(LEN= 1000)                :: lencorr
      character(LEN= 20)                  :: feature_file="tmp_GPqz_000000.dat"
      character(LEN=10000)                :: command_plus_arguments
      real(dl), dimension(CP%nb)          :: gpa !ora applico un gp ad a per poi convertirla in redshift
      integer :: status
      integer :: getpid
      integer :: system
      integer :: i,m,nlbins
      real(dl) :: redshift, wdetest, rhodetest, omegam, omegade,rhode0



      !allocating arrays
      if (allocated(binned_z) .eqv. .false.) allocate (binned_a(nsteps),binned_z(nsteps),binned_w(nsteps),binned_red(nsteps-1), rhodeint(nsteps), redint(nsteps))
      if (allocated(b1) .eqv. .false.) allocate (b1(nsteps), c1(nsteps), d1(nsteps))
      if (allocated(b2) .eqv. .false.) allocate (b2(nsteps), c2(nsteps), d2(nsteps))


      nlbins=(CP%nb)-1

      if (debugging) then
         if ((CP%mode.eq.theta_bin).or.(CP%mode.eq.smooth_bin)) then
            write(*,*) 'num_bins=',CP%nb
            do i=1,CP%nb
               write(*,*) 'redshift',i,'=',CP%zb(i)
               write(*,*) 'wde',i,'=',CP%wb(i)
            end do
         end if
      end if

!provo la ricostruzione in a(z), poi otterrò comunque dei binnedz ma meno fitti: binnedz = -1 + 1/binneda
!così posso usare la stessa funzione di correlazione della prior e non ho problemi di conversione in lunghezza di corr in redshift



      !Gaussian process interface
      if (CP%mode.eq.GP) then

         !Setting GP redshift to median redshift of each bin
         gpa(1) = (initial_a+CP%ab(1))/2
         do i=2,CP%nb
            gpa(i) = (CP%ab(i-1)+CP%ab(i))/2.
         end do
         !write(*,*) 'scale factors di input', CP%ab
         !write(*,*) 'scale factors mediani', gpa

      final_a   = gpa(CP%nb)

         !Creating command line

         !Generate tmp file name based on PID
!         ipid = getpid()
         write (feature_file(11:16), "(Z6.6)") getpid()
         !1. Prepare command and launch it!
         write(a_ini, "(E15.7)"      ) initial_a
         write(a_end, "(E15.7)"      ) final_a
         write(steps_de, "(I10)"     ) nsteps
         write(abins, "(10E15.7)"   ) (gpa(i),i=1,CP%nb)
         write(wbins, "(10f15.7)"     ) (CP%wb(i),i=1,CP%nb) !python parser struggles with scientific notation negatives: using floats here
         write(wn, "(10f15.7)"     ) CP%w0
         write(lencorr, "(10E15.7)"  ) CP%corrlen

         !write(*,*) 'scale factors mediani', abins

         if (CP%mode.eq.GP) then
            if (debugging) write(*,*) 'WORKING WITH GP'
            !here needs the call to script with no baseline

            if (debugging) write(*,*) 'WORKING WITH GP (with baseline)'
        
            command_plus_arguments = "python GP.py --inia "//trim(adjustl(a_ini))//&
            &" --enda "//trim(adjustl(a_end))//" --ODEsteps "//trim(adjustl(steps_de))// &
            & " --scalefactors "//trim(adjustl(abins))// " --eos "//trim(adjustl(wbins))//&
	    & " --eosn "//trim(adjustl(wn))//&
            & " --l "//trim(adjustl(lencorr))// " --outfile " // feature_file


            !calling script!!!
            if (debugging) then
               write(*,*) 'Calling Gaussian process script with command line:'
               write(*,*) trim(adjustl(command_plus_arguments))
            end if
            status = system(trim(adjustl(command_plus_arguments)))
            if (status/=0) then
               print *, "Failure in GP reconstruction of w(z) -- see above."
               call abort
            end if


         end if


        
         !Reading temporary file generated by GP script--------------
         open(unit=17, file=feature_file, action='read')
         do i=1,nsteps
            read(17, "(E15.8, 1X, E15.8)", iostat=status) binned_z(i), binned_w(i)
            if (status>0) then
               print *, "Error reading the tmp w(z) file."
               call abort
            end if
         end do
         if (debugging) then
            close(17)
         else
            close(17, status='delete')
         end if
         !-----------------------------------------------------------

         !Setting interpolation for GP arrays------------------------
         call newspline(binned_z,binned_w, b1, c1, d1, nsteps)
         !-----------------------------------------------------------

     else if (CP%mode.gt.3) then
         write(*,*) "THIS MODEL DOESN'T EXIST!!!!"
         stop
     end if

     if (debugging) write(*,*) 'w(z) computed'

     call get_integral_rhode(CP)

     if (debugging) write(*,*) 'done rho_de integral'

     if (debugging) then
         write(*,*) 'printing w(z)'
         open(40,file='printwde.dat')
         open(42,file='printomega.dat')
         do m=1,1001
            redshift=(m-1)*3._dl/1000
            call get_wofz(CP,1/(1+redshift), wdetest)
            call get_rhode(1/(1+redshift), rhodetest)
            rhode0=3._dl*((1000*CP%H0/c)**2.)*CP%omegav
            rhodetest = rhode0*rhodetest
            omegam = ((3*(1000*CP%H0/c)**2.*(1-CP%omegav)*(1+redshift)**3._dl)/(rhodetest+3*(1000*CP%H0/c)**2.*(1-CP%omegav)*(1+redshift)**3._dl))
            omegade = (rhodetest/(rhodetest+3*(1000*CP%H0/c)**2.*(1-CP%omegav)*(1+redshift)**3._dl))
            write(40,*) redshift, wdetest
            write(42,*) redshift, omegam, omegade 
         end do
         close(40)
         close(42)
         !stop
     end if

     end subroutine calc_w_de

!INTERPOLATION ROUTINES-----------------------------------------------
   subroutine newspline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
real(dl) x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
real(dl) h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine newspline

  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
  

end module binnedw

