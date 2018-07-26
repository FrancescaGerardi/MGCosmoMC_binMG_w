module binnedMG
use precision
use constants
use ModelParams

      implicit none
      logical, parameter :: debugging = .false.
      real(dl), dimension(:),allocatable :: derivm, derivs, derivg
      real(dl), dimension(:),allocatable :: binned_z, binned_mu, binned_sigma  !output arrays of GP reconstruction
      real(dl), dimension(:),allocatable :: g1, c1, d1                         !coefficients for interpolation
      real(dl), dimension(:),allocatable :: g2, c2, d2                         !coefficients for interpolation
      real(dl), dimension(:),allocatable :: g3, c3, d3                         !coefficients for interpolation
      real(dl), dimension(:),allocatable :: g4, c4, d4                         !coefficients for interpolation
      real(dl), dimension(:),allocatable :: g5, c5, d5                         !coefficients for interpolation
      real(dl)    :: multitheta !double theta function for binning
      integer  :: theta_bin=1, smooth_bin=2, GP=3
   
      real(dl), parameter :: initial_a = 1._dl
      real(dl) :: final_a
      integer, parameter  :: nsteps = 100

      real(dl), dimension(:),allocatable :: red !array to compute derivatives
      real(dl) :: red_interval
      real(dl) :: final_z
      integer :: nsteps_derivs = 200      


      contains
!Mu---------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_mu(CP, a, mu)
      Type(CAMBparams) CP
      real(dl), intent(in)  :: a
      real(dl)              :: z
      real(dl), intent(out) :: mu
      integer               :: i,j,k

      real(dl), parameter   :: eps=1.e-12 !avoids 1/0

      if (a.gt.0._dl) then
         z = -1+1._dl/a
      else
         z = -1+1._dl/(a+eps)
      end if

      if (CP%modemg.eq.theta_bin) then   
         write(*,*) 'STEP FUNCTION IS A BAD CHOICHE FOR DERIVATIVES! CHANGE IT' 
         stop    
         if (z.ge.CP%zbmg(CP%nbmg)) then
            mu = CP%mb(CP%nbmg)
         else
            mu = CP%mb(1)
            do i=1,CP%nbmg-1
               multitheta = (sign(1d0,(z-CP%zbmg(i)))+1)/2 - (sign(1d0,(z-CP%zbmg(i+1)))+1)/2
               mu = mu + (CP%mb(i+1)-CP%mb(1))*multitheta
            end do 
         
         end if

!OBS: the binned values of MG functions are given with respect to an inverse order of a, so the order in terms of the redshift have not to be inverted         
      else if (CP%modemg.eq.smooth_bin) then
            mu = CP%mb(1)
            do i=1,CP%nbmg-1
                  mu = mu + (CP%mb(i+1)-CP%mb(i))/2 * (1+tanh(CP%ms*(z-CP%zbmg(i+1))/(CP%zbmg(i+1)-CP%zbmg(i))) )
            end do


      else if (CP%modemg.eq.GP) then
         if ((z.ge.binned_z(1)).and.(z.le.binned_z(nsteps))) then
            mu = ispline(z, binned_z, binned_mu, g1, c1, d1, nsteps)
         else
            mu = binned_mu(nsteps)
         end if
!          mu = 0
      end if

      end subroutine get_mu

      subroutine get_mu_derivs(CP)
      Type(CAMBparams) CP
      real(dl)              :: muplus, muminus
      integer               :: i,j,k


!Binning for derivatives------------------
         final_z=CP%zbmg(CP%nbmg)

         do i=1,nsteps_derivs
            red(i)=(i-1)*final_z/(nsteps_derivs-1)
         end do
         red_interval=final_z/(nsteps_derivs-1)     

!-----------------------------------------

         call get_mu(CP, 1/(1+red(2)), muplus)
         call get_mu(CP, 1/(1+red(1)), muminus)
         derivm(1) = (muplus - muminus)/red_interval
         write(*,*) red(1), muminus, red(2), muplus
         do j=2, nsteps_derivs-1
               call get_mu(CP, 1/(1+red(j+1)), muplus)
               call get_mu(CP, 1/(1+red(j-1)), muminus)
               derivm(j) = (muplus - muminus)/(2*red_interval)
         end do 
         call get_mu(CP, 1/(1+red(nsteps_derivs)), muplus)
         call get_mu(CP, 1/(1+red(nsteps_derivs-1)), muminus)
         derivm(nsteps_derivs) = (muplus - muminus)/red_interval

!         if (debugging) write(*,*) derivm

         call newspline(red, derivm, g2, c2, d2, nsteps_derivs)

 
      end subroutine get_mu_derivs


      subroutine get_dotmu(a,dotmu)
      real(dl), intent(in)  :: a
      real(dl), intent(out) :: dotmu
      real(dl)              :: z
      real(dl), parameter   :: eps=1.e-12 !avoids 1/0
      real(dl), dimension(CP%nbmg+1) :: inter_red, inter_der
      integer                  :: i

      if (a.gt.0._dl) then
         z = -1+1._dl/a
      else
         z = -1+1._dl/(a+eps)
      end if


      if (CP%modemg .eq. theta_bin) then 
         write(*,*) 'STEP FUNCTION IS A BAD CHOICHE FOR DERIVATIVES! CHANGE IT'
      else if (CP%modemg .eq. smooth_bin) then
!computing derivatives at middle redshifts 
         inter_red(1) = CP%zbmg(1)/2
         inter_der(1) = 0               
         do i=1,CP%nbmg-1
           inter_red(i+1) = (CP%zbmg(i+1)+CP%zbmg(i))/2
           inter_der(i+1) = (CP%mb(i+1)-CP%mb(i))/(CP%zbmg(i+1)-CP%zbmg(i))
         end do  
         inter_red(CP%nbmg+1) = CP%zbmg(CP%nbmg)
         inter_der(CP%nbmg+1) = 0   
!computive the partial derivative respect to the redshift
           dotmu = inter_der(1)
           do i=1,CP%nbmg
                dotmu = dotmu + (inter_der(i+1)-inter_der(i))/2 * (1+tanh( CP%ms*(z-inter_red(i+1))/(inter_red(i+1)-inter_red(i))) )
!in equations_ppf.f90 computing the partial derivative respect to the scale factor, which is the one required
           end do      
      else
         if (z.le.red(nsteps_derivs)) then
            dotmu = ispline(z, red, derivm, g2, c2, d2, nsteps_derivs)
         else
            dotmu = 0
         end if
      end if

      end subroutine get_dotmu

      

!Sigma------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_sigma(CP, a, sigma)
      Type(CAMBparams) CP
      real(dl), intent(in)  :: a
      real(dl)              :: z
      real(dl), intent(out) :: sigma
      integer               :: i,j,k

      real(dl), parameter   :: eps=1.e-12 !avoids 1/0


      if (a.gt.0._dl) then
         z = -1+1._dl/a
      else
         z = -1+1._dl/(a+eps)
      end if

      if (CP%modemg.eq.theta_bin) then   
         write(*,*) 'STEP FUNCTION IS A BAD CHOICHE FOR DERIVATIVES! CHANGE IT'     
         if (z.ge.CP%zbmg(CP%nbmg)) then
            sigma = CP%sb(CP%nbmg)
         else
            sigma = CP%sb(1)
            do i=1,CP%nbmg-1
               multitheta = (sign(1d0,(z-CP%zbmg(i)))+1)/2 - (sign(1d0,(z-CP%zbmg(i+1)))+1)/2
               sigma = sigma + (CP%sb(i+1)-CP%sb(1))*multitheta
            end do 
         
         end if
         
      else if (CP%modemg.eq.smooth_bin) then
            sigma = CP%sb(1)
            do i=1,CP%nbmg-1
                  sigma = sigma + (CP%sb(i+1)-CP%sb(i))/2 * (1+tanh(CP%ss*(z-CP%zbmg(i+1))/(CP%zbmg(i+1)-CP%zbmg(i))) )
            end do

      else if (CP%modemg.eq.GP) then
         if ((z.ge.binned_z(1)).and.(z.le.binned_z(nsteps))) then
            sigma = ispline(z, binned_z, binned_sigma, g4, c4, d4, nsteps)
         else
            sigma = binned_sigma(nsteps)
         end if
      end if

      end subroutine get_sigma


      subroutine get_sigma_derivs(CP)
      Type(CAMBparams) CP
      real(dl)              :: sigplus, sigminus
      integer               :: i,j,k


!Binning for derivatives------------------
         final_z=CP%zbmg(CP%nbmg)

         do i=1,nsteps_derivs
            red(i)=(i-1)*final_z/(nsteps_derivs-1)
         end do
         red_interval=final_z/(nsteps_derivs-1)     
!-----------------------------------------

         call get_sigma(CP, 1/(1+red(2)), sigplus)
         call get_sigma(CP, 1/(1+red(1)), sigminus)
         derivs(1) = (sigplus - sigminus)/red_interval        
         do j=2, nsteps_derivs-1
               call get_sigma(CP, 1/(1+red(j+1)), sigplus)
               call get_sigma(CP, 1/(1+red(j-1)), sigminus)
               derivs(j) = (sigplus - sigminus)/(2*red_interval)
         end do
         call get_sigma(CP, 1/(1+red(nsteps_derivs)), sigplus)
         call get_sigma(CP, 1/(1+red(nsteps_derivs-1)), sigminus)
         derivs(nsteps_derivs) = (sigplus - sigminus)/red_interval

         call newspline(red, derivs, g3, c3, d3, nsteps_derivs)

      end subroutine get_sigma_derivs


      subroutine get_dotsigma(a,dotsigma)
      real(dl), intent(in)  :: a
      real(dl), intent(out) :: dotsigma
      real(dl)              :: z
      real(dl), parameter   :: eps=1.e-12 !avoids 1/0
      integer               :: i
      real(dl), dimension(CP%nbmg+1)    :: inter_red, inter_der

      if (a.gt.0._dl) then
        z = -1+1._dl/a
      else
         z = -1+1._dl/(a+eps)
      end if

      if (CP%modemg .eq. theta_bin) then 
         write(*,*) 'STEP FUNCTION IS A BAD CHOICHE FOR DERIVATIVES! CHANGE IT'
      else if (CP%modemg .eq. smooth_bin) then
!computing derivatives at middle redshifts 
         inter_red(1) = CP%zbmg(1)/2
         inter_der(1) = 0               
         do i=1,CP%nbmg-1
           inter_red(i+1) = (CP%zbmg(i+1)+CP%zbmg(i))/2
           inter_der(i+1) = (CP%sb(i+1)-CP%sb(i))/(CP%zbmg(i+1)-CP%zbmg(i))
         end do  
         inter_red(CP%nbmg+1) = CP%zbmg(CP%nbmg)
         inter_der(CP%nbmg+1) = 0   
!computive the partial derivative respect to the redshift
           dotsigma = inter_der(1)
           do i=1,CP%nbmg
                dotsigma = dotsigma + (inter_der(i+1)-inter_der(i))/2 * (1+tanh( CP%ss*(z-inter_red(i+1))/(inter_red(i+1)-inter_red(i))  ) )
           end do    
!in equations_ppf.f90 computing the partial derivative respect to the scale factor, which is the one required 
      else
         if (z.le.red(nsteps_derivs)) then
            dotsigma = ispline(z, red, derivs, g3, c3, d3, nsteps_derivs)
         else
            dotsigma = 0
         end if
      end if

      end subroutine get_dotsigma


!GAMMA--------------------------------------------------------------------------------------------------------

      subroutine get_gamma(CP, a, gamm)
      Type(CAMBparams) CP
      real(dl), intent(in)  :: a
      real(dl), intent(out) :: gamm
      integer               :: i,j,k
      real(dl)              :: mu, sigma

      call get_mu(CP, a, mu)
      call get_sigma(CP, a, sigma)

      gamm =  (2._dl*sigma/mu)-1._dl 

      end subroutine get_gamma
      


      subroutine get_gamma_derivs(CP)
      Type(CAMBparams) CP
      integer               :: i,j
      real(dl)              :: gamplus,gamminus
     
!Binning for derivatives------------------
         final_z=CP%zbmg(CP%nbmg)

         do i=1,nsteps_derivs
            red(i)=(i-1)*final_z/(nsteps_derivs-1)
         end do
         red_interval=final_z/(nsteps_derivs-1)     
!-----------------------------------------

         call get_gamma(CP, 1/(1+red(2)), gamplus)
         call get_gamma(CP, 1/(1+red(1)), gamminus)
         derivg(1) = (gamplus - gamminus)/red_interval        
         do j=2, nsteps_derivs-1
               call get_gamma(CP, 1/(1+red(j+1)), gamplus)
               call get_gamma(CP, 1/(1+red(j-1)), gamminus)
               derivg(j) = (gamplus - gamminus)/(2*red_interval)
         end do
         call get_gamma(CP, 1/(1+red(nsteps_derivs)), gamplus)
         call get_gamma(CP, 1/(1+red(nsteps_derivs-1)), gamminus)
         derivg(nsteps_derivs) = (gamplus - gamminus)/red_interval

         call newspline(red, derivg, g5, c5, d5, nsteps_derivs)

      end subroutine get_gamma_derivs  

      subroutine get_dotgamma(a,dotgamma)
      real(dl), intent(in)  :: a
      real(dl), intent(out) :: dotgamma
      real(dl)              :: z
      real(dl), parameter   :: eps=1.e-12 !avoids 1/0
      real(dl)              :: gammminus, gammplus
      integer               :: i
      real(dl), dimension(CP%nbmg+1)    :: inter_red, inter_der

      if (a.gt.0._dl) then
        z = -1+1._dl/a
      else
         z = -1+1._dl/(a+eps)
      end if

      if (CP%modemg .eq. theta_bin) then 
         write(*,*) 'STEP FUNCTION IS A BAD CHOICHE FOR DERIVATIVES! CHANGE IT'
      else if (CP%modemg .eq. smooth_bin) then
!computing derivatives at middle redshifts 
         inter_red(1) = CP%zbmg(1)/2
         inter_der(1) = 0               
         do i=1,CP%nbmg-1
           inter_red(i+1) = (CP%zbmg(i+1)+CP%zbmg(i))/2
           call get_gamma(CP, 1/(1+CP%zbmg(i)), gammminus)
           call get_gamma(CP, 1/(1+CP%zbmg(i+1)), gammplus)
           inter_der(i+1) = (gammplus-gammminus)/(CP%zbmg(i+1)-CP%zbmg(i))
         end do  
         inter_red(CP%nbmg+1) = CP%zbmg(CP%nbmg)
         inter_der(CP%nbmg+1) = 0   
!computive the partial derivative respect to the redshift
           dotgamma = inter_der(1)
           do i=1,CP%nbmg
                dotgamma = dotgamma + (inter_der(i+1)-inter_der(i))/2 * (1+tanh( CP%ss*(z-inter_red(i+1))/(inter_red(i+1)-inter_red(i))  ) )
           end do    
!in equations_ppf.f90 computing the partial derivative respect to the scale factor, which is the one required 
      else
         if (z.le.red(nsteps_derivs)) then
            dotgamma = ispline(z, red, derivg, g5, c5, d5, nsteps_derivs)
         else
            dotgamma = 0
         end if
      end if

      end subroutine get_dotgamma      


!-----------------------------------------------------------------------------------------------------------------

      subroutine calc_MGfunc(CP)

      Type(CAMBparams) CP

      !Interface with GP python script
      character(LEN= 1000)                :: abins
      character(LEN= 1000)                :: mbins, sbins
      character(LEN= 1000)                :: steps_de
      character(LEN= 1000)                :: a_ini
      character(LEN= 1000)                :: mn, sn
      character(LEN= 1000)                :: a_end
      character(LEN= 1000)                :: mcorr, scorr
      character(LEN= 20)                  :: feature_file_mu="tmp_GPmu_000000.dat"
      character(LEN= 20)                  :: feature_file_sig="tmp_GPsig_000000.dat"
      character(LEN=10000)                :: command_plus_arguments
      real(dl), dimension(CP%nbmg)          :: gpa !ora applico un gp ad a per poi convertirla in redshift
      integer :: status
      integer :: getpid
      integer :: system
      integer :: i,m,nlbins
      real(dl) :: redshift, mutest, sigmatest, gammatest, dotsigmatest, dotmutest, dotgammatest


      !allocating arrays
      if (allocated(binned_z) .eqv. .false.) allocate (binned_z(nsteps),binned_mu(nsteps),&
           & binned_sigma(nsteps),derivm(nsteps_derivs),derivs(nsteps_derivs),derivg(nsteps_derivs), red(nsteps_derivs))
      if (allocated(g1) .eqv. .false.) allocate (g1(nsteps), c1(nsteps), d1(nsteps))
      if (allocated(g2) .eqv. .false.) allocate (g2(nsteps_derivs), c2(nsteps_derivs), d2(nsteps_derivs))
      if (allocated(g3) .eqv. .false.) allocate (g3(nsteps_derivs), c3(nsteps_derivs), d3(nsteps_derivs))
      if (allocated(g4) .eqv. .false.) allocate (g4(nsteps), c4(nsteps), d4(nsteps))
      if (allocated(g5) .eqv. .false.) allocate (g5(nsteps_derivs), c5(nsteps_derivs), d5(nsteps_derivs))

      nlbins=(CP%nbmg)-1

      if (debugging) then
         if ((CP%modemg.eq.theta_bin).or.(CP%modemg.eq.smooth_bin)) then
            write(*,*) 'num_bins=',CP%nbmg
            do i=1,CP%nbmg
               write(*,*) 'redshift',i,'=',CP%zbmg(i)
               write(*,*) 'mu',i,'=',CP%mb(i)
               write(*,*) 'sigma', i, '=', CP%sb(i)
            end do
         end if
      end if


      !Gaussian process interface
      if (CP%modemg.eq.GP) then

         !Setting GP redshift to median redshift of each bin
         gpa(1) = (initial_a+CP%abmg(1))/2
         do i=2,CP%nbmg
            gpa(i) = (CP%abmg(i-1)+CP%abmg(i))/2.
         end do
         !write(*,*) 'scale factors di input', CP%abmg
         !write(*,*) 'scale factors mediani', gpa

         final_a   = gpa(CP%nbmg)

         !Creating command line

!RECONSTRUCTION OF MU--------------------------------------------------------------------------------------------------------------------------------------
         !Generate tmp file name based on PID
!         ipid = getpid()
         write (feature_file_mu(11:16), "(Z6.6)") getpid()
         !1. Prepare command and launch it!
         write(a_ini, "(E15.7)"      ) initial_a
         write(a_end, "(E15.7)"      ) final_a
         write(steps_de, "(I10)"     ) nsteps
         write(abins, "(10E15.7)"   ) (gpa(i),i=1,CP%nbmg)
         write(mbins, "(10f15.7)"     ) (CP%mb(i),i=1,CP%nbmg) !python parser struggles with scientific notation negatives: using floats here
         write(mn, "(10f15.7)"     ) CP%mu0
         write(mcorr, "(10E15.7)"  ) CP%mcorr



         if (CP%modemg.eq.GP) then
            if (debugging) write(*,*) 'WORKING WITH GP for mu(z)'
            !here needs the call to script with no baseline

            if (debugging) write(*,*) 'WORKING WITH GP (with baseline) for mu(z)'
        
            command_plus_arguments = "python GP_mu.py --inia "//trim(adjustl(a_ini))//&
            &" --enda "//trim(adjustl(a_end))//" --ODEsteps "//trim(adjustl(steps_de))// &
            & " --scalefactors "//trim(adjustl(abins))// " --mu "//trim(adjustl(mbins))//&
	    & " --mun "//trim(adjustl(mn))//&
            & " --l "//trim(adjustl(mcorr))// " --outfile " // feature_file_mu


            !calling script!!!
            if (debugging) then
               write(*,*) 'Calling Gaussian process script with command line:'
               write(*,*) trim(adjustl(command_plus_arguments))
            end if
            status = system(trim(adjustl(command_plus_arguments)))
            if (status/=0) then
               print *, "Failure in GP reconstruction of mu(z) -- see above."
               call abort
            end if


         end if


        
         !Reading temporary file generated by GP script--------------
         open(unit=17, file=feature_file_mu, action='read')
         do i=1,nsteps
            read(17, "(E15.8, 1X, E15.8, 1X, E15.8)", iostat=status) binned_z(i), binned_mu(i)
            if (status>0) then
               print *, "Error reading the tmp mu(z) file."
               call abort
            end if
         end do
         close(17, status='delete')


         !-----------------------------------------------------------

         !Setting interpolation for GP arrays------------------------
         call newspline(binned_z,binned_mu, g1, c1, d1, nsteps)
         !-----------------------------------------------------------


!RECONSTRUCTION OF SIGMA-----------------------------------------------------------------------------------------------------------------------------------------

         !Generate tmp file name based on PID
!         ipid = getpid()
         write (feature_file_sig(11:16), "(Z6.6)") getpid()
         !1. Prepare command and launch it!
         write(a_ini, "(E15.7)"      ) initial_a
         write(a_end, "(E15.7)"      ) final_a
         write(steps_de, "(I10)"     ) nsteps
         write(abins, "(10E15.7)"   ) (gpa(i),i=1,CP%nbmg)
         write(sbins, "(10f15.7)"     ) (CP%sb(i),i=1,CP%nbmg) !python parser struggles with scientific notation negatives: using floats here
         write(sn, "(10f15.7)"     ) CP%sig0
         write(scorr, "(10E15.7)"  ) CP%scorr

         if (CP%modemg.eq.GP) then
            if (debugging) write(*,*) 'WORKING WITH GP for sigma(z)'
            !here needs the call to script with no baseline

            if (debugging) write(*,*) 'WORKING WITH GP (with baseline) for sigma(z)'
        
            command_plus_arguments = "python GP_sigma.py --inia "//trim(adjustl(a_ini))//&
            &" --enda "//trim(adjustl(a_end))//" --ODEsteps "//trim(adjustl(steps_de))// &
            & " --scalefactors "//trim(adjustl(abins))// " --sigma "//trim(adjustl(sbins))//&
	    & " --sigman "//trim(adjustl(sn))//&
            & " --l "//trim(adjustl(scorr))// " --outfile " // feature_file_sig


            !calling script!!!
            if (debugging) then
               write(*,*) 'Calling Gaussian process script with command line:'
               write(*,*) trim(adjustl(command_plus_arguments))
            end if
            status = system(trim(adjustl(command_plus_arguments)))
            if (status/=0) then
               print *, "Failure in GP reconstruction of sigma(z) -- see above."
               call abort
            end if


         end if


        
         !Reading temporary file generated by GP script--------------
         open(unit=17, file=feature_file_sig, action='read')
         do i=1,nsteps
            read(17, "(E15.8, 1X,E15.8, 1X, E15.8)", iostat=status) binned_z(i), binned_sigma(i)
            if (status>0) then
               print *, "Error reading the tmp sigma(z) file."
               call abort
            end if
         end do
         close(17, status='delete')
         
         !-----------------------------------------------------------

         !Setting interpolation for GP arrays------------------------
         call newspline(binned_z,binned_sigma, g4, c4, d4, nsteps)
         !-----------------------------------------------------------





!------------------------------------------------------

     else if (CP%modemg.gt.3) then
         write(*,*) "THIS MODEL DOESN'T EXIST!!!!"
         stop
     end if

     if (debugging) write(*,*) 'mu(z), sigma(z) computed'

     if (CP%modemg .eq. GP) then
         call get_mu_derivs(CP)
         if (debugging) write(*,*) 'done derivs for mu(z)' 
         call get_sigma_derivs(CP)
         if (debugging) write(*,*) 'done derivs for sigma(z)'    
         call get_gamma_derivs(CP)
         if (debugging) write(*,*) 'done derivs for gamma(z)' 
     end if   


     if (debugging) then
         write(*,*) 'printing mu, sigma'
         write(*,*) 'sono qui'
         open(58,file='test_dotsigma.dat')
         write(*,*) 'sono qui'  
         open(57,file='test_dotmu.dat')
         open(56,file='test_dotgamma.dat')

         do m=1,1001
            redshift=(m-1)*3._dl/1000            
            call get_mu(CP,1/(1+redshift), mutest)
            call get_sigma(CP,1/(1+redshift), sigmatest)
            call get_gamma(CP,1/(1+redshift), gammatest)
            call get_dotmu(1/(1+redshift), dotmutest)
            call get_dotsigma(1/(1+redshift), dotsigmatest)
            call get_dotgamma(1/(1+redshift), dotgammatest)
            write(56,*) redshift, gammatest, dotgammatest
            write(57,*) redshift, mutest, dotmutest
            write(58,*) redshift, sigmatest, dotsigmatest
         end do
         !close(40)
         close(56)
         close(58)     
         close(57)        	
         !stop
     end if

     end subroutine calc_MGfunc

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
  

end module binnedMG
