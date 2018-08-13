    ! covariance matrix from paper (?)

    module priormg
    use CosmologyTypes    !da cui pesco CAMBParams
    use Likelihood_Cosmology	!da cui pesco LogLike..
    use MatrixUtils		!per avere l'inversione della matrice di covarianza
    use precision
    implicit none
    private

    logical :: debugging=.false.

    !likelihood variables
    type, extends(TCosmoCalcLikelihood) :: mgLikelihood
        real :: prior_n_mu, prior_xi_mu                    !correlation parameters for mu
        real :: prior_n_sigma, prior_xi_sigma              !correlation parameters for sigma
    contains

    procedure :: LogLikeDataParams => mg_LnLike
    end type mgLikelihood

    public mgLikelihood, mgLikelihood_Add
    contains

!...................................................................

    subroutine mgLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(mgLikelihood), pointer :: this
    
    if (Ini%Read_Logical('use_priormg',.false.)) then
       allocate(this)

       this%name= Ini%Read_String('priormg_name')
       !CORRELATION PARAMETERS
	   	this%prior_n_mu        = 1.74
	        this%prior_xi_mu       = 0.31
		this%prior_n_sigma     = 1.7
	        this%prior_xi_sigma    = 0.38

       !AUTOCORRELATION PARAMETERS
       this%needs_background_functions = .true.
       call LikeList%Add(this)   !added to the list of likelihoods
    end if

    end subroutine mgLikelihood_Add

!-------------------------------------------------------------------

    real(mcp) function mg_LnLike(this, CMB, DataParams)
    Class(mgLikelihood) :: this
    Class(CMBParams) CMB
    real(mcp) DataParams(:)                            

    integer :: i,j,d
    real :: chi2_mu, chi2_sigma
    real(dl), dimension(:), allocatable          :: diff_vecm,diff_vecs, mfid, sfid, mui, sigi, gpamg
    real(dl), dimension(:,:), allocatable        :: covmatm,inv_covmatm,covmats,inv_covmats
    real(dl), dimension(:), allocatable          :: autocorrm, autocorrs
    real(dl)                                     :: distancemg, autodistmg

write(*,*) 'sono in priormg'

    if(CMB%modemg.ne.3) then
        d=CMB%numbinsmg
        allocate (diff_vecm(d),diff_vecs(d), mfid(d), sfid(d), mui(d), sigi(d), gpamg(d))
        allocate (covmatm(d,d),inv_covmatm(d,d),autocorrm(d),covmats(d,d),inv_covmats(d,d),autocorrs(d))
        do i=1,CMB%numbinsmg
           mui(i)=CMB%binmu(i)
           sigi(i)=CMB%binsigma(i)         
        end do
        gpamg(1) = (1._dl + CMB%binamg(1))/2
         do i=2,CMB%numbinsmg
            gpamg(i) = (CMB%binamg(i-1)+CMB%binamg(i))/2.
         end do
    else
        d=CMB%numbinsmg+1
        allocate (diff_vecm(d),diff_vecs(d), mfid(d), sfid(d), mui(d), sigi(d), gpamg(d))
        allocate (covmatm(d,d),inv_covmatm(d,d),autocorrm(d),covmats(d,d),inv_covmats(d,d),autocorrs(d))
        mui(1)=CMB%binmu0
        sigi(1)=CMB%binsigma0
        do i=1,CMB%numbinsmg
           mui(i+1)=CMB%binmu(i)
           sigi(i+1)=CMB%binsigma(i)
        end do
        gpamg(1)=1._dl
         do i=1,CMB%numbinsmg
            gpamg(i+1) = (CMB%binamg(i-1)+CMB%binamg(i))/2.
         end do
    end if

    if (debugging) write(*,*) 'il vettore scale factors e', gpamg
    if (debugging) write(*,*) 'la dimensione dei vettori e', d
    if (debugging) write(*,*) 'il vettore v_mu', mui
    if (debugging) write(*,*) 'il vettore v_sigma', sigi
 

    !COMPUTING MEAN OF MU_I AND SIGMA_I VALUES-----------
    mfid(1) = (mui(1)+mui(2))/2
    sfid(1) = (sigi(1)+sigi(2))/2
    do i=2,d
       mfid(i)=(mui(i-1)+mui(i)+mui(i+1))/3
       sfid(i)=(sigi(i-1)+sigi(i)+sigi(i+1))/3
    end do
    mfid(d)= (mui(d-1)+mui(d)*2)/3
    sfid(d)= (sigi(d-1)+sigi(d)*2)/3
    if (debugging) write(*,*) 'il vettore v_mfid', mfid
    if (debugging) write(*,*) 'il vettore v_sfid', sfid

    !COMPUTING ARRAY OF mu_i-mu_fid AND sigma_i-sigma_fid
    do i=1,d
       diff_vecm(i) = mui(i)-mfid(i)
       diff_vecs(i) = sigi(i)-sfid(i)
    end do

    if (debugging) write(*,*) 'il vettore differenza di mu e', diff_vecm
    if (debugging) write(*,*) 'il vettore differenza di sigma e', diff_vecs

    !COMPUTING AUTOCORRELATION
!    ?????????????


    !COMPUTING COV MAT AND ITS INVERSE
       do i=1,d
          do j=1,d
                distancemg = abs(gpamg(i)-gpamg(j))
!            	covmat(i,j) = sqrt(autocorr(i)*autocorr(j))/(1+((distance/this%prior_xi)**this%prior_n))
                covmatm(i,j) = 1._dl/(1+((distancemg/this%prior_xi_mu)**this%prior_n_mu))
                covmats(i,j) = 1._dl/(1+((distancemg/this%prior_xi_sigma)**this%prior_n_sigma))
          end do
       end do

    if (debugging) write(*,*) 'la matrice di mu e'
    if (debugging) write(*,*) covmatm

    if (debugging) write(*,*) 'la matrice di sigma e'
    if (debugging) write(*,*) covmats

    inv_covmatm(:,:) = covmatm(:,:)
    call Matrix_Inverse(inv_covmatm)

    inv_covmats(:,:) = covmats(:,:)
    call Matrix_Inverse(inv_covmats)

    if (debugging) write(*,*) 'la matrice inversa di mu dopo e'
    if (debugging) write(*,*) inv_covmatm

    if (debugging) write(*,*) 'la matrice inversa di sigma dopo e'
    if (debugging) write(*,*) inv_covmats
    
    !COMPUTING CHI2
    chi2_mu = 0._dl
    chi2_sigma = 0._dl

    chi2_mu = dot_product( diff_vecm, MatMul(inv_covmatm,diff_vecm))
    chi2_sigma = dot_product( diff_vecs, MatMul(inv_covmats,diff_vecs))

    mg_Lnlike = chi2_mu/2._dl + chi2_sigma/2._dl  

    if (feedback.gt.0) write(*,*) 'MGPrior Like =', mg_Lnlike

    end function  mg_LnLike


    end module priormg
