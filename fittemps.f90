program template_fitting
  use healpix_types
  use pix_tools
  use fitstools
  implicit none

  ! This program will take maps for each band, given a commander run and will subtract off all 
  ! foregrounds excluding the foreground template we want to fit (will make optional). Once this
  ! has been done, we will use a chi-square minimization to find the best fit of the foreground
  ! template for a each band.

  ! Try to make it work for a single map first:

  integer(i4b)        :: i,j,k,l,total,nlheader,band,temp
  integer(i4b)        :: nside,ordering,nmaps,npix,fg,order_map,order_temp
  integer(i4b)        :: nside_fg,ordering_fg,nmaps_fg,npix_fg
  character(len=128)  :: band_file,residual_map,mask_file,temps,maps,offset_file,gain_file
  character(len=128)  :: no_fg_map, fg_map, dust_amp, foreground, version, output
  character(len=2)    :: number
  real(dp)            :: chisq,sum1,sum2,amp,count
  real(dp)            :: nullval
  real(dp)            :: missval = -1.6375d30
  logical(lgt)        :: anynull
  logical(lgt)        :: double_precision
  real(dp), allocatable, dimension(:,:,:) :: fitter
  real(dp), allocatable, dimension(:,:)   :: raw_map,mask,new_map
  real(dp), allocatable, dimension(:,:)   :: ame_temp,cmb_temp,ff_temp,synch_temp,dust_temp
  real(dp), allocatable, dimension(:,:)   :: hcn_temp,co100_temp,co217_temp,co353_temp,temp_map
  real(dp), allocatable, dimension(:)     :: gains,offsets
  character(len=65),  allocatable         :: bands(:)
  character(len=100), allocatable         :: templates(:,:)
  character(len=10),  dimension(9)        :: fgs
  character(len=80),  dimension(180)      :: header

  fgs(1) = 'cmb'
  fgs(2) = 'ame1'
  fgs(3) = 'ff'
  fgs(4) = 'synch'
  fgs(5) = 'dust'
  fgs(6) = 'hcn'
  fgs(7) = 'co-100'
  fgs(8) = 'co-217'
  fgs(9) = 'co-353'

  if  (iargc() /=5) then
     write(*,*) "Usage:"
     write(*,*) 'fittemps [Number of bands] [foreground number] [template band number (ex. "34")] [version tag (ex. "v1")]'
     write(*,*)
     write(*,*) 'Foregrounds: cmb = 1, ame = 2, ff = 3, synch = 4, dust = 5, hcn = 6'
     write(*,*) '             co-1 = 7, co-2 = 8, co-3 = 9 '
     stop
  endif

  call getarg(2, total)
  call getarg(3, fg)
  call getarg(4, band)
  call getarg(5,version)

  output = 'amplitudes/' // trim(version) // '/'

  call system('mkdir -p amplitudes/' // trim(version) // '/dust_amplitudes/')
  call system('mkdir -p amplitudes/' // trim(version) // '/maps/')
  call system('mkdir -p templates/'  // trim(version) // '/' )
  
  temps       = 'templates/' // trim(version) // '/'
  band_file   = 'bands.txt'
  mask_file   = 'masks/fg_mask.fits'
  maps        = 'maps/'
  offset_file = 'offset_'//trim(version)//'.dat'
  gain_file   = 'gains_'//trim(version)//'.dat'

  allocate(bands(total),gains(total),offsets(total))
  allocate(templates(9,total))


  ! Declare the names of the templates in an 9xband array
  do i=1,9
     write(number,10) i
     templates(1,i) = trim(temps) // trim(fgs(1)) // '_band0' // trim(number) // '_k00010.fits'
     templates(2,i) = trim(temps) // trim(fgs(2)) // '_band0' // trim(number) // '_k00010.fits'
     templates(3,i) = trim(temps) // trim(fgs(3)) // '_band0' // trim(number) // '_k00010.fits'
     templates(4,i) = trim(temps) // trim(fgs(4)) // '_band0' // trim(number) // '_k00010.fits'
     templates(5,i) = trim(temps) // trim(fgs(5)) // '_band0' // trim(number) // '_k00010.fits'
     templates(6,i) = trim(temps) // trim(fgs(6)) // '_band0' // trim(number) // '_k00010.fits'
     templates(7,i) = trim(temps) // trim(fgs(7)) // '_band0' // trim(number) // '_k00010.fits'
     templates(8,i) = trim(temps) // trim(fgs(8)) // '_band0' // trim(number) // '_k00010.fits'
     templates(9,i) = trim(temps) // trim(fgs(9)) // '_band0' // trim(number) // '_k00010.fits'
  end do

  do i=10,total
     write(number,11) i
     templates(1,i) = trim(temps) // trim(fgs(1)) // '_band' // trim(number) // '_k00010.fits'
     templates(2,i) = trim(temps) // trim(fgs(2)) // '_band' // trim(number) // '_k00010.fits'
     templates(3,i) = trim(temps) // trim(fgs(3)) // '_band' // trim(number) // '_k00010.fits'
     templates(4,i) = trim(temps) // trim(fgs(4)) // '_band' // trim(number) // '_k00010.fits'
     templates(5,i) = trim(temps) // trim(fgs(5)) // '_band' // trim(number) // '_k00010.fits'
     templates(6,i) = trim(temps) // trim(fgs(6)) // '_band' // trim(number) // '_k00010.fits'
     templates(7,i) = trim(temps) // trim(fgs(7)) // '_band' // trim(number) // '_k00010.fits'
     templates(8,i) = trim(temps) // trim(fgs(8)) // '_band' // trim(number) // '_k00010.fits'
     templates(9,i) = trim(temps) // trim(fgs(9)) // '_band' // trim(number) // '_k00010.fits'
  end do

10 format (I1)
11 format (I2)

  ! Read raw map names, gains, and offsets
  open(30,file=band_file)
  open(31,file=gain_file)
  open(32,file=offset_file)
  do i=1,total
     read(30,fmt='(a)') bands(i)
     read(31,fmt='(F8.6)') gains(i)
     read(32,fmt='(F10.6)') offsets(i)
  end do

  ! Read nside and nmaps from a mask file
  i=getsize_fits(mask_file,nside=nside,ordering=ordering,nmaps=nmaps)
  npix = nside2npix(nside)
  
  ! Allocate all necessary map/template arrays
  allocate(raw_map(0:npix-1,nmaps))
  allocate(new_map(0:npix-1,nmaps))
  allocate(mask(0:npix-1,nmaps))
  allocate(ame_temp(0:npix-1,nmaps))
  allocate(cmb_temp(0:npix-1,nmaps))
  allocate(ff_temp(0:npix-1,nmaps))
  allocate(synch_temp(0:npix-1,nmaps))
  allocate(dust_temp(0:npix-1,nmaps))
  allocate(hcn_temp(0:npix-1,nmaps))
  allocate(co100_temp(0:npix-1,nmaps))
  allocate(co217_temp(0:npix-1,nmaps))
  allocate(co353_temp(0:npix-1,nmaps))
  allocate(temp_map(0:npix-1,nmaps))

  allocate(fitter(9,0:npix-1,nmaps))

  write(*,*) nside
  write(*,*) npix

  nmaps=1

  ! Load the foreground template for fitting
  if (fg == 1) then
     foreground = 'cmb'
     call read_bintab(templates(fg,band), cmb_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else if (fg == 2) then
     foreground = 'ame'
     call read_bintab(templates(fg,band), ame_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else if (fg == 3) then
     foreground = 'ff'
     call read_bintab(templates(fg,band), ff_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else if (fg == 4) then
     foreground = 'synch'
     call read_bintab(templates(fg,band), synch_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else if (fg == 5) then
     foreground = 'dust'
     call read_bintab(templates(fg,band), dust_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else if (fg == 6) then
     foreground = 'hcn'
     call read_bintab(templates(fg,band), hcn_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else if (fg == 7) then
     foreground = 'co100'
     call read_bintab(templates(fg,band), co100_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else if (fg == 8) then
     foreground = 'co217'
     call read_bintab(templates(fg,band), co217_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else if (fg == 9) then
     foreground = 'co353'
     call read_bintab(templates(fg,band), co353_temp, npix, nmaps, nullval, anynull, header=header)
     i=getsize_fits(templates(fg,band),nside=nside,ordering=ordering,nmaps=nmaps)
  else 
     write(*,*) ' Choose a foreground (between 1 and 9) to fit a foreground template!'
     stop
  end if

  fitter(1,:,:) = cmb_temp
  fitter(2,:,:) = ame_temp
  fitter(3,:,:) = ff_temp
  fitter(4,:,:) = synch_temp
  fitter(5,:,:) = dust_temp
  fitter(6,:,:) = hcn_temp
  fitter(7,:,:) = co100_temp
  fitter(8,:,:) = co217_temp
  fitter(9,:,:) = co353_temp

  call read_bintab(mask_file, mask, npix, nmaps, nullval, anynull, header=header)
  nlheader = size(header)

  write(*,*) ' Let the fitting begin!'

  dust_amp    = trim(output) // 'dust_amplitudes/dust_amplitudes.dat'

  ! Begin the template fitting process
  open(35,file=dust_amp)
  do i = 1,9

     ! Initializing the input/output maps and foreground templates
     write(number,10) i
     no_fg_map    = trim(output) // 'maps/nofg_band0'// trim(number) //'.fits'
     fg_map       = trim(output) // 'maps/'//trim(foreground) //'_band0'// trim(number) //'.fits'
     residual_map = trim(output) // 'maps/no_' // trim(foreground) // '_band0'// trim(number) //'.fits'

     call read_bintab(trim(maps) // bands(i), raw_map, npix, nmaps, nullval, anynull, header=header)
     l=getsize_fits(trim(maps) // bands(i),nside=nside,ordering=order_map,nmaps=nmaps)
     nmaps = 1
     if (ordering /= order_map) then
        if (order_map == 1) then
           call convert_ring2nest(nside,raw_map(:,1))
        else
           call convert_nest2ring(nside,raw_map(:,1))
        end if
     end if
     do j = 1,9
        if (j .ne. fg) then
             if (j == 1) then
                call read_bintab(templates(j,i), cmb_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = cmb_temp
             else if (j == 2) then
                call read_bintab(templates(j,i), ame_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = ame_temp
             else if (j == 3) then
                call read_bintab(templates(j,i), ff_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = ff_temp
             else if (j == 4) then
                call read_bintab(templates(j,i), synch_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = synch_temp
             else if (j == 5) then
                call read_bintab(templates(j,i), dust_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = dust_temp
             else if (j == 6) then
                call read_bintab(templates(j,i), hcn_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = hcn_temp
             else if (j == 7) then
                call read_bintab(templates(j,i), co100_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = co100_temp
             else if (j == 8) then
                call read_bintab(templates(j,i), co217_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = co217_temp
             else if (j == 9) then
                call read_bintab(templates(j,i), co353_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = co353_temp
             end if
        end if
     end do
     ! Subtracting off all other foregrounds
     write(*,*) 'Fitting band 0'// trim(number)

     do k = 1, nmaps
        do j = 0, npix-1
           new_map(j,k) = gains(i)*raw_map(j,k)-offsets(i)
           do l=1,9
              if (l .ne. fg) then
                 new_map(j,k) = new_map(j,k) - fitter(l,j,k)
              end if
           end do
        end do
     end do

     call write_bintab(new_map, npix, nmaps, header, nlheader, no_fg_map)

     ! Computes average of dust_template weights (per pixel)
     ! -----------------------------------------------------
     ! Here we will compute a find the weight for the minimized chi-square distribution.
     ! The total chi-square takes the form:
     !
     !    chi^2 = sum ((map - amplitude*template)**2/rms**2)
     !
     ! So to find the minimized value (since there's only one free parameter), 
     ! we take the derivative wrt a (amplitude) which gives us:
     !
     !     sum((temp*map)/rms**2)
     ! a = ----------------------
     !     sum((temp**2)/rms**2)
     !
     ! But because rms is the same at each pixel, we don't need to calculate?
     
     sum1  = 0.d0
     sum2  = 0.d0
     do k=1,nmaps
        do j=0,npix-1
           sum1  = sum1 + dust_temp(j,k)*new_map(j,k) ! *mask(j,k)
           sum2  = sum2 + dust_temp(j,k)**2.d0 ! *mask(j,k)
        end do
     end do

     amp = sum1/sum2

     do k=1,nmaps
        ! Create mock dust map for a given band
        do j=0,npix-1
           temp_map(j,k) = amp*dust_temp(j,k)
        end do
     end do
     write(*,*) 'Band ', i, ' weight = ', amp
     
     write(35,'(I2,F20.8)') i, amp

     call write_bintab(temp_map, npix, nmaps, header, nlheader, fg_map)

     ! Compute chi-square of observed residual map vs mock dust map
     do k=1,nmaps
        chisq=0.d0
        do j=0,npix-1
           chisq = chisq + (new_map(j,k) - temp_map(j,k))**2.d0
           new_map(j,k) = new_map(j,k) - temp_map(j,k)
        end do
     end do
     write(*,*) 'Chi-square for band ', i, ' dust fit = ', chisq
     write(*,*) ''
     call write_bintab(new_map, npix, nmaps, header, nlheader, residual_map)
  end do
  
  do i = 10,total

     ! Initializing the input/output maps and foreground templates
     write(number,11) i
     no_fg_map    = trim(output) // 'maps/nofg_band'// trim(number) //'.fits'
     fg_map       = trim(output) // 'maps/'//trim(foreground) //'_band'// trim(number) //'.fits'
     residual_map = trim(output) // 'maps/no_' // trim(foreground) // '_band'// trim(number) //'.fits'

     call read_bintab(trim(maps) // bands(i), raw_map, npix, nmaps, nullval, anynull, header=header)
     l=getsize_fits(trim(maps) // bands(i),nside=nside,ordering=order_map,nmaps=nmaps)
     nmaps = 1
     if (ordering /= order_map) then
        if (order_map == 1) then
           call convert_ring2nest(nside,raw_map(:,1))
        else
           call convert_nest2ring(nside,raw_map(:,1))
        end if
     end if
     do j = 1,9
        if (j .ne. fg) then
             if (j == 1) then
                call read_bintab(templates(j,i), cmb_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = cmb_temp
             else if (j == 2) then
                call read_bintab(templates(j,i), ame_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = ame_temp
             else if (j == 3) then
                call read_bintab(templates(j,i), ff_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = ff_temp
             else if (j == 4) then
                call read_bintab(templates(j,i), synch_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = synch_temp
             else if (j == 5) then
                call read_bintab(templates(j,i), dust_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = dust_temp
             else if (j == 6) then
                call read_bintab(templates(j,i), hcn_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = hcn_temp
             else if (j == 7) then
                call read_bintab(templates(j,i), co100_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = co100_temp
             else if (j == 8) then
                call read_bintab(templates(j,i), co217_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = co217_temp
             else if (j == 9) then
                call read_bintab(templates(j,i), co353_temp, npix, nmaps, nullval, anynull, header=header)
                fitter(j,:,:) = co353_temp
             end if
        end if
     end do

     ! Subtracting off all other foregrounds
     write(*,*) 'Fitting band '// trim(number)

     do k = 1, nmaps
        do j = 0, npix-1
           new_map(j,k) = gains(i)*raw_map(j,k)-offsets(i)
           do l=1,9
              if (l .ne. fg) then
                 new_map(j,k) = new_map(j,k) - fitter(l,j,k)
              end if
           end do
        end do
     end do

     call write_bintab(new_map, npix, nmaps, header, nlheader, no_fg_map)

     ! Computes average of dust_template weights (per pixel)
     sum1  = 0.d0
     sum2  = 0.d0
     do k=1,nmaps
        do j=0,npix-1
           sum1  = sum1 + dust_temp(j,k)*new_map(j,k) ! *mask(j,k)
           sum2  = sum2 + dust_temp(j,k)**2.d0 ! *mask(j,k)
        end do
     end do

     amp = sum1/sum2

     ! Create mock dust map for a given band
     do k=1,nmaps
        do j=0,npix-1
           temp_map(j,k) = amp*dust_temp(j,k)
        end do
     end do
     write(*,*) 'Band ', i, ' weight = ', amp
     
     write(35,'(I2,F20.8)') i, amp

     call write_bintab(temp_map, npix, nmaps, header, nlheader, fg_map)

     ! Compute chi-square of observed residual map vs mock dust map
     do k=1,nmaps
        chisq=0.d0
        do j=0,npix-1
           chisq = chisq + (new_map(j,k) - temp_map(j,k))**2.d0
           new_map(j,k) = new_map(j,k) - temp_map(j,k)
        end do
     end do
     write(*,*) 'Chi-square for band ', i, ' dust fit = ', chisq
     write(*,*) ''
     call write_bintab(new_map, npix, nmaps, header, nlheader, residual_map)

  end do

  close(35)

end program template_fitting
