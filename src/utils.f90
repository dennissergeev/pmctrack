module utils

use datetime_module

use types, only : wp
use constants, only : missval, pi, deg2rad, fh_track, ra, rkilo
use params, only : proj, halo_r, nx1, nx2, ny1, ny2

implicit none

contains

  subroutine make_nc_file_name(nc_file_name, datadir, prefix, y, m, var_name)

    character(len=*)              , intent(inout) :: nc_file_name
    character(len=*)              , intent(in)    :: datadir
    character(len=*)              , intent(in)    :: prefix
    character(len=*)              , intent(in)    :: var_name
    integer                       , intent(in)    :: y
    integer                       , intent(in)    :: m
    
    write(nc_file_name, '(A,A,A,I4.4,A,I2.2,A,A,A)') trim(datadir), '/',      &
                                                   & trim(prefix),            &
                                                   & y, '.',                  &
                                                   & m, '.',                  &
                                                   & trim(var_name), '.nc'
  end subroutine make_nc_file_name


  subroutine write_vortrack(outdir, vortex, id1, id2, idt)
    
    character(len=*), intent(in) :: outdir
    integer         , intent(in) :: id1, id2
    type(datetime)  , intent(in) :: idt
    real(wp)        , intent(in) :: vortex(5)
    character(len=256)           :: fname_track


    write(fname_track, '(A,A,A,I4.4,A,I4.4,A)') trim(outdir), '/',        &
                                             & 'vortrack_', id1,    &
                                             & '_', id2, '.txt'
    open(unit=fh_track, file=fname_track, form='formatted',                 &
       & access='append', status='unknown')
    write(unit=fh_track, fmt='(3f12.5,A15,f15.5,I3)')                         &
      & vortex(1),                                                            &
      & vortex(2),                                                            &
      & vortex(3),                                                            &
      & trim(idt%strftime('%Y%m%d%H%M')),                                     &
      & vortex(4),                                                            &
      & int(vortex(5))
    close(unit=fh_track)
  end subroutine write_vortrack 


  subroutine apply_mask_2d(var, nx, ny, flag, lon, lat)
  
    integer    , intent (in)    :: nx, ny
    real   (wp), intent (inout) :: var (0:nx, 0:ny)
    real   (wp), intent (in)    :: flag(0:nx, 0:ny)
    real   (wp), intent(in)     :: lon (0:nx      )
    real   (wp), intent(in)     :: lat (      0:ny)
    
    integer                     :: i, j
    integer                     :: ii, jj
    integer                     :: halo_x
    integer                     :: halo_y
    real   (wp)                 :: dist
    real   (wp)                 :: lonin
    real   (wp)                 :: latin
    
  
    lonin = lon(2) - lon(1)
    latin = lat(2) - lat(1)

    if (proj==1) then
      halo_x = nint(halo_r / (ra * lonin * deg2rad * cosd(lat(ny1+ny2/2)) / rkilo))
      halo_y = nint(halo_r / (ra * latin * deg2rad / rkilo))
    elseif (proj==2) then
      halo_x = nint(halo_r / lonin / rkilo)
      halo_y = nint(halo_r / latin / rkilo)
    endif

    do j = 0, ny
      do i = 0, nx
        if (halo_x > 0 .and. halo_y > 0) then
          do jj = max(-halo_y, -j), min(halo_y, ny-j)
            do ii = max(-halo_x, -i), min(halo_x, nx-i)
              if (proj == 1) then
                dist = great_circle(lon(i), lon(i+ii), lat(j), lat(j+jj), ra)
              elseif(proj==2)then
                dist = sqrt((ii * lonin)**2 + (jj * latin)**2)
              endif

              if ((dist <= halo_r * rkilo) &
                & .and. (flag(i+ii, j+jj) == 1.)) then
                var(i, j) = missval
              endif
            enddo
          enddo
        else
          if (flag(i, j) == 1.) then
            var(i, j) = missval
          endif
        endif
      enddo
    enddo
  
  end subroutine apply_mask_2d
  
  
  function integral_p(var, p, nz)
    integer    , intent (in)  :: nz
    real   (wp), intent (in)  :: var(nz)
    real   (wp), intent (in)  :: p  (nz)
    real   (wp)               :: integral_p
    integer                   :: k
    
    integral_p = 0.
  
    do k = 1, nz-1
      integral_p = integral_p + 0.5 * (var(k) + var(k+1)) * (p(k) - p(k+1))
    end do
    integral_p = integral_p / (p(1) - p(nz))
  
  end function integral_p


  function sind(x)
    real(wp) :: x
    real(wp) :: sind
    sind = sin(deg2rad * x)
  end function sind 


  function cosd(x)
    real(wp) :: x
    real(wp) :: cosd
    cosd = cos(deg2rad * x)
  end function cosd 


  function great_circle(lon1, lon2, lat1, lat2, ra)
    ! TODO: Double precision version
    real(wp) :: lon1, lat1, lon2, lat2 ! in degrees
    real(wp) :: ra
    real(wp) :: ang_cos
    real(wp) :: great_circle

    ang_cos = cosd(lat1) * cosd(lat2) * cosd(lon2 - lon1)             &
          & + sind(lat1) * sind(lat2)

    if (abs(ang_cos) <= 1.0_wp) then
      great_circle = ra * acos(ang_cos)
    else 
      great_circle = 0.0_wp
      if ((lon1 /= lon2) .or. (lat1 /= lat2)) then
         write(*, *) 'ValueError: ang_cos is', ang_cos, & 
                   & '; Setting great_circle to 0', lon1, lon2, lat1, lat2
      endif
    endif
  end function


  subroutine check_exist( fname ) ! , ierr)
    ! Check if file exists
    implicit none
  
    character(len=*), intent(in)  :: fname
    ! integer         , intent(out) :: ierr
    logical                       :: f_exist

    inquire( file=trim(fname), exist=f_exist )
    if (.not. f_exist) then
      write(*, '(a)') "ERROR File '"//trim(fname)//"' not found"
      stop
    endif
  end subroutine check_exist


  subroutine makedirs_p(newdirpath, overwrite)
    ! Author:  Jess Vriesema
    ! Date:    Spring 2011
    ! Purpose: Creates a directory at ./newDirPath
  
    implicit none
  
    character(len=*), intent(in)  :: newdirpath
    logical, optional, intent(in) :: overwrite
    character(len=256)            :: mkdir_cmd
    logical                       :: dir_exists
    logical                       :: ovwrt
  
    ! Check if the directory exists first
    ! Works with gfortran, but not ifort
#ifdef gnu
    inquire(file=trim(newdirpath)//'/.', exist=dir_exists)
#else
    ! Works with ifort, but not gfortran
    inquire(directory=newdirpath, exist=dir_exists)
#endif
  
    if (present(overwrite)) then
      ovwrt = overwrite
    else
      ovwrt = .false.
    endif

    if (dir_exists) then
      write (*, '(a)') "Directory already exists: '"//trim(newdirpath)//"'"
      if (ovwrt) then
        write (*, '(a)') "Overwriting it"
        call system('rm -r '//trim(newdirpath))
        mkdir_cmd = 'mkdir -p '//trim(newdirpath)
        call system(mkdir_cmd)
      endif
    else
      mkdir_cmd = 'mkdir -p '//trim(newdirpath)
      write(*, '(a)') "Creating new directory: '"//trim(mkdir_cmd)//"'"
      call system(mkdir_cmd)
    endif
  end subroutine makedirs_p

end module utils