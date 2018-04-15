program main
  implicit none 
  integer,parameter  ::nx=116,ny=140
  real (4)::vor(0:nx,0:ny),vor_smth(0:nx,0:ny),vor_part(0:nx,0:ny)
  real (4)::steering_u(0:nx,0:ny),steering_v(0:nx,0:ny),psea(0:nx,0:ny)
  
  open (10,file="vor_out_0001.dat",form="unformatted",access='sequential')
  
  read(10)vor
  read(10)vor_smth
  read(10)vor_part
  read(10)steering_u
  read(10)steering_v
  read(10)psea

  stop
end program main
  
