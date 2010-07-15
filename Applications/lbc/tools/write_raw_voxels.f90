program load_voxel
implicit none
! read in binary voxel file created by voxelizer (Chris Sewell and Dan Morris)

   integer :: num_modifiers,num_columns,data_alloc_chunk,cc,int_ex
   integer :: ii,jj,kk
   integer :: unit1,out_unit,minx,miny,minz,maxx,maxy,maxz
   integer(4) :: header_size,num_objects,object_header_size,voxel_struct_size
   integer(4) :: num_voxels,voxel_resolution(3),resolution_product
   real(4) :: voxel_size(3),model_scale_factor,model_offset(3),zero_coordinate(3)

   integer(1) :: has_texture  !,exists
!   character(1) :: has_texture
   character(1) :: exists,namestr(260)
   character(200) :: filename,out_file
   
!   integer(1) :: 
   integer(2) :: namestr8(260)
   logical    :: fileexists

type datentyp
   integer(2) :: vcoords(3)
   integer(1) :: vhas_texture,vhas_normal,vhas_distance,vnum_mod,is_on_border
!   character(1) :: vhas_texture,vhas_normal,vhas_distance,vnum_mod,is_on_border
   real(4) :: v_uv(2),vnormal(3),vdistance_to_surf,vdistance_grad(3),vdistance_to_mod(5),vmod_grad(15) 
end type datentyp
   type(datentyp) :: daten
   integer(4) :: nx,ny,nz 
   integer(2) :: value
fileexists = .false.
minx =  10000
miny =  10000
minz =  10000
maxx = -10000  
maxy = -10000  
maxz = -10000  

nx=10
ny=10
nz=10


call getarg(1,filename)
filename=adjustl(filename)
out_file=trim(filename)//".raw"
unit1=84
out_unit=85

open(unit=out_unit,file=trim(out_file),access='stream',form='unformatted')
write(out_unit) nx,ny,nz 

do kk=1,nz
do jj=1,ny
do ii=1,nx
write(out_unit) value 
enddo
enddo
enddo


close(out_unit)


end program 
