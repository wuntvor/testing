program load_voxel
implicit none
! read in binary voxel file created by voxelizer (Chris Sewell and Dan Morris)

   integer :: num_modifiers,num_columns,data_alloc_chunk,ii,jj,cc,int_ex
   integer :: unit1,out_unit,minx,miny,minz,maxx,maxy,maxz
   integer(4) :: header_size,num_objects,object_header_size,voxel_struct_size
   integer(4) :: num_voxels,voxel_resolution(3),resolution_product
   real(4) :: voxel_size(3),model_scale_factor,model_offset(3),zero_coordinate(3)

   integer(1) :: has_texture  !,exists
!   character(1) :: has_texture
   character(1) :: exists,namestr(260)
   
!   integer(1) :: 
   integer(2) :: namestr8(260)

type datentyp
   integer(2) :: vcoords(3)
   integer(1) :: vhas_texture,vhas_normal,vhas_distance,vnum_mod,is_on_border
!   character(1) :: vhas_texture,vhas_normal,vhas_distance,vnum_mod,is_on_border
   real(4) :: v_uv(2),vnormal(3),vdistance_to_surf,vdistance_grad(3),vdistance_to_mod(5),vmod_grad(15) 
end type datentyp
   type(datentyp),allocatable :: daten(:)
!   type(datentyp) :: daten

minx =  10000
miny =  10000
minz =  10000
maxx = -10000  
maxy = -10000  
maxz = -10000  

unit1=84
out_unit=85
open(unit=out_unit,file='obstacles')
write(out_unit,*) "This is the header"


open(unit=unit1,file='data.vxl',access='stream',form='unformatted')
read(unit1) header_size,num_objects,object_header_size,voxel_struct_size
read(unit1) num_voxels,voxel_resolution,voxel_size,model_scale_factor,model_offset,zero_coordinate,has_texture 
read(unit1) namestr 
write(*,*) "header_size",header_size
write(*,*) "num_obj    ",num_objects
write(*,*) "obj_hd_size",object_header_size
write(*,*) "voxels_size",voxel_struct_size
write(*,*) num_voxels,voxel_resolution,voxel_size,model_scale_factor,model_offset,zero_coordinate
!write(*,*) "texture?: ",has_texture 
!write(*,*) "string  : ",namestr

!do ii=1,260
!if(namestr(ii) < 0) then
!   namestr8(ii) = 256+namestr(ii)
!else
!   namestr8(ii) = namestr(ii)
!endif
!enddo
num_modifiers = 5
num_columns   = 17+num_modifiers*4

resolution_product = voxel_resolution(1)*voxel_resolution(2)*voxel_resolution(3)
!resolution_product = 10 

!allocate(daten(resolution_product))

write(*,*) "resolution product",resolution_product
ii=0
do cc = 1,resolution_product
   read(unit1) exists
   if(exists=='2') then
      ii=ii+1
      read(unit1) daten(ii)%vcoords,daten(ii)%vhas_texture,daten(ii)%v_uv,daten(ii)%vhas_normal,daten(ii)%vnormal
      read(unit1) daten(ii)%vhas_distance,daten(ii)%vdistance_to_surf,daten(ii)%vdistance_grad,daten(ii)%vnum_mod
      do jj=1,num_modifiers
         read(unit1) daten(ii)%vdistance_to_mod(jj)
      enddo
      do jj=1,num_modifiers*3
         read(unit1) daten(ii)%vmod_grad(jj)
      enddo
      read(unit1) daten(ii)%is_on_border
      write(out_unit,*) daten(ii)%vcoords(1)+1,daten(ii)%vcoords(2)+1,daten(ii)%vcoords(3)+1
      
      ! Update boounding box
      if(daten(ii)%vcoords(1) > maxx) maxx=daten(ii)%vcoords(1)
      if(daten(ii)%vcoords(2) > maxy) maxy=daten(ii)%vcoords(2)
      if(daten(ii)%vcoords(3) > maxz) maxz=daten(ii)%vcoords(3)
      if(daten(ii)%vcoords(1) < minx) minx=daten(ii)%vcoords(1)
      if(daten(ii)%vcoords(2) < miny) miny=daten(ii)%vcoords(2)
      if(daten(ii)%vcoords(3) < minz) minz=daten(ii)%vcoords(3)

      ! Get status
      if(mod(real(resolution_product/100.)
   endif
enddo

write(*,*) "done."
write(*,*) "bounding box"
write(*,*) "x:  ",minx,maxx 
write(*,*) "y:  ",miny,maxy 
write(*,*) "z:  ",minz,maxz 


close(unit1)
close(out_unit)


end program load_voxel
