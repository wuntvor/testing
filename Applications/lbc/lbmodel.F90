module LBMODEL
   use nrt_lib
    
   ! Lattice Size
   ! Speed configuration 
   !      7   3   6
   !        \ | /
   !      4 - 1 - 2
   !        / | \
   !      8   5   9
    

#ifdef D2Q9
   integer, parameter                  :: nnod= 9

   real(R8B), dimension(9),parameter    :: t   = (/ 4._R8B/9._R8B, &
             1._R8B/9._R8B, 1._R8B/9._R8B, 1._R8B/9._R8B, 1._R8B/9._R8B, &
             1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B/) 
   integer, dimension(9,3) :: cx  
   integer, dimension(9) :: opp 

   data cx(:,1) / 0,1,0,-1,0,1,-1,-1,1/
   data cx(:,2) / 0,0,1,0,-1,1,1,-1,-1/
   data cx(:,3) / 0,0,0,0, 0,0,0, 0, 0/

   integer, parameter :: I__0= 1,I__E= 2,I__N= 3,I__W= 4,I__S= 5
   integer, parameter :: I_NE= 6,I_NW= 7,I_SW= 8,I_SE= 9

   ! 1  0     6  NE  
   ! 2  E     7  NW
   ! 3  N     8  SW
   ! 4  W     9  SE
   ! 5  S               
         

#ifdef SLIP_BB
   data opp /1,2,5,4,3,9,8,7,6 /
#else /* SLIP BB*/
   data opp /1,4,5,2,3,8,9,6,7 /
#endif /* SLIP BB */
#endif
#ifdef D3Q19
   integer, parameter                  :: nnod= 19
   real(R8B), dimension(19),parameter    :: t   = (/ 1._R8B/3._R8B, &
             1._R8B/18._R8B,1._R8B/18._R8B,1._R8B/18._R8B,1._R8B/18._R8B,1._R8B/18._R8B,1._R8B/18._R8B,& 
             1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B,&
             1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B,1._R8B/36._R8B/)
   integer, dimension(19,3)  :: cx 
   integer, dimension(19)    :: opp  


               ! 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
   data cx(:,1) /0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0/ 
   data cx(:,2) /0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1/
   data cx(:,3) /0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1/

   integer, parameter :: I__0= 1,I__E= 2,I__W= 3,I__N= 4,I__S= 5,I__T= 6,I__B= 7
   integer, parameter :: I_NE= 8,I_SW= 9,I_SE=10,I_NW=11,I_TE=12,I_BW=13,I_BE=14
   integer, parameter :: I_TW=15,I_TN=16,I_BS=17,I_BN=18,I_TS=19

   ! According to ILBDC ordering
   ! 1  0           8  NE y+x+         15  TW z+x-
   ! 2  E x+        9  SW y-x-         16  TN z+y+
   ! 3  W x-       10  SE y-x+         17  BS z-y-
   ! 4  N y+       11  NW y+x-         18  BN z-y+
   ! 5  S y-       12  TE z+x+         19  TS z+y-
   ! 6  T z+       13  BW z-x-      
   ! 7  B z-       14  BE z-x+      

   data opp /1,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18/

#endif
end module LBMODEL
!------------------------------------------------------------------------



