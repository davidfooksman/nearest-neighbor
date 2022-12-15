      program PCcluster
      implicit none
      
      integer PC_tot_num
      parameter (PC_tot_num=1174) ! total number of red cells
      integer PC_tot_num2
      parameter (PC_tot_num2=533) ! the number of green cells
      integer sample_num
      parameter (sample_num=200)
      integer iteration
      parameter (iteration=10) 
      integer NN_num
      parameter (NN_num=10) ! the number of neighbors you calculate the average distance from

      real*8 PC_x(PC_tot_num)
      real*8 PC_y(PC_tot_num)
      real*8 PC_z(PC_tot_num)

      real*8 PC_x2(PC_tot_num2)
      real*8 PC_y2(PC_tot_num2)
      real*8 PC_z2(PC_tot_num2)

      integer i,j,i2,j2,n_t
      real*8 min_dist(PC_tot_num)
      real*8 ave_NN_dist(PC_tot_num)
      real*8 ave_NN_dist2(PC_tot_num2)	  
      real*8 dist_ij
      integer cell_flag(PC_tot_num)
      integer pick_cell_num,pick_cell
      integer cell_flag2(PC_tot_num2)
      integer pick_cell_num2,pick_cell2
      integer NN_index
      real*8 NN_dist(PC_tot_num+PC_tot_num2)
      integer NN_index2
      real*8 NN_dist2(PC_tot_num+PC_tot_num2)	  
      real*8 temp
      integer index,index2
      real*8 sample1_aveNNdist(sample_num)
      real*8 sample2_aveNNdist(sample_num)

      real rand3
      double precision r3
 
      r3=5.0      

cccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>  read the input file

      open(unit=10,file='MoreCluster_RedCells.dat',status='old')
      do i=1,PC_tot_num
         read(10,*) PC_x(i),PC_y(i),PC_z(i)
      enddo
      close(10)

      open(unit=10,file='MoreCluster_GreenCells.dat',status='old')
      do i=1,PC_tot_num2
         read(10,*) PC_x2(i),PC_y2(i),PC_z2(i)
      enddo
      close(10)
    
cccccccccccccccccccccccccccccccccccccccccc

      do n_t=1,iteration  ! multiple iterations of random sampling


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>   randomly pick the same number of samples from red cells and green cells

         do i=1,PC_tot_num
            cell_flag(i)=0
         enddo
         pick_cell_num=0
         
 1978    continue
         
         pick_cell=int(real(PC_tot_num)*rand3(r3))+1
         
         if((pick_cell.gt.0).and.(pick_cell.le.PC_tot_num))then
            if(cell_flag(pick_cell).eq.0)then
               cell_flag(pick_cell)=1
               pick_cell_num=pick_cell_num+1
               if(pick_cell_num.eq.sample_num)then
                  goto 1977
               endif
            endif
         endif
         
         goto 1978
         
 1977    continue

         do i=1,PC_tot_num2
            cell_flag2(i)=0
         enddo
         pick_cell_num2=0
         
 1979    continue
         
         pick_cell2=int(real(PC_tot_num2)*rand3(r3))+1
         
         if((pick_cell2.gt.0).and.(pick_cell2.le.PC_tot_num2))then
            if(cell_flag2(pick_cell2).eq.0)then
               cell_flag2(pick_cell2)=1
               pick_cell_num2=pick_cell_num2+1
               if(pick_cell_num2.eq.sample_num)then
                  goto 1980
               endif
            endif
         endif
         
         goto 1979
         
 1980    continue
         
ccccccccccccccccccccccccccccccccccccccccccccccccc
         
         do i=1,PC_tot_num
            ave_NN_dist(i)=9999.0
            if(cell_flag(i).eq.1)then
               do j=1,PC_tot_num+PC_tot_num2
                  NN_dist(j)=0
               enddo
               NN_index=0
               do j=1,PC_tot_num
                  if(i.ne.j)then
                     NN_index=NN_index+1
                     dist_ij=sqrt((PC_x(i)-PC_x(j))**2
     &                    +(PC_y(i)-PC_y(j))**2
     &                    +(PC_z(i)-PC_z(j))**2)
                     NN_dist(NN_index)=dist_ij
                  endif
               enddo
               do j=1,PC_tot_num2
                  NN_index=NN_index+1
                  dist_ij=sqrt((PC_x(i)-PC_x2(j))**2
     &                 +(PC_y(i)-PC_y2(j))**2
     &                 +(PC_z(i)-PC_z2(j))**2)
                  NN_dist(NN_index)=dist_ij
               enddo
c>>   rank
               do i2=1,NN_index-1
                  do j2=i2+1,NN_index
                     if(NN_dist(j2).lt.NN_dist(i2))then
                        temp=NN_dist(i2)
                        NN_dist(i2)=NN_dist(j2)
                        NN_dist(j2)=temp
                     endif
                  enddo
               enddo
               ave_NN_dist(i)=0
               do i2=1,NN_num
                  ave_NN_dist(i)=ave_NN_dist(i)+NN_dist(i2)
               enddo
               ave_NN_dist(i)=ave_NN_dist(i)/real(NN_num)
            endif
         enddo

cccccccccccccccccccccccccccc

         do i=1,PC_tot_num2
            ave_NN_dist2(i)=9999.0
            if(cell_flag2(i).eq.1)then
               do j=1,PC_tot_num+PC_tot_num2
                  NN_dist2(j)=0
               enddo
               NN_index2=0
               do j=1,PC_tot_num2
                  if(i.ne.j)then
                     NN_index2=NN_index2+1
                     dist_ij=sqrt((PC_x2(i)-PC_x2(j))**2
     &                    +(PC_y2(i)-PC_y2(j))**2
     &                    +(PC_z2(i)-PC_z2(j))**2)
                     NN_dist2(NN_index2)=dist_ij
                  endif
               enddo
               do j=1,PC_tot_num
                  NN_index2=NN_index2+1
                  dist_ij=sqrt((PC_x2(i)-PC_x(j))**2
     &                 +(PC_y2(i)-PC_y(j))**2
     &                 +(PC_z2(i)-PC_z(j))**2)
                  NN_dist2(NN_index2)=dist_ij
               enddo
c>>   rank
               do i2=1,NN_index2-1
                  do j2=i2+1,NN_index2
                     if(NN_dist2(j2).lt.NN_dist2(i2))then
                        temp=NN_dist2(i2)
                        NN_dist2(i2)=NN_dist2(j2)
                        NN_dist2(j2)=temp
                     endif
                  enddo
               enddo
               ave_NN_dist2(i)=0
               do i2=1,NN_num
                  ave_NN_dist2(i)=ave_NN_dist2(i)+NN_dist2(i2)
               enddo
               ave_NN_dist2(i)=ave_NN_dist2(i)/real(NN_num)
            endif
         enddo

         
c>>>  output the calculated distance for the cells you picked

         index=0
         do i=1,PC_tot_num
            if(cell_flag(i).eq.1)then
               index=index+1
               sample1_aveNNdist(index)=ave_NN_dist(i)
            endif
         enddo
         index2=0
         do i=1,PC_tot_num2
            if(cell_flag2(i).eq.1)then
               index2=index2+1
               sample2_aveNNdist(index2)=ave_NN_dist2(i)
            endif
         enddo         
         
         do i=1,sample_num
            print*,n_t,i,sample1_aveNNdist(i),sample2_aveNNdist(i)
         enddo
         
      enddo
      
      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
c>>   random number generator
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc        
