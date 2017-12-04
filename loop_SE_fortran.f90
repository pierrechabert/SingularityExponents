  
      subroutine sh_treat(im,lst5)
!im=sst
      implicit none
        
!f2py intent(in) sst
!f2py intent(inout) lst5
            
      real(kind=8),dimension(:,:),intent(in) :: im

      integer :: i,j,k,l,t1,t2,som,indx
      real(kind=8)::n,m
      real(kind=8),dimension(size(im,1),size(im,2)),intent(inout)::lst5

   
      real(kind=8),dimension(size(im,1),size(im,2),5) :: lst1
      real(kind=8),dimension(size(im,1),size(im,2),5)::lst2
      real(kind=8),dimension(size(im,1),size(im,2),2):: lst3
      real(kind=8),dimension(size(im,1),size(im,2)) :: lst4
      real(kind=8),dimension(size(im,1),size(im,2)) :: lst41
      complex(kind=8) :: compl0 = (0., 1.),compl,zer=(0.,0.)&
           ,un=(1.,0.)

      complex(kind=8),dimension(5) :: dxhat,dyhat,dxhatinv,dyhatinv,&
           si,pi,p0
      complex(kind=8),dimension(5,1)::r
      complex(kind=8):: Ax = (0.,0.)
      complex(kind=8):: Ay = (0.,0.)
      real(kind=8),dimension(1,1)::gx,gy,r0,rox,roy,htemp
      real(kind=8)::gx0,gy0,LCSM,LCSM0,s,LCSM01
      real(kind=8)::gx1,gy1
      complex(kind=8)::rox1,roy1
      
      compl = compl0*sqrt(3.0)

      dxhat = [zer,compl,-compl,zer,zer]
      dyhat = [zer,zer,zer,compl,-compl]
      dxhatinv = [zer,-un/compl,un/compl,zer,zer]
      dyhatinv = [zer,zer,zer,-un/compl,un/compl]
      
      n = size(im,1)
      m = size(im,2)
      zer = 0.
      un = 1.

      indx = 0
      do i = 1,n-1
         do j = 1,n-1
            si = [im(i,j),im(i+1,j),im(i-1,j),im(i,j+1),im(i,j-1)]
            lst1(i,j,:) = si
            s = (1.0/3.0)*sum(si)
            pi = [im(i,j)+s,im(i+1,j)-s,im(i-1,j)-s,im(i,j+1)-s,im(i,j-1)-s]
            lst2(i,j,:) = pi
            
            p0 = pi
            gx = imagpart(matmul(reshape(pi,(/1,5/)),reshape(dxhat,(/5,1/))))
            gx1 = gx(1,1)
            gy = imagpart(matmul(reshape(pi,(/1,5/)),reshape(dyhat,(/5,1/))))
            gy1 = gy(1,1)
            lst3(i,j,:) = cmplx(gx1,gy1)

            r = matmul(reshape(dxhatinv,(/5,1/)),gx) + &
                  matmul(reshape(dyhatinv,(/5,1/)),gy)
            rox = matmul(reshape(r,(/1,5/)),reshape(dxhat,(/5,1/)))
            roy = matmul(reshape(r,(/1,5/)),reshape(dyhat,(/5,1/)))
            rox1 = rox(1,1)
            roy1 = roy(1,1)
            LCSM = real(sqrt((Ax-rox1)**2+(Ay-roy1)**2))
            lst4(i,j) = LCSM
         end do
      end do


      lst41(:,:) = 0
      som = 0
      do t1 = 1,n-1
         do t2 = 1,n-1
            if (lst4(t1,t2)>-10 .and. lst4(t1,t2)<10) then
            lst41(t1,t2) = lst4(t1,t2)
            som = som + 1
            end if
         end do
      end do
    
      LCSM01 = sum(lst41)/som       
      r0 = 1/sqrt(n*m)
      print*, 'LCSM0', LCSM0, '- r0',r0,'LCSM011',LCSM01
      do k = 1,n-1
         do l = 1,n-1
            htemp = log10(lst4(k,l)/LCSM01)/log10(r0)
            lst5(k,l) = htemp(1,1)
         end do
      end do
      end subroutine sh_treat
